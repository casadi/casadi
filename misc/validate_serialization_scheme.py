"""Check native debug serialization against the committed decoding layouts.

This test utility validates structure; it neither reconstructs CasADi objects nor
imports a reader package. Only values used by layout expressions are retained.
Opaque payloads are checked for extent and skipped without decoding a second copy.
"""
import json
from pathlib import Path
import re
import struct

SCHEME_PATH = Path(__file__).with_name('serialization_scheme.json')


class SchemeValidator:
    def __init__(self, scheme=None):
        self.scheme = json.loads(SCHEME_PATH.read_text()) if scheme is None else scheme

    def validate(self, text):
        return _Validation(self.scheme, text).run()


class _Validation:
    def __init__(self, scheme, text):
        self.scheme, self.text = scheme, text
        self.pos, self.fields, self.depth, self.skipped = 0, 0, 0, 0
        self.shared, self.context = 0, []
        self.debug = False
        self.referenced_fields = set()
        self.dynamic_fields = False
        def references(value):
            if isinstance(value, list):
                if len(value)==2 and value[0]=='field':
                    if isinstance(value[1], str): self.referenced_fields.add(value[1])
                    else: self.dynamic_fields = True
                for child in value: references(child)
            elif isinstance(value, dict):
                for child in value.values(): references(child)
        references(scheme['reader'])
        if len(text) % 2: self.fail('Odd-length serialization')
        self.length = len(text)//2

    def fail(self, message):
        context = ' / '.join(self.context)
        raise ValueError(f'Serialization scheme mismatch at byte {self.pos} ({context}): {message}')

    def take(self, n):
        if not isinstance(n, int) or n < 0 or n > self.length-self.pos:
            self.fail('Truncated or invalid payload length')
        start = self.pos
        self.pos += n
        return start

    def raw(self, n):
        start = self.take(n)
        out = bytearray(n)
        for i in range(n):
            a, b = (ord(c)-97 for c in self.text[2*(start+i):2*(start+i)+2])
            if not (0 <= a < 16 and 0 <= b < 16): self.fail('Invalid nibble encoding')
            out[i] = a | (b << 4)
        return out

    def tag(self, tag):
        if self.debug and self.raw(1)[0] != ord(tag): self.fail('Expected decoration '+tag)

    def number(self, kind):
        tag, fmt = {'int': ('i', 'i'), 'unsignedint': ('u', 'I'),
                    'casadi_int': ('J', 'q'), 'size_t': ('K', 'Q'),
                    'double': ('d', 'd')}[kind]
        self.tag(tag)
        return struct.unpack('='+fmt, self.raw(struct.calcsize('='+fmt)))[0]

    def size(self, n):
        # Every element must consume at least a byte in debug encoding.
        if not isinstance(n, int) or not 0 <= n <= self.length-self.pos:
            self.fail('Invalid collection length')
        return n

    def string(self, retain=True):
        self.tag('s')
        n = self.number('int')
        if retain:
            try: return self.raw(n).decode('utf8')
            except UnicodeDecodeError: self.fail('Expected UTF-8 metadata')
        self.take(n)
        self.skipped += n

    def name(self, expected):
        actual = self.string()
        if actual != expected: self.fail(f'Expected field {expected!r}, got {actual!r}')
        self.fields += 1

    def expression(self, operand, scope):
        op = operand[0]
        if op == 'literal': return operand[1]
        if op in ('field', 'parameter'):
            name = operand[1] if isinstance(operand[1], str) else self.expression(operand[1], scope)
            if name not in scope: self.fail('Missing scheme field or parameter: '+name)
            return scope[name]
        if op == 'length': return len(self.expression(operand[1], scope))
        a = self.expression(operand[1], scope)
        if op == 'not': return not a
        b = self.expression(operand[2], scope)
        if op == 'or': return bool(a or b)
        if op == 'and': return bool(a and b)
        if op == 'equal': return a == b
        if op == 'not_equal': return a != b
        if op == 'concat': return str(a)+str(b)
        self.fail('Unknown scheme operand: '+op)

    def program(self, steps, scope):
        self.depth += 1
        if self.depth > 256: self.fail('Layout nesting limit exceeded')
        try:
            for step in steps:
                op = step['op']
                if op == 'field':
                    name = step.get('name')
                    if name is None: name = self.expression(step['name_expression'], scope)
                    self.context.append(name)
                    self.name(name)
                    kind = step.get('type')
                    if not kind: self.fail('Unresolved field type')
                    retain = self.dynamic_fields or name in self.referenced_fields
                    value = self.value(kind, retain)
                    if retain: scope[name] = value
                    self.context.pop()
                elif op == 'version':
                    name = step.get('name')
                    if name is None: name = self.expression(step['name_expression'], scope)
                    self.name(name+'::serialization::version')
                    version = self.number('int')
                    if version != step['value']:
                        self.fail(f"{name}: expected version {step['value']}, got {version}")
                elif op == 'call':
                    name = step['layout']
                    self.context.append(name)
                    layout = self.scheme['reader']['layouts'].get(name)
                    if layout is None: self.fail('Missing layout '+name)
                    scope.update(step.get('params', {}))
                    self.program(layout, scope)
                    self.context.pop()
                elif op == 'if':
                    self.program(step['body'] if self.expression(step['condition'], scope)
                                 else step.get('else', []), scope)
                elif op == 'repeat':
                    for _ in range(self.size(self.expression(step['count'], scope))):
                        self.program(step['body'], scope)
                elif op == 'select':
                    tag = self.expression(step['value'], scope)
                    key = str(tag).lower() if isinstance(tag, bool) else str(tag)
                    body = step['cases'].get(key)
                    if body is None: self.fail('Missing discriminator case '+key)
                    self.program(body, scope)
                else:
                    self.fail('Unsupported layout instruction: '+step.get('reason', op))
        finally:
            self.depth -= 1

    def value(self, kind, retain=False):
        if isinstance(kind, dict):
            shape = kind['name']
            if shape == 'std::vector':
                self.tag('V'); n = self.size(self.number('casadi_int'))
                values = [] if retain else None
                for _ in range(n):
                    value = self.value(kind['element'], retain)
                    if retain: values.append(value)
                return values
            if shape == 'std::map':
                self.tag('D')
                for _ in range(self.size(self.number('casadi_int'))):
                    self.value(kind['key']); self.value(kind['value'])
                return
            if shape == 'std::pair':
                self.tag('p')
                return [self.value(kind['first'], retain), self.value(kind['second'], retain)]
            self.fail('Unknown container kind: '+shape)

        if kind in ('int', 'unsignedint', 'casadi_int', 'size_t', 'double'):
            return self.number(kind)
        if kind == 'char': return self.raw(1)[0]
        if kind == 'bool':
            self.tag('b')
            value = self.raw(1)[0]
            if value > 1: self.fail('Invalid boolean')
            return bool(value)
        if kind == 'std::string': return self.string(retain)
        if kind in ('std::istream', 'std::stringstream'):
            self.tag('B')
            n = self.number('size_t')
            self.take(n)
            self.skipped += n
            return
        if kind == 'Dict': return self.value({'name':'std::map', 'key':'std::string', 'value':'GenericType'}, retain)
        definition = self.scheme['reader']['types'].get(kind)
        if definition is None: self.fail('Missing type '+kind)
        if definition.get('decoration'): self.tag(definition['decoration'])
        if definition.get('shared'):
            self.name('Shared::flag')
            flag = self.raw(1)[0]
            if flag == ord('r'):
                self.name('Shared::reference')
                index = self.number('casadi_int')
                # The native map is keyed by pointer; null objects can be shared
                # across different wrapper types. Validate the index, not its type.
                if not 0 <= index < self.shared:
                    self.fail('Invalid shared reference')
                return
            if flag != ord('d'): self.fail('Invalid shared definition')
        self.program(definition['body'], {})
        if definition.get('shared'): self.shared += 1

    def run(self):
        if self.number('casadi_int') != self.scheme['wire']['magic']: self.fail('Invalid magic')
        if self.number('casadi_int') != self.scheme['wire']['protocol']: self.fail('Invalid protocol')
        if self.raw(1)[0] != 1: self.fail('Validation requires debug serialization')
        self.debug = True
        # Function.serialize() writes the body directly, without FileSerializer
        # framing or the decoration/shared marker of a nested Function value.
        self.program(self.scheme['reader']['types']['Function']['body'], {})
        if self.pos != self.length: self.fail('Trailing serialized data')
        return {'fields': self.fields, 'shared_objects': self.shared,
                'skipped_payload_bytes': self.skipped}
