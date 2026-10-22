#!/usr/bin/env python3
"""Extract release serialization metadata from C++ without building CasADi.

The output describes decoding, not the C++ source used to derive it.
Source bodies, locations and extraction records stay internal to this generator.
Unresolved types are explicit; consumers must reject layouts they cannot read.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess

TOKEN = re.compile(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[{}()]')
METHOD = re.compile(r'\b(?P<owner>[A-Za-z_]\w*(?:<[^;{}]*?>)?(?:::[A-Za-z_]\w*)*)\s*::\s*'
                    r'(?P<method>serialize_body|serialize_type|serialize|delayed_serialize_members)\s*'
                    r'\(\s*SerializingStream\s*&\s*(?P<stream>\w+)\s*\)\s*(?:const\s*)?\{')


def mask_comments(text):
    def replace(m):
        value = m.group()
        return re.sub(r'[^\n]', ' ', value) if value.startswith(('//', '/*')) else value
    return mask_alternatives(TOKEN.sub(replace, text))


def mask_alternatives(text):
    """Keep the first branch of preprocessor conditionals; a serializer writes one of them."""
    lines, depth, masked = text.split('\n'), [], []
    for line in lines:
        directive = re.match(r'\s*#\s*(if|ifdef|ifndef|elif|else|endif)\b', line)
        if directive:
            kind = directive[1]
            if kind.startswith('if'): depth.append(False)
            elif kind == 'endif': depth.pop()
            else: depth[-1] = True
            line = ''
        elif any(depth):
            line = ''
        masked.append(line)
    return '\n'.join(masked)


def balanced(text, start, opening='{', closing='}'):
    depth = 0
    for match in TOKEN.finditer(text, start):
        token = match.group()
        if token == opening:
            depth += 1
        elif token == closing:
            depth -= 1
            if depth == 0:
                return match.start()
    raise ValueError('Unbalanced C++ at offset ' + str(start))


def split_args(text):
    result, start, depth = [], 0, 0
    # Braced initializers and template types also contain commas.
    tokens = re.finditer(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[(){}<>,]', text)
    for match in tokens:
        token = match.group()
        if token in ('(', '{', '<'):
            depth += 1
        elif token in (')', '}', '>'):
            depth -= 1
        elif token == ',' and depth == 0:
            result.append(text[start:match.start()].strip())
            start = match.end()
    result.append(text[start:].strip())
    return result


def declarations(text):
    """Collect simple member/local declarations, retaining ambiguous candidates."""
    result = {}
    pattern = re.compile(r'(?:^|[;{}\n])\s*(?:mutable\s+|const\s+)?'
                         r'((?:std::)?[A-Za-z_]\w*(?:\s+[A-Za-z_]\w*)?(?:\s*<[^;{}\n]+>)?)'
                         r'\s+([A-Za-z_]\w*(?:\s*,\s*[A-Za-z_]\w*)*)\s*(?:=[^;\n]*)?;')
    for match in pattern.finditer(text):
        cpp_type = re.sub(r'\s+', ' ', match[1]).strip()
        if cpp_type in ('return', 'using', 'typedef') or cpp_type.startswith(('return ', 'using ', 'typedef ')):
            continue
        for declaration in split_args(text[match.start(2):match.end()-1]):
            name = re.match(r'\s*([A-Za-z_]\w*)', declaration)
            if name:
                result.setdefault(name[1], set()).add(cpp_type)
    return result


def infer_type(expression, decls):
    cpp_type = infer_cpp_type(expression, decls)
    # Normalize serialized scalar aliases, including those inside containers.
    aliases = {'T1': 'double', 'libmad_int': 'casadi_int', 'fmi2ValueReference': 'unsigned int',
               'fmi2Real': 'double', 'fmi2Integer': 'int', 'fmi2Boolean': 'int',
               'fmi3ValueReference': 'unsigned int', 'fmi3Float64': 'double',
               'fmi3Int32': 'int', 'fmi3Boolean': 'bool'}
    return re.sub(r'\b\w+\b', lambda m: aliases.get(m[0], m[0]), cpp_type) if cpp_type is not None else None


STRUCTS = {}


def struct_member_type(expression, decls):
    chain = re.fullmatch(r'([A-Za-z_]\w*)((?:(?:\.|->)[A-Za-z_]\w*)+)', expression)
    if not chain:
        return None
    cpp_type = infer_cpp_type(chain[1], decls)
    for member in re.findall(r'\w+', chain[2]):
        owner = re.sub(r'<.*>|\b(?:struct|const)\b|[&\s]', '', cpp_type or '')
        candidates = STRUCTS.get(owner, {}).get(member, set())
        cpp_type = next(iter(candidates)) if len(candidates) == 1 else None
    return cpp_type


def infer_cpp_type(expression, decls):
    member_type = struct_member_type(expression, decls)
    if member_type:
        return member_type
    cast = re.match(r'static_cast<([^>]+)>', expression)
    if cast:
        return cast[1]
    if expression in ('true', 'false'):
        return 'bool'
    if re.fullmatch(r"'(?:\\.|[^'\\])'", expression):
        return 'char'
    if expression.startswith('"'):
        return 'std::string'
    if expression.endswith('.size()'):
        return 'size_t'
    if expression in ('class_name()', 'serialize_base_function()'):
        return 'std::string'
    if expression.startswith('std::string('):
        return 'std::string'
    constructed = re.match(r'(std::\w+<.*>)[{(]', expression)
    if constructed:
        return constructed[1]
    member = re.search(r'(?:\.|->)([A-Za-z_]\w*)$', expression)
    variable = member[1] if member else expression.removesuffix('()')
    candidates = decls.get(variable, set())
    if len(candidates) == 1:
        return next(iter(candidates))
    for container in ('std::vector', 'std::map', 'std::pair'):
        if candidates and all(t.startswith(container+'<') for t in candidates):
            return container+'<type-dependent>'
    return None


def operations(text):
    match = re.search(r'enum Operation\s*\{', text)
    body = text[match.end():balanced(text, match.end()-1)]
    value, result = -1, {}
    for item in split_args(body):
        if not item:
            continue
        name, *assignment = item.split('=')
        value = int(assignment[0], 0) if assignment else value + 1
        result[name.strip()] = value
    return result


def scheme_operand(text):
    """Compile the scheme's small expression vocabulary to explicit operands."""
    text = text.strip()
    if text in ('true', 'false'):
        return ['literal', text == 'true']
    if re.fullmatch(r'"(?:[^"\\]|\\.)*"', text):
        return ['literal', json.loads(text)]
    if re.fullmatch(r'-?\d+', text):
        return ['literal', int(text)]
    for token, op in [('||', 'or'), ('&&', 'and'), ('==', 'equal'), ('!=', 'not_equal'),
                      ('>=', 'greater_equal'), ('<=', 'less_equal'), ('>', 'greater'), ('<', 'less'),
                      ('+', 'concat')]:
        # Do not split operators inside string literals.
        quoted = escaped = False
        for i, char in enumerate(text):
            if escaped: escaped = False; continue
            if char == '\\' and quoted: escaped = True; continue
            if char == '"': quoted = not quoted
            if quoted or not text.startswith(token, i): continue
            # A lone < or > is neither part of <=, >=, -> nor of a two-character operator.
            if token in ('<', '>') and (text[i+1:i+2] == '=' or text[i-1:i] in '-<>='): continue
            return [op, scheme_operand(text[:i]), scheme_operand(text[i+len(token):])]
    if text.startswith('!'):
        return ['not', scheme_operand(text[1:])]
    if re.fullmatch(r'[\w.]+(?:\(\))?', text):
        return ['binding', text]
    raise ValueError('Unsupported scheme expression: '+text)


def explicit_scheme_types(data):
    """Instantiate serializer templates and describe containers structurally."""
    layouts = data['layouts']
    for name in list(layouts):
        if name.startswith('XFunction<MXFunction,MX,MXNode>') or name.startswith('XFunction<SXFunction,SX,SXNode>'):
            matrix = 'MX' if name.startswith('XFunction<MXFunction') else 'SX'
            generic = 'XFunction<DerivedType,MatType,NodeType>::'+name.rsplit('::', 1)[1]
            layouts[name] = json.loads(re.sub(r'\bMatType\b', matrix, json.dumps(layouts[generic])))
    matrix = layouts.pop('Matrix<Scalar>::serialize', None)
    if matrix:
        for scalar in ('double', 'SXElem', 'casadi_int'):
            layouts['Matrix<'+scalar+'>::serialize'] = json.loads(
                re.sub(r'\bScalar\b', scalar, json.dumps(matrix)))
    # These generic implementation helpers are superseded by concrete dispatch programs.
    for name in list(layouts):
        if name.startswith('XFunction<DerivedType,MatType,NodeType>') or name == 'GenericTypeInternal::serialize':
            del layouts[name]

    def descriptor(text):
        kind = re.sub(r'\s+', '', text).replace('std::size_t', 'size_t')
        if re.search(r'\b(Scalar|MatType|T)\b', kind):
            raise ValueError('Unresolved serializer type: '+kind)
        if '<' not in kind: return kind
        base, rest = kind.split('<', 1)
        args = split_args(rest[:-1])
        fields = {'std::vector': ('element',), 'std::map': ('key', 'value'), 'std::pair': ('first', 'second')}
        if base not in fields or len(args) != len(fields[base]):
            raise ValueError('Unknown serialized container: '+kind)
        return {'name': base, **{key: descriptor(arg) for key, arg in zip(fields[base], args)}}

    def visit(value):
        if isinstance(value, list):
            for child in value: visit(child)
        elif isinstance(value, dict):
            if value.get('op') == 'field' and value.get('type'):
                value['type'] = descriptor(value['type'])
            if 'params' in value:
                value['params'] = {k:v for k,v in value['params'].items() if k not in ('MatType','Scalar','T')}
                if not value['params']: del value['params']
            for child in value.values(): visit(child)
    visit(data)
    for key in ('file_types', 'file_prefixes'):
        data[key] = {tag: descriptor(kind) for tag, kind in data.get(key, {}).items()}
    return data


def lower_scheme_operands(data):
    """Resolve member references to serialized fields before publishing the scheme."""
    import copy
    raw = copy.deepcopy(explicit_scheme_types(data))
    candidates, parameters = {}, set()

    def inventory(value):
        if isinstance(value, list):
            for child in value: inventory(child)
        elif isinstance(value, dict):
            if value.get('op') == 'field' and value.get('bind') and value.get('name'):
                candidates.setdefault(value['bind'], set()).add(value['name'])
            parameters.update(value.get('params', {}))
            for child in value.values(): inventory(child)
    inventory(raw)
    initial = {key: ['field', next(iter(names))] for key, names in candidates.items() if len(names) == 1}
    initial.update({key: ['parameter', key] for key in parameters})

    def resolve(operand, scope):
        op = operand[0]
        if op == 'literal': return operand
        if op == 'binding':
            key = operand[1]
            if key in scope: return scope[key]
            if key.endswith('.size()') and key[:-7] in scope:
                return ['length', scope[key[:-7]]]
            raise ValueError('No serialized field for scheme reference: '+key)
        return [op, *(resolve(child, scope) for child in operand[1:])]

    def field_operand(step, scope):
        name = step.get('name')
        if name is not None: return ['field', name]
        return ['field', resolve(scheme_operand(step['name_expression']), scope)]

    def collect(steps, scope, seen=()):
        for step in steps:
            op = step['op']
            if op == 'field' and step.get('bind'):
                scope[step['bind']] = field_operand(step, scope)
            elif op == 'version' and step.get('bind'):
                scope[step['bind']] = ['version', step['name']]
            elif op == 'call' and step['layout'] not in seen:
                collect(raw['layouts'].get(step['layout'], []), scope, (*seen, step['layout']))
            elif op in ('if', 'repeat'):
                collect(step['body'], scope, seen)
                collect(step.get('else', []), scope, seen)

    def program(steps, scope):
        result = []
        for original in steps:
            step = copy.deepcopy(original)
            op = step['op']
            for key in ('condition', 'count', 'name_expression'):
                if isinstance(step.get(key), str):
                    try:
                        step[key] = resolve(scheme_operand(step[key]), scope)
                    except ValueError:
                        # Unreachable standalone serializer helpers may lack their caller's context.
                        if key != 'condition': raise
                        step = {'op': 'unsupported', 'reason': 'Predicate has no serialized field'}
                        break
            if step['op'] == 'unsupported':
                result.append(step)
                continue
            if step.get('name_expression', [None])[0] == 'literal':
                step['name'] = step.pop('name_expression')[1]
            if op == 'field':
                bind = step.pop('bind', None)
                if bind: scope[bind] = field_operand(original, scope)
            elif op == 'version':
                bind = step.pop('bind', None)
                if bind: scope[bind] = ['version', step['name']]
            elif op == 'assign':
                # A default the deserializer sets in place of an absent field, read by a later condition.
                target = scope.get(step.pop('bind'))
                if not target or target[0] != 'field' or not isinstance(target[1], str): continue
                step['name'] = target[1]
            elif op == 'select':
                step['value'] = resolve(scheme_operand(step.pop('bind')), scope)
                step['cases'] = {key: program(body, dict(scope)) for key,body in original['cases'].items()}
            elif op == 'call':
                collect(raw['layouts'].get(step['layout'], []), scope)
            if 'body' in original: step['body'] = program(original['body'], dict(scope))
            if 'else' in original: step['else'] = program(original['else'], dict(scope))
            result.append(step)
        return result

    data['layouts'] = {name: tidy(program(steps, dict(initial))) for name, steps in raw['layouts'].items()}
    if 'writer_layouts' in raw:
        data['writer_layouts'] = {name: tidy(program(steps, dict(initial))) for name, steps in raw['writer_layouts'].items()}
    for name, definition in data['types'].items():
        definition['body'] = program(raw['types'][name]['body'], dict(initial))
    data['version'] = 1
    return data


def unbound(value):
    """A layout without its source bindings, to compare what two layouts read."""
    if isinstance(value, list): return [unbound(v) for v in value]
    if isinstance(value, dict): return {k: unbound(v) for k, v in value.items() if k != 'bind'}
    return value


def tidy(steps):
    """Drop branches left without decoding steps."""
    result = []
    for step in steps:
        if step['op'] in ('if', 'repeat'):
            step = {**step, 'body': tidy(step['body'])}
            if 'else' in step:
                step['else'] = tidy(step['else'])
                if not step['else']: del step['else']
            if step['op'] == 'if' and not step['body'] and not step.get('else'): continue
        elif step['op'] == 'select':
            step = {**step, 'cases': {k: tidy(v) for k, v in step['cases'].items()}}
        result.append(step)
    return result


def generate(root):
    source = root/'casadi/core'
    if (root/'.git').exists():
        # Ignore untracked build products and merge-tool backup sources.
        tracked = subprocess.check_output(
            ['git', '-C', str(root), 'ls-files', '-z', '--', 'casadi'])
        files = sorted(root/Path(name.decode()) for name in tracked.split(b'\0')
                       if name and Path(name.decode()).suffix in ('.cpp', '.hpp'))
    else:
        files = sorted([*(root/'casadi').rglob('*.cpp'), *(root/'casadi').rglob('*.hpp')])
    texts = {p: mask_comments(p.read_text()) for p in files}
    headers = '\n'.join(text for p, text in texts.items() if p.suffix == '.hpp')
    global_decls = declarations(headers)
    for m in re.finditer(r'\b(?:struct|class)\s+(?:CASADI(?:_\w+)?_EXPORT\s+)?(\w+)(?:\s*:[^{;]+)?\s*\{', headers):
        STRUCTS[m[1]] = declarations(headers[m.end():balanced(headers, m.end()-1)])
    layouts, versions, readable_versions = {}, {}, {}
    contracts = {}
    for path, text in texts.items():
        for call in re.finditer(r'\b\w+\.version\(\s*"([^"]+)"\s*,\s*(\d+)(?:\s*,\s*(\d+))?\s*\)', text):
            first, last = int(call[2]), int(call[3] or call[2])
            readable_versions.setdefault(call[1], set()).update(range(first, last+1))
        local = declarations(text)
        stem = path.stem.removesuffix('_impl')
        own_header = (texts.get(path.with_name(stem+'.hpp'), '')
                      + texts.get(path.with_name(stem+'_impl.hpp'), ''))
        decls = {**global_decls, **local, **declarations(own_header)}
        # Collect named-field type hints, including local deserializer declarations.
        named = re.compile(r'(?<![\w])(?:\w+\.)?(?:pack|unpack)\s*\(')
        for call in named.finditer(text):
            finish = balanced(text, call.end()-1, '(', ')')
            args = split_args(text[call.end():finish])
            if len(args) != 2:
                continue
            literals = re.findall(r'"(?:\\.|[^"\\])*"', args[0])
            if not literals:
                continue
            literal = ''.join(json.loads(v) for v in literals)
            if '::' not in literal:
                continue
            pattern = literal if args[0].strip() == literals[0] else '*' + literal
            entry = contracts.setdefault(pattern, {'cpp_types': set()})
            cpp_type = 'int' if pattern.endswith('::serialization::version') else infer_type(args[1], {
                **decls, **declarations(text[text.rfind('{', 0, call.start()):call.start()])})
            if pattern == 'Shared::reference':
                # Protocol-3 shared references use casadi_int, independently of object type.
                cpp_type = 'casadi_int'
            if pattern == 'Matrix::sparsity':
                cpp_type = 'Sparsity'
            if pattern == 'Matrix::nonzeros':
                cpp_type = 'std::vector<Scalar>'
            if cpp_type:
                entry['cpp_types'].add(cpp_type)
        for match in METHOD.finditer(text):
            end = balanced(text, match.end()-1)
            body = text[match.end():end].strip()
            owner, method, stream = match['owner'], match['method'], match['stream']
            owner = re.sub(r'\s+', ' ', owner)
            # Members of the owning class shadow same-named members of its neighbours.
            cls = re.search(r'\bclass\s+(?:CASADI(?:_\w+)?_EXPORT\s+)?'+re.escape(re.sub(r'<.*>', '', owner).split('::')[-1])
                            + r'\b(?:\s*:[^{;]+)?\s*\{', own_header)
            scoped = {**decls, **declarations(own_header[cls.end():balanced(own_header, cls.end()-1)])} if cls else decls
            fields, calls = [], []
            pack = re.compile(r'\b'+re.escape(stream)+r'\.(pack|version)\s*\(')
            for call in pack.finditer(body):
                finish = balanced(body, call.end()-1, '(', ')')
                args = split_args(body[call.end():finish])
                if len(args) != 2:
                    continue
                label = json.loads(args[0]) if args[0].startswith('"') else None
                entry = {'operation': call[1], 'name': label, 'expression': args[1],
                         'offset': call.start(), 'cpp_type': infer_type(args[1], scoped)}
                if call[1] == 'version' and label and args[1].isdigit():
                    versions.setdefault(label, set()).add(int(args[1]))
                fields.append(entry)
            for call in re.finditer(r'([\w:<>, ]+)::(serialize_body|serialize_type|delayed_serialize_members)\s*\('+re.escape(stream)+r'\)', body):
                calls.append({'owner': call[1].strip(), 'method': call[2], 'offset': call.start()})
            key = owner+'::'+method
            if key in layouts:
                raise ValueError('Duplicate serializer '+key)
            layouts[key] = {'fields': fields, 'calls': calls, 'cpp_body': body}
    stream = texts[source/'serializing_stream.cpp']
    cmake = (root/'CMakeLists.txt').read_text()
    release = '.'.join(re.search(r'set\(CASADI_'+part+r'_VERSION (\d+)\)', cmake)[1]
                       for part in ('MAJOR', 'MINOR', 'PATCH'))
    reader = generate_layouts(root, texts, layouts, contracts,
                              operations(texts[source/'calculus.hpp']),
                              {name: (min(v), max(v)) for name, v in readable_versions.items()})
    reader = lower_scheme_operands(reader)
    check_readers(reader, reader.pop('writer_layouts'), versions)
    return {
        'reader': reader,
        'format': 'casadi_serialization_scheme', 'version': 1, 'casadi_version': release,
        'wire': {'protocol': int(re.search(r'serialization_protocol_version\s*=\s*(\d+)', stream)[1]),
                 'magic': int(re.search(r'serialization_check\s*=\s*(\d+)', stream)[1]),
                 'encoding': 'two a-p characters per byte, low nibble first',
                 'byte_order': 'writer-native',
                 'blob': {'length': 'uint64', 'unit': 'bytes', 'payload_encoding': 'same as stream',
                          'can_defer_payload': True}},
        'class_versions': {name: sorted(v) for name, v in sorted(versions.items())},
        'readable_class_versions': {name: sorted(v) for name, v in sorted(readable_versions.items())},
        'field_versions': field_versions(reader, {name: sorted(v) for name, v in readable_versions.items()}),
        'operations': operations(texts[source/'calculus.hpp']),
    }


def generate_layouts(root, texts, indexed, contracts, ops, readable):
    """Lower serializers into decoding rules; unsupported constructs fail explicitly."""
    layouts, parents, records = {}, {}, {}
    def make_record(body, decls, stream):
        return extract_record(body, decls, stream)
    def normalize(owner):
        owner = owner.split(' void ')[-1]
        return re.sub(r'\s+', '', owner)
    for key, value in indexed.items():
        records[normalize(key)] = value
    # Capture inline serializers with the declarations of their enclosing class.
    for path, text in texts.items():
        classes = []
        for m in re.finditer(r'\bclass\s+(?:CASADI(?:_\w+)?_EXPORT\s+)?(\w+)(?:\s*:[^{;]+)?\s*\{', text):
            end = balanced(text, m.end()-1)
            classes.append((m.start(), end, m[1], text[m.end():end]))
            parent = re.search(r':\s*public\s+([\w:<>, ]+)', m[0])
            if parent:
                parents[m[1]] = normalize(parent[1].rstrip('{ '))
        for m in re.finditer(r'\bvoid\s+(serialize_node|serialize_type|serialize_body|serialize)\s*\(SerializingStream\s*&\s*(\w+)\)\s*const\s*(?:override\s*)?\{', text):
            enclosing = [c for c in classes if c[0] < m.start() < c[1]]
            if not enclosing: continue
            cls = max(enclosing, key=lambda c:c[0])
            body = text[m.end():balanced(text,m.end()-1)].strip()
            records[cls[2]+'::'+m[1]] = make_record(body, declarations(cls[3]), m[2])
        # Helpers which serialize data on behalf of a class, optionally under a name prefix.
        for m in re.finditer(r'\bvoid\s+([\w:]+)\s*\(\s*SerializingStream\s*&\s*(\w+)\s*,([^{;]*)\)\s*(?:const\s*)?\{', text):
            body=text[m.end():balanced(text,m.end()-1)].strip()
            params, prefix = {}, None
            for arg in split_args(m[3]):
                arg = re.fullmatch(r'\s*(?:const\s+)?(.+?)\s*&?\s*(\w+)\s*', arg)
                if not arg: continue
                params.setdefault(arg[2], set()).add(arg[1])
                if arg[1] == 'std::string': prefix = arg[2]
            record=make_record(body,{**params,**declarations(body)},m[2])
            record['helper']=prefix
            records[normalize(m[1])]=record

    def resolve_type(field):
        t=field['cpp_type']
        if t is None:
            types=contracts.get(field.get('name'),{}).get('cpp_types',[])
            if len(types)==1: t=next(iter(types))
        return t

    def compile_record(key, record):
        body=record['cpp_body']
        owner=key.rsplit('::',1)[0]
        events=[]
        for f in record['fields']:
            if f['operation']=='version':
                step={'op':'version','name':f['name'],'value':int(f['expression'])}
                if f['name'] is None: step['name_expression']=f['name_expression']
            else:
                step={'op':'field','name':f['name'],'type':resolve_type(f),'bind':f['expression']}
                if f.get('name_expression'): step['name_expression']=f['name_expression']
            events.append((f['offset'],step))
        for c in record.get('calls',[]):
            events.append((c['offset'],{'op':'call','layout':normalize(c['owner']+'::'+c['method'])}))
        for m in re.finditer(r'\b([\w:]+)\s*\(\s*s\s*,',body):
            callee=normalize(m[1])
            if callee not in records and owner+'::'+callee in records: callee=owner+'::'+callee
            if 'helper' not in records.get(callee,{}): continue
            step={'op':'call','layout':callee}
            if records[callee]['helper']:
                opening=body.index('(',m.start())
                args=split_args(body[opening+1:balanced(body,opening,'(',')')])
                if not args[1].startswith('"'): continue
                step['params']={records[callee]['helper']:json.loads(args[1])}
            events.append((m.start(),step))
        return lower_body(body, unlowered(body, events, 's'))

    def unlowered(body, events, stream):
        for m in re.finditer(r'\b(?!(?:if|for|while|switch)\b)[\w:]+(?:\.\w+|->\w+)?\s*\([^;]*?\b'+stream+r'\b[^;]*?\)\s*;', body):
            if not any(body.rfind(';', 0, m.start())+1 <= position < m.end() for position, _ in events):
                events.append((m.start(), {'op':'unsupported','reason':'unlowered serialization call'}))
        return events

    def lower_body(body, events):
        ranges=[]
        alternatives={}
        for m in re.finditer(r'\b(if|for|while|switch)\s*\(',body):
            end=balanced(body,m.end()-1,'(',')'); condition=body[m.end():end].strip()
            start=end+1
            while start<len(body) and body[start].isspace():start+=1
            braced = start<len(body) and body[start]=='{'
            finish = balanced(body,start) if braced else body.find(';',start)
            if finish<0:continue
            ranges.append((m.start(),start,finish,m[1],condition,braced))
            after = re.match(r'\s*else\s*\{', body[finish+1:])
            if after:
                other = finish+1+after.end()-1
                alternatives[m.start()] = (other, balanced(body, other))
        def region(start,end):
            children=[r for r in ranges if start<=r[0]<end and not any(
                start<=p[0]<r[0]<alternatives.get(p[0], (0, p[2]))[1]<=end for p in ranges)]
            ordered=[(p,s) for p,s in events if start<=p<end and not any(r[0]<=p<=alternatives.get(r[0], (0, r[2]))[1] for r in children)]
            for p,a,b,kind,condition,braced in children:
                inner=region(a+1 if braced else a,b)
                other = alternatives.get(p)
                alternate = region(other[0]+1,other[1]) if other else []
                if not inner and not alternate: continue  # Computation without serialized output.
                if kind=='if' and alternate and unbound(inner)==unbound(alternate):
                    # Both branches read the same: the condition does not shape the stream.
                    ordered.extend((p, x) for x in inner); continue
                if kind=='if':
                    step={'op':'if','condition':condition,'body':inner}
                    if alternate: step['else']=alternate
                elif kind=='for' and ':' in condition:
                    count=condition.split(':',1)[1].strip()+'.size()'
                    # Reuse serialized counts asserted equal to the container extent.
                    alias=re.search(r'casadi_assert_dev\(\s*'+re.escape(count)
                                    +r'\s*==\s*(\w+)\s*\)', body[:p])
                    if alias: count=alias[1]
                    # FunctionInternal reconstructs these counts from its sparsities.
                    count={'n_in_':'sparsity_in_.size()',
                           'n_out_':'sparsity_out_.size()'}.get(count,count)
                    step={'op':'repeat','count':count,'body':inner}
                elif kind=='for' and re.fullmatch(r'[^;]*;\s*\w+\s*<\s*[^;]+;\s*(?:\+\+\w+|\w+\+\+)\s*', condition):
                    # Counted loops in deserializers: the bound is a field read before the loop.
                    step={'op':'repeat','count':re.search(r';\s*\w+\s*<\s*([^;]+?)\s*;', condition)[1],'body':inner}
                else:step={'op':'unsupported','reason':'unlowered '+kind+' loop or branch'}
                ordered.append((p,step))
            return [s for _,s in sorted(ordered,key=lambda x:x[0])]
        result=region(0,len(body))
        # An unhandled else branch must not be read unconditionally.
        if len(re.findall(r'\belse\b',body)) != len(alternatives):
            return [{'op':'unsupported','reason':'unlowered else branch'}]
        return result

    def deserializer(cls, body, init, stream, decls):
        """Events of a deserializing constructor, named after the writer layouts they mirror."""
        events=[]
        for m in re.finditer(r'([\w:]+(?:<[^(){};]*>)?)\s*\(\s*'+stream+r'\s*\)', init or ''):
            events.append((0, {'op':'call','layout':normalize(m[1])+'::serialize_body'}))
        for m in re.finditer(r'\b'+stream+r'\.(unpack|version)\s*\(',body):
            args=split_args(body[m.end():balanced(body,m.end()-1,'(',')')])
            if not re.fullmatch(r'"[^"\\]*"', args[0]):
                events.append((m.start(), {'op':'unsupported','reason':'computed field name'})); continue
            name=json.loads(args[0])
            if m[1]=='version':
                step={'op':'version','name':name,'value':int(args[-1]),'minimum':int(args[1])}
                local=re.search(r'\b(\w+)\s*=\s*$', body[:m.start()])
                if local: step['bind']=local[1]
            else:
                scope={**decls, **declarations(body[:m.start()])}
                step={'op':'field','name':name,'type':resolve_type({'name':name,'cpp_type':infer_type(args[1], scope)}),
                      'bind':args[1]}
            events.append((m.start(), step))
        for m in re.finditer(r'([\w:<>, ]+?)::delayed_deserialize_members\s*\(\s*'+stream+r'\s*\)', body):
            events.append((m.start(), {'op':'call','layout':normalize(m[1])+'::delayed_serialize_members'}))
        # Helpers reading on behalf of the class mirror a writer helper: deserialize/unpack -> serialize/pack.
        for m in re.finditer(r'\b((?:[\w]+::)*)(deserialize|unpack)(\w*)\s*\(\s*'+stream+r'\s*,', body):
            helper=normalize(m[1]+{'deserialize':'serialize','unpack':'pack'}[m[2]]+m[3])
            if helper not in records and normalize(cls+'::'+helper) in records: helper=normalize(cls+'::'+helper)
            if 'helper' not in records.get(helper,{}): continue
            step={'op':'call','layout':helper}
            if records[helper]['helper']:
                opening=body.index('(',m.start())
                args=split_args(body[opening+1:balanced(body,opening,'(',')')])
                if not args[1].startswith('"'): continue
                step['params']={records[helper]['helper']:json.loads(args[1])}
            events.append((m.start(), step))
        for m in re.finditer(r'\bcasadi_error\s*\(\s*("(?:\\.|[^"\\])*")', body):
            events.append((m.start(), {'op':'unsupported','reason':json.loads(m[1])}))
        # Defaults set in place of fields that older versions lack; kept if a condition reads them.
        for m in re.finditer(r'(?<![\w.>])([A-Za-z_]\w*)\s*=\s*("(?:\\.|[^"\\])*"|-?\d+|true|false)\s*;', body):
            events.append((m.start(), {'op':'assign','bind':m[1],'value':json.loads(m[2])}))
        return unlowered(body, events, stream)

    multi={name for name,(first,last) in readable.items() if first<last}
    derived={}
    for path, text in texts.items():
        for m in re.finditer(r'\b([A-Za-z_]\w*(?:<[^;{}()]*?>)?)\s*::\s*([A-Za-z_]\w*)\s*\(\s*DeserializingStream\s*&\s*(\w+)\s*\)\s*(:[^{;]*)?\{', text):
            cls=m[2]
            if re.sub(r'<.*>','',m[1])!=cls or cls not in multi: continue
            body=text[m.end():balanced(text,m.end()-1)]
            if not re.search(r'\b'+m[3]+r'\.version\(\s*"'+cls+'"', body): continue
            stem=path.stem.removesuffix('_impl')
            header=texts.get(path.with_name(stem+'.hpp'),'')+texts.get(path.with_name(stem+'_impl.hpp'),'')
            owner=re.search(r'\bclass\s+(?:CASADI(?:_\w+)?_EXPORT\s+)?'+cls+r'\b(?:\s*:[^{;]+)?\s*\{', header)
            decls={**declarations(header), **(declarations(header[owner.end():balanced(header,owner.end()-1)]) if owner else {})}
            derived[cls+'::serialize_body']=lower_body(body, deserializer(cls, body, m[4], m[3], decls))

    writer={}
    for key, record in records.items():
        layouts[key]=compile_record(key, record)
        if key in derived:
            writer[key]=layouts[key]
            layouts[key]=derived[key]
    functions={}
    text=texts[root/'casadi/core/function_internal.cpp']
    for name,cls in re.findall(r'\{"([^"]+)",\s*(\w+)::deserialize\}',text):functions[name]=cls
    mx={str(ops[op]) if op in ops else op:normalize(cls) for op,cls in re.findall(r'\{(OP_\w+|-1),\s*([\w<>:, ]+)::deserialize\}',texts[root/'casadi/core/mx_node.cpp'])}
    # The native dispatcher uses these byte-layout families for arithmetic nodes.
    math=texts[root/'casadi/core/calculus.hpp']
    def family(kind):
        macro=re.search(r'#define CASADI_MATH_'+kind+r'_BUILTIN[^\n]*(?:\n[^\n]*\\)*\n[^\n]*',math)[0]
        method=re.search(r'bool casadi_math<T>::is_'+kind.lower()+r'\([^)]*\)\s*\{',math)
        body=math[method.end():balanced(math,method.end()-1)]
        return set(re.findall(r'case (OP_\w+)',macro+body))
    unary,binary=family('UNARY'),family('BINARY')
    sx={str(ops['OP_PARAMETER']):'SymbolicSX',str(ops['OP_CALL']):'CallSX','-1':'OutputSX'}
    for name in unary:
        mx[str(ops[name])]='UnaryMX';sx[str(ops[name])]='UnarySX'
    for name in binary:
        mx[str(ops[name])]='BinaryMX';sx[str(ops[name])]='BinarySX'
    generic_text=texts[root/'casadi/core/generic_type.hpp']
    names=re.findall(r'\bOT_\w+',generic_text[generic_text.index('enum TypeID'):generic_text.index('};',generic_text.index('enum TypeID'))])
    generic={}
    for m in re.finditer(r'typedef\s+GenericTypeInternal<(OT_\w+),\s*([\s\S]*?)>\s+\w+Type;',texts[root/'casadi/core/generic_type.cpp']):
        generic[str(names.index(m[1]))]=normalize(m[2])
    plugins = {}
    for text in texts.values():
        creator = re.search(r'plugin->creator\s*=\s*(\w+)::creator', text)
        name = re.search(r'plugin->name\s*=\s*"([^"]+)"', text)
        if creator and name:
            cls = creator[1]
            base = cls
            seen = set()
            while base not in seen:
                seen.add(base)
                for function, owner in {**functions, 'Linsol':'LinsolInternal'}.items():
                    if owner == base:
                        plugins.setdefault(function, {})[name[1]] = cls
                base = parents.get(base, '').split('<')[0]
                if not base: break
    programs = type_programs(layouts, parents, functions, mx, sx, generic, ops, plugins)
    stream = texts[root/'casadi/core/serializing_stream.cpp']
    for match in re.finditer(r'void SerializingStream::pack\(const (\w+)& \w+\)\s*\{', stream):
        body=stream[match.end():balanced(stream,match.end()-1)]
        decoration=re.search(r"decorate\('([^']+)'\)",body)
        if match[1] in programs:
            programs[match[1]]['shared']='shared_pack' in body
            if decoration:programs[match[1]]['decoration']=decoration[1]
    materialize_layouts(layouts, parents, programs)
    result = {'writer_layouts': writer, 'file_types': {'0':'Sparsity','2':'DM','4':'Linsol','5':'Function','6':'GenericType','7':'casadi_int','8':'double','9':'std::string','10':'std::vector<Sparsity>','12':'std::vector<DM>','15':'std::vector<Function>','16':'std::vector<GenericType>','17':'std::vector<casadi_int>','18':'std::vector<double>','19':'std::vector<std::string>'}, 'types': programs, 'version':1,'layouts':dict(sorted(layouts.items()))}
    # Expression files store ordered dependencies before the expression root(s).
    for tag, kind, prefix in [(1, 'MX', 'Function'), (3, 'SX', 'Function'),
                              (11, 'std::vector<MX>', 'Function'),
                              (13, 'std::vector<SX>', 'Function'),
                              (20, 'MX', 'std::vector<MX>'),
                              (21, 'SX', 'std::vector<SX>'),
                              (22, 'std::vector<MX>', 'std::vector<MX>'),
                              (23, 'std::vector<SX>', 'std::vector<SX>')]:
        result['file_types'][str(tag)] = kind
        result.setdefault('file_prefixes', {})[str(tag)] = prefix
    return compact_bindings(result)



def materialize_layouts(layouts, parents, programs):
    """Make inherited and template serializer calls explicit at generation time."""
    def calls(value):
        if isinstance(value, dict):
            if value.get('op') == 'call': yield value['layout']
            for item in value.values(): yield from calls(item)
        elif isinstance(value, list):
            for item in value: yield from calls(item)

    implementations = tuple(layouts)

    def resolve(name, seen):
        if name in layouts: return True
        if name in seen: raise ValueError('Cyclic serializer inheritance: '+name)
        seen = seen | {name}
        cls, method = name.rsplit('::', 1)
        base = re.sub(r'<.*>', '', cls)
        candidates = [key for key in implementations
                      if key.startswith(base+'<') and key.endswith('::'+method)]
        if len(candidates) == 1:
            target = candidates[0]
        elif base in parents:
            target = parents[base]+'::'+method
            if not resolve(target, seen): return False
        else:
            return False
        layouts[name] = [{'op': 'call', 'layout': target}]
        return True

    for name in sorted(set(calls(layouts)) | set(calls(programs))):
        resolve(name, set())


def compact_bindings(reader):
    """Retain field variables only when a decoding expression can reference them.

    Calls share scope, so collect uses across all layouts conservatively. This
    intentionally keeps ambiguous uses rather than changing decoding semantics.
    """
    expressions = []
    fields = []
    def visit(value):
        if isinstance(value, list):
            for item in value: visit(item)
        elif isinstance(value, dict):
            if value.get('op') == 'field': fields.append(value)
            for key, item in value.items():
                if key in ('condition', 'count', 'name_expression') or (
                        key == 'bind' and value.get('op') == 'select'):
                    if isinstance(item, str): expressions.append(re.sub(r'"(?:[^"\\]|\\.)*"', '', item))
                if key == 'params':
                    expressions.extend(v for v in item.values() if isinstance(v, str))
                visit(item)
    visit(reader)
    for field in fields:
        bind = field.get('bind')
        if bind and not any(bind in expression for expression in expressions):
            del field['bind']
    def defaults(value):
        if isinstance(value, list):
            value[:] = [v for v in value if not (isinstance(v, dict) and v.get('op') == 'assign'
                        and not any(re.search(r'\b'+re.escape(v['bind'])+r'\b', e) for e in expressions))]
            for item in value: defaults(item)
        elif isinstance(value, dict):
            for item in value.values(): defaults(item)
    defaults(reader)
    return reader


def at_version(steps, versions):
    """Specialize a layout to fixed class versions, dropping defaults and untaken branches."""
    def value(operand):
        op = operand[0]
        if op == 'literal': return operand[1]
        if op == 'version': return versions.get(operand[1])
        if op in ('field', 'parameter', 'length', 'concat', 'next'): return None
        args = [value(a) for a in operand[1:]]
        if op == 'not': return None if args[0] is None else not args[0]
        if op == 'and' and False in args or op == 'or' and True in args: return op == 'or'
        if None in args: return None
        a, b = args
        return {'and': lambda: a and b, 'or': lambda: a or b, 'equal': lambda: a == b,
                'not_equal': lambda: a != b, 'less': lambda: a < b, 'less_equal': lambda: a <= b,
                'greater': lambda: a > b, 'greater_equal': lambda: a >= b}[op]()
    result = []
    for step in steps:
        op = step['op']
        if op == 'assign': continue
        if op == 'version':
            result.append({k: v for k, v in step.items() if k != 'minimum'})
        elif op == 'if':
            taken = value(step['condition'])
            if taken is None:
                result.append({**step, 'body': at_version(step['body'], versions),
                               **({'else': at_version(step['else'], versions)} if 'else' in step else {})})
            else:
                result += at_version(step['body'] if taken else step.get('else', []), versions)
        elif op == 'repeat':
            result.append({**step, 'body': at_version(step['body'], versions)})
        elif op == 'select':
            result.append({**step, 'cases': {k: at_version(v, versions) for k, v in step['cases'].items()}})
        else:
            result.append(step)
    return tidy(result)


def check_readers(reader, writers, class_versions):
    """A layout derived from a deserializer must read exactly what the writer writes today."""
    current = {name: max(v) for name, v in class_versions.items()}
    mismatches = [name for name, steps in writers.items()
                  if at_version(reader['layouts'][name], current) != steps]
    if mismatches:
        raise ValueError('Deserializer and serializer disagree at the current version: '+', '.join(sorted(mismatches)))
    # A layout may only branch on the versions it reads itself.
    def visit(value, own, layout):
        if isinstance(value, dict):
            if value.get('op') == 'version': own.add(value['name'])
            for item in value.values(): visit(item, own, layout)
        elif isinstance(value, list):
            if len(value) == 2 and value[0] == 'version' and value[1] not in own:
                raise ValueError(layout+' branches on the version of '+value[1])
            for item in value: visit(item, own, layout)
    for name, steps in reader['layouts'].items():
        visit(steps, set(), name)


def field_versions(reader, readable):
    """Fields of versioned layouts that only some readable versions carry, for tools; readers ignore it."""
    def fields(value, out):
        if isinstance(value, dict):
            if value.get('op') == 'field' and value.get('name'): out.add(value['name'])
            for item in value.values(): fields(item, out)
        elif isinstance(value, list):
            for item in value: fields(item, out)
        return out
    result = {}
    for cls, versions in readable.items():
        steps = reader['layouts'].get(cls+'::serialize_body')
        if len(versions) < 2 or steps is None: continue
        present = {}
        for v in versions:
            specialized = at_version(steps, {cls: v})
            # A version the deserializer refuses carries no fields.
            if any(step['op'] == 'unsupported' for step in specialized): continue
            present[v] = fields(specialized, set())
        readable_versions = sorted(present)
        partial = {name: [v for v in readable_versions if name in present[v]]
                   for name in set().union(*present.values())}
        partial = {name: vs for name, vs in partial.items() if vs != readable_versions}
        if partial: result[cls] = dict(sorted(partial.items()))
    return result


def extract_record(body, decls, stream):
    fields=[];calls=[]
    for m in re.finditer(r'\b'+stream+r'\.(pack|version)\s*\(',body):
        end=balanced(body,m.end()-1,'(',')');args=split_args(body[m.end():end])
        if len(args)!=2:continue
        name=json.loads(args[0]) if re.fullmatch(r'"[^"\\]*"',args[0]) else None
        fields.append({'operation':m[1],'name':name,'name_expression':args[0],
                       'expression':args[1],'cpp_type':infer_type(args[1],decls),'offset':m.start()})
    for m in re.finditer(r'([\w:<>, ]+)::(serialize_body|serialize_type|delayed_serialize_members)\s*\('+stream+r'\)',body):
        calls.append({'owner':m[1].strip(),'method':m[2],'offset':m.start()})
    return {'cpp_body':body,'fields':fields,'calls':calls}


def type_programs(layouts, parents, functions, mx, sx, generic, ops, plugins):
    def field(name, t, bind=None):
        return {'op':'field','name':name,'type':t,'bind':bind or name}
    def call(name, **params):
        return {'op':'call','layout':name, **({'params':params} if params else {})}
    def select(bind, cases):
        return {'op':'select','bind':bind,'cases':cases}
    def body(cls):
        # Concrete template substitutions remain layout metadata, not runtime logic.
        if cls in ('MXFunction','SXFunction'):
            mat='MX' if cls=='MXFunction' else 'SX'
            return [call(cls+'::serialize_body',MatType=mat)]
        return [call(cls+'::serialize_body')]
    function_cases={name:body(cls) for name,cls in functions.items()}
    for name, cls in functions.items():
        prefix = layouts.get(cls+'::serialize_type', [])
        prefix = [s for s in prefix if s != {'op':'call','layout':'FunctionInternal::serialize_type'}]
        if prefix and all(s['op']=='field' for s in prefix):
            selector = next((s for s in prefix if s.get('bind')=='class_name()'), None)
            if selector:
                candidates={cls}
                for child in parents:
                    parent=child;seen=set()
                    while parent and parent not in seen:
                        seen.add(parent)
                        if parent==cls:candidates.add(child)
                        parent=parents.get(parent, '').split('<')[0]
                function_cases[name]=prefix+[select('class_name()', {c:body(c) for c in sorted(candidates)})]
            else:
                function_cases[name]=prefix+body(cls)

    for family, registrations in plugins.items():
        if family not in functions: continue
        function_cases[family]=[field('PluginInterface::plugin_name','std::string','plugin'),
                                select('plugin',{name:body(cls) for name,cls in registrations.items()})]
    if 'Interpolant' in function_cases and 'linear' in plugins.get('Interpolant', {}):
        function_cases['Interpolant'][1]['cases']['linear'] = [
            {'op':'version','name':'LinearInterpolant','value':1},
            field('LinearInterpolant::type','char','linear_type'),
            select('linear_type', {str(ord('f')):body('LinearInterpolant'),
                                   str(ord('j')):body('LinearInterpolantJac')})]
    if 'External' in function_cases:
        function_cases['External'] = [
            {'op':'version','name':'GenericExternal','value':1},
            field('GenericExternal::type','char','subtype'),
            select('subtype', {str(ord('g')):body('External')})]
    types={
        'Function':{'shared':True,'body':[field('Function::null','bool','null'),
            {'op':'if','condition':'!null','body':[
                field('FunctionInternal::base_function','std::string','base'),
                select('base',function_cases)]}]},
        'Sparsity':{'shared':True,'body':[field('SparsityInternal::compressed','std::vector<casadi_int>')]},
        'Slice':{'body':[call('Slice::serialize')]},
        'DM':{'body':[field('Matrix::sparsity','Sparsity'),field('Matrix::nonzeros','std::vector<double>')]},
        'IM':{'body':[field('Matrix::sparsity','Sparsity'),field('Matrix::nonzeros','std::vector<casadi_int>')]},
        'SX':{'body':[field('Matrix::sparsity','Sparsity'),field('Matrix::nonzeros','std::vector<SXElem>')]},
        'GenericType':{'body':[field('GenericType::type','int','type'),select('type',{
            k:[field('GenericType::d',v)] for k,v in generic.items()})]},
    }
    # Expand dispatch families to their serialized layouts, without evaluating nodes.
    mx_cases={k:body(v) for k,v in mx.items()}
    for k,v in mx.items():
        if v=='BinaryMX':
            mx_cases[k]=[field('BinaryMX::scalar_flags','char')]+body('BinaryMX<ScX,ScY>')
    nonzeros = {'GetNonzeros': ['Vector', 'Slice', 'Slice2'],
                'SetNonzeros': ['Vector<Add>', 'Slice<Add>', 'Slice2<Add>'],
                'SetNonzerosParam': ['ParamVector<Add>', 'ParamSlice<Add>', 'SliceParam<Add>', 'ParamParam<Add>']}
    for name,tag in [('OP_GETNONZEROS','GetNonzeros'),('OP_SETNONZEROS','SetNonzeros'),('OP_ADDNONZEROS','SetNonzeros'),
                     ('OP_SETNONZEROS_PARAM','SetNonzerosParam'),('OP_ADDNONZEROS_PARAM','SetNonzerosParam')]:
        stem = tag.removesuffix('Param')
        mx_cases[str(ops[name])]=[field(tag+'::type','char','subtype'),select('subtype',{
            str(ord('a')+i):body(stem+variant) for i,variant in enumerate(nonzeros[tag])
        })]
    mx_cases[str(ops['OP_SOLVE'])]=[field('Solve::Tr','bool')]+body('LinsolCall<Tr>')
    # These subtype tags are literal fields in serialize_type. Derive their
    # cases from those writers, including subclasses which inherit the body.
    for op, name in [('OP_TRANSPOSE','Transpose::dense'),
                     ('OP_MTIMES','Multiplication::kind'), ('OP_BSPLINE','BSpline::type'),
                     ('OP_PROJECT','Project::type'), ('OP_KRON','Kron::kind'),
                     ('OP_KRON_CONTRACT','KronContract::kind'),
                     ('OP_GETNONZEROS_PARAM','GetNonzerosParam::type')]:
        cases, field_type = {}, None
        for layout, steps in layouts.items():
            if not layout.endswith('::serialize_type'): continue
            for step in steps:
                if step.get('name') != name: continue
                literal = step['bind']
                if literal in ('true','false'): key = literal
                elif re.fullmatch(r"'[^']'", literal): key = str(ord(literal[1]))
                else:
                    match = re.fullmatch(r'std::string\((".*")\)', literal)
                    if not match: raise ValueError('Nonliteral subtype tag '+name)
                    key = json.loads(match[1])
                field_type = step['type']
                cases[key] = body(layout.removesuffix('::serialize_type'))
        if not cases: raise ValueError('No serializer variants for '+op)
        mx_cases[str(ops[op])] = [field(name,field_type,'subtype'), select('subtype',cases)]
    # Before 3.8, Multiplication wrote a bool and no BLAS field; Multiplication::deserialize
    # tells the formats apart by the upcoming field (debug) or a bool-valued next byte.
    mx_cases[str(ops['OP_MTIMES'])] = [{'op':'if',
        'condition':['next',['literal','Multiplication::dense'],['literal','bool']],
        'body':[field('Multiplication::dense','bool'), call('MXNode::serialize_body')],
        'else':mx_cases[str(ops['OP_MTIMES'])]}]
    constants={str(ord(c)):[call('MXNode::serialize_body')] for c in ['0','1','m']}
    constants[str(ord('z'))]=[]
    for char,typ in [('D','double'),('I','casadi_int')]:
        constants[str(ord(char))]=[field('Constant::value',typ),call('MXNode::serialize_body')]
    for char,cls in [('a','ConstantDM'),('f','ConstantFile'),('p','ConstantPool')]:
        constants[str(ord(char))]=body(cls)
    mx_cases[str(ops['OP_CONST'])]=[field('ConstantMX::type','char','subtype'),select('subtype',constants)]
    sx_cases={k:[call(v+'::serialize_node')] for k,v in sx.items()}
    constant_cases={}
    for cls in ['RealtypeSX','IntegerSX','ZeroSX','OneSX','MinusOneSX','InfSX','MinusInfSX','NanSX']:
        key=cls+'::serialize_node'
        if key not in layouts: continue
        prefix=layouts[key][0]
        literal=prefix.get('bind','')
        if len(literal)==3 and literal[0]=="'":
            constant_cases[str(ord(literal[1]))]=layouts[key][1:]
    sx_cases[str(ops['OP_CONST'])]=[field('ConstantSX::type','char','subtype'),select('subtype',constant_cases)]
    types['MX']={'shared':True,'body':[field('MXNode::op','int','op'),select('op',mx_cases)]}
    types['SXElem']={'shared':True,'body':[field('SXNode::op','casadi_int','op'),select('op',sx_cases)]}
    resources={cls:body(cls) for cls in ['ZipMemResource','ZipResource','DirResource']}
    types['Linsol']={'shared':True,'body':[field('PluginInterface::plugin_name','std::string','plugin'),
        select('plugin', {name:body(cls) for name,cls in plugins.get('Linsol',{}).items()})]}
    types['Fmu']={'shared':True,'decoration':'F','body':[
        field('FmuInternal::type','std::string','fmu_type'),
        select('fmu_type', {name:body(name) for name in ['Fmu2','Fmu3']})]}
    types['Importer']={'shared':True,'decoration':'M','body':[
        field('ImporterInternal::type','std::string','importer_type'),
        select('importer_type', {'DllLibrary':body('ImporterInternal')})]}
    types['Resource']={'shared':True,'body':[{'op':'version','name':'ResourceInternal','value':1},
        field('ResourceInternal::type','std::string','resource_type'),select('resource_type',resources)]}
    for name, tag in {'Function':'F','Sparsity':'S','Slice':'S','SXElem':'E','MX':'X','GenericType':'G','Resource':'R','Linsol':'L'}.items():
        types[name]['decoration'] = tag
    return types


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check', action='store_true', help='Fail if the committed scheme is stale')
    args = parser.parse_args()
    scheme = generate(args.root.resolve())
    args.output = args.output or args.root/'misc/serialization_scheme.json'
    rendered = json.dumps(scheme, indent=2, sort_keys=True)+'\n'
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            parser.exit(1, 'Serialization scheme is stale; regenerate it.\n')
    else:
        args.output.write_text(rendered)
    print(f"{'Checked' if args.check else 'Wrote'} {args.output}: "
          f"{len(scheme['reader']['layouts'])} reader layouts")


if __name__ == '__main__':
    main()
