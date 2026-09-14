"""Validate named fields and primitive wire tags of a debug-decorated stream.

Native deserialization still checks object layout and semantics. This checker
independently verifies the generated field contract, class versions and payload
boundaries. It does not execute C++ conditions from the source-derived scheme.
"""
import fnmatch
import json
from pathlib import Path
import struct

SCHEME_PATH = Path(__file__).with_name('serialization_scheme.json')


def codec(cpp_type):
    cpp_type = cpp_type.replace(' ', '')
    if cpp_type.startswith('std::vector<'): return 'V'
    if cpp_type.startswith('std::map<') or cpp_type == 'Dict': return 'D'
    if cpp_type.startswith('std::pair<'): return 'p'
    return {'bool': 'b', 'int': 'i', 'unsignedint': 'u', 'casadi_int': 'J',
            'size_t': 'K', 'double': 'd', 'std::string': 's', 'char': 'char',
            'Sparsity': 'S', 'MX': 'X', 'SXElem': 'E', 'Function': 'F', 'Fmu': 'F',
            'Resource': 'R', 'Importer': 'M', 'Linsol': 'L', 'GenericType': 'G',
            'std::stringstream': 'B', 'std::istream': 'B',
            'Slice': 'inline', 'DM': 'inline', 'SX': 'inline', 'IM': 'inline'}.get(cpp_type)


class SchemeValidator:
    def __init__(self, scheme=None):
        self.scheme = scheme or json.loads(SCHEME_PATH.read_text())
        self.fields = self.scheme['field_contracts']
        self.patterns = [(p, v) for p, v in self.fields.items() if '*' in p]
        self.resolved = {}

    def validate(self, text):
        if len(text) % 2:
            raise ValueError('Invalid serialization encoding')
        byte_length = len(text)//2
        pos, fields, unresolved = 0, [], set()

        def skip(n):
            nonlocal pos
            if n < 0 or pos+n > byte_length: raise ValueError('Truncated serialization')
            start = pos
            pos += n
            return start

        def raw(n):
            start = skip(n)
            out = bytearray(n)
            for j in range(n):
                i = 2*(start+j)
                a, b = ord(text[i])-97, ord(text[i+1])-97
                if not (0 <= a <= 15 and 0 <= b <= 15):
                    raise ValueError('Invalid serialization encoding')
                out[j] = a + (b << 4)
            return bytes(out)

        def number(tag):
            return struct.unpack({'i': '<i', 'u': '<I', 'J': '<q', 'K': '<Q', 'd': '<d'}[tag],
                                 raw(4 if tag in ('i', 'u') else 8))[0]

        def string():
            if raw(1) != b'i': raise ValueError('Missing string-length decoration')
            data = raw(number('i'))
            try:
                return data.decode('utf8')
            except UnicodeDecodeError:
                return data

        def value(tag=None):
            tag = tag or raw(1).decode('ascii')
            if tag in ('i', 'u', 'J', 'K', 'd'): return tag, number(tag)
            if tag == 'b':
                v = raw(1)[0]
                if v > 1: raise ValueError('Invalid serialized boolean')
                return tag, bool(v)
            if tag == 's': return tag, string()
            if tag == 'B':
                if raw(1) != b'K': raise ValueError('Missing blob-length decoration')
                size = number('K')
                skip(size)
                return tag, {'bytes': size}
            if tag in 'VDpSXEFRMLG': return tag, None
            raise ValueError('Unknown wire decoration '+repr(tag)+' at '+str(pos-1))

        if number('J') != self.scheme['wire']['magic']: raise ValueError('Invalid serialization magic')
        if number('J') != self.scheme['wire']['protocol']: raise ValueError('Unsupported serialization protocol')
        if raw(1) != b'\x01': raise ValueError('Scheme validation requires debug serialization')
        # Function.serialize() has no FileSerializer type byte or root shared-object marker.
        while pos < byte_length:
            tag = raw(1).decode('ascii')
            if tag != 's':
                value(tag)
                continue
            name = string()
            if not isinstance(name, str) or '::' not in name:
                continue
            contract = self.fields.get(name) or self.resolved.get(name)
            if contract is None:
                matches = [v for pattern, v in self.patterns
                           if fnmatch.fnmatchcase(name, pattern)]
                if not matches: raise ValueError('Field absent from serialization scheme: '+name)
                contract = self.resolved[name] = matches[0]
            codes = {codec(t) for t in contract['cpp_types']} - {None}
            if codes == {'inline'}:
                tag, payload = 'inline', None
            elif codes == {'char'}:
                tag, payload = 'char', raw(1)[0]
            else:
                tag, payload = value()
            if codes and tag not in codes:
                raise ValueError(f'{name}: scheme expects {sorted(codes)}, wire has {tag}')
            if not codes: unresolved.add(name)
            if name.endswith('::serialization::version'):
                cls = name.removesuffix('::serialization::version')
                versions = self.scheme.get('readable_class_versions', self.scheme['class_versions']).get(cls)
                if versions is not None and payload not in versions:
                    raise ValueError(f'{cls}: version {payload} absent from scheme {versions}')
            fields.append(name)
        return {'fields': len(fields), 'unresolved_types': sorted(unresolved)}
