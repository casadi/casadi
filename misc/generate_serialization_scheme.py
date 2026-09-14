#!/usr/bin/env python3
"""Extract release serialization metadata from C++ without building CasADi.

The output is a source-derived scheme, not a universal executable C++ grammar.
Every serializer retains its body and every pack its expression and context.
Unresolved types are explicit; consumers must reject layouts they cannot read.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

TOKEN = re.compile(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[{}()]')
METHOD = re.compile(r'\b(?P<owner>[A-Za-z_]\w*(?:<[^;{}]*?>)?(?:::[A-Za-z_]\w*)*)\s*::\s*'
                    r'(?P<method>serialize_body|serialize_type|serialize|delayed_serialize_members)\s*'
                    r'\(\s*SerializingStream\s*&\s*(?P<stream>\w+)\s*\)\s*(?:const\s*)?\{')


def mask_comments(text):
    def replace(m):
        value = m.group()
        return re.sub(r'[^\n]', ' ', value) if value.startswith(('//', '/*')) else value
    return TOKEN.sub(replace, text)


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


def generate(root):
    source = root/'casadi/core'
    files = sorted([*(root/'casadi').rglob('*.cpp'), *(root/'casadi').rglob('*.hpp')])
    texts = {p: mask_comments(p.read_text()) for p in files}
    headers = '\n'.join(text for p, text in texts.items() if p.suffix == '.hpp')
    global_decls = declarations(headers)
    layouts, versions, readable_versions = {}, {}, {}
    contracts = {}
    for path, text in texts.items():
        for call in re.finditer(r'\b\w+\.version\(\s*"([^"]+)"\s*,\s*(\d+)(?:\s*,\s*(\d+))?\s*\)', text):
            first, last = int(call[2]), int(call[3] or call[2])
            readable_versions.setdefault(call[1], set()).update(range(first, last+1))
        local = declarations(text)
        own_header = texts.get(path.with_suffix('.hpp'), '')
        decls = {**global_decls, **local, **declarations(own_header)}
        # Include inline/template/plugin serializers in the named-field contract.
        named = re.compile(r'(?<![\w])(?:\w+\.)?pack\s*\(')
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
            entry = contracts.setdefault(pattern, {'cpp_types': set(), 'sources': []})
            cpp_type = 'int' if pattern.endswith('::serialization::version') else infer_type(args[1], decls)
            if pattern == 'Shared::reference':
                # Protocol-3 shared references use casadi_int, independently of object type.
                cpp_type = 'casadi_int'
            if pattern == 'Matrix::sparsity':
                cpp_type = 'Sparsity'
            if pattern == 'Matrix::nonzeros':
                cpp_type = 'std::vector<Scalar>'
            if cpp_type:
                entry['cpp_types'].add(cpp_type)
            entry['sources'].append({'path': str(path.relative_to(root)),
                                    'line': text.count('\n', 0, call.start())+1,
                                    'expression': args[1]})
        for match in METHOD.finditer(text):
            end = balanced(text, match.end()-1)
            body = text[match.end():end].strip()
            owner, method, stream = match['owner'], match['method'], match['stream']
            owner = re.sub(r'\s+', ' ', owner)
            fields, calls, controls = [], [], []
            pack = re.compile(r'\b'+re.escape(stream)+r'\.(pack|version)\s*\(')
            for call in pack.finditer(body):
                finish = balanced(body, call.end()-1, '(', ')')
                args = split_args(body[call.end():finish])
                if len(args) != 2:
                    continue
                label = json.loads(args[0]) if args[0].startswith('"') else None
                entry = {'operation': call[1], 'name': label, 'expression': args[1],
                         'offset': call.start(), 'cpp_type': infer_type(args[1], decls)}
                if call[1] == 'version' and label and args[1].isdigit():
                    versions.setdefault(label, set()).add(int(args[1]))
                fields.append(entry)
            for call in re.finditer(r'([\w:<>, ]+)::(serialize_body|serialize_type|delayed_serialize_members)\s*\('+re.escape(stream)+r'\)', body):
                calls.append({'owner': call[1].strip(), 'method': call[2], 'offset': call.start()})
            for control in re.finditer(r'\b(if|for|while|switch)\s*\(', body):
                finish = balanced(body, control.end()-1, '(', ')')
                controls.append({'kind': control[1], 'condition': body[control.end():finish], 'offset': control.start()})
            key = owner+'::'+method
            if key in layouts:
                raise ValueError('Duplicate serializer '+key)
            layouts[key] = {'source': str(path.relative_to(root)),
                            'line': text.count('\n', 0, match.start())+1,
                            'fields': fields, 'calls': calls, 'control_flow': controls,
                            'cpp_body': body, 'sha256': hashlib.sha256(body.encode()).hexdigest()}
    stream = texts[source/'serializing_stream.cpp']
    cmake = (root/'CMakeLists.txt').read_text()
    release = '.'.join(re.search(r'set\(CASADI_'+part+r'_VERSION (\d+)\)', cmake)[1]
                       for part in ('MAJOR', 'MINOR', 'PATCH'))
    return {
        'format': 'casadi_serialization_scheme', 'version': 1, 'casadi_version': release,
        'coverage': {'kind': 'source-derived-index',
                     'scope': 'Out-of-line core/plugin serializers plus named pack fields from all core/plugin headers',
                     'limitations': ['Inline serializer bodies are not indexed; their named fields are included in field_contracts',
                                    'C++ control flow/expressions require reader support; null cpp_type is unresolved',
                                    'This is not a complete executable deserialization grammar']},
        'wire': {'protocol': int(re.search(r'serialization_protocol_version\s*=\s*(\d+)', stream)[1]),
                 'magic': int(re.search(r'serialization_check\s*=\s*(\d+)', stream)[1]),
                 'encoding': 'two a-p characters per byte, low nibble first',
                 'byte_order': 'writer-native',
                 'blob': {'length': 'uint64', 'unit': 'bytes', 'payload_encoding': 'same as stream',
                          'source': 'casadi/core/serializing_stream.cpp:SerializingStream::pack(std::istream&)',
                          'can_defer_payload': True}},
        'class_versions': {name: sorted(v) for name, v in sorted(versions.items())},
        'readable_class_versions': {name: sorted(v) for name, v in sorted(readable_versions.items())},
        'field_contracts': {name: {**entry, 'cpp_types': sorted(entry['cpp_types'])}
                            for name, entry in sorted(contracts.items())},
        'operations': operations(texts[source/'calculus.hpp']),
        'serializers': dict(sorted(layouts.items())),
    }


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
    unresolved = sum(f['cpp_type'] is None for s in scheme['serializers'].values()
                     for f in s['fields'] if f['operation'] == 'pack')
    print(f"{'Checked' if args.check else 'Wrote'} {args.output}: {len(scheme['serializers'])} serializers, "
          f"{unresolved} unresolved field types (explicitly marked)")


if __name__ == '__main__':
    main()
