"""Check the generated serialization scheme against the CasADi sources."""
import importlib.util
from pathlib import Path
import functools
import json
import copy
import struct
import casadi as ca
import unittest
from helpers import casadiTestCase, memory_heavy

root = Path(__file__).resolve().parents[2]


def load_module(name):
  spec = importlib.util.spec_from_file_location(name, root/'misc'/(name+'.py'))
  assert spec is not None and spec.loader is not None
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


@functools.lru_cache(maxsize=None)
def validator():
  return load_module('validate_serialization_scheme').SchemeValidator()


def validate_serialization(text):
  return validator().validate(text)


class SerializationSchemeTests(casadiTestCase):
  def test_native_serialization(self):
    for T in [ca.SX, ca.MX]:
      x = T.sym('x', 3)
      f = ca.Function('f', [x], [ca.sin(x)+x*x])
      self.check_serialize(f, [ca.DM([1, 2, 3])])
      for g in [f.map(2), f.jacobian(), f.forward(2), f.reverse(2)]:
        self.assertGreater(validate_serialization(g.serialize({'debug': True}))['fields'], 0)
    sx = ca.SX.sym('s', 3)
    inner = ca.Function('inner', [sx], [ca.sin(sx)])
    outer = ca.Function('outer', [x], [inner(x)])
    self.check_serialize(outer, [ca.DM([1, 2, 3])])

  def test_constant_serialization(self):
    for value in [-1, 0, 1, 2, 0.5]:
      with self.subTest(value=value):
        f = ca.Function('constant', [], [ca.MX(value)])
        self.assertGreater(validate_serialization(f.serialize({'debug': True}))['fields'], 0)
        self.check_serialize(f, [])

  def test_layout_mutations(self):
    x = ca.MX.sym('x')
    data = ca.Function('f', [x], [ca.sin(x)]).serialize({'debug': True})
    validate_serialization(data)
    for mutation in ['name', 'type', 'version', 'order', 'missing', 'branch', 'count']:
      scheme = copy.deepcopy(validator().scheme)
      layouts = scheme['reader']['layouts']
      fields = layouts['ProtoFunction::serialize_body']
      name = next(f for f in fields if f.get('name') == 'ProtoFunction::name')
      if mutation == 'name': name['name'] = 'ProtoFunction::wrong'
      elif mutation == 'type': name['type'] = 'bool'
      elif mutation == 'version': next(f for f in fields if f['op'] == 'version')['value'] = 999
      elif mutation == 'order': fields.reverse()
      elif mutation == 'missing': del layouts['MXFunction::serialize_body']
      elif mutation == 'branch': scheme['reader']['types']['Function']['body'][1]['condition'] = ['literal', False]
      else:
        next(f for f in layouts['MXFunction::serialize_body'] if f['op'] == 'repeat')['count'] = ['literal', 0]
      with self.subTest(mutation=mutation), self.assertRaises(ValueError):
        type(validator())(scheme).validate(data)
    with self.assertRaises(ValueError): validate_serialization(data[:-2])
    with self.assertRaises(ValueError): validate_serialization(data+'aa')

  def test_payload_extent(self):
    # Large opaque payloads must be checked without materializing decoded bytes.
    x = ca.SX.sym('x')
    validate_serialization(ca.Function('f', [x], [x]).serialize({'debug': True}))
    scheme = copy.deepcopy(validator().scheme)
    scheme['reader']['types']['Function']['body'] = [
      {'op': 'field', 'name': 'Test::payload', 'type': 'std::istream'}]
    name = b'Test::payload'
    size = 2*1024*1024
    header = struct.pack('=qqB', scheme['wire']['magic'], scheme['wire']['protocol'], 1)
    prefix = header+b'si'+struct.pack('=i', len(name))+name+b'BK'+struct.pack('=Q', size)
    text = ''.join(chr(97+(b & 15))+chr(97+(b >> 4)) for b in prefix)+'aa'*size
    mutated = type(validator())(scheme)
    self.assertEqual(mutated.validate(text)['skipped_payload_bytes'], size)
    with self.assertRaises(ValueError): mutated.validate(text[:-2])

  def test_scheme_contains_decoding_rules_only(self):
    scheme = json.loads((root/'misc/serialization_scheme.json').read_text())
    self.assertNotIn('serializers', scheme)
    self.assertNotIn('dispatch', scheme['reader'])
    self.assertNotIn('parents', scheme['reader'])
    self.assertEqual(scheme['reader']['layouts']['Horzcat::serialize_body'],
                     [{'op': 'call', 'layout': 'Concat::serialize_body'}])
    self.assertEqual(scheme['reader']['layouts']['Concat::serialize_body'],
                     [{'op': 'call', 'layout': 'MXNode::serialize_body'}])
    self.assertNotRegex(json.dumps(scheme), r'\bT1\b')
    self.assertNotIn('field_contracts', scheme)
    self.assertEqual(scheme['reader']['version'], 1)
    fields = scheme['reader']['layouts']['Ipqp::serialize_body']
    self.assertEqual(next(f for f in fields if f.get('name') == 'Ipqp::pr_tol')['type'],
                     'double')
    def check(value):
      if isinstance(value, dict):
        for key, item in value.items():
          if key in ('condition', 'count', 'name_expression', 'type_expression') or (key == 'bind' and value.get('op') == 'select'):
            self.assertIsInstance(item, list)
          self.assertNotIn(key, ('cpp_body', 'expression', 'offset', 'sources', 'sha256', 'bind', 'type_expression'))
          check(item)
      elif isinstance(value, list):
        if value: self.assertNotEqual(value[0], 'binding')
        for item in value: check(item)
    check(scheme)
    fields = scheme['reader']['layouts']['XFunction<MXFunction,MX,MXNode>::serialize_body']
    self.assertEqual(next(f for f in fields if f.get('name') == 'XFunction::in')['type'],
                     {'name': 'std::vector', 'element': 'MX'})
    self.assertNotRegex(json.dumps(scheme), r'\b(?:Scalar|MatType|DerivedType)\b')
    fields = scheme['reader']['layouts']['FunctionInternal::serialize_body']
    condition = next(f['condition'] for f in fields if f.get('op') == 'if')
    self.assertIn(['equal', ['field', 'FunctionInternal::jit_serialize'], ['literal', 'link']], condition)
    fields = scheme['reader']['layouts']['SundialsInterface::serialize_body']
    self.assertEqual(next(f for f in fields if f.get('name') == 'SundialsInterface::abstol'),
                     {'op': 'field', 'name': 'SundialsInterface::abstol', 'type': 'double'})

  @memory_heavy()
  def test_generated_scheme_is_current(self):
    self.message("committed scheme matches the generator and current C++ sources")
    module = load_module('generate_serialization_scheme')
    self.assertEqual(module.infer_type('value', {'value': {'std::vector<T1>'}}),
                     'std::vector<double>')
    self.assertEqual(module.infer_type('value', {'value': {'T10'}}), 'T10')
    self.assertEqual(module.infer_type('value', {'value': {'std::vector<libmad_int>'}}),
                     'std::vector<casadi_int>')
    self.assertEqual(module.generate(root), json.loads((root/"misc/serialization_scheme.json").read_text()))

if __name__ == '__main__':
  unittest.main()
