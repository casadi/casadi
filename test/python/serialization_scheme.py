"""Shared serialization-scheme check, loaded only by check_serialize."""
import importlib.util
from pathlib import Path

_validator = None


def validate_serialization(text):
  global _validator
  if _validator is None:
    path = Path(__file__).resolve().parents[2]/"misc/validate_serialization_scheme.py"
    spec = importlib.util.spec_from_file_location("casadi_scheme_validator", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _validator = module.SchemeValidator()
  return _validator.validate(text)


import copy
import json
import unittest
import casadi as ca
from helpers import casadiTestCase, memory_heavy

class SerializationSchemeTests(casadiTestCase):
  def test_roundtrip_contract(self):
    self.message("check_serialize validates scheme fields for MX and SX")
    for T in [ca.MX, ca.SX]:
      x = T.sym("x", 4)
      f = ca.Function("f", [x], [ca.sin(x)+x*x])
      self.check_serialize(f, [ca.DM([1, 2, 3, 4])])
      for g in [f.map(2), f.jacobian(), f.forward(2), f.reverse(2)]:
        report = validate_serialization(g.serialize({"debug": True}))
        self.assertGreater(report['fields'], 0)
        self.assertEqual(report['unresolved_types'], [])

  def test_scheme_mutations(self):
    self.message("wrong field names, codecs and class versions fail independently")
    x = ca.MX.sym("x")
    data = ca.Function("f", [x], [ca.sin(x)]).serialize({"debug": True})
    validate_serialization(data)
    for mutation in ['field', 'codec', 'version']:
      scheme = copy.deepcopy(_validator.scheme)
      if mutation == 'field':
        del scheme['field_contracts']['ProtoFunction::name']
      elif mutation == 'codec':
        scheme['field_contracts']['ProtoFunction::name']['cpp_types'] = ['bool']
      else:
        scheme['readable_class_versions']['ProtoFunction'] = [999]
      with self.assertRaises(ValueError):
        type(_validator)(scheme).validate(data)
    with self.assertRaises(ValueError):
      validate_serialization(data[:-2])

  @memory_heavy()
  def test_generated_scheme_is_current(self):
    self.message("committed scheme matches the generator and current C++ sources")
    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location("generate_scheme", root/"misc/generate_serialization_scheme.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    self.assertEqual(module.generate(root), json.loads((root/"misc/serialization_scheme.json").read_text()))

if __name__ == '__main__':
  unittest.main()
