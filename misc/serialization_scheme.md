# Serialization scheme

`serialization_scheme.json` is generated from the C++ sources and committed with
them. Update it after changing serialization:

```sh
python misc/generate_serialization_scheme.py
python test/python/serialization_scheme.py
```

No build, CasADi Python installation or third-party Python package is needed to
run the generator. The validation tests require the matching CasADi build.

The scheme records protocol constants, operation IDs, written/readable class
versions, serializer fields and their source expressions/types, base calls,
control-flow expressions, source locations and body hashes. It includes a named
field contract for inline and plugin serializers as well. Generation is
deterministic; the tests compare it with the committed file to catch stale output.

Every `casadiTestCase.check_serialize` validates its existing debug serialization
against the committed field contract before native deserialization and numerical
comparison. Unknown field names, known wire-type mismatches, missing supported
class versions, malformed encodings and truncated payload extents fail. Tests
also deliberately mutate the scheme to demonstrate that validation is effective.
Archive payloads are skipped by length, avoiding an additional decoded archive
copy just for validation.

This first scheme is **not a complete executable grammar**. Serializer order and
C++ branches/loops are preserved in the indexed source, but the independent
validator does not execute them or replace the native deserializer's structural
checks. Unresolved C++ types are explicit and reported by the validator; those
fields receive name/decoration validation, not inferred-type validation. Container
validation checks the outer decoration, not every element's declared C++ type.
This boundary is intentional and must not be confused with complete format
coverage. Improving type resolution and compiling more class layouts into an
executable schema can be done incrementally with the same regression tests.

`casadi-reader` vendors a copy. Its current specialized MX/Resource readers use
the protocol constants, operation IDs and class-version inventory, while still
implementing their supported positional layouts explicitly. Regenerating this
file does not automatically teach an older reader a new node layout.

The `reader` section contains executable structural layouts for casadi-reader.
`serialization_reader_layouts.py` lowers indexed serializers, inline SX serializers,
base-class calls, conditions, repeated fields and serialization helpers such as
`pack_tensors`. Shared-object flags and decorations are extracted from the wire
codec; plugin registrations and operation dispatch families come from native
source. None of these layouts evaluates an operation or reconstructs a viewer graph.

Lowering is intentionally explicit about its limits: unresolved field types,
unlowered control flow/calls, and missing dispatch cases cannot be decoded by
the reader. Native plain/debug fixture pairs in casadi-reader validate the same
layouts against both encodings, including independent names and wire tags in
the debug stream. The source index remains broader than the tested executable
layout coverage.
