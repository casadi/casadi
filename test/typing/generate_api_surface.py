#!/usr/bin/env python3
"""Generate strict CasADi stubs and pyright access checks from a live module."""

from __future__ import annotations

import argparse
import ast
import difflib
import importlib
import inspect
import keyword
import re
import textwrap
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType


# These dunder methods are part of the usable Python API even though they are not
# "public" in the conventional underscore sense.
SPECIAL_METHODS = frozenset(
    {
        "__abs__",
        "__add__",
        "__array__",
        "__bool__",
        "__call__",
        "__contains__",
        "__eq__",
        "__float__",
        "__ge__",
        "__getitem__",
        "__gt__",
        "__init__",
        "__int__",
        "__iter__",
        "__le__",
        "__len__",
        "__lt__",
        "__matmul__",
        "__mul__",
        "__ne__",
        "__neg__",
        "__next__",
        "__pos__",
        "__pow__",
        "__radd__",
        "__rmatmul__",
        "__rmul__",
        "__rpow__",
        "__rsub__",
        "__rtruediv__",
        "__setitem__",
        "__sub__",
        "__truediv__",
    }
)

SIGNATURE_LINE_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*\(.*\)(?: -> .+)?$")

HEADER = [
    "from __future__ import annotations",
    "",
    "from collections.abc import Callable, Iterator, Mapping, Sequence",
    "from types import ModuleType, TracebackType",
    "from typing import IO, Never, ParamSpec, Self, TypeAlias, TypeVar, overload",
    "",
    "import numpy as np",
    "from numpy.typing import NDArray",
    "from scipy.sparse import csc_matrix",
    "",
    "_P = ParamSpec(\"_P\")",
    "_R = TypeVar(\"_R\")",
    "_T = TypeVar(\"_T\")",
    "",
    "_ScalarLike: TypeAlias = bool | int | float",
    "_ArrayKey: TypeAlias = int | slice | tuple[int | slice, int | slice]",
    "_ArrayData: TypeAlias = (",
    "    _ScalarLike",
    "    | Sequence[_ScalarLike]",
    "    | Sequence[Sequence[_ScalarLike]]",
    "    | NDArray[np.bool_]",
    "    | NDArray[np.int_]",
    "    | NDArray[np.float64]",
    ")",
    "_FunctionOutput: TypeAlias = \"DM\" | \"SX\" | \"MX\"",
    "_FunctionCallResult: TypeAlias = _FunctionOutput",
    "_SerializedValue: TypeAlias = (",
    "    bool | int | float | str | bytes | Sequence[\"_SerializedValue\"] | Mapping[str, \"_SerializedValue\"]",
    ")",
    "_PackedElement: TypeAlias = (",
    "    _SerializedValue | \"Sparsity\" | \"DM\" | \"SX\" | \"MX\" | \"Function\" | \"GenericType\" | \"Linsol\"",
    ")",
    "_PackedValue: TypeAlias = _PackedElement | Sequence[\"_PackedElement\"]",
    "_TypeNameValue: TypeAlias = _SerializedValue | tuple[\"_TypeNameValue\", ...] | NDArray[np.bool_] | NDArray[np.int_] | NDArray[np.float64]",
    "_ReturnAnnotation: TypeAlias = None | type[int] | type[\"DM\"] | type[\"SX\"] | type[\"MX\"] | type[\"Sparsity\"]",
    "",
]

RAW_TYPE_MAP = {
    "bool": "bool",
    "casadi_int": "int",
    "char const *": "str",
    "char*": "str",
    "double": "float",
    "false": "bool",
    "float": "float",
    "IM": "DM",
    "int": "int",
    "integer": "int",
    "PyObject *": "_PackedValue",
    "ConstraintType": "int",
    "SerializerBase::SerializationType": "int",
    "ptrdiff_t": "int",
    "size_t": "int",
    "std::string": "str",
    "std::size_t": "int",
    "str": "str",
    "true": "bool",
    "void": "None",
    "casadi::SerializerBase::SerializationType": "int",
    "casadi::VariableType": "int",
}

PROPERTY_TYPE_OVERRIDES = {
    ("DM", "T"): "DM",
    ("DM", "nz"): "NZproxy",
    ("DM", "shape"): "tuple[int, int]",
    ("MX", "T"): "MX",
    ("MX", "nz"): "NZproxy",
    ("MX", "shape"): "tuple[int, int]",
    ("Opti", "debug"): "OptiAdvanced",
    ("Opti", "advanced"): "OptiAdvanced",
    ("Opti", "casadi_solver"): "Function",
    ("Opti", "f"): "MX",
    ("Opti", "f_linear_scale"): "float",
    ("Opti", "g"): "MX",
    ("Opti", "g_linear_scale"): "DM",
    ("Opti", "lbg"): "MX",
    ("Opti", "lam_g"): "MX",
    ("Opti", "ng"): "int",
    ("Opti", "np"): "int",
    ("Opti", "nx"): "int",
    ("Opti", "p"): "MX",
    ("Opti", "ubg"): "MX",
    ("Opti", "x"): "MX",
    ("Opti", "x_linear_scale"): "DM",
    ("Opti", "x_linear_scale_offset"): "DM",
    ("OptiAdvanced", "advanced"): "OptiAdvanced",
    ("OptiAdvanced", "casadi_solver"): "Function",
    ("OptiAdvanced", "debug"): "OptiAdvanced",
    ("OptiAdvanced", "f"): "MX",
    ("OptiAdvanced", "f_linear_scale"): "float",
    ("OptiAdvanced", "g"): "MX",
    ("OptiAdvanced", "g_linear_scale"): "DM",
    ("OptiAdvanced", "lam_g"): "MX",
    ("OptiAdvanced", "lbg"): "MX",
    ("OptiAdvanced", "ng"): "int",
    ("OptiAdvanced", "np"): "int",
    ("OptiAdvanced", "nx"): "int",
    ("OptiAdvanced", "p"): "MX",
    ("OptiAdvanced", "ubg"): "MX",
    ("OptiAdvanced", "x"): "MX",
    ("OptiAdvanced", "x_linear_scale"): "DM",
    ("OptiAdvanced", "x_linear_scale_offset"): "DM",
    ("OptiSol", "debug"): "OptiAdvanced",
    ("OptiSol", "opti"): "Opti",
    ("SX", "T"): "SX",
    ("SX", "nz"): "NZproxy",
    ("SX", "shape"): "tuple[int, int]",
    ("Sparsity", "T"): "Sparsity",
    ("Sparsity", "shape"): "tuple[int, int]",
}

METHOD_OVERRIDES = {
    ("Callback", "__call__"): [
        "def __call__(self, *args: \"DM | SX | MX\") -> _FunctionOutput: ...",
        "def __call__(self, **kwargs: \"DM | SX | MX\") -> Mapping[str, _FunctionOutput]: ...",
    ],
    ("Callback", "cache"): ["def cache(self) -> Mapping[str, GenericType]: ..."],
    ("Callback", "repr"): ["def repr(self) -> str: ..."],
    ("DM", "__abs__"): ["def __abs__(self) -> DM: ..."],
    ("DM", "__array__"): ["def __array__(self, *args: str | tuple[np.ufunc, tuple[NDArray[np.float64], ...], int] | tuple[np.ufunc, tuple[NDArray[np.float64], ...], int, NDArray[np.float64]]) -> NDArray[np.float64]: ..."],
    ("DM", "__bool__"): ["def __bool__(self) -> bool: ..."],
    ("DM", "__getitem__"): ["def __getitem__(self, s: _ArrayKey) -> DM: ..."],
    ("DM", "__iter__"): ["def __iter__(self) -> Iterator[Never]: ..."],
    ("DM", "__setitem__"): ["def __setitem__(self, s: _ArrayKey, val: _ArrayData | DM) -> None: ..."],
    ("DM", "repr"): ["def repr(self) -> str: ..."],
    ("DM", "sym"): [
        "@staticmethod",
        "def sym(name: str) -> DM: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int = ...) -> DM: ...",
        "@staticmethod",
        "def sym(name: str, rc: tuple[int, int]) -> DM: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity) -> DM: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity, p: int) -> Sequence[DM]: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int, p: int) -> Sequence[DM]: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity, p: int, r: int) -> Sequence[Sequence[DM]]: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int, p: int, r: int) -> Sequence[Sequence[DM]]: ...",
    ],
    ("DM", "toarray"): ["def toarray(self, simplify: bool = ...) -> float | NDArray[np.float64]: ..."],
    ("DM", "tocsc"): ["def tocsc(self) -> csc_matrix: ..."],
    ("DeserializerBase", "__init__"): ["def __init__(self) -> None: ..."],
    ("DeserializerBase", "unpack"): [
        "def unpack(self) -> _PackedValue: ..."
    ],
    ("Function", "__call__"): [
        "def __call__(self, *args: \"DM | SX | MX\") -> _FunctionOutput: ...",
        "def __call__(self, **kwargs: \"DM | SX | MX\") -> Mapping[str, _FunctionOutput]: ...",
    ],
    ("Function", "bspline"): [
        "@staticmethod",
        "def bspline(name: str, knots: Sequence[Sequence[float]], coeffs: Sequence[float], degree: Sequence[int], m: int, opts: Mapping[str, GenericType]) -> Function: ...",
    ],
    ("Function", "conditional"): [
        "@staticmethod",
        "def conditional(name: str, f: Function, opts: Mapping[str, GenericType]) -> Function: ...",
        "@staticmethod",
        "def conditional(name: str, f: Sequence[Function], f_def: Function, opts: Mapping[str, GenericType]) -> Function: ...",
    ],
    ("Function", "__init__"): [
        "def __init__(self) -> None: ...",
        "def __init__(self, fname: str) -> None: ...",
        "def __init__(self, other: Function) -> None: ...",
        "def __init__(self, name: str, ex_in: Sequence[SX], ex_out: Sequence[SX], opts: Mapping[str, GenericType] = ...) -> None: ...",
        "def __init__(self, name: str, ex_in: Sequence[MX], ex_out: Sequence[MX], opts: Mapping[str, GenericType] = ...) -> None: ...",
        "def __init__(self, name: str, dict: Mapping[str, SX], name_in: Sequence[str], name_out: Sequence[str], opts: Mapping[str, GenericType] = ...) -> None: ...",
        "def __init__(self, name: str, dict: Mapping[str, MX], name_in: Sequence[str], name_out: Sequence[str], opts: Mapping[str, GenericType] = ...) -> None: ...",
        "def __init__(self, name: str, ex_in: Sequence[SX], ex_out: Sequence[SX], name_in: Sequence[str], name_out: Sequence[str], opts: Mapping[str, GenericType] = ...) -> None: ...",
        "def __init__(self, name: str, ex_in: Sequence[MX], ex_out: Sequence[MX], name_in: Sequence[str], name_out: Sequence[str], opts: Mapping[str, GenericType] = ...) -> None: ...",
    ],
    ("Function", "buffer"): ["def buffer(self) -> tuple[FunctionBuffer, Callable[..., None]]: ..."],
    ("Function", "cache"): ["def cache(self) -> Mapping[str, GenericType]: ..."],
    ("Function", "repr"): ["def repr(self) -> str: ..."],
    ("MX", "sym"): [
        "@staticmethod",
        "def sym(name: str) -> MX: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int = ...) -> MX: ...",
        "@staticmethod",
        "def sym(name: str, rc: tuple[int, int]) -> MX: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity) -> MX: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity, p: int) -> Sequence[MX]: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int, p: int) -> Sequence[MX]: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity, p: int, r: int) -> Sequence[Sequence[MX]]: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int, p: int, r: int) -> Sequence[Sequence[MX]]: ...",
    ],
    ("NZproxy", "__getitem__"): ["def __getitem__(self, s: _ArrayKey) -> DM | SX | MX: ..."],
    ("NZproxy", "__init__"): ["def __init__(self, matrix: DM | SX | MX) -> None: ..."],
    ("NZproxy", "__iter__"): ["def __iter__(self) -> Iterator[DM | SX | MX]: ..."],
    ("NZproxy", "__len__"): ["def __len__(self) -> int: ..."],
    ("NZproxy", "__setitem__"): ["def __setitem__(self, s: _ArrayKey, val: _ArrayData | DM | SX | MX) -> None: ..."],
    ("Opti", "callback"): ["def callback(self, fh: Callable[[Function], None] | None = ...) -> None: ..."],
    ("Opti", "variable"): [
        "def variable(self) -> MX: ...",
        "def variable(self, n: int, m: int = ..., attribute: str = ...) -> MX: ...",
        "def variable(self, sp: Sparsity, attribute: str = ...) -> MX: ...",
        "def variable(self, symbol: MX, attribute: str = ...) -> MX: ...",
    ],
    ("OptiAdvanced", "callback"): ["def callback(self, fh: Callable[[Function], None] | None = ...) -> None: ..."],
    ("OptiAdvanced", "variable"): [
        "def variable(self) -> MX: ...",
        "def variable(self, n: int, m: int = ..., attribute: str = ...) -> MX: ...",
        "def variable(self, sp: Sparsity, attribute: str = ...) -> MX: ...",
        "def variable(self, symbol: MX, attribute: str = ...) -> MX: ...",
    ],
    ("SerializerBase", "__init__"): ["def __init__(self) -> None: ..."],
    ("SX", "sym"): [
        "@staticmethod",
        "def sym(name: str) -> SX: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int = ...) -> SX: ...",
        "@staticmethod",
        "def sym(name: str, rc: tuple[int, int]) -> SX: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity) -> SX: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity, p: int) -> Sequence[SX]: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int, p: int) -> Sequence[SX]: ...",
        "@staticmethod",
        "def sym(name: str, sp: Sparsity, p: int, r: int) -> Sequence[Sequence[SX]]: ...",
        "@staticmethod",
        "def sym(name: str, nrow: int, ncol: int, p: int, r: int) -> Sequence[Sequence[SX]]: ...",
    ],
    ("Sparsity", "__array__"): ["def __array__(self, *args: str | tuple[np.ufunc, tuple[NDArray[np.float64], ...], int] | tuple[np.ufunc, tuple[NDArray[np.float64], ...], int, NDArray[np.float64]]) -> NDArray[np.float64]: ..."],
    ("Sparsity", "__getitem__"): ["def __getitem__(self, s: _ArrayKey) -> Sparsity: ..."],
    ("Sparsity", "repr"): ["def repr(self) -> str: ..."],
    ("SwigPyIterator", "__init__"): ["def __init__(self, *args: Never, **kwargs: Never) -> None: ..."],
    ("SwigPyIterator", "__iter__"): ["def __iter__(self) -> Self: ..."],
}

FUNCTION_OVERRIDES = {
    "Callback_swigregister": ["def Callback_swigregister(cls: type[Callback]) -> None: ..."],
    "CasadiMeta_swigregister": ["def CasadiMeta_swigregister(cls: type[CasadiMeta]) -> None: ..."],
    "CodeGenerator_swigregister": ["def CodeGenerator_swigregister(cls: type[CodeGenerator]) -> None: ..."],
    "DM_from_array": ["def DM_from_array(m: _ArrayData, check_only: bool = ...) -> DM: ..."],
    "DM_from_csc": ["def DM_from_csc(m: csc_matrix, check_only: bool = ...) -> DM: ..."],
    "DM_swigregister": ["def DM_swigregister(cls: type[DM]) -> None: ..."],
    "DaeBuilder_swigregister": ["def DaeBuilder_swigregister(cls: type[DaeBuilder]) -> None: ..."],
    "DeserializerBase_swigregister": ["def DeserializerBase_swigregister(cls: type[DeserializerBase]) -> None: ..."],
    "FileDeserializer_swigregister": ["def FileDeserializer_swigregister(cls: type[FileDeserializer]) -> None: ..."],
    "FileSerializer_swigregister": ["def FileSerializer_swigregister(cls: type[FileSerializer]) -> None: ..."],
    "FunctionBuffer_swigregister": ["def FunctionBuffer_swigregister(cls: type[FunctionBuffer]) -> None: ..."],
    "Function_bspline": [
        "def Function_bspline(name: str, knots: Sequence[Sequence[float]], coeffs: Sequence[float], degree: Sequence[int], m: int, opts: Mapping[str, GenericType]) -> Function: ...",
    ],
    "Function_conditional": [
        "def Function_conditional(name: str, f: Function, opts: Mapping[str, GenericType]) -> Function: ...",
        "def Function_conditional(name: str, f: Sequence[Function], f_def: Function, opts: Mapping[str, GenericType]) -> Function: ...",
    ],
    "Function_swigregister": ["def Function_swigregister(cls: type[Function]) -> None: ..."],
    "GenDM_swigregister": ["def GenDM_swigregister(cls: type[GenDM]) -> None: ..."],
    "GenMX_swigregister": ["def GenMX_swigregister(cls: type[GenMX]) -> None: ..."],
    "GenSX_swigregister": ["def GenSX_swigregister(cls: type[GenSX]) -> None: ..."],
    "GenSharedObject_swigregister": ["def GenSharedObject_swigregister(cls: type[GenSharedObject]) -> None: ..."],
    "GenWeakRef_swigregister": ["def GenWeakRef_swigregister(cls: type[GenWeakRef]) -> None: ..."],
    "GenericExpressionCommon_swigregister": ["def GenericExpressionCommon_swigregister(cls: type[GenericExpressionCommon]) -> None: ..."],
    "GenericMatrixCommon_swigregister": ["def GenericMatrixCommon_swigregister(cls: type[GenericMatrixCommon]) -> None: ..."],
    "GenericType_swigregister": ["def GenericType_swigregister(cls: type[GenericType]) -> None: ..."],
    "GlobalOptions_swigregister": ["def GlobalOptions_swigregister(cls: type[GlobalOptions]) -> None: ..."],
    "IM_from_array": ["def IM_from_array(m: _ArrayData, check_only: bool = ...) -> DM: ..."],
    "Importer_swigregister": ["def Importer_swigregister(cls: type[Importer]) -> None: ..."],
    "IndexAbstraction_swigregister": ["def IndexAbstraction_swigregister(cls: type[IndexAbstraction]) -> None: ..."],
    "Linsol_swigregister": ["def Linsol_swigregister(cls: type[Linsol]) -> None: ..."],
    "MX_swigregister": ["def MX_swigregister(cls: type[MX]) -> None: ..."],
    "MatrixCommon_swigregister": ["def MatrixCommon_swigregister(cls: type[MatrixCommon]) -> None: ..."],
    "MetaCon_swigregister": ["def MetaCon_swigregister(cls: type[MetaCon]) -> None: ..."],
    "MetaVar_swigregister": ["def MetaVar_swigregister(cls: type[MetaVar]) -> None: ..."],
    "NlpBuilder_swigregister": ["def NlpBuilder_swigregister(cls: type[NlpBuilder]) -> None: ..."],
    "Opti_swigregister": ["def Opti_swigregister(cls: type[Opti]) -> None: ..."],
    "OptiAdvanced_swigregister": ["def OptiAdvanced_swigregister(cls: type[OptiAdvanced]) -> None: ..."],
    "OptiCallback_swigregister": ["def OptiCallback_swigregister(cls: type[OptiCallback]) -> None: ..."],
    "OptiSol_swigregister": ["def OptiSol_swigregister(cls: type[OptiSol]) -> None: ..."],
    "Options_swigregister": ["def Options_swigregister(cls: type[Options]) -> None: ..."],
    "PrintableCommon_swigregister": ["def PrintableCommon_swigregister(cls: type[PrintableCommon]) -> None: ..."],
    "Resource_swigregister": ["def Resource_swigregister(cls: type[Resource]) -> None: ..."],
    "SX_from_array": ["def SX_from_array(m: _ArrayData, check_only: bool = ...) -> SX: ..."],
    "SXElem_swigregister": ["def SXElem_swigregister(cls: type[SXElem]) -> None: ..."],
    "SX_swigregister": ["def SX_swigregister(cls: type[SX]) -> None: ..."],
    "SerializerBase_swigregister": ["def SerializerBase_swigregister(cls: type[SerializerBase]) -> None: ..."],
    "SharedObject_swigregister": ["def SharedObject_swigregister(cls: type[SharedObject]) -> None: ..."],
    "Slice_swigregister": ["def Slice_swigregister(cls: type[Slice]) -> None: ..."],
    "Sparsity_swigregister": ["def Sparsity_swigregister(cls: type[Sparsity]) -> None: ..."],
    "SparsityInterfaceCommon_swigregister": ["def SparsityInterfaceCommon_swigregister(cls: type[SparsityInterfaceCommon]) -> None: ..."],
    "StreamStateGuard_swigregister": ["def StreamStateGuard_swigregister(cls: type[StreamStateGuard]) -> None: ..."],
    "StringDeserializer_swigregister": ["def StringDeserializer_swigregister(cls: type[StringDeserializer]) -> None: ..."],
    "StringSerializer_swigregister": ["def StringSerializer_swigregister(cls: type[StringSerializer]) -> None: ..."],
    "SwigPyIterator_swigregister": ["def SwigPyIterator_swigregister(cls: type[SwigPyIterator]) -> None: ..."],
    "WeakRef_swigregister": ["def WeakRef_swigregister(cls: type[WeakRef]) -> None: ..."],
    "XmlFile_swigregister": ["def XmlFile_swigregister(cls: type[XmlFile]) -> None: ..."],
    "attach_return_type": ["def attach_return_type(f: Callable[_P, _R], t: _ReturnAnnotation) -> Callable[_P, _R]: ..."],
    "dcat": [
        "def dcat(args: Sequence[_ScalarLike | DM]) -> DM: ...",
        "def dcat(args: Sequence[_ScalarLike | SX]) -> SX: ...",
        "def dcat(args: Sequence[_ScalarLike | MX]) -> MX: ...",
        "def dcat(args: Sequence[_ScalarLike | Sparsity]) -> Sparsity: ...",
    ],
    "diagcat": [
        "def diagcat(*args: _ScalarLike | DM) -> DM: ...",
        "def diagcat(*args: _ScalarLike | SX) -> SX: ...",
        "def diagcat(*args: _ScalarLike | MX) -> MX: ...",
        "def diagcat(*args: _ScalarLike | Sparsity) -> Sparsity: ...",
    ],
    "global_pickle_context": [],
    "global_unpickle_context": [],
    "hcat": [
        "def hcat(args: Sequence[_ScalarLike | DM]) -> DM: ...",
        "def hcat(args: Sequence[_ScalarLike | SX]) -> SX: ...",
        "def hcat(args: Sequence[_ScalarLike | MX]) -> MX: ...",
        "def hcat(args: Sequence[_ScalarLike | Sparsity]) -> Sparsity: ...",
    ],
    "horzcat": [
        "def horzcat(*args: _ScalarLike | DM) -> DM: ...",
        "def horzcat(*args: _ScalarLike | SX) -> SX: ...",
        "def horzcat(*args: _ScalarLike | MX) -> MX: ...",
        "def horzcat(*args: _ScalarLike | Sparsity) -> Sparsity: ...",
    ],
    "pycallback": ["def pycallback(f: Callable[_P, _R]) -> Callable[_P, _R]: ..."],
    "pyevaluate": ["def pyevaluate(f: Callable[_P, _R]) -> Callable[_P, _R]: ..."],
    "swig_typename_convertor_python2cpp": ["def swig_typename_convertor_python2cpp(a: _TypeNameValue) -> str: ..."],
    "swigtypeconvertor": ["def swigtypeconvertor(*args: _TypeNameValue) -> str: ..."],
    "vcat": [
        "def vcat(args: Sequence[_ScalarLike | DM]) -> DM: ...",
        "def vcat(args: Sequence[_ScalarLike | SX]) -> SX: ...",
        "def vcat(args: Sequence[_ScalarLike | MX]) -> MX: ...",
        "def vcat(args: Sequence[_ScalarLike | Sparsity]) -> Sparsity: ...",
    ],
    "veccat": [
        "def veccat(*args: _ScalarLike | DM) -> DM: ...",
        "def veccat(*args: _ScalarLike | SX) -> SX: ...",
        "def veccat(*args: _ScalarLike | MX) -> MX: ...",
        "def veccat(*args: _ScalarLike | Sparsity) -> Sparsity: ...",
    ],
    "vertcat": [
        "def vertcat(*args: _ScalarLike | DM) -> DM: ...",
        "def vertcat(*args: _ScalarLike | SX) -> SX: ...",
        "def vertcat(*args: _ScalarLike | MX) -> MX: ...",
        "def vertcat(*args: _ScalarLike | Sparsity) -> Sparsity: ...",
    ],
    "vvcat": [
        "def vvcat(args: Sequence[_ScalarLike | DM]) -> DM: ...",
        "def vvcat(args: Sequence[_ScalarLike | SX]) -> SX: ...",
        "def vvcat(args: Sequence[_ScalarLike | MX]) -> MX: ...",
        "def vvcat(args: Sequence[_ScalarLike | Sparsity]) -> Sparsity: ...",
    ],
    "weakref_proxy": ["def weakref_proxy(value: _T, callback: Callable[[WeakRef], None] | None = ...) -> _T: ..."],
    "wrapper": ["def wrapper(f: Callable[_P, _R], warning: str, error: bool = ...) -> Callable[_P, _R]: ..."],
}

CLASS_OVERRIDES = {
    "backup_object": [
        "class backup_object:",
        "    ...",
    ],
    "Deprecate": [
        "class Deprecate:",
        "    def __new__(cls, o: _T, warning: str, error: bool = ...) -> _T: ...",
    ],
    "NZproxy": [
        "class NZproxy:",
        "    def __init__(self, matrix: DM | SX | MX) -> None: ...",
        "    def __getitem__(self, s: _ArrayKey) -> DM | SX | MX: ...",
        "    def __setitem__(self, s: _ArrayKey, val: _ArrayData | DM | SX | MX) -> None: ...",
        "    def __len__(self) -> int: ...",
        "    def __iter__(self) -> Iterator[DM | SX | MX]: ...",
    ],
    "PyFunction": [
        "class PyFunction:",
        "    def __call__(",
        "        self,",
        "        name: str,",
        "        obj: SupportsEvaluate,",
        "        inputs: Sequence[DM | SX | MX],",
        "        outputs: Sequence[DM | SX | MX],",
        "        opts: Mapping[str, GenericType] = ...,",
        "    ) -> Function: ...",
    ],
    "global_pickle_context": [
        "class global_pickle_context:",
        "    ctx: StringSerializer",
        "    def __enter__(self) -> StringSerializer: ...",
        "    def __exit__(",
        "        self,",
        "        exc_type: type[BaseException] | None,",
        "        exc: BaseException | None,",
        "        tb: TracebackType | None,",
        "    ) -> None: ...",
    ],
    "global_unpickle_context": [
        "class global_unpickle_context:",
        "    ctx: StringDeserializer",
        "    def __enter__(self) -> StringDeserializer: ...",
        "    def __exit__(",
        "        self,",
        "        exc_type: type[BaseException] | None,",
        "        exc: BaseException | None,",
        "        tb: TracebackType | None,",
        "    ) -> None: ...",
    ],
}

UNARY_SELF_METHODS = {
    "__abs__",
    "__neg__",
    "__pos__",
    "arccos",
    "arccosh",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctanh",
    "ceil",
    "cos",
    "cosh",
    "erf",
    "exp",
    "expm1",
    "fabs",
    "floor",
    "log",
    "log10",
    "log1p",
    "sign",
    "sin",
    "sinh",
    "sqrt",
    "tan",
    "tanh",
}

BINARY_SELF_METHODS = {
    "__add__",
    "__eq__",
    "__ge__",
    "__getitem__",
    "__gt__",
    "__le__",
    "__lt__",
    "__matmul__",
    "__mul__",
    "__ne__",
    "__pow__",
    "__radd__",
    "__rmatmul__",
    "__rmul__",
    "__rpow__",
    "__rsub__",
    "__rtruediv__",
    "__setitem__",
    "__sub__",
    "__truediv__",
    "constpow",
    "copysign",
    "fmax",
    "fmin",
    "fmod",
    "hypot",
    "logic_and",
    "logic_or",
    "rconstpow",
    "rcopysign",
    "remainder",
}


@dataclass(frozen=True)
class Parameter:
    name: str
    type_expr: str | None
    default: bool = False
    vararg: bool = False
    kwarg: bool = False

    def render(self) -> str:
        prefix = "**" if self.kwarg else "*" if self.vararg else ""
        if self.type_expr is None:
            return f"{prefix}{self.name}"
        suffix = " = ..." if self.default else ""
        return f"{prefix}{self.name}: {self.type_expr}{suffix}"


@dataclass(frozen=True)
class Signature:
    parameters: tuple[Parameter, ...]
    return_type: str
    decorators: tuple[str, ...] = ()

    def render(self, name: str, indent: str = "") -> list[str]:
        decorators = tuple(
            decorator
            for decorator in ("@overload", "@staticmethod")
            if decorator in self.decorators
        ) + tuple(
            decorator
            for decorator in self.decorators
            if decorator not in {"@overload", "@staticmethod"}
        )
        lines: list[str] = []
        if len(self.parameters) > 4:
            lines.extend(f"{indent}{decorator}" for decorator in decorators)
            lines.append(f"{indent}def {name}(")
            for parameter in self.parameters:
                lines.append(f"{indent}    {parameter.render()},")
            lines.append(f"{indent}) -> {self.return_type}: ...")
            return lines

        lines.extend(f"{indent}{decorator}" for decorator in decorators)
        params = ", ".join(parameter.render() for parameter in self.parameters)
        lines.append(f"{indent}def {name}({params}) -> {self.return_type}: ...")
        return lines


def exported_names(value: ModuleType | type) -> list[str]:
    return sorted(
        name
        for name in dir(value)
        if not name.startswith("_")
        or (
            name in SPECIAL_METHODS
            and not isinstance(value, type)
            or name in SPECIAL_METHODS
            and isinstance(value, type)
            and any(name in base.__dict__ for base in value.__mro__ if base is not object)
        )
    )


def split_top_level(text: str, delimiter: str = ",") -> list[str]:
    if not text.strip():
        return []

    parts: list[str] = []
    current: list[str] = []
    round_depth = 0
    square_depth = 0
    angle_depth = 0

    for char in text:
        if char == "(":
            round_depth += 1
        elif char == ")":
            round_depth -= 1
        elif char == "[":
            square_depth += 1
        elif char == "]":
            square_depth -= 1
        elif char == "<":
            angle_depth += 1
        elif char == ">":
            angle_depth -= 1

        if (
            char == delimiter
            and round_depth == 0
            and square_depth == 0
            and angle_depth == 0
        ):
            parts.append("".join(current).strip())
            current = []
            continue

        current.append(char)

    parts.append("".join(current).strip())
    return [part for part in parts if part]


def sanitize_name(name: str) -> str:
    sanitized = name.replace("-", "_").replace(":", "_")
    if keyword.iskeyword(sanitized):
        sanitized += "_"
    return sanitized


def strip_namespace(type_text: str) -> str:
    text = type_text.strip()
    text = text.replace("casadi::", "")
    text = text.replace("std::allocator< ", "")
    text = text.replace("std::allocator<", "")
    text = text.replace(" > >", ">>")
    text = text.replace(" > >", ">>")
    text = text.replace(" const &", "")
    text = text.replace(" &", "")
    text = text.replace(" const *", "*")
    return text


def normalize_type(type_text: str, known_types: set[str]) -> str:
    text = strip_namespace(type_text)
    text = re.sub(r"\s+", " ", text).strip()

    if text in RAW_TYPE_MAP:
        return RAW_TYPE_MAP[text]

    if text.endswith(" OUTPUT") or re.search(r" OUTPUT\d+$", text):
        return normalize_type(re.sub(r" OUTPUT\d*$", "", text), known_types)

    if text.endswith(" INOUT") or re.search(r" INOUT\d+$", text):
        return normalize_type(re.sub(r" INOUT\d*$", "", text), known_types)

    if text.endswith("]") and text.startswith("[["):
        return f"Sequence[{normalize_type(text[1:-1], known_types)}]"

    if text.endswith("]") and text.startswith("["):
        return f"Sequence[{normalize_type(text[1:-1], known_types)}]"

    if text == "(int,int)":
        return "tuple[int, int]"

    if text.startswith("(") and text.endswith(")"):
        inner = text[1:-1].strip()
        parts = [normalize_type(part, known_types) for part in split_top_level(inner)]
        return f"tuple[{', '.join(parts)}]"

    if text.startswith("dict:"):
        return f"Mapping[str, {normalize_type(text[5:], known_types)}]"

    if text == "dict":
        return "Mapping[str, GenericType]"

    if text.startswith("std::pair<") and text.endswith(">"):
        inner = text[len("std::pair<") : -1]
        left, right = split_top_level(inner)
        return f"tuple[{normalize_type(left, known_types)}, {normalize_type(right, known_types)}]"

    if text.startswith("std::vector<") and text.endswith(">"):
        inner = text[len("std::vector<") : -1]
        element = split_top_level(inner)[0]
        return f"Sequence[{normalize_type(element, known_types)}]"

    if text == "bytes_or_buffer[":
        return "bytes | bytearray | memoryview"

    if text == "encoding[" or text == "errors]]":
        return "str"

    if text.startswith("memoryview"):
        return "memoryview"

    if text in {"double const **", "double **", "double**"}:
        return "Sequence[Sequence[float]]"

    if text in {"double *", "double*"}:
        return "Sequence[float]"

    if text in {"istream *", "istream", "std::istream", "std::ostream"}:
        return "IO[str]"

    if text == "DeserializingStream":
        return "DeserializerBase"

    if text == "SerializingStream":
        return "SerializerBase"

    if text == "casadi_int *":
        return "Sequence[int]"

    if text in {"SharedObjectInternal *", "SharedObjectInternal*"}:
        return "SharedObject"

    if text == "object=''" or text == "object":
        return "str"

    if text.startswith("base="):
        return "int"

    if text in known_types:
        return text

    if text.endswith("Type"):
        return "int"

    if re.fullmatch(r"-?\d+", text):
        return "int"

    raise ValueError(f"unknown type syntax: {type_text!r}")


def split_typed_name(token: str) -> tuple[str, str]:
    round_depth = 0
    square_depth = 0
    angle_depth = 0
    for index in range(len(token) - 1, -1, -1):
        char = token[index]
        if char == ")":
            round_depth += 1
        elif char == "(":
            round_depth -= 1
        elif char == "]":
            square_depth += 1
        elif char == "[":
            square_depth -= 1
        elif char == ">":
            angle_depth += 1
        elif char == "<":
            angle_depth -= 1
        elif char.isspace() and round_depth == 0 and square_depth == 0 and angle_depth == 0:
            return token[:index].strip(), token[index + 1 :].strip()
    return "", token.strip()


def parse_parameter(token: str, known_types: set[str], index: int) -> Parameter:
    if token == "self":
        return Parameter("self", None)

    type_text, name_text = split_typed_name(token)
    if not type_text:
        if re.fullmatch(r"-?\d+", token):
            return Parameter(f"arg{index}", "int")
        if token in {"true", "false"}:
            return Parameter(f"arg{index}", "bool")
        if token == "*args":
            return Parameter("args", "Never", vararg=True)
        if token == "**kwargs":
            return Parameter("kwargs", "Never", kwarg=True)
        raise ValueError(f"could not split parameter token: {token!r}")

    if "=" in name_text:
        name_text = name_text.split("=", 1)[0].strip()
        default = True
    else:
        default = False

    return Parameter(
        sanitize_name(name_text),
        normalize_type(type_text, known_types),
        default=default,
    )


def parse_signature_line(
    signature_text: str,
    target_name: str,
    known_types: set[str],
    *,
    constructor: bool,
    method: bool,
    static: bool,
) -> Signature:
    if "->" in signature_text:
        call_text, return_text = signature_text.split("->", 1)
        return_type = normalize_type(return_text.strip(), known_types)
    else:
        call_text = signature_text
        return_type = "None"

    params_text = call_text[call_text.index("(") + 1 : call_text.rindex(")")]
    parameters = [
        parse_parameter(token, known_types, index)
        for index, token in enumerate(split_top_level(params_text))
    ]
    if not method and parameters and parameters[0].name == "self":
        parameters = parameters[1:]
    if constructor and (not parameters or parameters[0].name != "self"):
        parameters.insert(0, Parameter("self", None))
    if method and not static and not constructor and (not parameters or parameters[0].name != "self"):
        parameters.insert(0, Parameter("self", None))

    decorators: tuple[str, ...] = ("@staticmethod",) if static else ()
    return Signature(tuple(parameters), return_type, decorators)


def extract_doc_signatures(
    target_name: str,
    doc: str,
    known_types: set[str],
    *,
    constructor: bool,
    method: bool,
    static: bool,
) -> list[Signature]:
    signatures: list[Signature] = []
    seen: set[str] = set()
    for raw_line in doc.splitlines():
        line = raw_line.strip()
        if not line.startswith(f"{target_name}("):
            continue
        if not SIGNATURE_LINE_RE.fullmatch(line):
            continue
        if line.count("->") > 1:
            continue
        if "->" in line:
            tail = line.split("->", 1)[1].strip()
            if "(" in tail and not tail.startswith("("):
                continue
        if line in seen:
            continue
        seen.add(line)
        try:
            signatures.append(
                parse_signature_line(
                    line,
                    target_name,
                    known_types,
                    constructor=constructor,
                    method=method,
                    static=static,
                )
            )
        except ValueError:
            continue
    return collapse_signatures(signatures)


def collapse_signatures(signatures: list[Signature]) -> list[Signature]:
    by_parameters: dict[tuple[Parameter, ...], Signature] = {}
    ordered: list[tuple[Parameter, ...]] = []
    for signature in signatures:
        key = signature.parameters
        current = by_parameters.get(key)
        if current is None:
            by_parameters[key] = signature
            ordered.append(key)
            continue
        if current.return_type == "None" and signature.return_type != "None":
            by_parameters[key] = signature
    return [by_parameters[key] for key in ordered]


def alias_target_from_source(source: str) -> tuple[str | None, list[str]]:
    try:
        node = ast.parse(textwrap.dedent(source))
    except SyntaxError:
        return None, []
    function_node: ast.FunctionDef | ast.AsyncFunctionDef | None = None
    lambda_node: ast.Lambda | None = None

    if len(node.body) == 1 and isinstance(node.body[0], ast.Assign):
        value = node.body[0].value
        if isinstance(value, ast.Lambda):
            lambda_node = value
    elif len(node.body) == 1 and isinstance(node.body[0], (ast.FunctionDef, ast.AsyncFunctionDef)):
        function_node = node.body[0]

    def parse_call(call: ast.Call) -> tuple[str | None, list[str]]:
        if isinstance(call.func, ast.Attribute):
            if isinstance(call.func.value, ast.Name):
                target = f"{call.func.value.id}.{call.func.attr}"
            else:
                return None, []
        elif isinstance(call.func, ast.Name):
            target = call.func.id
        else:
            return None, []

        arguments: list[str] = []
        for arg in call.args:
            if isinstance(arg, ast.Name):
                arguments.append(arg.id)
            elif isinstance(arg, ast.Starred) and isinstance(arg.value, ast.Name):
                arguments.append(f"*{arg.value.id}")
            else:
                arguments.append(ast.unparse(arg))
        return target, arguments

    if lambda_node is not None and isinstance(lambda_node.body, ast.Call):
        return parse_call(lambda_node.body)

    if function_node is None:
        return None, []

    for statement in function_node.body:
        if isinstance(statement, ast.Return) and isinstance(statement.value, ast.Call):
            return parse_call(statement.value)
        if isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Call):
            targets = [target.id for target in statement.targets if isinstance(target, ast.Name)]
            if targets == ["ret"]:
                return parse_call(statement.value)
    return None, []


def apply_call_argument_order(
    target: Signature,
    wrapper_param_names: list[str],
    call_arguments: list[str],
    *,
    method: bool,
) -> Signature:
    if any(argument.startswith("*") for argument in call_arguments):
        return target

    target_parameters = list(target.parameters)
    if len(target_parameters) != len(call_arguments):
        return target

    type_by_wrapper_name: dict[str, str] = {}
    for parameter, argument in zip(target_parameters, call_arguments):
        if argument == "self":
            continue
        type_by_wrapper_name[argument] = parameter.type_expr or "Never"

    ordered_parameters: list[Parameter] = []
    for index, wrapper_name in enumerate(wrapper_param_names):
        if method and index == 0:
            ordered_parameters.append(Parameter("self", None))
            continue
        param_type = type_by_wrapper_name.get(wrapper_name)
        if param_type is None:
            return target
        ordered_parameters.append(Parameter(sanitize_name(wrapper_name), param_type))

    return Signature(tuple(ordered_parameters), target.return_type)


def signature_from_alias_target(
    module: ModuleType,
    owner_name: str | None,
    member_name: str,
    source: str,
    known_types: set[str],
    resolving: set[tuple[str | None, str]],
) -> list[Signature]:
    target_name, call_arguments = alias_target_from_source(source)
    if target_name is None:
        return []

    if (
        owner_name is not None
        and target_name.startswith("_casadi.")
        and member_name in UNARY_SELF_METHODS | BINARY_SELF_METHODS
    ):
        return []

    if target_name.startswith("_casadi."):
        resolved_owner = None
        resolved_name = target_name.split(".", 1)[1]
    elif target_name.startswith("self."):
        if owner_name is None:
            return []
        resolved_owner = owner_name
        resolved_name = target_name.split(".", 1)[1]
    else:
        resolved_owner = None
        resolved_name = target_name

    if resolved_owner is None:
        target_obj = getattr(module, resolved_name, None)
        if target_obj is None:
            return []
        signatures = render_callable_signatures(
            module,
            None,
            resolved_name,
            target_obj,
            known_types,
            resolving,
        )
    else:
        owner = getattr(module, resolved_owner)
        target_obj = getattr(owner, resolved_name, None)
        if target_obj is None:
            return []
        signatures = render_callable_signatures(
            module,
            resolved_owner,
            resolved_name,
            target_obj,
            known_types,
            resolving,
        )

    try:
        wrapper_signature = inspect.signature(getattr(getattr(module, owner_name), member_name) if owner_name else getattr(module, member_name))
    except (TypeError, ValueError, AttributeError):
        return signatures

    wrapper_param_names = list(wrapper_signature.parameters)
    method = owner_name is not None
    return [
        apply_call_argument_order(signature, wrapper_param_names, call_arguments, method=method)
        for signature in signatures
    ]


def render_callable_signatures(
    module: ModuleType,
    owner_name: str | None,
    member_name: str,
    callable_value,
    known_types: set[str],
    resolving: set[tuple[str | None, str]] | None = None,
) -> list[Signature]:
    if resolving is None:
        resolving = set()
    key = (owner_name, member_name)
    if key in resolving:
        return []
    resolving = set(resolving)
    resolving.add(key)

    static = False
    constructor = member_name == "__init__"
    if owner_name is not None:
        descriptor = inspect.getattr_static(getattr(module, owner_name), member_name)
        static = isinstance(descriptor, staticmethod)

    override_key = (owner_name, member_name) if owner_name is not None else member_name
    if owner_name is not None and override_key in METHOD_OVERRIDES:
        return parse_rendered_signatures(
            METHOD_OVERRIDES[override_key],
            member_name,
            known_types,
        )
    if owner_name is None and override_key in FUNCTION_OVERRIDES and FUNCTION_OVERRIDES[override_key]:
        return parse_rendered_signatures(
            FUNCTION_OVERRIDES[override_key],
            member_name,
            known_types,
        )

    doc = getattr(callable_value, "__doc__", "") or ""
    signatures = extract_doc_signatures(
        owner_name if constructor and owner_name is not None else member_name,
        doc,
        known_types,
        constructor=constructor,
        method=owner_name is not None,
        static=static,
    )
    if not signatures and owner_name is None and "_" in member_name:
        signatures = extract_doc_signatures(
            member_name.split("_", 1)[1],
            doc,
            known_types,
            constructor=constructor,
            method=False,
            static=static,
        )
    if signatures:
        return signatures

    try:
        source = inspect.getsource(callable_value)
    except (OSError, TypeError):
        source = ""

    if owner_name is not None:
        owner = getattr(module, owner_name)
        for base in owner.__mro__[1:]:
            base_name = getattr(base, "__name__", "")
            if base_name not in known_types or not hasattr(base, member_name):
                continue
            base_attr = getattr(base, member_name)
            try:
                return render_callable_signatures(
                    module,
                    base_name,
                    member_name,
                    base_attr,
                    known_types,
                    resolving,
                )
            except ValueError:
                continue

    if source:
        alias_signatures = signature_from_alias_target(
            module,
            owner_name,
            member_name,
            source,
            known_types,
            resolving,
        )
        if alias_signatures:
            return alias_signatures

    if owner_name is not None:
        if member_name == "repr":
            return [Signature((Parameter("self", None),), "str")]
        if member_name in UNARY_SELF_METHODS:
            return [Signature((Parameter("self", None),), "Self")]
        if member_name in {"__bool__"}:
            return [Signature((Parameter("self", None),), "bool")]
        if member_name in {"__iter__"}:
            return [Signature((Parameter("self", None),), "Iterator[Self]")]
        if member_name in {"__len__", "__int__", "__index__"}:
            return [Signature((Parameter("self", None),), "int")]
        if member_name in {"__float__"}:
            return [Signature((Parameter("self", None),), "float")]
        if member_name in {"__contains__"}:
            return [
                Signature(
                    (
                        Parameter("self", None),
                        Parameter("value", "bool | int | float | Self"),
                    ),
                    "bool",
                )
            ]
        if member_name in {"__array__"}:
            return [
                Signature(
                    (
                        Parameter("self", None),
                        Parameter("args", "str", vararg=True),
                        Parameter("kwargs", "str", kwarg=True),
                    ),
                    "NDArray[np.float64]",
                )
            ]
        if member_name in BINARY_SELF_METHODS:
            return [
                Signature(
                    (
                        Parameter("self", None),
                        Parameter("other", "bool | int | float | Self"),
                    ),
                    "bool" if member_name in {"__eq__", "__ge__", "__gt__", "__le__", "__lt__", "__ne__"} else "Self",
                )
            ]

    raise ValueError(
        f"no precise signature strategy for "
        f"{owner_name + '.' if owner_name else ''}{member_name}"
    )


def parse_rendered_signatures(
    rendered_lines: list[str], name: str, known_types: set[str]
) -> list[Signature]:
    signatures: list[Signature] = []
    current_decorators: list[str] = []
    for line in rendered_lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("@"):
            current_decorators.append(stripped)
            continue
        if not stripped.startswith("def "):
            continue
        signatures.append(parse_stub_signature(stripped, tuple(current_decorators), known_types))
        current_decorators = []
    if not signatures:
        raise ValueError(f"override for {name} does not contain a def")
    return signatures


def parse_stub_signature(
    line: str, decorators: tuple[str, ...], known_types: set[str]
) -> Signature:
    match = re.fullmatch(r"def ([A-Za-z_][A-Za-z0-9_]*)\((.*)\) -> (.*): \.\.\.", line)
    if match is None:
        raise ValueError(f"cannot parse stub signature: {line!r}")
    params_text = match.group(2)
    parameters: list[Parameter] = []
    for token in split_top_level(params_text):
        if token == "self":
            parameters.append(Parameter("self", None))
            continue
        if token.startswith("**"):
            name_text, type_text = token[2:].split(":", 1)
            parameters.append(Parameter(sanitize_name(name_text.strip()), type_text.strip(), kwarg=True))
            continue
        if token.startswith("*"):
            name_text, type_text = token[1:].split(":", 1)
            parameters.append(Parameter(sanitize_name(name_text.strip()), type_text.strip(), vararg=True))
            continue
        name_text, type_text = token.split(":", 1)
        parameters.append(
            Parameter(
                sanitize_name(name_text.strip()),
                type_text.replace(" = ...", "").strip(),
                default=token.endswith(" = ..."),
            )
        )
    return Signature(tuple(parameters), match.group(3).strip(), decorators)


def render_property_type(
    module: ModuleType, owner_name: str, member_name: str, known_types: set[str]
) -> str:
    override = PROPERTY_TYPE_OVERRIDES.get((owner_name, member_name))
    if override is not None:
        return override

    owner = getattr(module, owner_name)
    prop = inspect.getattr_static(owner, member_name)
    fget = getattr(prop, "fget", None)
    doc = getattr(fget, "__doc__", "") or ""
    getter_name = f"{owner_name}_{member_name}_get"
    signatures = extract_doc_signatures(
        getter_name,
        doc,
        known_types,
        constructor=False,
        method=True,
        static=False,
    )
    if signatures:
        return signatures[0].return_type

    value = getattr(module, owner_name)
    sample_factories = {
        "DM": lambda: getattr(module, "DM")(1),
        "MX": lambda: getattr(module, "MX").sym("x"),
        "Opti": lambda: getattr(module, "Opti")(),
        "SX": lambda: getattr(module, "SX").sym("x"),
        "Sparsity": lambda: getattr(module, "Sparsity").dense(1, 1),
        "NlpBuilder": lambda: getattr(module, "NlpBuilder")(),
    }
    factory = sample_factories.get(owner_name)
    if factory is not None:
        sample = factory()
        runtime_value = getattr(sample, member_name)
        return infer_runtime_type(runtime_value)

    try:
        sample = getattr(module, owner_name)()
    except Exception as exc:
        raise ValueError(f"no precise property type for {owner_name}.{member_name}") from exc
    return infer_runtime_type(getattr(sample, member_name))


def infer_runtime_type(value) -> str:
    class_name = value.__class__.__name__
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    if isinstance(value, list):
        if not value:
            raise ValueError("cannot infer element type from empty list")
        inner = infer_runtime_type(value[0])
        return f"list[{inner}]"
    if isinstance(value, tuple):
        return "tuple[" + ", ".join(infer_runtime_type(item) for item in value) + "]"
    if isinstance(value, ModuleType):
        return "ModuleType"
    if class_name in {"DM", "SX", "MX", "Sparsity", "OptiAdvanced", "NZproxy"}:
        return class_name
    raise ValueError(f"cannot infer runtime type for {value!r}")


def render_class(name: str, cls, module: ModuleType, known_types: set[str]) -> list[str]:
    override = CLASS_OVERRIDES.get(name)
    if override is not None:
        return override

    lines = [f"class {name}:"]
    members = exported_names(cls)
    if not members:
        lines.append("    ...")
        return lines

    for member in members:
        if inspect.isdatadescriptor(inspect.getattr_static(cls, member)):
            lines.append(f"    {member}: {render_property_type(module, name, member, known_types)}")
            continue

        attr = getattr(cls, member)
        if callable(attr):
            signatures = render_callable_signatures(module, name, member, attr, known_types)
            if len(signatures) > 1:
                signatures = [
                    Signature(
                        signature.parameters,
                        signature.return_type,
                        signature.decorators
                        if "@overload" in signature.decorators
                        else signature.decorators + ("@overload",),
                    )
                    for signature in signatures
                ]
            for signature in signatures:
                lines.extend(signature.render(member, indent="    "))
        else:
            lines.append(f"    {member}: {infer_runtime_type(attr)}")

    return lines


def render_support_protocols() -> list[str]:
    return [
        "class SupportsEvaluate:",
        "    def evaluate(self, inputs: Sequence[DM | SX | MX], outputs: list[DM | SX | MX]) -> None: ...",
        "",
    ]


def render_special_functions() -> list[str]:
    return [
        "@overload",
        "def PyFunction(",
        "    name: str,",
        "    obj: SupportsEvaluate,",
        "    inputs: Sequence[DM | SX | MX],",
        "    outputs: Sequence[DM | SX | MX],",
        "    opts: Mapping[str, GenericType] = ...,",
        ") -> Function: ...",
        "",
        "def pyfunction(",
        "    inputs: Sequence[DM | SX | MX],",
        "    outputs: Sequence[DM | SX | MX],",
        ") -> Callable[[Callable[[Sequence[DM | SX | MX]], DM | SX | MX | list[DM | SX | MX]]], Function]: ...",
        "",
    ]


def render_stub(module_name: str) -> str:
    module = importlib.import_module(module_name)
    known_types = {
        name
        for name in dir(module)
        if isinstance(getattr(module, name, None), type) and getattr(module, name, None) is not object
    }

    lines = HEADER.copy()

    for line in render_support_protocols():
        lines.append(line)

    for name in exported_names(module):
        obj = getattr(module, name)
        if isinstance(obj, ModuleType):
            continue
        if name == "object" and obj is object:
            lines.append("object: type[backup_object]")
            lines.append("")
            continue
        if isinstance(obj, type):
            lines.extend(render_class(name, obj, module, known_types))
        elif callable(obj):
            if name in {"PyFunction", "pyfunction"}:
                pass
            else:
                signatures = render_callable_signatures(module, None, name, obj, known_types)
                if len(signatures) > 1:
                    signatures = [
                        Signature(
                            signature.parameters,
                            signature.return_type,
                            signature.decorators
                            if "@overload" in signature.decorators
                            else signature.decorators + ("@overload",),
                        )
                        for signature in signatures
                    ]
                for signature in signatures:
                    lines.extend(signature.render(name))
        else:
            lines.append(f"{name}: {infer_runtime_type(obj)}")
        lines.append("")

    for line in render_special_functions():
        lines.append(line)

    return "\n".join(lines).rstrip() + "\n"


def render_access_check(module_name: str) -> str:
    module = importlib.import_module(module_name)
    lines = [
        "from __future__ import annotations",
        "",
        f"import {module_name} as ca",
        "",
        "def touch(value: object) -> None: ...",
        "",
    ]

    for name in exported_names(module):
        obj = getattr(module, name)
        if isinstance(obj, ModuleType):
            continue
        lines.append(f"touch(ca.{name})")
        if isinstance(obj, type):
            for member in exported_names(obj):
                lines.append(f"touch(ca.{name}.{member})")

    return "\n".join(lines).rstrip() + "\n"


def update_file(content: str, path: Path, check: bool) -> int:
    if check:
        current = path.read_text() if path.exists() else ""
        if current == content:
            return 0

        diff = difflib.unified_diff(
            current.splitlines(),
            content.splitlines(),
            fromfile=str(path),
            tofile=f"{path} (generated)",
            lineterm="",
        )
        print("\n".join(diff))
        return 1

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["stub", "access-check"])
    parser.add_argument("--module", default="casadi")
    parser.add_argument("--output", required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.mode == "stub":
        content = render_stub(args.module)
    else:
        content = render_access_check(args.module)

    return update_file(content, Path(args.output), args.check)


if __name__ == "__main__":
    raise SystemExit(main())
