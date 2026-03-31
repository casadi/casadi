from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from typing import Never, Self, TypeAlias, overload

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csc_matrix

from ._api_surface import *
from ._api_surface import Callback as _SurfaceCallback
from ._api_surface import DM as _SurfaceDM
from ._api_surface import Function as _SurfaceFunction
from ._api_surface import Linsol as _SurfaceLinsol
from ._api_surface import MX as _SurfaceMX
from ._api_surface import Opti as _SurfaceOpti
from ._api_surface import OptiSol as _SurfaceOptiSol
from ._api_surface import Sparsity as _SurfaceSparsity
from ._api_surface import SX as _SurfaceSX

__version__: str
__all__: list[str]

inf: float
nan: float
pi: float

OP_CONST: int
OP_INPUT: int
OP_OUTPUT: int
OP_ADD: int
OP_MUL: int
OP_MTIMES: int
OP_TWICE: int

_ArrayKey: TypeAlias = int | slice | tuple[int | slice, int | slice]
_MatrixScalar: TypeAlias = bool | int | float
_MatrixData: TypeAlias = (
    _MatrixScalar
    | Sequence[_MatrixScalar]
    | Sequence[Sequence[_MatrixScalar]]
    | NDArray[np.bool_]
    | NDArray[np.int_]
    | NDArray[np.float64]
)
_FunctionOutput: TypeAlias = "DM" | "SX" | "MX"
_FunctionCallResult: TypeAlias = _FunctionOutput


class _CasadiObject:
    ...


class Sparsity(_SurfaceSparsity, _CasadiObject):
    shape: tuple[int, int]
    T: Sparsity
    def __getitem__(self, key: _ArrayKey) -> Sparsity: ...
    def __add__(self, other: bool | int | float | Sparsity) -> Sparsity: ...
    def __radd__(self, other: bool | int | float | Sparsity) -> Sparsity: ...


class _CasadiMatrix(_CasadiObject):
    shape: tuple[int, int]
    T: Self
    nz: NZproxy
    def size1(self) -> int: ...
    def size2(self) -> int: ...
    def sparsity(self) -> Sparsity: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Never]: ...
    def __getitem__(self, key: _ArrayKey) -> Self: ...
    def __setitem__(self, key: _ArrayKey, value: _MatrixData | Self) -> None: ...
    def __add__(self, other: bool | int | float | Self) -> Self: ...
    def __radd__(self, other: bool | int | float | Self) -> Self: ...
    def __sub__(self, other: bool | int | float | Self) -> Self: ...
    def __rsub__(self, other: bool | int | float | Self) -> Self: ...
    def __mul__(self, other: bool | int | float | Self) -> Self: ...
    def __rmul__(self, other: bool | int | float | Self) -> Self: ...
    def __truediv__(self, other: bool | int | float | Self) -> Self: ...
    def __rtruediv__(self, other: bool | int | float | Self) -> Self: ...
    def __pow__(self, other: bool | int | float | Self) -> Self: ...
    def __rpow__(self, other: bool | int | float | Self) -> Self: ...
    def __matmul__(self, other: bool | int | float | Self) -> Self: ...
    def __rmatmul__(self, other: bool | int | float | Self) -> Self: ...
    def __neg__(self) -> Self: ...
    def __lt__(self, other: bool | int | float | Self) -> Self: ...
    def __le__(self, other: bool | int | float | Self) -> Self: ...
    def __gt__(self, other: bool | int | float | Self) -> Self: ...
    def __ge__(self, other: bool | int | float | Self) -> Self: ...


class DM(_SurfaceDM, _CasadiMatrix):
    shape: tuple[int, int]
    T: DM
    nz: NZproxy
    @staticmethod
    def eye(n: int) -> DM: ...
    def full(self) -> NDArray[np.float64]: ...
    def sparse(self) -> csc_matrix: ...


class SX(_SurfaceSX, _CasadiMatrix):
    shape: tuple[int, int]
    T: SX
    nz: NZproxy
    @overload
    @staticmethod
    def sym(name: str) -> SX: ...
    @overload
    @staticmethod
    def sym(name: str, nrow: int, ncol: int = ...) -> SX: ...
    @overload
    @staticmethod
    def sym(name: str, rc: tuple[int, int]) -> SX: ...
    @overload
    @staticmethod
    def sym(name: str, sp: Sparsity) -> SX: ...


class MX(_SurfaceMX, _CasadiMatrix):
    shape: tuple[int, int]
    T: MX
    nz: NZproxy
    @overload
    @staticmethod
    def sym(name: str) -> MX: ...
    @overload
    @staticmethod
    def sym(name: str, nrow: int, ncol: int = ...) -> MX: ...
    @overload
    @staticmethod
    def sym(name: str, rc: tuple[int, int]) -> MX: ...
    @overload
    @staticmethod
    def sym(name: str, sp: Sparsity) -> MX: ...


class Function(_SurfaceFunction, _CasadiObject):
    @overload
    def __call__(self, *args: DM | SX | MX) -> _FunctionOutput: ...
    @overload
    def __call__(self, **kwargs: DM | SX | MX) -> Mapping[str, _FunctionOutput]: ...


class OptiSol(_SurfaceOptiSol, _CasadiObject):
    debug: OptiAdvanced
    opti: Opti
    def value(self, x: DM | SX | MX, values: Sequence[MX] = ...) -> float: ...


class Opti(_SurfaceOpti, _CasadiObject):
    advanced: OptiAdvanced
    debug: OptiAdvanced
    f: MX
    f_linear_scale: float
    g: MX
    g_linear_scale: DM
    lam_g: MX
    lbg: MX
    ng: int
    np: int
    nx: int
    p: MX
    ubg: MX
    x: MX
    x_linear_scale: DM
    x_linear_scale_offset: DM
    def __init__(self) -> None: ...
    @overload
    def variable(self) -> MX: ...
    @overload
    def variable(self, n: int, m: int = ..., attribute: str = ...) -> MX: ...
    @overload
    def variable(self, sp: Sparsity, attribute: str = ...) -> MX: ...
    @overload
    def variable(self, symbol: MX, attribute: str = ...) -> MX: ...
    def solver(
        self,
        solver: str,
        plugin_options: Mapping[str, GenericType] = ...,
        solver_options: Mapping[str, GenericType] = ...,
    ) -> None: ...
    def solve(self) -> OptiSol: ...
    def value(self, x: DM | SX | MX, values: Sequence[MX] = ...) -> float: ...


class Callback(_SurfaceCallback, _CasadiObject):
    @overload
    def __call__(self, *args: DM | SX | MX) -> _FunctionOutput: ...
    @overload
    def __call__(self, **kwargs: DM | SX | MX) -> Mapping[str, _FunctionOutput]: ...


class Linsol(_SurfaceLinsol, Function):
    ...


@overload
def sin(x: float) -> float: ...
@overload
def sin(x: DM) -> DM: ...
@overload
def sin(x: SX) -> SX: ...
@overload
def sin(x: MX) -> MX: ...


@overload
def dot(x: DM, y: DM) -> DM: ...
@overload
def dot(x: SX, y: SX) -> SX: ...
@overload
def dot(x: MX, y: MX) -> MX: ...


@overload
def gradient(ex: DM, arg: DM, opts: Mapping[str, GenericType] = ...) -> DM: ...
@overload
def gradient(ex: SX, arg: SX, opts: Mapping[str, GenericType] = ...) -> SX: ...
@overload
def gradient(ex: MX, arg: MX, opts: Mapping[str, GenericType] = ...) -> MX: ...


@overload
def mtimes(args: Sequence[DM]) -> DM: ...
@overload
def mtimes(args: Sequence[SX]) -> SX: ...
@overload
def mtimes(args: Sequence[MX]) -> MX: ...
@overload
def mtimes(x: DM, y: DM) -> DM: ...
@overload
def mtimes(x: SX, y: SX) -> SX: ...
@overload
def mtimes(x: MX, y: MX) -> MX: ...


@overload
def vertcat(*args: _MatrixScalar | DM) -> DM: ...
@overload
def vertcat(*args: _MatrixScalar | SX) -> SX: ...
@overload
def vertcat(*args: _MatrixScalar | MX) -> MX: ...
@overload
def vertcat(*args: _MatrixScalar | Sparsity) -> Sparsity: ...


@overload
def nlpsol(
    name: str,
    solver: str,
    nlp: Mapping[str, SX],
    opts: Mapping[str, GenericType] = ...,
) -> Function: ...
@overload
def nlpsol(
    name: str,
    solver: str,
    nlp: Mapping[str, MX],
    opts: Mapping[str, GenericType] = ...,
) -> Function: ...
