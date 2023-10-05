from typing import ClassVar, List, Union, overload

import numpy as np
import numpy.typing as npt

from fastscapelib.grid import ProfileGrid, RasterGrid, TriMesh

Grid = Union[ProfileGrid, RasterGrid, TriMesh]

class FlowGraphImpl:
    @property
    def single_flow(self) -> bool: ...
    @property
    def receivers(self) -> npt.NDArray[np.uint64]: ...
    @property
    def receivers_count(self) -> npt.NDArray[np.uint64]: ...
    @property
    def receivers_distance(self) -> npt.NDArray[np.float64]: ...
    @property
    def receivers_weight(self) -> npt.NDArray[np.float64]: ...
    @property
    def donors(self) -> npt.NDArray[np.uint64]: ...
    @property
    def donors_count(self) -> npt.NDArray[np.uint64]: ...
    @property
    def dfs_indices(self) -> npt.NDArray[np.uint64]: ...
    @property
    def basins(self) -> npt.NDArray[np.uint64]: ...

class FlowDirection:
    __members__: ClassVar[dict] = ...  # read-only
    UNDEFINED: ClassVar[FlowDirection] = ...
    SINGLE: ClassVar[FlowDirection] = ...
    MULTI: ClassVar[FlowDirection] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class FlowOperator:
    @property
    def name(self) -> str: ...

class SingleFlowRouter(FlowOperator):
    def __init__(self) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    def __repr__(self) -> str: ...

class MultiFlowRouter(FlowOperator):
    def __init__(self, slope_exp: float = 1.0) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    slope_exp: float
    def __repr__(self) -> str: ...

class PFloodSinkResolver(FlowOperator):
    def __init__(self) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    def __repr__(self) -> str: ...

class MSTMethod:
    __members__: ClassVar[dict] = ...  # read-only
    KRUSKAL: ClassVar[MSTMethod] = ...
    BORUVKA: ClassVar[MSTMethod] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class MSTRouteMethod:
    __members__: ClassVar[dict] = ...  # read-only
    BASIC: ClassVar[MSTRouteMethod] = ...
    CARVE: ClassVar[MSTRouteMethod] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class MSTSinkResolver(FlowOperator):
    def __init__(
        self,
        basin_method: MSTMethod = MSTMethod.KRUSKAL,
        route_method: MSTRouteMethod = MSTRouteMethod.CARVE,
    ) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    basin_method: MSTMethod
    route_method: MSTRouteMethod
    def __repr__(self) -> str: ...

class FlowSnapshot(FlowOperator):
    def __init__(
        self, snapshot_name: str, save_graph: bool = True, save_elevation: bool = False
    ) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    @property
    def snapshot_name(self) -> str: ...
    @property
    def save_graph(self) -> bool: ...
    @property
    def save_elevation(self) -> bool: ...
    def __repr__(self) -> str: ...

class FlowGraph:
    def __init__(self, grid: Grid, operators: List[FlowOperator]) -> None: ...
    @property
    def operators(self) -> List[FlowOperator]: ...
    @property
    def single_flow(self) -> bool: ...
    def impl(self) -> FlowGraphImpl: ...
    @property
    def graph_snapshot_keys(self) -> List[str]: ...
    def graph_snapshot(self, name: str) -> FlowGraph: ...
    @property
    def elevation_snapshot_keys(self) -> List[str]: ...
    def elevation_snapshot(self, name: str) -> npt.NDArray[np.float64]: ...
    def update_routes(
        self, elevation: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]: ...
    @property
    def base_levels(self) -> list[int]: ...
    @base_levels.setter
    def base_levels(self, value: list[int]) -> None: ...
    @property
    def mask(self) -> npt.NDArray[np.bool_]: ...
    @mask.setter
    def mask(self, value: npt.NDArray[np.bool_]) -> None: ...
    @overload
    def accumulate(
        self, acc: npt.NDArray[np.float64], src: npt.NDArray[np.float64]
    ) -> None: ...
    @overload
    def accumulate(self, acc: npt.NDArray[np.float64], src: float) -> None: ...
    @overload
    def accumulate(self, src: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]: ...
    @overload
    def accumulate(self, src: float) -> npt.NDArray[np.float64]: ...
    def basins(self) -> npt.NDArray[np.uint64]: ...
    def __repr__(self) -> str: ...
