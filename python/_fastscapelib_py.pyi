from typing import ClassVar, List, Tuple, Type, overload

import numpy as np
import numpy.typing as npt

#
# Grid types
#

class NodeStatus:
    __members__: ClassVar[dict] = ...  # read-only
    CORE: ClassVar[NodeStatus] = ...
    FIXED_VALUE_BOUNDARY: ClassVar[NodeStatus] = ...
    FIXED_GRADIENT_BOUNDARY: ClassVar[NodeStatus] = ...
    LOOPED_BOUNDARY: ClassVar[NodeStatus] = ...
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

class Node:
    def __init__(self, idx: int, status: NodeStatus) -> None: ...
    idx: int
    status: NodeStatus

class Neighbor:
    def __init__(self, idx: int, distance: float, status: NodeStatus) -> None: ...
    def __eq__(self, other: Neighbor) -> bool: ...
    idx: int
    distance: float
    status: NodeStatus

class ProfileBoundaryStatus:
    @overload
    def __init__(self, status: NodeStatus) -> None: ...
    @overload
    def __init__(self, left_status: NodeStatus, right_status: NodeStatus) -> None: ...
    @overload
    def __init__(self, status: List[NodeStatus]) -> None: ...

class ProfileGrid:
    @overload
    def __init__(
        self,
        size: int,
        spacing: float,
        status_at_bounds: ProfileBoundaryStatus,
        status_at_nodes: List[Node],
    ) -> None: ...
    @overload
    def __init__(
        self,
        size: int,
        spacing: float,
        status_at_bounds: List[NodeStatus],
        status_at_nodes: List[Tuple[int, NodeStatus]],
    ) -> None: ...
    @classmethod
    def from_length(
        cls: Type[ProfileGrid],
        size: int,
        length: float,
        status_at_bounds: ProfileBoundaryStatus,
        status_at_nodes: List[Node],
    ) -> ProfileGrid: ...
    is_structured: bool
    is_uniform: bool
    n_neighbors_max: int
    @property
    def size(self) -> int: ...
    @property
    def shape(self) -> List[int]: ...
    @property
    def status_at_nodes(self) -> npt.NDArray[np.uint8]: ...
    def neighbors_count(self, idx: int) -> int: ...
    def neighbors_indices(self, idx: int) -> npt.NDArray[np.uint64]: ...
    def neighbors_distances(self, idx: int) -> npt.NDArray[np.float64]: ...
    def neighbors(self, idx: int) -> List[Neighbor]: ...

class UnstructuredMesh:
    def __init__(
        self,
        points: npt.NDArray[np.float64],
        neighbors_indices_ptr: npt.NDArray[np.uint64],
        neighbors_indices: npt.NDArray[np.uint64],
        convex_hull_indices: npt.NDArray[np.uint64],
        areas: npt.NDArray[np.float64],
        status_at_nodes: List[Node],
    ) -> None: ...
    is_structured: bool
    is_uniform: bool
    n_neighbors_max: int
    @property
    def size(self) -> int: ...
    @property
    def shape(self) -> List[int]: ...
    @property
    def status_at_nodes(self) -> npt.NDArray[np.uint8]: ...
    def neighbors_count(self, idx: int) -> int: ...
    def neighbors_indices(self, idx: int) -> npt.NDArray[np.uint64]: ...
    def neighbors_distances(self, idx: int) -> npt.NDArray[np.float64]: ...
    def neighbors(self, idx: int) -> List[Neighbor]: ...
