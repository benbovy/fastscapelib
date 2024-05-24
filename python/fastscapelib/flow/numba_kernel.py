import time
from contextlib import contextmanager
from textwrap import dedent, indent

import numba as nb
import numpy as np

from fastscapelib.flow import Kernel, KernelApplicationOrder, KernelData


@contextmanager
def timer(msg: str, do_print: bool):
    """Times and prints a code block."""

    start_time = time.time()
    yield
    end_time = time.time()
    elapsed_time = end_time - start_time
    if do_print:
        print(f"Time spent in {msg}: {elapsed_time:.1f} seconds")


class NumbaFlowKernel:
    def __init__(
        self,
        flow_graph,
        kernel_func,
        spec,
        application_order,
        outputs=(),
        max_receivers: int = -1,
        n_threads: int = 1,
        print_generated_code: bool = False,
        print_stats: bool = False,
    ):
        with timer("flow kernel init", print_stats):
            self._flow_graph = flow_graph
            self._py_flow_kernel = kernel_func
            self._outputs = outputs
            self._max_receivers = max_receivers
            self._spec = spec
            self._print_generated_code = print_generated_code
            self._print_stats = print_stats
            self._n_threads = n_threads
            self._application_order = application_order

            self._build_fs_kernel()

    def _build_fs_kernel(self):
        """Builds the fastscapelib flow kernel.

        The fastscapelib flow kernel is a collection of data and
        function pointers to be used by the `apply_kernel` runtime
        to call the flow kernel on all the grid.

        The flow kernel must be thread-safe. It can be called sequentially
        or in parallel depending on the caller implementation.
        """

        self.kernel = kernel = Kernel()
        self.kernel_data = KernelData()

        with timer("build data classes", self._print_stats):
            self._build_and_set_data_classes()
        with timer("build node data create/free", self._print_stats):
            self._build_and_set_node_data_create_free()
        with timer("build node data getter/setter", self._print_stats):
            self._build_and_set_node_data_getter()
            self._build_and_set_node_data_setter()
        with timer("build flow kernel", self._print_stats):
            self._build_and_set_flow_kernel_ptr()
        self.kernel.n_threads = self._n_threads
        self.kernel.application_order = self._application_order

    def _build_and_set_flow_kernel_ptr(self):
        """Builds and sets the flow kernel jitted function.

        The flow kernel function is called by the computation thread with
        a node data and the integration time step.
        """
        self._jitted_flow_kernel = jitted_func = nb.njit(inline="always")(
            self._py_flow_kernel
        )

        compiled_func = jitted_func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (self._node_data_jitclass.class_type.instance_type,),
                None,
            )
        )
        if self._print_stats:
            print("flow kernel compilation stats", compiled_func.metadata["timers"])

        self.kernel.func = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _build_and_set_node_data_create_free(self):
        """Builds the node data create and free functions."""

        node_data_cls = self._node_data_jitclass

        scalar_data = []

        @nb.njit
        def node_data_create():
            return node_data_cls()

        self._set_node_data_create(self._node_data_jitclass, node_data_create)
        self._build_and_set_node_data_init()
        self._set_node_data_free()

    def _build_and_set_node_data_init(self):
        """Sets the pointer to the node data init function.

        The node data init function is called after instantiation/creation
        to set scalar variables. Scalar variables are shared by all the nodes
        and as such are considered immutable. They can be set once before 
        applying the kernel on each node.

        This function pointer may be set to a null pointer if there is no
        scalar variable to be set:
          - in case there is no scalar variable at all
          - if the scalar variable value is given in the kernel specification
            (in which case the constant scalar value is intialized inline in the 
            node data constructor)
        """

        func_tmpl = dedent(
            """
        def node_data_init(node_data, data):
        {content}
        """
        )

        content = "\n".join(
            [
                f"node_data.{name} = data.{name}"
                for name, ty in self._constants.items()
                if issubclass(ty.__class__, nb.core.types.Type)
            ]
        )
        if content == "":
            self.node_data_init = None
            return

        init_source = func_tmpl.format(content=indent(content, self._indent4))

        if self._print_generated_code:
            print(f"Node data init source code:\n{indent(init_source, self._indent4)}")

        glbls = {}
        exec(init_source, glbls)
        func = glbls["node_data_init"]

        self.node_data_init = func = nb.njit(inline="always", boundscheck=False)(func)
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (
                    self._node_data_jitclass.class_type.instance_type,
                    self._data_jitclass.class_type.instance_type,
                ),
                None,
            )
        )
        if self._print_stats:
            print(
                "Node data init compilation stats", compiled_func.metadata["timers"]
            )

        self.kernel.node_data_init = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _build_and_set_node_data_getter(self):
        """Builds the jitted node data getter function.

        The node data getter is called prior to the flow kernel to
        copy the global data for one specific node and its receivers
        to a node data instance.
        """

        self.node_data_getter = func = self._build_node_data_getter(self._max_receivers)
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (
                    nb.uint64,
                    self._data_jitclass.class_type.instance_type,
                    self._node_data_jitclass.class_type.instance_type,
                ),
                None,
            )
        )
        if self._print_stats:
            print(
                "Node data getter compilation stats", compiled_func.metadata["timers"]
            )

        self.kernel.node_data_getter = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _build_and_set_node_data_setter(self):
        """Builds the jitted node data setter function.

        The node data setter is called after the flow kernel to
        copy back the node data in the global data instance.
        """

        self.node_data_setter = func = self._build_node_data_setter()
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (
                    nb.uint64,
                    self._node_data_jitclass.class_type.instance_type,
                    self._data_jitclass.class_type.instance_type,
                ),
                None,
            )
        )
        if self._print_stats:
            print(
                "Node data setter compilation stats", compiled_func.metadata["timers"]
            )

        self.kernel.node_data_setter = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _set_node_data_create(self, cls, func):
        """Sets the pointer to the node data create function.

        The node data create function is called by a compute thread
        to create a new node data then used on each node of the grid.
        """

        self.node_data_create = func
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                cls.class_type.instance_type,
                (),
                None,
            )
        )

        self.kernel.node_data_create = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _set_node_data_free(self):
        """Sets the pointer to the node data free function.

        The node data free function is called by a compute thread
        to delete a node data (it calls the jitclass destructor and
        deallocates related NRT meminfo and data pointers).
        """

        self.kernel.node_data_free = (
            nb.core.runtime.nrt.rtsys.library.get_pointer_to_function("NRT_decref")
        )

    def _build_and_set_data_classes(self):
        """Builds data and node data jitclasses."""

        flow_graph = self._flow_graph

        # FIXME
        self._grid_data = {name: ty for name, ty in self._spec if ty == nb.float64[::1]}
        self._constants = {name: ty for name, ty in self._spec if ty != nb.float64[::1]}

        # for name, value in self._grid_data.items():
        #     if value.size != flow_graph.size:
        #         raise ValueError("Invalid size")
        #     if value.shape != (flow_graph.size,):
        #         raise ValueError("Invalid shape")

        self._build_node_data_jitclass()
        self._build_data_jitclass()

        flow_graph_data = {
            "donors_idx": flow_graph.impl().donors.view(),
            "donors_count": flow_graph.impl().donors_count.view(),
            "receivers_idx": flow_graph.impl().receivers.view(),
            "receivers_count": flow_graph.impl().receivers_count.view(),
            "receivers_distance": flow_graph.impl().receivers_distance.view(),
            "receivers_weight": flow_graph.impl().receivers_weight.view(),
        }

        from numba.experimental.jitclass import _box

        self._data = self._data_jitclass(**flow_graph_data)
        self.kernel_data.data.meminfo = _box.box_get_meminfoptr(self._data)
        self.kernel_data.data.data = _box.box_get_dataptr(self._data)

    @staticmethod
    def _generate_jitclass(name, spec, init_source, glbls={}):
        exec(init_source, glbls)
        ctor = glbls["generated_init"]
        return nb.experimental.jitclass(spec)(type(name, (), {"__init__": ctor}))

    def _build_node_data_jitclass(self):
        """Builds a node data jitclass.

        A node data instance contains all the data related to a node and its receivers.
        It is used to call the flow kernel function with this specific node data.

        The grid data specified by the user are aggregated with the constants and
        also with the same data at the receivers. The count, distance and weight to
        the current node's receivers are also provided.

        We need to workaround some `numba` limitations:
        - `setattr` is not implemented -> use a template source code to be executed
        """

        grid_data = self._grid_data
        constants = self._constants
        max_receivers = self._max_receivers

        receivers_flow_graph_spec = [
            ("distance", nb.float64[::1]),
            ("weight", nb.float64[::1]),
        ]
        receivers_grid_data_spec = [(name, ty) for name, ty in grid_data.items()]
        receivers_internal_spec = [
            ("count", nb.uint64),
        ]

        receivers_spec = (
            receivers_flow_graph_spec
            + receivers_grid_data_spec
            + receivers_internal_spec
            + [
                ("_" + name, ty)
                for name, ty in receivers_flow_graph_spec + receivers_grid_data_spec
            ]
        )

        @nb.experimental.jitclass(receivers_spec)
        class ReceiversData(object):
            def __init__(self):
                pass

        def get_type(value):
            if issubclass(value.__class__, nb.core.types.Type):
                return value
            return nb.typeof(value)

        base_spec = [("receivers", ReceiversData.class_type.instance_type)]
        grid_data_spec = [(name, ty.dtype) for name, ty in grid_data.items()]
        constants_spec = [(name, get_type(ty)) for name, ty in constants.items()]

        __init___template = dedent(
            """
        def generated_init(self):
            {constants_content}

            self.receivers = ReceiversData()
            self.receivers._distance = np.ones({default_size})
            self.receivers._weight = np.ones({default_size})
            {receivers_content_init}

            self.receivers.distance = self.receivers._distance[:]
            self.receivers.weight = self.receivers._weight[:]
            {receivers_content_view}

            self.receivers.count = 0
        """
        )

        default_size = max_receivers if max_receivers > 0 else 0
        constants_content = "\n    ".join(
            [
                f"self.{name} = {ty}"
                for name, ty in constants.items()
                if not issubclass(ty.__class__, nb.core.types.Type)
            ]
        )
        receivers_content_init = "\n    ".join(
            [
                f"self.receivers._{name} = np.ones({default_size}, dtype=np.{value.dtype})"
                for name, value in grid_data.items()
            ]
        )
        receivers_content_view = "\n    ".join(
            [f"self.receivers.{name} = self.receivers._{name}[:]" for name in grid_data]
        )
        init_source = __init___template.format(
            constants_content=constants_content,
            receivers_content_init=receivers_content_init,
            receivers_content_view=receivers_content_view,
            default_size=default_size,
        )

        spec = base_spec + grid_data_spec + constants_spec

        if self._print_generated_code:
            print(f"Node data jitclass constructor source code:\n{indent(init_source, self._indent4)}")

        self._node_data_jitclass = NumbaFlowKernel._generate_jitclass(
            "FlowKernelNodeData",
            base_spec + grid_data_spec + constants_spec,
            init_source,
            {"ReceiversData": ReceiversData, "np": np},
        )

    def _build_data_jitclass(self):
        """Builds a data jitclass.

        A data instance contains all the data maped on the grid. It contains the
        initial and final data after applying the kernel.

        The grid data specified by the user are aggregated with the flow graph ones
        to give access to the receivers (indices, count, distance, weight).
        The flow graph donors are also exposed but currently not used.
        """

        grid_data = self._grid_data
        scalars = self._constants

        base_spec = [
            ("donors_idx", nb.uint64[:, ::1]),
            ("donors_count", nb.uint64[::1]),
            ("receivers_idx", nb.uint64[:, ::1]),
            ("receivers_count", nb.uint64[::1]),
            ("receivers_distance", nb.float64[:, ::1]),
            ("receivers_weight", nb.float64[:, ::1]),
        ]
        grid_data_spec = [(name, ty) for name, ty in grid_data.items()]
        scalars_spec = [
            (name, ty)
            for name, ty in scalars.items()
            if issubclass(ty.__class__, nb.core.types.Type)
        ]
        spec = base_spec + grid_data_spec + scalars_spec

        __init___template = dedent(
            """
        def generated_init(self, {args}):
            {content}
        """
        )

        content = "\n    ".join([f"self.{name} = {name}" for name, _ in base_spec])
        args = ", ".join([name for name, _ in base_spec])
        init_source = __init___template.format(content=content, args=args)

        if self._print_generated_code:
            print(f"Data jitclass constructor source code:\n{indent(init_source, self._indent4)}")

        self._data_jitclass = NumbaFlowKernel._generate_jitclass(
            "FlowKernelData",
            spec,
            init_source,
        )

    _node_data_getter_tmpl = dedent(
        """
        def node_data_getter(index, data, node_data):
            receivers_count = data.receivers_count[index]
            receivers = node_data.receivers

        {node_content}
            {resize_content}

            receivers.count = receivers_count

            for i in range(receivers_count):
                receiver_idx = data.receivers_idx[index, i]
                receivers._distance[i] = data.receivers_distance[index, i]
                receivers._weight[i] = data.receivers_weight[index, i]
        {receivers_set_content}

            return 0
        """
    )

    _node_data_getter_fixed_resize_tmpl = dedent(
        """
        if {max_receivers} < receivers_count:
            return 1

        if receivers_count != receivers.count:
        {set_views}
        """
    ).rstrip("\n")

    _node_data_getter_dynamic_resize_tmpl = dedent(
        """
        if receivers_count != receivers.count:
            if receivers_count > receivers.count:
        {receivers_resize_source}
        {set_views}
        """
    ).rstrip("\n")

    _set_view_tmpl = dedent(
        """
        set_view(
            (
        {view_data},
            ),
            receivers_count
        )"""
    ).lstrip("\n")

    _indent4 = " " * 4
    _indent8 = _indent4 * 2
    _indent12 = _indent4 * 3

    def _build_node_data_getter(self, max_receivers):
        """Builds a node data getter from the global data

        The node data getter is called prior to the flow kernel to
        copy the global data for one specific node and its receivers
        to a node data instance.

        We need to workaround some `numba` limitations:
        - `setattr` is not implemented -> use a template source code to be executed
        - passing a tuple of various array types doesn't work -> use `set_view` with
          pack of consistent data

        Note: instead of using `set_view` multiple times, an inlining of the views has
        been tested but generates very poor performances (not understood in details)
        """

        data_dtypes = {}
        for name, value in self._grid_data.items():
            if value.dtype in data_dtypes:
                data_dtypes[value.dtype].append(name)
            else:
                data_dtypes[value.dtype] = [name]

        node_content = "\n".join(
            [f"node_data.{name} = data.{name}[index]" for name in self._grid_data]
        )

        receivers_view_data = [
            f",\n".join([f"(receivers.{name}, receivers._{name})" for name in names])
            for names in data_dtypes.values()
        ]

        receivers_resize_source = "\n".join(
            [
                f"receivers._{name} = np.empty(receivers_count, dtype=np.{value.dtype})"
                for name, value in self._grid_data.items()
            ]
        )
        receivers_set_content = "\n".join(
            [
                f"receivers._{name}[i] = data.{name}[receiver_idx]"
                for name in self._grid_data
            ]
        )

        if max_receivers > 0:
            resize_tmpl = self._node_data_getter_fixed_resize_tmpl
        else:
            resize_tmpl = self._node_data_getter_dynamic_resize_tmpl

        resize_source = indent(
            resize_tmpl.format(
                set_views=indent(
                    "\n".join(
                        self._set_view_tmpl.format(
                            view_data=indent(content, self._indent8)
                        )
                        for content in receivers_view_data
                    ),
                    self._indent4,
                ),
                receivers_resize_source=indent(receivers_resize_source, self._indent8),
                max_receivers=max_receivers,
            ),
            self._indent4,
        )

        getter_source = self._node_data_getter_tmpl.format(
            node_content=indent(node_content, self._indent4),
            resize_content=resize_source,
            receivers_set_content=indent(receivers_set_content, self._indent8),
        )
        if self._print_generated_code:
            print(
                f"Node data getter source code:\n{indent(getter_source, self._indent4)}"
            )

        @nb.njit(inline="always")
        def set_view(data, size):
            for view, source in data:
                view = source[:size]

        glbls = {"np": np, "set_view": set_view}
        exec(getter_source, glbls)
        getter_fun = glbls["node_data_getter"]

        return nb.njit(inline="always", boundscheck=False)(getter_fun)

    node_data_setter_tmpl = dedent(
        """
        def node_data_setter(index, node_data, data):
        {content}

            return 0
        """
    )

    def _build_node_data_setter(self):
        """Builds a node data setter from the global data

        The node data setter is called after the flow kernel to
        copy back the node data in the global data instance.

        We need to workaround some `numba` limitations:
        - `setattr` is not implemented -> use a template source code to be executed
        """

        invalid_outputs = [
            name for name in self._outputs if name not in self._grid_data
        ]
        if invalid_outputs:
            raise KeyError(
                f"Output name{'s' if len(invalid_outputs)>1 else ''} {invalid_outputs} not defined as 'grid_data'"
            )

        content = "\n".join(
            [
                f"data.{name}[index] = node_data.{name}"
                for name in self._grid_data
                if name in self._outputs
            ]
        )

        setter_source = self.node_data_setter_tmpl.format(
            content=indent(content, self._indent4)
        )
        if self._print_generated_code:
            print(f"Node data setter source code:\n{indent(setter_source, self._indent4)}")

        glbls = {}
        exec(setter_source, glbls)
        setter_fun = glbls["node_data_setter"]
        return nb.njit(inline="always", boundscheck=False)(setter_fun)

    def bind_data(self, **kwargs):
        for name, value in kwargs.items():
            setattr(self._data, name, value)


@nb.njit
def py_apply_kernel_impl(
    indices,
    func,
    data,
    node_data,
    node_data_getter,
    node_data_setter,
):
    """Applies a kernel on a grid.

    This Python implementation allows numba to inline
    the calls to the node data getter, the flow kernel,
    and node data setter.

    In case of sequential execution, it may give better
    performances than using the C++ fastscapelib `apply_kernel`
    implementation.
    """

    for i in indices:
        node_data_getter(i, data, node_data)
        func(node_data)
        node_data_setter(i, node_data, data)


def py_apply_kernel(nb_kernel, data):
    """Applies a kernel on a grid.

    This wrapper function calls the Python jitted implementation
    of a sequential call of the flow kernel on the grid nodes.
    """

    kernel = nb_kernel.kernel
    node_data = nb_kernel.node_data_create()
    if nb_kernel.node_data_init:
        nb_kernel.node_data_init(node_data, nb_kernel._data)

    if kernel.application_order == KernelApplicationOrder.ANY:
        indices = np.arange(0, nb_kernel._flow_graph.size, 1)
    elif kernel.application_order == KernelApplicationOrder.BREADTH_UPSTREAM:
        indices = nb_kernel._flow_graph.impl().bfs_indices
    else:
        raise RuntimeError("Unsupported kernel application order")

    py_apply_kernel_impl(
        indices,
        nb_kernel._jitted_flow_kernel,
        nb_kernel._data,
        node_data,
        nb_kernel.node_data_getter,
        nb_kernel.node_data_setter,
    )

    return 0
