---
myst:
  substitutions:
    condalogo: |-
      ```{image} _static/conda_logo.svg
      :class: no-scaled-link
      :width: 20%
      ```
---

(install)=

# Install Fastscapelib

Choose the way to install Fastscapelib that best suits your needs. In general we
recommend installing either the C++ or Python library using [conda] (or
[mamba]).

::::{grid} 2
:gutter: 3

:::{grid-item-card} C++ / {{ condalogo }}
:img-top: _static/cpp_logo.svg
:link: install-cpp-conda
:link-type: ref
:text-align: center

Install the C++ (header-only) library with conda_

(stable versions)
:::
:::{grid-item-card} C++ / from source
:img-top: _static/cpp_logo.svg
:link: install-cpp-source
:link-type: ref
:text-align: center

Download and install the C++ (header-only) library with CMake

(stable & development versions)
:::
:::{grid-item-card} Python / {{ condalogo }}
:img-top: _static/python_logo.svg
:link: install-python-conda
:link-type: ref
:text-align: center

Install the pre-compiled Python library with conda_

(stable versions)
:::
:::{grid-item-card} Python / from source
:img-top: _static/python_logo.svg
:link: install-python-source
:link-type: ref
:text-align: center

Download, build and install the Python library from source using Pip

(stable & development versions)
:::
::::

(download-fastscapelib)=

## Download the Fastscapelib Source

You can either download the last stable version as an archive
[here](https://github.com/fastscape-lem/fastscapelib/releases/latest) or you can
clone the source repository using git (by default it will clone the active
development `main` branch):

```bash
$ git clone https://github.com/fastscape-lem/fastscapelib
```

(install-cpp)=

## Install the C++ Library

Fastscapelib C++ is header only and requires a recent compiler that supports the
C++17 standards. It also depends on:

- [xtensor] (C++ array computing)
- [cmake] (build and configuration)

See Section {doc}`build_options` for more details on how to include Fastscapelib
in a (CMake) project.

(install-cpp-conda)=

### Using Conda

Fastscapelib [conda] packages are available for all stable versions via the
conda-forge channel. You can install Fastscapelib's C++ headers all dependencies
using the following command (alternatively you can use [mamba]):

```bash
$ conda install fastscapelib -c conda-forge
```

(install-cpp-source)=

### From Source Using CMake

In addition to a C++ compiler supporting C++17, you need [xtensor] and [cmake]
that you can install, e.g., using [conda]:

```bash
$ conda install xtensor cmake -c conda-forge
```

After {ref}`downloading the Fastscapelib source <download-fastscapelib>`, run
the commands below from the source root directory to install the fastscapelib's
header files using CMake:

```bash
$ cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
$ cmake --build build
$ cmake --install build
```

Where `/path/to/prefix` is the path where the header files will be installed
(skip this option if you want to install Fastscapelib in a default location).

See Section {doc}`build_options` for more information on the available build
options.

(install-python)=

## Install the Python Library

(install-python-conda)=

### Using Conda

Fastscapelib's Python bindings are available for all stable versions as binary
[conda] packages for Linux, MacOS and Windows via the conda-forge channel. You
can install it using the following command (alternatively you can use [mamba]):

```bash
$ conda install fastscapelib-python -c conda-forge
```

(install-python-source)=

### From Source Using Pip

Fastscapelib's Python bindings require Python (3.8+), pip, numpy, [pybind11],
[xtensor-python] and [scikit-build-core] which are all available on conda-forge:

```bash
$ conda install python pip numpy pybind11 xtensor-python scikit-build-core -c conda-forge
```

After {ref}`downloading the Fastscapelib source <download-fastscapelib>`, you
can build and install the Python package using `pip`. Run the following commands
from the source root directory:

```bash
$ python -m pip install . --no-build-isolation
```

The ``--no-build-isolation`` option is required since xtensor-python cannot yet
be installed using pip. You can pass extra options like in the example below,
which builds the Python extension in ``Debug`` mode in a given directory (useful
for cached builds):

```bash
$ python -m pip install . \
$   --no-build-isolation \
$   --config-settings=cmake.build-type=Debug \
$   --config-settings=build-dir=build/skbuild
```

See [scikit-build-core]'s documentation for more available options.

[cmake]: https://cmake.org/
[conda]: https://conda.io/docs/
[mamba]: https://mamba.readthedocs.io/en/latest/
[pybind11]: https://github.com/pybind/pybind11
[xtensor]: https://xtensor.readthedocs.io
[xtensor-python]: https://xtensor-python.readthedocs.io
[scikit-build-core]: https://scikit-build-core.readthedocs.io