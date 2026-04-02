# Installation

`smolpack` can be installed from PyPI, GitHub, or built from source.

---

## Prerequisites

- **Python 3.10+**
- **NumPy** (installed automatically as a dependency)

For source builds you additionally need:

- A C compiler (`gcc` or `clang`)
- `meson` and `meson-python` build system
- `numpy` (for `f2py` compilation)

## [PyPI](https://pypi.org/project/smolpack)

For using the PyPI package in your project, add the following to your
configuration file:

=== "pyproject.toml"

    ```toml
    [project]
    dependencies = [
        "smolpack"
    ]
    ```

=== "requirements.txt"

    ```text
    smolpack
    ```

### pip

```bash
pip install --upgrade smolpack
```

### uv

```bash
# Add to a uv project
uv add smolpack

# Or install into the current environment
uv pip install smolpack
```

### pipenv

```bash
pipenv install smolpack
```

### poetry

```bash
poetry add smolpack
```

### pdm

```bash
pdm add smolpack
```

### hatch

```bash
hatch add smolpack
```

## [git](https://github.com/eggzec/smolpack)

Install the latest development version directly from the repository:

```bash
pip install --upgrade "git+https://github.com/eggzec/smolpack.git#egg=smolpack"
```

### Building locally

Clone and build from source to modify the C code or test local changes:

```bash
git clone https://github.com/eggzec/smolpack.git
cd smolpack
uv pip install .
```

This invokes the `meson` build system to compile the C sources via `f2py`
and install the resulting extension module.

!!! warning "C compiler required"
    Source builds require a working C compiler.

    ```bash
    # Debian/Ubuntu
    sudo apt install gcc

    # Fedora
    sudo dnf install gcc

    # macOS (Clang via Xcode Command Line Tools)
    xcode-select --install
    ```

    On Windows, install MSVC via Visual Studio Build Tools. MinGW is not supported.

## Verifying the installation

After installation, verify that the package loads correctly:

```python
import numpy as np
import smolpack

result = smolpack.int_smolyak(lambda d, x: np.exp(np.sum(x)), dim=3, qq=5)
print("smolpack is working! Integral:", result)
```

## Dependencies

- Python >=3.10
- [numpy](https://pypi.org/project/numpy)
