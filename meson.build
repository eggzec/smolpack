# ===========================================================================
# Root Meson build for SMOLPACK
# ===========================================================================

project('smolpack', 'c',
    version: '0.1.0',
    meson_version: '>= 1.1.0',
    default_options: [
        'c_std=c11',
        'warning_level=2',
    ],
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# --- NumPy / f2py include directories ---
# Use compile_args instead of include_directories() because the venv
# may live inside the source tree, which meson disallows for inc dirs.
incdir_numpy = run_command(py,
    ['-c', 'import numpy; print(numpy.get_include())'],
    check: true,
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check: true,
).stdout().strip()

np_compile_args = ['-I' + incdir_numpy, '-I' + incdir_f2py]

# f2py's fortranobject.c (compiled into every f2py extension)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# --- Subdirectories ---
subdir('src')
subdir('smolpack')
