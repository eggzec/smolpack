
# Windows static link args logic (local scope)
if host_machine.system() == 'windows'
    c_link_args = [
        '-static-libgcc',
        '-Wl,-Bstatic',
        '-lwinpthread',
        '-Wl,-Bdynamic',
    ]
else
    c_link_args = []
endif

_smolpack_module = custom_target('_smolpackmodule',
    input: '_smolpack.pyf',
    output: ['_smolpackmodule.c'],
    command: [py, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@'],
)

py.extension_module('_smolpack',
    [_smolpack_module, '_smolpack_helpers.c', fortranobject_c],
    include_directories: smolpack_inc,
    c_args: np_compile_args,
    link_args: c_link_args,
    link_with: smolpack_lib,
    dependencies: [py_dep, m_dep],
    install: true,
    subdir: 'smolpack',
)

py.install_sources(
    '__init__.py',
    subdir: 'smolpack',
)
