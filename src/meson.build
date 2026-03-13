# ===========================================================================
# Build SMOLPACK C library as a static archive
# ===========================================================================

smolpack_sources = files(
    'smolpack.c',
    'smolyak_init.c',
    'smolyak_weights.c',
    'smolyak_alg.c',
    'cc_init.c',
    'cc_weights.c',
    'cc_alg.c',
)

smolpack_inc = include_directories('.')

smolpack_lib = static_library('smolpack',
    smolpack_sources,
    include_directories: smolpack_inc,
    dependencies: m_dep,
)
