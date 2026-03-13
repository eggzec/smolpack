/**
 * @file _smolpack_helpers.c
 * @brief Minimal helper for f2py Python bindings.
 *
 * Exposes the global function-evaluation counter as a callable.
 */

#include "smolpack.h"

int get_count(void) {
    return count;
}
