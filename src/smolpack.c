/**
 * @file smolpack.c
 * @brief Shared global state for the SMOLPACK library.
 *
 * Provides the single definition of the function-evaluation counter
 * that is shared between both Smolyak algorithm implementations and
 * user-supplied integrand callbacks.
 */

#include "smolpack.h"

/**
 * @brief Global function-evaluation counter.
 *
 * Reset to zero at the start of each integration call.
 * User callbacks should increment this on every evaluation.
 */
int count = 0;
