/**
 * @file
 * @brief This file provides C bindings for some ArgoDSM synchronization primitives
 * @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
 */

#ifndef argo_synchronization_h
#define argo_synchronization_h argo_synchronization_h

#include <stddef.h>

/**
 * @brief a barrier for threads
 * @param threadcount number of threads on each ArgoDSM node
 * @details this barrier waits until threadcount threads have reached this
 *          barrier call on EACH ArgoDSM node.
 * @todo better explanation
 */
void argo_barrier(size_t threadcount);

/** @fortran{argo_barrier} */
void argo_barrier_(size_t* threadcount);

#endif /* argo_synchronization_h */
