/**
 * @file
 * @brief This file provides the full C interface to ArgoDSM
 * @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
 */

#ifndef argo_argo_h
#define argo_argo_h argo_argo_h

#include <stddef.h>

#include "allocators/allocators.h"
#include "synchronization/synchronization.h"

/**
 * @brief initialize ArgoDSM system
 * @param size the desired size of the global memory in bytes
 */
void argo_init(size_t size);

/** @fortran{argo_init} */
void argo_init_(size_t* size);

/**
 * @brief shut down ArgoDSM system
 */
void argo_finalize();

/** @fortran{argo_finalize} */
void argo_finalize_();

/**
 * @brief reset ArgoDSM system
 * @warning this is mainly intended for testing, not for public use
 */
void argo_reset();

/** @fortran{argo_reset} */
void argo_reset_();

/**
 * @brief A unique ArgoDSM node identifier. Counting starts from 0.
 * @return The node id
 */
int argo_node_id();

/** @fortran{argo_node_id} */
int argo_node_id_();

/**
 * @brief Number of ArgoDSM nodes being run
 * @return The total number of ArgoDSM nodes
 */
int argo_number_of_nodes();

/** @fortran{argo_number_of_nodes} */
int argo_number_of_nodes_();

#endif /* argo_argo_h */
