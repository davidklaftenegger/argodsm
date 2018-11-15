/**
 * @file
 * @brief This file provides the collective allocator for ArgoDSM
 * @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
 */

#ifndef argo_collective_allocators_h
#define argo_collective_allocators_h argo_collective_allocators_h

/**
 * @brief basic collective allocation function, for C interface
 * @param size number of bytes to allocate
 * @return pointer to the newly allocated memory
 */
void* collective_alloc(size_t size);

/** @fortran{collective_alloc} */
void* collective_alloc_(size_t* size);

/**
 * @brief basic free function for collective allocations, for C interface
 * @param ptr pointer to free
 */
void collective_free(void* ptr);

/** @fortran{collective_free} */
void collective_free_(void* ptr);

#endif /* argo_collective_allocators_h */
