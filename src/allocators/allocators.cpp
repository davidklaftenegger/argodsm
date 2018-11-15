/**
 * @file
 * @brief This file provides the default allocators
 * @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
 */

#include "collective_allocator.hpp"
#include "dynamic_allocator.hpp"

namespace mem = argo::mempools;
namespace alloc = argo::allocators;

/* default allocators */
alloc::global_allocator<char> alloc::default_global_allocator;
alloc::default_dynamic_allocator_t alloc::default_dynamic_allocator;
alloc::collective_allocator alloc::default_collective_allocator;

extern "C"
void* collective_alloc(size_t size) {
	using namespace argo::allocators;
	/** @bug this is wrong: either it should not be done at all, or also when using the C++ interface */
	argo::backend::barrier();
	return static_cast<void*>(default_collective_allocator.allocate(size));
}

extern "C"
void* collective_alloc_(size_t* size) {
	return collective_alloc(*size);
}

extern "C"
void collective_free(void* ptr) {
	using namespace argo::allocators;
	using atype = decltype(default_collective_allocator)::value_type;
	if (ptr == NULL)
		return;
	default_collective_allocator.free(static_cast<atype*>(ptr));
}

extern "C"
void collective_free_(void* ptr) {
	return collective_free(ptr);
}

extern "C"
void* dynamic_alloc(size_t size) {
	using namespace argo::allocators;
	return static_cast<void*>(default_dynamic_allocator.allocate(size));
}

extern "C"
void* dynamic_alloc_(size_t* size) {
	dynamic_alloc(*size);
}

extern "C"
void dynamic_free(void* ptr) {
	using namespace argo::allocators;
	using atype = decltype(default_dynamic_allocator)::value_type;
	if (ptr == NULL)
		return;
	default_dynamic_allocator.free(static_cast<atype*>(ptr));
}

extern "C"
void dynamic_free_(void* ptr) {
	dynamic_free(ptr);
}
