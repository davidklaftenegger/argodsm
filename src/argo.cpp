/**
 * @file
 * @brief This file implements some of the basic ArgoDSM calls
 * @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
 */

#include "argo.hpp"

#include "allocators/collective_allocator.hpp"
#include "allocators/dynamic_allocator.hpp"

namespace mem = argo::mempools;
namespace alloc = argo::allocators;

/* some memory pools for default use */
/** @todo should be static? */
mem::global_memory_pool<>* default_global_mempool;
mem::dynamic_memory_pool<alloc::global_allocator, mem::NODE_ZERO_ONLY> collective_prepool(&alloc::default_global_allocator);
mem::dynamic_memory_pool<alloc::global_allocator, mem::ALWAYS> dynamic_prepool(&alloc::default_global_allocator);

namespace argo {
	void init(size_t size) {
		default_global_mempool = new mem::global_memory_pool<>(size);
		argo_reset();
	}

	void finalize() {
		delete default_global_mempool;
	}

	int node_id() {
		return static_cast<int>(argo::backend::node_id());
	}

	int number_of_nodes() {
		return static_cast<int>(argo::backend::number_of_nodes());
	}

} // namespace argo

extern "C" {
	void argo_init(size_t size) {
		argo::init(size);
	}
	void argo_init_(size_t* size) {
		argo_init(*size);
	}

	void argo_finalize() {
		argo::finalize();
	}
	void argo_finalize_() {
		argo_finalize();
	}

	void argo_reset() {
		default_global_mempool->reset();
		using namespace alloc;
		using namespace mem;
		collective_prepool = dynamic_memory_pool<global_allocator, NODE_ZERO_ONLY>(&default_global_allocator);
		dynamic_prepool = dynamic_memory_pool<global_allocator, ALWAYS>(&default_global_allocator);
		default_global_allocator = global_allocator<char>();
		default_dynamic_allocator = default_dynamic_allocator_t();
		default_collective_allocator = collective_allocator();
		default_global_allocator.set_mempool(default_global_mempool);
		default_dynamic_allocator.set_mempool(&dynamic_prepool);
		default_collective_allocator.set_mempool(&collective_prepool);
	}

	void argo_reset_() {
		argo_reset();
	}

	int argo_node_id() {
		return argo::node_id();
	}

	int argo_node_id_() {
		return argo_node_id();
	}

	int argo_number_of_nodes() {
		return argo::number_of_nodes();
	}

	int argo_number_of_nodes_() {
		return argo_number_of_nodes();
	}
}
