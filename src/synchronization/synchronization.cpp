/**
 * @file
 * @brief This file implements synchronization primitives for ArgoDSM
 * @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
 */

#include "../backend/backend.hpp"

namespace argo {
	void barrier(std::size_t threadcount) {
		backend::barrier(threadcount);
	}
} // namespace argo

extern "C" {
	void argo_barrier(size_t threadcount) {
		argo::barrier(threadcount);
	}

	void argo_barrier_(size_t* threadcount) {
		argo_barrier(*threadcount);
	}
}
