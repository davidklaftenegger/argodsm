/**
 * @file
 * @brief This file provides a tas lock for the ArgoDSM system based on the TAS lock made by David Klaftenegger
 * @copyright Eta Scale AB. Licensed under the Eta Scale Open Source License. See the LICENSE file for details.
 */

#ifndef argo_global_tas_lock_hpp
#define argo_global_tas_lock_hpp argo_global_lock_hpp

#include "../backend/backend.hpp"
#include "../data_distribution/data_distribution.hpp"
#include <chrono>
#include <thread>

namespace argo {
	namespace globallock {
		/** @brief a global test-and-set lock */
		class global_tas_lock {
			private:
				/** @brief constant signifying lock is in an initial state and free */
				static const std::size_t init = -2;
				/** @brief constant signifying lock is taken */
				static const std::size_t locked = -1;

				/** @brief import global_ptr */
				using global_size_t = typename argo::data_distribution::global_ptr<std::size_t>;

				/**
				 * @brief pointer to lock field
				 * @todo should be replaced with an ArgoDSM-specific atomic type
				 *       to allow efficient synchronization over more backends
				 */
				global_size_t lastuser;

			public:
				/**
				 * @brief construct global tas lock from existing memory in global address space
				 * @param f pointer to global field for storing lock state
				 */
				global_tas_lock(std::size_t* f) : lastuser(global_size_t(f)) {
					*lastuser = init;
				};

				/**
				 * @brief try to lock
				 * @return true if lock was successfully taken,
				 *         false otherwise
				 */
				bool try_lock() {
					auto old = backend::atomic::exchange(lastuser, locked, atomic::memory_order::relaxed);
					if(old != locked) {
						std::size_t self = backend::node_id();
						if(old == self || old == init) {
							/* note: this here is NOT a release.
							 *       instead, it should be a backend-provided
							 *       call to enforce ordering for writebacks.
							 *       A trivial implementation would call a
							 *       self-downgrade (as release() does), but a
							 *       better implementation could be thought of
							 *       if a better write-buffer is also implemented.
							 * why: semantically, we acquire at the beginning of
							 *      a lock. any OTHER node can therefore rely on
							 *      us having completed any previous work before
							 *      entering the lock if they see any writes from
							 *      within this lock's critical section. As write
							 *      buffers can be cleared at any time without
							 *      ordering guarantees, this must be prevented.
							 */
							backend::release();
							/* note: on non-TSO architectures, here a node-local
							 *       acquire synchronization is needed. To ensure
							 *       this happens we might be forced to put a
							 *       local lock somewhere in this global lock.
							 */
						} else {
							backend::acquire();
						}
						return true;
					}
					else {
						return false;
					}
				}

				/**
				 * @brief release the lock
				 */
				void unlock() {
					std::size_t self = backend::node_id();
					backend::release();
					backend::atomic::store(lastuser, self);
				}

				/**
				 * @brief take the lock
				 */
				void lock() {
					while(!try_lock())
						std::this_thread::yield();
				}
		};
	} // namespace globallock
} // namespace argo

#endif /* argo_global_tas_lock_hpp */
