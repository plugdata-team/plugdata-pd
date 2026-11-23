/* Copyright (c) 2026 Timothy Schoen
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */

/* Windows plugins run on host threads (often 1 MB of stack) and macOS audio
threads get 512 KB, so there recursion also stops when the stack runs low.
Linux threads have 8mb stack by default, so we should be in the clear there.
*/

#if defined(__APPLE__) || (PDTHREADS && defined(PDINSTANCE) && defined(_WIN32))

#include <stdint.h>
#ifdef _WIN32
#include <windows.h>
#define STACK_GUARD_THREAD PERTHREAD
#else
#include <pthread.h>
#define STACK_GUARD_THREAD __thread
#endif

#ifdef _MSC_VER
#include <intrin.h>
#define STACK_GUARD_SP ((uintptr_t)_AddressOfReturnAddress())
#else
#define STACK_GUARD_SP ((uintptr_t)__builtin_frame_address(0))
#endif

    /* shallower messages skip the check, deeper ones check every 4th hop */
#define STACK_GUARD_DEPTH 32
    /* stack left for four more hops and for reporting the overflow */
#define STACK_GUARD_MARGIN (128 * 1024)

    /* starts above any stack, so the first check looks up the real limit */
static STACK_GUARD_THREAD uintptr_t stack_guard_limit = UINTPTR_MAX;

    /* lowest address of the stack that sp is on, or 0 if unknown */
static uintptr_t stack_guard_findbase(uintptr_t sp)
{
#ifdef _WIN32
    MEMORY_BASIC_INFORMATION info;
    return VirtualQuery((void *)sp, &info, sizeof(info)) ? (uintptr_t)info.AllocationBase : 0;
#else
    pthread_t self = pthread_self();
    uintptr_t top = (uintptr_t)pthread_get_stackaddr_np(self);
    size_t size = pthread_get_stacksize_np(self);
    return (sp <= top && top - sp < size) ? top - size : 0;
#endif
}

static NOINLINE int stack_guard_check(void)
{
    uintptr_t sp = STACK_GUARD_SP;
    if (sp >= stack_guard_limit)
        return 0;
        /* first check on this thread, or sp moved to another stack (e.g. a coroutine) */
    if (stack_guard_limit == UINTPTR_MAX || sp < stack_guard_limit - STACK_GUARD_MARGIN)
    {
        uintptr_t base = stack_guard_findbase(sp);
        stack_guard_limit = base ? base + STACK_GUARD_MARGIN : 0;
    }
    return sp < stack_guard_limit;
}

#define stack_guard_exhausted(depth) \
    ((depth) >= STACK_GUARD_DEPTH && !((depth) & 3) && stack_guard_check())

#else
#define stack_guard_exhausted(depth) 0
#endif
