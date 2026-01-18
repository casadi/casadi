//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//


// FILTER-MACROS OFF

#define CASADI_THREAD_TYPE_NONE 0
#define CASADI_THREAD_TYPE_POSIX 1
#define CASADI_THREAD_TYPE_C11 2
#define CASADI_THREAD_TYPE_OMP 3
#define CASADI_THREAD_TYPE_WINDOWS 4

#ifndef CASADI_THREAD_TYPE
  #define CASADI_THREAD_TYPE CASADI_THREAD_TYPE_POSIX
#endif

#if CASADI_THREAD_TYPE == CASADI_THREAD_TYPE_NONE
  #ifndef CASADI_MUTEX_USE_STATIC_INIT
    #define CASADI_MUTEX_USE_STATIC_INIT 0
  #endif
  #ifndef CASADI_MUTEX_TYPE
    #define CASADI_MUTEX_TYPE int
  #endif
  #ifndef CASADI_MUTEX_LOCK
    #define CASADI_MUTEX_LOCK(m) ((void)(m))
  #endif
  #ifndef CASADI_MUTEX_UNLOCK
    #define CASADI_MUTEX_UNLOCK(m) ((void)(m))
  #endif
  #ifndef CASADI_MUTEX_INIT
    #define CASADI_MUTEX_INIT(m) ((void)(m))
  #endif
  #ifndef CASADI_MUTEX_DESTROY
    #define CASADI_MUTEX_DESTROY(m) ((void)(m))
  #endif
  #ifndef CASADI_MUTEX_STATIC_INIT
    #define CASADI_MUTEX_STATIC_INIT 0
  #endif
#elif CASADI_THREAD_TYPE == CASADI_THREAD_TYPE_POSIX
  #ifndef CASADI_MUTEX_USE_STATIC_INIT
    #define CASADI_MUTEX_USE_STATIC_INIT 1
  #endif

  #include <pthread.h>
  #ifndef CASADI_MUTEX_TYPE
    #define CASADI_MUTEX_TYPE pthread_mutex_t
  #endif
  #ifndef CASADI_MUTEX_LOCK
    #define CASADI_MUTEX_LOCK(m) pthread_mutex_lock(m)
  #endif
  #ifndef CASADI_MUTEX_UNLOCK
    #define CASADI_MUTEX_UNLOCK(m) pthread_mutex_unlock(m)
  #endif
  #ifndef CASADI_MUTEX_INIT
    #define CASADI_MUTEX_INIT(m) pthread_mutex_init(m, NULL)
  #endif
  #ifndef CASADI_MUTEX_DESTROY
    #define CASADI_MUTEX_DESTROY(m) pthread_mutex_destroy(m)
  #endif
  #ifndef CASADI_MUTEX_STATIC_INIT
    #define CASADI_MUTEX_STATIC_INIT PTHREAD_MUTEX_INITIALIZER
  #endif

  #ifndef CASADI_THREAD_HANDLE
    #define CASADI_THREAD_HANDLE pthread_t
  #endif
  #ifndef CASADI_THREAD_CREATE
    #define CASADI_THREAD_CREATE(t, f, a) pthread_create(&(t), NULL, f, a)
  #endif
  #ifndef CASADI_THREAD_JOIN
    #define CASADI_THREAD_JOIN(t) pthread_join(t, NULL)
  #endif
  #ifndef CASADI_THREAD_WORKER_RETURN
    #define CASADI_THREAD_WORKER_RETURN void*
  #endif
  #ifndef CASADI_THREAD_WORKER_ARG
    #define CASADI_THREAD_WORKER_ARG void*
  #endif
  #ifndef CASADI_THREAD_RETURN_VALUE
    #define CASADI_THREAD_RETURN_VALUE NULL
  #endif

#elif CASADI_THREAD_TYPE == CASADI_THREAD_TYPE_C11
  #ifndef CASADI_MUTEX_USE_STATIC_INIT
    #define CASADI_MUTEX_USE_STATIC_INIT 0
  #endif

  #include <threads.h>
  #ifndef CASADI_MUTEX_TYPE
    #define CASADI_MUTEX_TYPE mtx_t
  #endif
  #ifndef CASADI_MUTEX_LOCK
    #define CASADI_MUTEX_LOCK(m) mtx_lock(m)
  #endif
  #ifndef CASADI_MUTEX_UNLOCK
    #define CASADI_MUTEX_UNLOCK(m) mtx_unlock(m)
  #endif
  #ifndef CASADI_MUTEX_INIT
    #define CASADI_MUTEX_INIT(m) mtx_init(m, mtx_plain)
  #endif
  #ifndef CASADI_MUTEX_DESTROY
    #define CASADI_MUTEX_DESTROY(m) mtx_destroy(m)
  #endif

  #ifndef CASADI_THREAD_HANDLE
    #define CASADI_THREAD_HANDLE thrd_t
  #endif
  #ifndef CASADI_THREAD_CREATE
    #define CASADI_THREAD_CREATE(t, f, a) thrd_create(&(t), f, a)
  #endif
  #ifndef CASADI_THREAD_JOIN
    #define CASADI_THREAD_JOIN(t) thrd_join(t, NULL)
  #endif
  #ifndef CASADI_THREAD_WORKER_RETURN
    #define CASADI_THREAD_WORKER_RETURN int
  #endif
  #ifndef CASADI_THREAD_WORKER_ARG
    #define CASADI_THREAD_WORKER_ARG void*
  #endif
  #ifndef CASADI_THREAD_RETURN_VALUE
    #define CASADI_THREAD_RETURN_VALUE 0
  #endif

#elif CASADI_THREAD_TYPE == CASADI_THREAD_TYPE_OMP
  #ifndef CASADI_MUTEX_USE_STATIC_INIT
    #define CASADI_MUTEX_USE_STATIC_INIT 0
  #endif

  #include <omp.h>
  #ifndef CASADI_MUTEX_TYPE
    #define CASADI_MUTEX_TYPE omp_lock_t
  #endif
  #ifndef CASADI_MUTEX_LOCK
    #define CASADI_MUTEX_LOCK(m) omp_set_lock(m)
  #endif
  #ifndef CASADI_MUTEX_UNLOCK
    #define CASADI_MUTEX_UNLOCK(m) omp_unset_lock(m)
  #endif
  #ifndef CASADI_MUTEX_INIT
    #define CASADI_MUTEX_INIT(m) omp_init_lock(m)
  #endif
  #ifndef CASADI_MUTEX_DESTROY
    #define CASADI_MUTEX_DESTROY(m) omp_destroy_lock(m)
  #endif

#elif CASADI_THREAD_TYPE == CASADI_THREAD_TYPE_WINDOWS
  #ifndef CASADI_MUTEX_USE_STATIC_INIT
    #define CASADI_MUTEX_USE_STATIC_INIT 1
  #endif

  #include <windows.h>
  #ifdef CASADI_MUTEX_USE_STATIC_INIT
    #ifndef CASADI_MUTEX_TYPE
      #define CASADI_MUTEX_TYPE SRWLOCK
    #endif
    #ifndef CASADI_MUTEX_LOCK
      #define CASADI_MUTEX_LOCK(m) AcquireSRWLockExclusive(m)
    #endif
    #ifndef CASADI_MUTEX_UNLOCK
      #define CASADI_MUTEX_UNLOCK(m) ReleaseSRWLockExclusive(m)
    #endif
    #ifndef CASADI_MUTEX_INIT
      #define CASADI_MUTEX_INIT(m) InitializeSRWLock(m)
    #endif
    #ifndef CASADI_MUTEX_DESTROY
      #define CASADI_MUTEX_DESTROY(m) ((void)(m))
    #endif
    #ifndef CASADI_MUTEX_STATIC_INIT
      #define CASADI_MUTEX_STATIC_INIT SRWLOCK_INIT
    #endif
  #else
    #ifndef CASADI_MUTEX_TYPE
      #define CASADI_MUTEX_TYPE CRITICAL_SECTION
    #endif
    #ifndef CASADI_MUTEX_LOCK
      #define CASADI_MUTEX_LOCK(m) EnterCriticalSection(m)
    #endif
    #ifndef CASADI_MUTEX_UNLOCK
      #define CASADI_MUTEX_UNLOCK(m) LeaveCriticalSection(m)
    #endif
    #ifndef CASADI_MUTEX_INIT
      #define CASADI_MUTEX_INIT(m) InitializeCriticalSection(m)
    #endif
    #ifndef CASADI_MUTEX_DESTROY
      #define CASADI_MUTEX_DESTROY(m) DeleteCriticalSection(m)
    #endif
  #endif

  #ifndef CASADI_THREAD_HANDLE
    #define CASADI_THREAD_HANDLE HANDLE
  #endif
  #ifndef CASADI_THREAD_CREATE
    #define CASADI_THREAD_CREATE(t, f, a) ((t) = CreateThread(NULL, 0, f, a, 0, NULL))
  #endif
  #ifndef CASADI_THREAD_JOIN
    #define CASADI_THREAD_JOIN(t) (WaitForSingleObject(t, INFINITE), CloseHandle(t))
  #endif
  #ifndef CASADI_THREAD_WORKER_RETURN
    #define CASADI_THREAD_WORKER_RETURN DWORD WINAPI
  #endif
  #ifndef CASADI_THREAD_WORKER_ARG
    #define CASADI_THREAD_WORKER_ARG LPVOID
  #endif
  #ifndef CASADI_THREAD_RETURN_VALUE
    #define CASADI_THREAD_RETURN_VALUE 0
  #endif

#else
  #error "Unknown CASADI_THREAD_TYPE. Use: CASADI_THREAD_TYPE_NONE, CASADI_THREAD_TYPE_POSIX, CASADI_THREAD_TYPE_C11, CASADI_THREAD_TYPE_OMP, or CASADI_THREAD_TYPE_WINDOWS" // NOLINT(whitespace/line_length)
#endif

// FILTER-MACROS ON
