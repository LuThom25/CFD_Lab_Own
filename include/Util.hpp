#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _OPENMP
#define OMP_CALL(_ARGSTR) _Pragma(_ARGSTR)
#else
#define OMP_CALL(_ARGSTR) /* no-op */
#endif

#ifdef __GNUG__
#define HOT [[gnu::hot]]

#define INLINE [[gnu::always_inline]]

#define NOTHROW [[gnu::nothrow]]
#else
#define HOT     /* no-op */
#define INLINE  /* no-op */
#define NOTHROW /* no-op */
#endif
