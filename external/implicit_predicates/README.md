Vendored from Marco Attene's "Indirect predicates for geometric constructions"
(Computer-Aided Design, 2020) -- exact predicates over points that may be
explicit or an implicit combination of other points (e.g. a Steiner point
introduced by an edge split).

Single include in your code:
    #include "indirect_predicates.h"
That one include pulls the entire chain. Just keep all 8 files in the same directory.

One runtime call needed before using any predicates (once per thread):
    initFPU();
This sets up the floating-point environment (defined in numerics.hpp).

Compiler flags required:
    - GCC/Clang: -frounding-math (needed for interval arithmetic's FPU rounding switches to work correctly)
    - MSVC: /fp:strict

## Local modification and the geometry-facing interface

`implicitPoint3D_LNC` (the linear-combination point type -- the only implicit
point type this repo uses; `LPI`/`TPI`/`SSI` and `explicitPoint2D` below are
unmodified upstream code) has been changed to store two point ids plus a
pointer to an external `s_dynarray`, instead of two direct `explicitPoint3D`
references. `P()` / `Q()` resolve the ids against that dynarray on every
call, which lets the point registry live in ordinary realloc-based storage
instead of a pointer-stable container.

`implicit_predicate_registry.cpp` in this directory is the sole compiled translation
unit for this otherwise header-only library. It owns a caller-allocated
`implicit_predicate_reg` (see the file for the struct layout) and exposes it
through the pure-C API declared in `geometry/include/gtests_implicit.h` --
no global state, so independent callers (e.g. one per thread) can each hold
their own registry.
