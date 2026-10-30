# Computational Geometry & Spatial Data Indexer (3D)

A dependency-free C++17 MVP for a few useful 3D geometry and spatial-query building blocks. It is intentionally small: it is a reference implementation and learning project, not a numerically certified geometry kernel.

## Build and run

Requirements: CMake 3.16+ and a C++17 compiler. No third-party libraries are used.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/geom3d_demo
```

The tests can be disabled with `-DGEOM3D_BUILD_TESTS=OFF`. Public API is in `include/geom3d/geometry.hpp`; implementations are in `src/geometry.cpp`.

## Included

- `Vec3`: basic arithmetic, indexing, dot/cross products and squared lengths/distances.
- `AABB`: incremental bounds, point-to-box squared distance, and a ray/slab interval test (including axis-parallel rays).
- `orient3d`: a determinant-sign predicate for four points. It first evaluates a double-precision determinant with a conservative error filter and uses a `long double` recomputation when the filtered result is uncertain.
- `KDTree3`: median-built 3D point index, single nearest-neighbor query and radius query. Radius results are ordered by distance and then original input index; nearest ties choose the lowest input index.
- `TriangleBVH`: median-split bounding-volume hierarchy with closest ray/triangle hit queries, two-sided Möller–Trumbore triangle tests and a configurable leaf size.
- `geom3d_demo`: a small CLI demonstrating each implemented feature.

## What is tested

The dependency-free CTest executable covers vector arithmetic; AABB distance and ray intervals (including parallel slab misses); positive/negative/coplanar orientation cases; empty, duplicate, tie and radius-boundary KD-tree cases; 60 deterministic randomized KD nearest/radius comparisons against brute-force references; and BVH closest hit, miss, interval clipping, edge hit and empty-input behavior.

## Limitations and numerical notes

- `orient3d` is **not formally exact**. Its `long double` fallback may have no extra precision on some platforms and remains subject to floating-point rounding/overflow/underflow. Use an adaptive exact-predicate library if topology decisions must be certified.
- Inputs to geometry indexes and predicates must be finite. The code rejects non-finite points/rays and invalid query intervals. Ray/AABB and ray/triangle intersections use widened or scaled intermediates to avoid common overflow/underflow failures, but this is not a general robustness guarantee: `long double` may equal `double` on some platforms, and operations such as KD-tree squared distances can still overflow at extreme magnitudes.
- The triangle routine uses a scale-relative floating-point parallelism threshold and exact range comparisons for barycentric coordinates. Near-degenerate triangles, nearly parallel rays and boundary cases can be sensitive to roundoff. Degenerate triangles are not treated as meaningful surfaces.
- Ray parameter `t` is measured in the supplied direction's units; directions are not normalized. Default ray interval is `[0, 1e300]`.
- The BVH is static after construction, uses median splits, and reports one closest hit; there are no mesh mutation, refit, all-hit, or self-intersection APIs.
- The KD tree is static, stores one point per node, and provides no k-nearest or dynamic update API.
- No convex hull routine is included in this MVP. In particular, there is no placeholder pretending to compute a hull.

## Roadmap

1. Add a clearly named educational convex-hull reference routine with explicit handling/documentation for duplicate, collinear, coplanar and lower-dimensional inputs.
2. Add k-nearest queries and larger randomized/property-style tests for the BVH against brute force.
3. Consider robust/adaptive exact predicates as an optional dependency and document their guarantees separately from this floating-point implementation.
4. Add serialization or point-cloud ingestion only after defining a stable data format and scale/precision expectations.

## License

No license has been selected or included. Add one only after choosing the intended terms.
