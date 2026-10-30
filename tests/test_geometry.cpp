#include "geom3d/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
int checks = 0;

void check(bool condition, const std::string& message) {
    ++checks;
    if (!condition) throw std::runtime_error("CHECK failed: " + message);
}

bool close(double a, double b, double epsilon = 1e-10) {
    return std::abs(a - b) <= epsilon;
}

void test_vec_and_aabb() {
    using namespace geom3d;
    check(close(dot({1, 2, 3}, {4, -1, 2}), 8.0), "dot product");
    const Vec3 c = cross({1, 0, 0}, {0, 1, 0});
    check(close(c.z, 1.0), "cross product");
    AABB box;
    box.expand({-1, -2, -3});
    box.expand({1, 2, 3});
    check(close(box.squared_distance({0, 0, 0}), 0.0), "point inside box");
    check(close(box.squared_distance({4, 2, 3}), 9.0), "point-box distance");
    const auto interval = box.intersect_ray({{0, 0, 5}, {0, 0, -1}});
    check(interval && close(interval->first, 2.0) && close(interval->second, 8.0), "AABB ray interval");
    check(!box.intersect_ray({{5, 0, 0}, {0, 1, 0}}), "parallel ray outside slab misses");
}

void test_orientation() {
    using namespace geom3d;
    check(orient3d({0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}) == Orientation::negative,
          "orient3d sign for canonical tetrahedron");
    check(orient3d({0,0,0}, {1,0,0}, {0,1,0}, {0.2,0.3,0}) == Orientation::coplanar,
          "coplanar points");
    check(orient3d({0,0,0}, {0,1,0}, {1,0,0}, {0,0,1}) == Orientation::positive,
          "orientation changes sign on swapping vertices");
}

void test_kdtree() {
    using namespace geom3d;
    KDTree3 empty({});
    check(!empty.nearest({0,0,0}), "empty tree nearest");
    check(empty.radius({0,0,0}, 1.0).empty(), "empty tree radius");

    KDTree3 simple({{1,0,0}, {-1,0,0}, {0,0,0}, {0,0,0}});
    const auto tie = simple.nearest({0,0,0});
    check(tie && tie->index == 2 && close(tie->distance, 0.0), "nearest tie selects lowest original index");
    const auto around = simple.radius({0,0,0}, 1.0);
    check(around.size() == 4, "radius includes boundary and duplicates");
    check(around[0].index == 2 && around[1].index == 3, "radius deterministic tie ordering");

    std::mt19937_64 rng(20261028);
    std::uniform_real_distribution<double> dist(-50.0, 50.0);
    std::vector<Vec3> data;
    for (int i = 0; i < 250; ++i) data.push_back({dist(rng), dist(rng), dist(rng)});
    KDTree3 tree(data);
    for (int q = 0; q < 60; ++q) {
        const Vec3 query{dist(rng), dist(rng), dist(rng)};
        std::size_t expected = 0;
        double best2 = squared_distance(query, data[0]);
        for (std::size_t i = 1; i < data.size(); ++i) {
            const double candidate = squared_distance(query, data[i]);
            if (candidate < best2 || (candidate == best2 && i < expected)) { best2 = candidate; expected = i; }
        }
        const auto actual = tree.nearest(query);
        check(actual && actual->index == expected, "KD nearest agrees with brute-force reference");
        const double radius = 15.0;
        const auto found = tree.radius(query, radius);
        std::vector<std::size_t> expected_indices;
        for (std::size_t i = 0; i < data.size(); ++i)
            if (squared_distance(query, data[i]) <= radius * radius) expected_indices.push_back(i);
        check(found.size() == expected_indices.size(), "KD radius agrees with brute-force count");
    }
}

void test_bvh() {
    using namespace geom3d;
    const Triangle lower{{-1,-1,0}, {1,-1,0}, {0,1,0}};
    const Triangle upper{{-1,-1,2}, {1,-1,2}, {0,1,2}};
    TriangleBVH bvh({lower, upper}, 1);
    const auto hit = bvh.closest_hit({{0,0,3}, {0,0,-1}});
    check(hit && hit->triangle_index == 1 && close(hit->t, 1.0), "BVH returns closest triangle");
    check(hit && close(hit->position.z, 2.0), "BVH hit position");
    check(!bvh.closest_hit({{3,3,3}, {0,0,-1}}), "BVH miss");
    check(!bvh.closest_hit({{0,0,3}, {0,0,-1},}, 0.0, 0.5), "ray interval clips hits");
    const auto edge_hit = bvh.closest_hit({{1,-1,3}, {0,0,-1}});
    check(edge_hit.has_value(), "triangle edge hit is accepted");
    TriangleBVH empty({}, 4);
    check(!empty.closest_hit({{0,0,1}, {0,0,-1}}), "empty BVH");
}

void test_ray_scale_edges() {
    using namespace geom3d;

    AABB unit_box({-1, -1, -1}, {1, 1, 1});
    bool zero_direction_rejected = false;
    try {
        (void)unit_box.intersect_ray({{0, 0, 0}, {0, 0, 0}});
    } catch (const std::invalid_argument&) {
        zero_direction_rejected = true;
    }
    check(zero_direction_rejected, "AABB rejects a zero-direction ray");

    const Triangle unit_triangle{{-1, -1, 0}, {1, -1, 0}, {0, 1, 0}};
    TriangleBVH unit_bvh({unit_triangle});
    const auto tiny_direction_hit = unit_bvh.closest_hit({{0, 0, 2e-200}, {0, 0, -1e-200}});
    check(tiny_direction_hit && close(tiny_direction_hit->t, 2.0),
          "ray direction remains valid when its squared length underflows");

    const double tiny = 1e-200;
    const Triangle tiny_triangle{{0, 0, 0}, {tiny, 0, 0}, {0, tiny, 0}};
    TriangleBVH tiny_bvh({tiny_triangle});
    const auto tiny_geometry_hit = tiny_bvh.closest_hit({{0.2*tiny, 0.2*tiny, 2*tiny}, {0, 0, -tiny}});
    check(tiny_geometry_hit && close(tiny_geometry_hit->t, 2.0),
          "triangle intersection remains stable when squared edge scales underflow");

    const double large = 8e307;
    const double larger = 9e307;
    const Triangle large_triangle{{large, large, large}, {larger, large, large}, {large, larger, large}};
    TriangleBVH large_bvh({large_triangle});
    const auto large_hit = large_bvh.closest_hit({{8.2e307, 8.2e307, larger}, {0, 0, -1e307}});
    check(large_hit && close(large_hit->t, 1.0),
          "BVH handles large finite coordinates without overflowing centroid or determinant");
    check(large_hit && close(large_hit->position.z / large, 1.0),
          "large-coordinate hit position remains finite and on the triangle");

    AABB extreme_slab({-1, -1, -1.6e308}, {1, 1, -1.5e308});
    const auto extreme_interval = extreme_slab.intersect_ray({{0, 0, 1.6e308}, {0, 0, -1.6e308}});
    check(extreme_interval && close(extreme_interval->first, 1.9375) &&
          close(extreme_interval->second, 2.0),
          "AABB ray slab handles a finite coordinate difference beyond double range");

    AABB large_box({large, large, large}, {larger, larger, larger});
    const Vec3 large_center = large_box.centroid();
    check(close(large_center.x / 8.5e307, 1.0) && close(large_center.y / 8.5e307, 1.0),
          "AABB centroid avoids overflow for same-sign large bounds");
}

} // namespace

int main() {
    try {
        test_vec_and_aabb();
        test_orientation();
        test_kdtree();
        test_bvh();
        test_ray_scale_edges();
        std::cout << "All " << checks << " checks passed.\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
