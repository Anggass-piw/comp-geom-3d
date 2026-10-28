#include "geom3d/geometry.hpp"

#include <iostream>
#include <vector>

int main() {
    using namespace geom3d;

    const auto side = orient3d({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    std::cout << "orient3d sign: " << static_cast<int>(side) << " (-1=negative, 0=coplanar, 1=positive)\n";

    KDTree3 points({{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2}});
    const Vec3 query{0.8, 0.1, 0.1};
    if (const auto nearest = points.nearest(query)) {
        std::cout << "nearest point index: " << nearest->index
                  << ", distance: " << nearest->distance << '\n';
    }
    std::cout << "points within radius 1.0: " << points.radius(query, 1.0).size() << '\n';

    const AABB box({-1, -1, -1}, {1, 1, 1});
    const Ray ray{{0, 0, 3}, {0, 0, -1}};
    if (const auto interval = box.intersect_ray(ray)) {
        std::cout << "AABB ray interval: [" << interval->first << ", " << interval->second << "]\n";
    }

    TriangleBVH bvh({{{-1, -1, 0}, {1, -1, 0}, {0, 1, 0}},
                     {{-1, -1, 2}, {1, -1, 2}, {0, 1, 2}}});
    if (const auto hit = bvh.closest_hit(ray)) {
        std::cout << "closest triangle hit: index " << hit->triangle_index
                  << " at t=" << hit->t << " (z=" << hit->position.z << ")\n";
    } else {
        std::cout << "no triangle hit\n";
    }
    return 0;
}
