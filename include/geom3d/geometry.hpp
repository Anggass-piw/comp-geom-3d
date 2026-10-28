#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace geom3d {

struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};

    double& operator[](std::size_t axis);
    double operator[](std::size_t axis) const;

    Vec3& operator+=(const Vec3& rhs);
    Vec3& operator-=(const Vec3& rhs);
};

Vec3 operator+(Vec3 lhs, const Vec3& rhs);
Vec3 operator-(Vec3 lhs, const Vec3& rhs);
Vec3 operator*(const Vec3& v, double scalar);
Vec3 operator*(double scalar, const Vec3& v);
Vec3 operator/(const Vec3& v, double scalar);
double dot(const Vec3& a, const Vec3& b);
Vec3 cross(const Vec3& a, const Vec3& b);
double squared_length(const Vec3& v);
double squared_distance(const Vec3& a, const Vec3& b);

struct Ray {
    Vec3 origin;
    Vec3 direction;
};

class AABB {
public:
    AABB();
    AABB(const Vec3& min_corner, const Vec3& max_corner);

    const Vec3& min() const noexcept { return min_; }
    const Vec3& max() const noexcept { return max_; }
    bool empty() const noexcept { return empty_; }
    void expand(const Vec3& point);
    void expand(const AABB& box);
    Vec3 centroid() const;
    double squared_distance(const Vec3& point) const;

    // Slab test. Returns the clipped [t_min, t_max] interval, or nullopt on a miss.
    std::optional<std::pair<double, double>> intersect_ray(
        const Ray& ray, double t_min = 0.0, double t_max = 1.0e300) const;

private:
    Vec3 min_{};
    Vec3 max_{};
    bool empty_{true};
};

// Sign of orient3d(a,b,c,d) = det(a-d, b-d, c-d).
// Positive and negative describe the determinant sign; coplanar is returned when
// both the double filter and long-double fallback evaluate to zero. This is a
// filtered floating-point predicate, NOT a formally exact/adaptive predicate.
enum class Orientation { negative = -1, coplanar = 0, positive = 1 };
Orientation orient3d(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

struct Neighbor {
    std::size_t index{};
    Vec3 point{};
    double distance{};
};

class KDTree3 {
public:
    explicit KDTree3(std::vector<Vec3> points);
    std::size_t size() const noexcept { return points_.size(); }

    std::optional<Neighbor> nearest(const Vec3& query) const;
    // Results are ordered by increasing distance, then by original point index.
    std::vector<Neighbor> radius(const Vec3& query, double radius) const;

private:
    struct Node {
        std::size_t point_index{};
        std::size_t left{static_cast<std::size_t>(-1)};
        std::size_t right{static_cast<std::size_t>(-1)};
        std::size_t axis{};
        AABB bounds;
    };

    std::size_t build(std::size_t begin, std::size_t end, std::size_t depth);
    void nearest_impl(std::size_t node, const Vec3& query, std::size_t& best,
                      double& best_distance2) const;
    void radius_impl(std::size_t node, const Vec3& query, double radius2,
                     std::vector<std::size_t>& found) const;

    std::vector<Vec3> points_;
    std::vector<std::size_t> order_;
    std::vector<Node> nodes_;
};

struct Triangle {
    Vec3 a;
    Vec3 b;
    Vec3 c;
};

struct TriangleHit {
    std::size_t triangle_index{};
    double t{};
    double u{};
    double v{};
    Vec3 position{};
};

class TriangleBVH {
public:
    explicit TriangleBVH(std::vector<Triangle> triangles, std::size_t leaf_size = 4);
    std::size_t size() const noexcept { return triangles_.size(); }

    // Two-sided Moller-Trumbore intersection; returns the closest hit in the ray interval.
    std::optional<TriangleHit> closest_hit(
        const Ray& ray, double t_min = 0.0, double t_max = 1.0e300) const;

private:
    struct Node {
        AABB bounds;
        std::size_t begin{};
        std::size_t end{};
        std::size_t left{static_cast<std::size_t>(-1)};
        std::size_t right{static_cast<std::size_t>(-1)};
        bool leaf{true};
    };

    std::size_t build(std::size_t begin, std::size_t end);
    void hit_impl(std::size_t node, const Ray& ray, double t_min,
                  double& best_t, std::optional<TriangleHit>& best) const;

    std::vector<Triangle> triangles_;
    std::vector<std::size_t> order_;
    std::vector<Node> nodes_;
    std::size_t leaf_size_;
};

} // namespace geom3d
