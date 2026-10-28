#include "geom3d/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>

namespace geom3d {
namespace {

constexpr std::size_t kNone = static_cast<std::size_t>(-1);

bool finite(const Vec3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

void require_finite(const Vec3& v, const char* what) {
    if (!finite(v)) throw std::invalid_argument(std::string(what) + " must be finite");
}

void require_ray(const Ray& ray, double t_min, double t_max) {
    require_finite(ray.origin, "ray origin");
    require_finite(ray.direction, "ray direction");
    if (squared_length(ray.direction) == 0.0)
        throw std::invalid_argument("ray direction must be nonzero");
    if (!std::isfinite(t_min) || !std::isfinite(t_max) || t_min > t_max)
        throw std::invalid_argument("ray interval must be finite and ordered");
}

AABB triangle_bounds(const Triangle& tri) {
    AABB box;
    box.expand(tri.a);
    box.expand(tri.b);
    box.expand(tri.c);
    return box;
}

Vec3 triangle_centroid(const Triangle& tri) {
    return {(tri.a.x + tri.b.x + tri.c.x) / 3.0,
            (tri.a.y + tri.b.y + tri.c.y) / 3.0,
            (tri.a.z + tri.b.z + tri.c.z) / 3.0};
}

std::optional<TriangleHit> intersect_triangle(const Triangle& tri, std::size_t index,
                                               const Ray& ray, double t_min,
                                               double t_max) {
    const Vec3 edge1 = tri.b - tri.a;
    const Vec3 edge2 = tri.c - tri.a;
    const Vec3 p = cross(ray.direction, edge2);
    const double determinant = dot(edge1, p);
    const double scale = std::sqrt(squared_length(edge1) * squared_length(edge2) *
                                   squared_length(ray.direction));
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon() * scale;
    if (std::abs(determinant) <= tolerance) return std::nullopt;

    const double inverse = 1.0 / determinant;
    const Vec3 from_a = ray.origin - tri.a;
    const double u = dot(from_a, p) * inverse;
    if (u < 0.0 || u > 1.0) return std::nullopt;
    const Vec3 q = cross(from_a, edge1);
    const double v = dot(ray.direction, q) * inverse;
    if (v < 0.0 || u + v > 1.0) return std::nullopt;
    const double t = dot(edge2, q) * inverse;
    if (t < t_min || t > t_max) return std::nullopt;
    return TriangleHit{index, t, u, v, ray.origin + ray.direction * t};
}

} // namespace

 double& Vec3::operator[](std::size_t axis) {
    if (axis > 2) throw std::out_of_range("Vec3 axis must be 0, 1, or 2");
    return axis == 0 ? x : axis == 1 ? y : z;
}

double Vec3::operator[](std::size_t axis) const {
    if (axis > 2) throw std::out_of_range("Vec3 axis must be 0, 1, or 2");
    return axis == 0 ? x : axis == 1 ? y : z;
}

Vec3& Vec3::operator+=(const Vec3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
Vec3& Vec3::operator-=(const Vec3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
Vec3 operator+(Vec3 lhs, const Vec3& rhs) { return lhs += rhs; }
Vec3 operator-(Vec3 lhs, const Vec3& rhs) { return lhs -= rhs; }
Vec3 operator*(const Vec3& v, double scalar) { return {v.x * scalar, v.y * scalar, v.z * scalar}; }
Vec3 operator*(double scalar, const Vec3& v) { return v * scalar; }
Vec3 operator/(const Vec3& v, double scalar) { return {v.x / scalar, v.y / scalar, v.z / scalar}; }
double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
double squared_length(const Vec3& v) { return dot(v, v); }
double squared_distance(const Vec3& a, const Vec3& b) { return squared_length(a - b); }

AABB::AABB() = default;
AABB::AABB(const Vec3& min_corner, const Vec3& max_corner) {
    require_finite(min_corner, "AABB minimum");
    require_finite(max_corner, "AABB maximum");
    if (min_corner.x > max_corner.x || min_corner.y > max_corner.y || min_corner.z > max_corner.z)
        throw std::invalid_argument("AABB minimum must not exceed maximum");
    min_ = min_corner;
    max_ = max_corner;
    empty_ = false;
}

void AABB::expand(const Vec3& point) {
    require_finite(point, "AABB point");
    if (empty_) { min_ = max_ = point; empty_ = false; return; }
    min_.x = std::min(min_.x, point.x); min_.y = std::min(min_.y, point.y); min_.z = std::min(min_.z, point.z);
    max_.x = std::max(max_.x, point.x); max_.y = std::max(max_.y, point.y); max_.z = std::max(max_.z, point.z);
}

void AABB::expand(const AABB& box) {
    if (box.empty_) return;
    expand(box.min_);
    expand(box.max_);
}

Vec3 AABB::centroid() const {
    if (empty_) throw std::logic_error("empty AABB has no centroid");
    return {(min_.x + max_.x) * 0.5, (min_.y + max_.y) * 0.5, (min_.z + max_.z) * 0.5};
}

double AABB::squared_distance(const Vec3& point) const {
    if (empty_) return std::numeric_limits<double>::infinity();
    double result = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double delta = point[axis] < min_[axis] ? min_[axis] - point[axis]
                           : point[axis] > max_[axis] ? point[axis] - max_[axis] : 0.0;
        result += delta * delta;
    }
    return result;
}

std::optional<std::pair<double, double>> AABB::intersect_ray(
    const Ray& ray, double t_min, double t_max) const {
    require_ray(ray, t_min, t_max);
    if (empty_) return std::nullopt;
    double near_t = t_min;
    double far_t = t_max;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double direction = ray.direction[axis];
        if (direction == 0.0) {
            if (ray.origin[axis] < min_[axis] || ray.origin[axis] > max_[axis]) return std::nullopt;
            continue;
        }
        double t0 = (min_[axis] - ray.origin[axis]) / direction;
        double t1 = (max_[axis] - ray.origin[axis]) / direction;
        if (t0 > t1) std::swap(t0, t1);
        near_t = std::max(near_t, t0);
        far_t = std::min(far_t, t1);
        if (near_t > far_t) return std::nullopt;
    }
    return std::make_pair(near_t, far_t);
}

Orientation orient3d(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    require_finite(a, "orientation point"); require_finite(b, "orientation point");
    require_finite(c, "orientation point"); require_finite(d, "orientation point");
    const double adx = a.x-d.x, ady = a.y-d.y, adz = a.z-d.z;
    const double bdx = b.x-d.x, bdy = b.y-d.y, bdz = b.z-d.z;
    const double cdx = c.x-d.x, cdy = c.y-d.y, cdz = c.z-d.z;
    const double m1 = bdy*cdz - bdz*cdy;
    const double m2 = bdx*cdz - bdz*cdx;
    const double m3 = bdx*cdy - bdy*cdx;
    const double det = adx*m1 - ady*m2 + adz*m3;
    const double permanent = (std::abs(bdy*cdz) + std::abs(bdz*cdy))*std::abs(adx)
                           + (std::abs(bdx*cdz) + std::abs(bdz*cdx))*std::abs(ady)
                           + (std::abs(bdx*cdy) + std::abs(bdy*cdx))*std::abs(adz);
    const double bound = 16.0 * std::numeric_limits<double>::epsilon() * permanent;
    if (std::isfinite(det) && std::isfinite(bound) && std::abs(det) > bound)
        return det > 0.0 ? Orientation::positive : Orientation::negative;

    // Recompute differences and products with wider intermediates when available.
    // long double can equal double on some platforms and is not an exact expansion.
    const long double ladx = static_cast<long double>(a.x) - d.x;
    const long double lady = static_cast<long double>(a.y) - d.y;
    const long double ladz = static_cast<long double>(a.z) - d.z;
    const long double lbdx = static_cast<long double>(b.x) - d.x;
    const long double lbdy = static_cast<long double>(b.y) - d.y;
    const long double lbdz = static_cast<long double>(b.z) - d.z;
    const long double lcdx = static_cast<long double>(c.x) - d.x;
    const long double lcdy = static_cast<long double>(c.y) - d.y;
    const long double lcdz = static_cast<long double>(c.z) - d.z;
    const long double wide_det = ladx*(lbdy*lcdz-lbdz*lcdy)
                               - lady*(lbdx*lcdz-lbdz*lcdx)
                               + ladz*(lbdx*lcdy-lbdy*lcdx);
    if (wide_det > 0.0L) return Orientation::positive;
    if (wide_det < 0.0L) return Orientation::negative;
    return Orientation::coplanar;
}

KDTree3::KDTree3(std::vector<Vec3> points) : points_(std::move(points)) {
    for (const Vec3& point : points_) require_finite(point, "KD-tree point");
    order_.resize(points_.size());
    std::iota(order_.begin(), order_.end(), 0);
    nodes_.reserve(points_.size());
    if (!points_.empty()) build(0, points_.size(), 0);
}

std::size_t KDTree3::build(std::size_t begin, std::size_t end, std::size_t depth) {
    const std::size_t axis = depth % 3;
    const std::size_t mid = begin + (end - begin) / 2;
    std::nth_element(order_.begin() + static_cast<std::ptrdiff_t>(begin),
                     order_.begin() + static_cast<std::ptrdiff_t>(mid),
                     order_.begin() + static_cast<std::ptrdiff_t>(end),
                     [&](std::size_t lhs, std::size_t rhs) {
                         const double lv = points_[lhs][axis], rv = points_[rhs][axis];
                         return lv < rv || (lv == rv && lhs < rhs);
                     });
    Node node;
    node.point_index = order_[mid];
    node.axis = axis;
    node.bounds.expand(points_[node.point_index]);
    const std::size_t node_index = nodes_.size();
    nodes_.push_back(node);
    const std::size_t left = mid > begin ? build(begin, mid, depth + 1) : kNone;
    const std::size_t right = mid + 1 < end ? build(mid + 1, end, depth + 1) : kNone;
    nodes_[node_index].left = left;
    nodes_[node_index].right = right;
    if (left != kNone) nodes_[node_index].bounds.expand(nodes_[left].bounds);
    if (right != kNone) nodes_[node_index].bounds.expand(nodes_[right].bounds);
    return node_index;
}

std::optional<Neighbor> KDTree3::nearest(const Vec3& query) const {
    require_finite(query, "KD-tree query");
    if (points_.empty()) return std::nullopt;
    std::size_t best = kNone;
    double best_d2 = std::numeric_limits<double>::infinity();
    nearest_impl(0, query, best, best_d2);
    return Neighbor{best, points_[best], std::sqrt(best_d2)};
}

void KDTree3::nearest_impl(std::size_t node_index, const Vec3& query, std::size_t& best,
                           double& best_distance2) const {
    if (node_index == kNone || nodes_[node_index].bounds.squared_distance(query) > best_distance2) return;
    const Node& node = nodes_[node_index];
    const double d2 = geom3d::squared_distance(points_[node.point_index], query);
    if (d2 < best_distance2 || (d2 == best_distance2 && node.point_index < best)) {
        best = node.point_index;
        best_distance2 = d2;
    }
    const std::size_t first = query[node.axis] <= points_[node.point_index][node.axis] ? node.left : node.right;
    const std::size_t second = first == node.left ? node.right : node.left;
    if (first != kNone) nearest_impl(first, query, best, best_distance2);
    if (second != kNone) nearest_impl(second, query, best, best_distance2);
}

std::vector<Neighbor> KDTree3::radius(const Vec3& query, double radius_value) const {
    require_finite(query, "KD-tree query");
    if (!std::isfinite(radius_value) || radius_value < 0.0)
        throw std::invalid_argument("radius must be finite and nonnegative");
    std::vector<Neighbor> result;
    if (points_.empty()) return result;
    std::vector<std::size_t> found;
    radius_impl(0, query, radius_value * radius_value, found);
    result.reserve(found.size());
    for (std::size_t index : found)
        result.push_back({index, points_[index], std::sqrt(geom3d::squared_distance(points_[index], query))});
    std::sort(result.begin(), result.end(), [](const Neighbor& lhs, const Neighbor& rhs) {
        return lhs.distance < rhs.distance || (lhs.distance == rhs.distance && lhs.index < rhs.index);
    });
    return result;
}

void KDTree3::radius_impl(std::size_t node_index, const Vec3& query, double radius2,
                          std::vector<std::size_t>& found) const {
    if (node_index == kNone || nodes_[node_index].bounds.squared_distance(query) > radius2) return;
    const Node& node = nodes_[node_index];
    if (geom3d::squared_distance(points_[node.point_index], query) <= radius2) found.push_back(node.point_index);
    radius_impl(node.left, query, radius2, found);
    radius_impl(node.right, query, radius2, found);
}

TriangleBVH::TriangleBVH(std::vector<Triangle> triangles, std::size_t leaf_size)
    : triangles_(std::move(triangles)), leaf_size_(std::max<std::size_t>(1, leaf_size)) {
    for (const Triangle& tri : triangles_) {
        require_finite(tri.a, "triangle vertex"); require_finite(tri.b, "triangle vertex");
        require_finite(tri.c, "triangle vertex");
    }
    order_.resize(triangles_.size());
    std::iota(order_.begin(), order_.end(), 0);
    nodes_.reserve(triangles_.empty() ? 0 : 2 * triangles_.size());
    if (!triangles_.empty()) build(0, triangles_.size());
}

std::size_t TriangleBVH::build(std::size_t begin, std::size_t end) {
    Node node;
    node.begin = begin;
    node.end = end;
    AABB centroid_bounds;
    for (std::size_t i = begin; i < end; ++i) {
        node.bounds.expand(triangle_bounds(triangles_[order_[i]]));
        centroid_bounds.expand(triangle_centroid(triangles_[order_[i]]));
    }
    const std::size_t node_index = nodes_.size();
    nodes_.push_back(node);
    if (end - begin <= leaf_size_) return node_index;

    const Vec3 extent = centroid_bounds.max() - centroid_bounds.min();
    std::size_t axis = 0;
    if (extent.y > extent.x) axis = 1;
    if (extent.z > extent[axis]) axis = 2;
    const std::size_t mid = begin + (end - begin) / 2;
    std::nth_element(order_.begin() + static_cast<std::ptrdiff_t>(begin),
                     order_.begin() + static_cast<std::ptrdiff_t>(mid),
                     order_.begin() + static_cast<std::ptrdiff_t>(end),
                     [&](std::size_t lhs, std::size_t rhs) {
                         const double lv = triangle_centroid(triangles_[lhs])[axis];
                         const double rv = triangle_centroid(triangles_[rhs])[axis];
                         return lv < rv || (lv == rv && lhs < rhs);
                     });
    const std::size_t left = build(begin, mid);
    const std::size_t right = build(mid, end);
    nodes_[node_index].left = left;
    nodes_[node_index].right = right;
    nodes_[node_index].leaf = false;
    return node_index;
}

std::optional<TriangleHit> TriangleBVH::closest_hit(const Ray& ray, double t_min, double t_max) const {
    require_ray(ray, t_min, t_max);
    if (triangles_.empty()) return std::nullopt;
    std::optional<TriangleHit> best;
    double best_t = t_max;
    hit_impl(0, ray, t_min, best_t, best);
    return best;
}

void TriangleBVH::hit_impl(std::size_t node_index, const Ray& ray, double t_min,
                           double& best_t, std::optional<TriangleHit>& best) const {
    const Node& node = nodes_[node_index];
    const auto interval = node.bounds.intersect_ray(ray, t_min, best_t);
    if (!interval) return;
    if (node.leaf) {
        for (std::size_t i = node.begin; i < node.end; ++i) {
            const std::size_t index = order_[i];
            auto hit = intersect_triangle(triangles_[index], index, ray, t_min, best_t);
            if (hit && (!best || hit->t < best_t || (hit->t == best_t && hit->triangle_index < best->triangle_index))) {
                best_t = hit->t;
                best = hit;
            }
        }
        return;
    }
    const auto left_interval = nodes_[node.left].bounds.intersect_ray(ray, t_min, best_t);
    const auto right_interval = nodes_[node.right].bounds.intersect_ray(ray, t_min, best_t);
    if (left_interval && right_interval) {
        const std::size_t first = left_interval->first <= right_interval->first ? node.left : node.right;
        const std::size_t second = first == node.left ? node.right : node.left;
        hit_impl(first, ray, t_min, best_t, best);
        hit_impl(second, ray, t_min, best_t, best);
    } else if (left_interval) hit_impl(node.left, ray, t_min, best_t, best);
    else if (right_interval) hit_impl(node.right, ray, t_min, best_t, best);
}

} // namespace geom3d
