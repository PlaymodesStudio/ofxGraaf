// delaunay_triangulation.h
#pragma once

#include <graaflib/graph.h>
#include <array>
#include <vector>

namespace graaf::algorithm {

template <typename T>
concept Point2D = requires(T p) {
    { p.x } -> std::convertible_to<double>;
    { p.y } -> std::convertible_to<double>;
};

namespace detail {
    template <typename V>
        requires Point2D<V>
    struct Triangle {
        std::array<vertex_id_t, 3> vertices;
        std::array<V, 3> points;
        
        [[nodiscard]] bool contains_vertex(vertex_id_t v) const {
            return std::find(vertices.begin(), vertices.end(), v) != vertices.end();
        }
        
        [[nodiscard]] bool circumcircle_contains(const V& p) const;
        
        [[nodiscard]] bool equals(const Triangle& other) const {
            return vertices == other.vertices &&
                   points[0].x == other.points[0].x &&
                   points[0].y == other.points[0].y &&
                   points[1].x == other.points[1].x &&
                   points[1].y == other.points[1].y &&
                   points[2].x == other.points[2].x &&
                   points[2].y == other.points[2].y;
        }
    };

    template <typename V>
        requires Point2D<V>
    struct Edge {
        vertex_id_t v1;
        vertex_id_t v2;
        
        bool operator==(const Edge& other) const {
            return (v1 == other.v1 && v2 == other.v2) ||
                   (v1 == other.v2 && v2 == other.v1);
        }
    };

    template <typename V>
        requires Point2D<V>
    [[nodiscard]] double cross_product(const V& a, const V& b, const V& c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    template <typename V>
        requires Point2D<V>
    [[nodiscard]] Triangle<V> create_super_triangle(const std::vector<V>& points);
    
} // namespace detail

template <typename V, typename E = int>
    requires Point2D<V>
[[nodiscard]] undirected_graph<V, E> delaunay_triangulation(
    const std::vector<V>& points,
    std::function<E(const V&, const V&)> edge_weight_calculator =
        [](const V&, const V&) { return E{1}; });

} // namespace graaf::algorithm

// Implementation in same file to avoid linker issues
namespace graaf::algorithm {
namespace detail {

template <typename V>
    requires Point2D<V>
bool Triangle<V>::circumcircle_contains(const V& p) const {
    const auto& a = points[0];
    const auto& b = points[1];
    const auto& c = points[2];
    
    const double ax = static_cast<double>(a.x) - static_cast<double>(p.x);
    const double ay = static_cast<double>(a.y) - static_cast<double>(p.y);
    const double bx = static_cast<double>(b.x) - static_cast<double>(p.x);
    const double by = static_cast<double>(b.y) - static_cast<double>(p.y);
    const double cx = static_cast<double>(c.x) - static_cast<double>(p.x);
    const double cy = static_cast<double>(c.y) - static_cast<double>(p.y);
    
    const double det = (ax * ax + ay * ay) * (bx * cy - cx * by) -
                      (bx * bx + by * by) * (ax * cy - cx * ay) +
                      (cx * cx + cy * cy) * (ax * by - bx * ay);
                      
    return det > 0;
}

template <typename V>
    requires Point2D<V>
Triangle<V> create_super_triangle(const std::vector<V>& points) {
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    
    for (const auto& p : points) {
        min_x = std::min(min_x, static_cast<double>(p.x));
        min_y = std::min(min_y, static_cast<double>(p.y));
        max_x = std::max(max_x, static_cast<double>(p.x));
        max_y = std::max(max_y, static_cast<double>(p.y));
    }
    
    const double dx = (max_x - min_x) * 10;
    const double dy = (max_y - min_y) * 10;
    
    V p1, p2, p3;
    p1.x = static_cast<typename std::remove_reference<decltype(p1.x)>::type>(min_x - dx);
    p1.y = static_cast<typename std::remove_reference<decltype(p1.y)>::type>(min_y - dy);
    p2.x = static_cast<typename std::remove_reference<decltype(p2.x)>::type>(max_x + dx);
    p2.y = static_cast<typename std::remove_reference<decltype(p2.y)>::type>(min_y - dy);
    p3.x = static_cast<typename std::remove_reference<decltype(p3.x)>::type>((min_x + max_x) / 2);
    p3.y = static_cast<typename std::remove_reference<decltype(p3.y)>::type>(max_y + dy);

    return Triangle<V>{
        {0, 1, 2},
        {p1, p2, p3}
    };
}

} // namespace detail

template <typename V, typename E>
requires Point2D<V>
undirected_graph<V, E> delaunay_triangulation(
    const std::vector<V>& points,
    std::function<E(const V&, const V&)> edge_weight_calculator) {

    if (points.size() < 3) {
        throw std::invalid_argument("Need at least 3 points for triangulation");
    }

    undirected_graph<V, E> result;
    std::vector<detail::Triangle<V>> triangles;

    // Step 1: Create super triangle
    auto super_triangle = detail::create_super_triangle(points);

    // Add super triangle vertices to the graph
    for (size_t i = 0; i < 3; ++i) {
        result.add_vertex(super_triangle.points[i], i);
    }
    
    triangles.push_back(super_triangle);

    // Step 2: Add each input point and perform triangulation
    vertex_id_t next_id = 3; // ID for original points, starting after super-triangle vertices
    for (const auto& point : points) {
        std::vector<detail::Triangle<V>> bad_triangles;
        std::vector<detail::Edge<V>> polygon;

        // Find all triangles where the point lies within the circumcircle
        for (const auto& triangle : triangles) {
            if (triangle.circumcircle_contains(point)) {
                bad_triangles.push_back(triangle);
            }
        }

        // Find the boundary of the polygon hole created by bad triangles
        for (const auto& triangle : bad_triangles) {
            for (size_t i = 0; i < 3; ++i) {
                detail::Edge<V> edge{triangle.vertices[i], triangle.vertices[(i + 1) % 3]};

                bool is_shared = false;
                for (const auto& other : bad_triangles) {
                    if (&triangle != &other &&
                        other.contains_vertex(edge.v1) &&
                        other.contains_vertex(edge.v2)) {
                        is_shared = true;
                        break;
                    }
                }

                if (!is_shared) {
                    polygon.push_back(edge);
                }
            }
        }

        // Remove bad triangles
        triangles.erase(
            std::remove_if(triangles.begin(), triangles.end(),
                [&bad_triangles](const detail::Triangle<V>& t) {
                    return std::find_if(bad_triangles.begin(), bad_triangles.end(),
                        [&t](const detail::Triangle<V>& bad) {
                            return t.equals(bad);
                        }) != bad_triangles.end();
                }), triangles.end());

        // Add the new point to the graph
        vertex_id_t new_vertex_id = result.add_vertex(point, next_id++);

        // Re-triangulate the polygon hole with the new point
        for (const auto& edge : polygon) {
            detail::Triangle<V> new_triangle{
                {edge.v1, edge.v2, new_vertex_id},
                {
                    result.get_vertex(edge.v1),
                    result.get_vertex(edge.v2),
                    point
                }
            };
            triangles.push_back(new_triangle);
        }
    }

    // Step 3: Remove triangles that involve super-triangle vertices
    triangles.erase(
        std::remove_if(triangles.begin(), triangles.end(),
            [](const detail::Triangle<V>& t) {
                return t.contains_vertex(0) || t.contains_vertex(1) || t.contains_vertex(2);
            }), triangles.end());

    // Remove super-triangle vertices from the graph
    result.remove_vertex(0);
    result.remove_vertex(1);
    result.remove_vertex(2);

    // Step 4: Add edges from remaining triangles to the result graph
    for (const auto& triangle : triangles) {
        for (size_t i = 0; i < 3; ++i) {
            const auto v1 = triangle.vertices[i];
            const auto v2 = triangle.vertices[(i + 1) % 3];
            if (!result.has_edge(v1, v2)) {
                result.add_edge(v1, v2,
                    edge_weight_calculator(
                        result.get_vertex(v1),
                        result.get_vertex(v2)));
            }
        }
    }

    return result;
}
} // namespace graaf::algorithm
