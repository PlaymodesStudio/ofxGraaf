#pragma once

#include <graaflib/properties/vertex_properties.h>

namespace graaf {

template <typename VERTEX_T, typename EDGE_T, graph_type GRAPH_TYPE_V>
tree<VERTEX_T, EDGE_T, GRAPH_TYPE_V>::tree(graph_t&& graph)
    : graph_{std::move(graph)} {
  std::size_t highest_degree{0};
  for (const auto& [vertex_id, _] : graph_.get_vertices()) {
    const auto degree{properties::vertex_degree(graph_, vertex_id)};

    if (degree > highest_degree) {
      highest_degree = degree;
      root_ = vertex_id;
    }

    if (degree == 1) {
      leaves_.insert(vertex_id);
    }
  }
}

template <typename VERTEX_T, typename EDGE_T, graph_type GRAPH_TYPE_V>
tree<VERTEX_T, EDGE_T, GRAPH_TYPE_V>::edges_t
tree<VERTEX_T, EDGE_T, GRAPH_TYPE_V>::get_edges() const {
  edges_t result{};

  for (const auto& [id, edge] : graph_.get_edges()) {
    result.push_back(id);
  }

  return result;
}

template <typename V, typename E, graph_type T>
tree<V, E, T> tree_from_graph(const graph<V, E, T>& input_graph,
                              const std::vector<edge_id_t>& tree_edges) {
  graph<V, E, T> result{};

  const auto add_vertex_if_missing{[&](vertex_id_t id) {
    if (!result.has_vertex(id)) {
      result.add_vertex(input_graph.get_vertex(id), id);
    }
  }};

  for (const auto& [lhs_id, rhs_id] : tree_edges) {
    add_vertex_if_missing(lhs_id);
    add_vertex_if_missing(rhs_id);

    result.add_edge(lhs_id, rhs_id, input_graph.get_edge(lhs_id, rhs_id));
  }

  return tree{std::move(result)};
}

template <typename VERTEX_T, typename EDGE_T, graph_type GRAPH_TYPE_V>
std::vector<vertex_id_t> tree<VERTEX_T, EDGE_T, GRAPH_TYPE_V>::get_vertex_ids() const {
    std::vector<vertex_id_t> vertex_ids;
    vertex_ids.reserve(graph_.vertex_count());  // Reserve space for efficiency

    for (const auto& [vertex_id, vertex] : graph_.get_vertices()) {
        vertex_ids.push_back(vertex_id);
    }

    return vertex_ids;
}

}  // namespace graaf
