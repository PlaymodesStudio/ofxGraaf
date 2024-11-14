#pragma once

#include <graaflib/graph.h>
#include <graaflib/types.h>

#include <vector>

namespace graaf {

template <typename VERTEX_T, typename EDGE_T, graph_type GRAPH_TYPE_V>
class tree {
  using graph_t = graph<VERTEX_T, EDGE_T, GRAPH_TYPE_V>;
  using edges_t = std::vector<edge_id_t>;

 public:
  tree(graph_t&& graph);

  [[nodiscard]] vertex_id_t get_root() const { return root_; }

  [[nodiscard]] const graph_t::vertices_t& get_leaves() const { return leaves_; }

  [[nodiscard]] edges_t get_edges() const;

  [[nodiscard]] bool empty() const { return graph_.vertex_count() == 0; }

  [[nodiscard]] graph_t::vertex_t& get_vertex(vertex_id_t vertex_id) {
    return graph_.get_vertex(vertex_id);
  }
  [[nodiscard]] const graph_t::vertex_t& get_vertex(
      vertex_id_t vertex_id) const {
    return graph_.get_vertex(vertex_id);
  }

  [[nodiscard]] graph_t::edge_t& get_edge(vertex_id_t vertex_id_lhs,
                                          vertex_id_t vertex_id_rhs) {
    return graph_.get_edge(vertex_id_lhs, vertex_id_rhs);
  }
  [[nodiscard]] const graph_t::edge_t& get_edge(
      vertex_id_t vertex_id_lhs, vertex_id_t vertex_id_rhs) const {
    return graph_.get_edge(vertex_id_lhs, vertex_id_rhs);
  }
  [[nodiscard]] graph_t::edge_t& get_edge(edge_id_t edge_id) {
    return graph_.get_edge(edge_id);
  }
  [[nodiscard]] const graph_t::edge_t& get_edge(edge_id_t edge_id) const {
    return graph_.get_edge(edge_id);
  }

  [[nodiscard]] graph_t::vertices_t get_neighbors(vertex_id_t vertex_id) const {
    return graph_.get_neighbors(vertex_id);
  }

  std::vector<vertex_id_t> get_vertex_ids() const;
    
 private:
  graph_t graph_{};
  vertex_id_t root_{};
  graph_t::vertices_t leaves_{};
};

template <typename V, typename E, graph_type T>
[[nodiscard]] tree<V, E, T> tree_from_graph(
    const graph<V, E, T>& input_graph,
    const std::vector<edge_id_t>& tree_edges);

}  // namespace graaf

#include "tree.tpp"
