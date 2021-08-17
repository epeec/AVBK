#ifndef LINE_DETECTION_HPP

#define LINE_DETECTION_HPP

#include <vector>

/**
 * @brief Detect the lines for re-numbering
 * @tparam T_Integer
 * @param start_from_patch_id
 * @param edge_vertices
 * @param x
 * @param y
 * @param boundary_vertices
 * @param boundary_offset
 * @param lines
 */
template<typename T_Integer>
void detect_lines_2d( const T_Integer start_from_patch_id,
                      const std::vector<T_Integer> &edge_vertices,
                      const std::vector<double> &x,
                      const std::vector<double> &y,
                      const std::vector<T_Integer> &boundary_vertices,
                      const std::vector<T_Integer> &boundary_offset,
                      std::map<int, std::list<T_Integer>> &lines ) {
  std::vector<bool> is_edge_occupied(edge_vertices.size() / 2, false);
  // Create the vertex-to-edge graph
  std::vector<std::set<int>> vertex_edges;
  // Start from the nodes along the boundary and
  // march in the `dir` coordiante to obtain the
  // lines
  auto num_bnodes = boundary_offset[start_from_patch_id + 1]
                    - boundary_offset[start_from_patch_id];
  // Loop over each nodes of this patch and keep adding the
  // edges along the direction specified by the user
  for( int i=0; i < num_bnodes; ++i ) {

  }
}

#endif