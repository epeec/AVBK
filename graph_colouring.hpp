/** @file graph_colouring
 *  @author Pavanakumar Mohanamuraly
 *  @date 2/22/20
 *  @brief File part of TreePart but included here for quick and
 *         dirty colouring interface
 *  @copyright CERFACS
 */
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifndef TREEPART_GRAPH_COLOURING_HPP
#define TREEPART_GRAPH_COLOURING_HPP
#define NOT_ASSIGNED -1

/**
 *
 * @tparam IntegerT
 * @param in
 * @return
 */
template <typename IntegerT> static inline IntegerT down_by_one(IntegerT in) {
  return in - 1;
}

/**
 *
 * @tparam IntegerT
 * @param i_LeftVertexCount
 * @return
 */
template <typename IntegerT>
std::vector<IntegerT> RowNaturalOrdering(const IntegerT i_LeftVertexCount) {

  std::vector<IntegerT> m_vi_OrderedVertices(i_LeftVertexCount);
  for (IntegerT i = 0; i < i_LeftVertexCount; i++)
    m_vi_OrderedVertices[i] = i;
  return m_vi_OrderedVertices;
}

/**
 *
 * @tparam IntegerT
 * @param i_LeftVertexCount
 * @param i_RightVertexCount
 * @return
 */
template <typename IntegerT>
std::vector<IntegerT> ColumnNaturalOrdering(const IntegerT i_LeftVertexCount,
                                            const IntegerT i_RightVertexCount) {

  std::vector<IntegerT> m_vi_OrderedVertices;
  m_vi_OrderedVertices.reserve(i_RightVertexCount);
  for (int i = 0; i < i_RightVertexCount; i++)
    m_vi_OrderedVertices.push_back(i + i_LeftVertexCount);
  return m_vi_OrderedVertices;
}

/**
 *
 * @tparam IntegerT
 * @param m_vi_LeftVertices
 * @param m_vi_RightVertices
 * @param m_vi_Edges
 * @return
 */
template <typename IntegerT>
std::vector<IntegerT>
ColumnLargestFirstOrdering(const std::vector<IntegerT> &m_vi_LeftVertices,
                           const std::vector<IntegerT> &m_vi_RightVertices,
                           const std::vector<IntegerT> &m_vi_Edges) {

  IntegerT i_LeftVertexCount = down_by_one(m_vi_LeftVertices.size());
  IntegerT i_RightVertexCount = down_by_one(m_vi_RightVertices.size());
  IntegerT i_VertexCount = i_RightVertexCount;
  IntegerT m_i_MaximumVertexDegree = 0;
  IntegerT m_i_MinimumVertexDegree = i_VertexCount;
  std::vector<std::vector<IntegerT>> vvi_GroupedVertexDegree(i_VertexCount);
  std::vector<IntegerT> vi_Visited(i_VertexCount, NOT_ASSIGNED);

  IntegerT i_DegreeCount = 0;
  for (IntegerT i = 0; i < i_VertexCount; ++i) {
    // reset the degree count
    i_DegreeCount = 0;
    // let's loop from mvi_RightVertices[i] to mvi_RightVertices[i+1] for the
    // i'th column
    for (auto j = m_vi_RightVertices[i]; j < m_vi_RightVertices[i + 1]; ++j) {
      auto i_Current = m_vi_Edges[j];
      for (auto k = m_vi_LeftVertices[i_Current];
           k < m_vi_LeftVertices[i_Current + 1]; ++k) {
        if (m_vi_Edges[k] != i && vi_Visited[m_vi_Edges[k]] != i) {
          ++i_DegreeCount;
          vi_Visited[m_vi_Edges[k]] = i;
        }
      }
    }
    vvi_GroupedVertexDegree[i_DegreeCount].push_back(i);
    if (m_i_MaximumVertexDegree < i_DegreeCount) {
      m_i_MaximumVertexDegree = i_DegreeCount;
    } else if (m_i_MinimumVertexDegree > i_DegreeCount) {
      m_i_MinimumVertexDegree = i_DegreeCount;
    }
  }

  if (i_VertexCount < 2)
    m_i_MinimumVertexDegree = i_DegreeCount;

  // take the bucket and place it in the vertex order
  std::vector<IntegerT> m_vi_OrderedVertices;
  for (auto i = m_i_MaximumVertexDegree; i >= m_i_MinimumVertexDegree; --i) {
    // j = size of the bucket
    // get the size of the i-th bucket
    auto j = vvi_GroupedVertexDegree[i].size();
    // place it into vertex ordering
    for (auto k = 0; k < j; ++k) {
      m_vi_OrderedVertices.push_back(vvi_GroupedVertexDegree[i][k] +
                                     i_LeftVertexCount);
    }
  }
  return m_vi_OrderedVertices;
}

/**
 *
 * @tparam IntegerT
 * @param m_vi_LeftVertices
 * @param m_vi_RightVertices
 * @param m_vi_Edges
 * @param m_vi_LeftVertexColors
 * @return
 */
template <typename IntegerT>
bool CheckPartialDistanceTwoRowColoring(
    const std::vector<IntegerT> &m_vi_LeftVertices,
    const std::vector<IntegerT> &m_vi_RightVertices,
    const std::vector<IntegerT> &m_vi_Edges,
    const std::vector<IntegerT> &m_vi_LeftVertexColors) {

  for (IntegerT i = 0; i < down_by_one(m_vi_LeftVertices.size()); i++) {
    // for each of left vertices, find its D1 neighbour (right vertices)
    for (auto j = m_vi_LeftVertices[i]; j < m_vi_LeftVertices[i + 1]; j++) {
      for (auto k = m_vi_RightVertices[m_vi_Edges[j]];
           k < m_vi_RightVertices[m_vi_Edges[j] + 1]; k++) {
        // for each of the right vertices, find its D1 neighbour (left vertices
        // exclude the original left)
        if (m_vi_Edges[k] == i)
          continue;
        if (m_vi_LeftVertexColors[m_vi_Edges[k]] == m_vi_LeftVertexColors[i]) {
          std::cout << "Left vertices " << i + 1 << " and " << m_vi_Edges[k] + 1
                    << " (connected by right vertex " << m_vi_Edges[j] + 1
                    << ") have the same color (" << m_vi_LeftVertexColors[i]
                    << ")" << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}

/**
 *
 * @tparam IntegerT
 * @param x
 * @return
 */
template <typename IntegerT> IntegerT foobar(const IntegerT x) { return x; }

/**
 *
 * @tparam IntegerT
 * @param m_vi_LeftVertices
 * @param m_vi_RightVertices
 * @param m_vi_Edges
 * @return
 */
template <typename IntegerT>
std::vector<IntegerT> ParallelPartialDistanceTwoRowColoring(
    const std::vector<IntegerT> &m_vi_LeftVertices,
    const std::vector<IntegerT> &m_vi_RightVertices,
    const std::vector<IntegerT> &m_vi_Edges) {

  IntegerT i_LeftVertexCount = down_by_one(m_vi_LeftVertices.size());
  std::vector<IntegerT> vi_forbiddenColors(i_LeftVertexCount, NOT_ASSIGNED);
  std::vector<IntegerT> vi_VerticesToBeColored(i_LeftVertexCount);
  IntegerT i_NumOfVerticesToBeColored = vi_VerticesToBeColored.size();
  std::vector<IntegerT> m_vi_LeftVertexColors(i_LeftVertexCount, NOT_ASSIGNED);
  std::vector<IntegerT> vi_verticesNeedNewColor;

  vi_verticesNeedNewColor.reserve(i_LeftVertexCount);
  auto m_i_LeftVertexColorCount = 0;
  // Algo 4 - Line 3: while U != 0 ; do
  while (i_NumOfVerticesToBeColored != 0) {
    // Phase 1: tentative coloring
    // Algo 4 - Line 4: for each right vertex v in U (in parallel) do
#pragma omp parallel for default(none) schedule(dynamic)                       \
    shared(i_NumOfVerticesToBeColored, vi_VerticesToBeColored)                 \
        firstprivate(vi_forbiddenColors)
    for (IntegerT i = 0; i < i_NumOfVerticesToBeColored; i++) {
      auto myVertex = vi_VerticesToBeColored[i];
      // Algo 4 - Line 5: for each left vertex w in adj (v) do
      for (auto w = m_vi_LeftVertices[myVertex];
           w < m_vi_LeftVertices[myVertex + 1]; w++) {
        // Algo 4 - Line 6: mark color [w] as forbidden to vertex v. NOTE: !!!
        // Not needed Algo 4 - Line 7: for each right vertex x in adj (w) and x
        // != v do
        for (auto x = m_vi_RightVertices[m_vi_Edges[w]];
             x < m_vi_RightVertices[m_vi_Edges[w] + 1]; x++) {
          // Algo 4 - Line 8: mark color [x] as forbidden to vertex v
          if (m_vi_LeftVertexColors[m_vi_Edges[x]] != NOT_ASSIGNED) {
            // !!! each thread should have its private vi_forbiddenColors[]
            // vector to ensure we don't override each other
            vi_forbiddenColors[m_vi_LeftVertexColors[m_vi_Edges[x]]] = myVertex;
          }
        }
      }
      // Algo 4 - Line 9: Pick a permissible color c for vertex v using some
      // strategy
      IntegerT i_cadidateColor = 0;
      // First fit
      while (vi_forbiddenColors[i_cadidateColor] == myVertex)
        i_cadidateColor++;
      m_vi_LeftVertexColors[myVertex] = i_cadidateColor;
      if (m_i_LeftVertexColorCount < i_cadidateColor)
        m_i_LeftVertexColorCount = i_cadidateColor;
    }
    // Algo 4 - Line 10: R.clear()   ; R denotes the set of vertices to be
    // recolored
    vi_verticesNeedNewColor.clear();
    // Phase 2: conflict detection. For each vertex v in U, check and see if v
    // need to be recolored Algo 4 - Line 11: for each vertex v in U (in
    // parallel) do
#pragma omp parallel for default(none) schedule(dynamic)                       \
    shared(i_NumOfVerticesToBeColored, vi_VerticesToBeColored,                 \
           vi_verticesNeedNewColor)
    for (IntegerT i = 0; i < i_NumOfVerticesToBeColored; i++) {
      // Algo 4 - Line 12: cont  <- true ; cont is used to break from the outer
      // loop below
      auto continueColouring = true;
      auto myVertex = vi_VerticesToBeColored[i];
      // Algo 4 - Line 13: for each vertex w in adj (v) and cont = true do
      for (auto w = m_vi_LeftVertices[myVertex];
           (w < m_vi_LeftVertices[myVertex + 1]) && continueColouring; w++) {
        // Algo 4 - Line 14: if color [v] = color [w] and f (v) > f (w) then .
        // NOTE: !!! Not needed Algo 4 - Line 15: add [v] to R ; break . NOTE:
        // !!! Not needed Algo 4 - Line 16: for each vertex x in adj (w) and v
        // != x do
        for (auto x = m_vi_RightVertices[m_vi_Edges[w]];
             x < m_vi_RightVertices[m_vi_Edges[w] + 1]; x++) {
          // Algo 4 - Line 17: if color [v] = color [x] and f (v) > f (x) then
          if (m_vi_LeftVertexColors[m_vi_Edges[x]] ==
                  m_vi_LeftVertexColors[myVertex] &&
              foobar(myVertex) > foobar(m_vi_Edges[x])) {
            // Algo 4 - Line 18: add [v] to R ; cont <- false; break
#pragma omp critical
            { vi_verticesNeedNewColor.push_back(myVertex); }
            continueColouring = false;
            break;
          }
        }
      }
    }
    // Algo 4 - Line 19: U <- R , i.e., vi_VerticesToBeColored <-
    // vi_verticesNeedNewColor
    vi_VerticesToBeColored.swap(vi_verticesNeedNewColor);
  }
  // Note that m_i_LeftVertexColorCount has not been updated yet
  return m_vi_LeftVertexColors;
}

/**
 * @brief Do the partial-distance-2 row colouring of the Bipartite graph
 * @tparam IntegerT
 * @param m_vi_LeftVertices
 * @param m_vi_RightVertices
 * @param m_vi_Edges
 * @return
 */
template <typename IntegerT>
std::vector<IntegerT>
PartialDistanceTwoRowColoring(const std::vector<IntegerT> &m_vi_LeftVertices,
                              const std::vector<IntegerT> &m_vi_RightVertices,
                              const std::vector<IntegerT> &m_vi_Edges) {

  auto i_LeftVertexCount = down_by_one(m_vi_LeftVertices.size());
  std::vector<IntegerT> m_vi_LeftVertexColors(i_LeftVertexCount, NOT_ASSIGNED);
  std::vector<IntegerT> vi_forbiddenColors(i_LeftVertexCount, NOT_ASSIGNED);
  auto m_vi_OrderedVertices = std::move(RowNaturalOrdering(i_LeftVertexCount));

  auto m_i_LeftVertexColorCount = 0;
  for (auto i = 0; i < i_LeftVertexCount; ++i) {
    auto i_CurrentVertex = m_vi_OrderedVertices[i];
    for (auto w = m_vi_LeftVertices[i_CurrentVertex];
         w < m_vi_LeftVertices[i_CurrentVertex + 1]; ++w) {
      for (auto x = m_vi_RightVertices[m_vi_Edges[w]];
           x < m_vi_RightVertices[m_vi_Edges[w] + 1]; ++x) {
        if (m_vi_LeftVertexColors[m_vi_Edges[x]] != NOT_ASSIGNED) {
          vi_forbiddenColors[m_vi_LeftVertexColors[m_vi_Edges[x]]] =
              i_CurrentVertex;
        }
      }
    }
    // do color[vi] <-min {c>0:forbiddenColors[c]=/=vi
    for (auto c = 0; c < i_LeftVertexCount; ++c) {
      if (vi_forbiddenColors[c] != i_CurrentVertex) {
        m_vi_LeftVertexColors[i_CurrentVertex] = c;
        if (m_i_LeftVertexColorCount < c) {
          m_i_LeftVertexColorCount = c;
        }
        break;
      }
    }
  }
  return m_vi_LeftVertexColors;
}

/**
 * @brief Do the partial-distance-2 column colouring of the Bipartite graph
 * @tparam IntegerT
 * @param m_vi_LeftVertices
 * @param m_vi_RightVertices
 * @param m_vi_Edges
 * @return
 */
template <typename IntegerT>
std::vector<IntegerT> PartialDistanceTwoColumnColoring(
    const std::vector<IntegerT> &m_vi_LeftVertices,
    const std::vector<IntegerT> &m_vi_RightVertices,
    const std::vector<IntegerT> &m_vi_Edges) {

  IntegerT i_LeftVertexCount = down_by_one(m_vi_LeftVertices.size());
  IntegerT i_RightVertexCount = down_by_one(m_vi_RightVertices.size());
  std::vector<IntegerT> m_vi_RightVertexColors(i_RightVertexCount,
                                               NOT_ASSIGNED);
  std::vector<IntegerT> vi_forbiddenColors(i_RightVertexCount, NOT_ASSIGNED);
  // auto m_vi_OrderedVertices =
  // std::move(ColumnNaturalOrdering(i_LeftVertexCount, i_RightVertexCount));
  auto m_vi_OrderedVertices = std::move(ColumnLargestFirstOrdering(
      m_vi_LeftVertices, m_vi_RightVertices, m_vi_Edges));

  for (auto i = 0; i < i_RightVertexCount; ++i) {
    auto i_CurrentVertex = m_vi_OrderedVertices[i] - i_LeftVertexCount;
    for (auto w = m_vi_RightVertices[i_CurrentVertex];
         w < m_vi_RightVertices[i_CurrentVertex + 1]; ++w) {
      for (auto x = m_vi_LeftVertices[m_vi_Edges[w]];
           x < m_vi_LeftVertices[m_vi_Edges[w] + 1]; ++x) {
        if (m_vi_RightVertexColors[m_vi_Edges[x]] != NOT_ASSIGNED) {
          vi_forbiddenColors[m_vi_RightVertexColors[m_vi_Edges[x]]] =
              i_CurrentVertex;
        }
      }
    }
    // do color[vi] <-min {c>0:forbiddenColors[c]=/=vi
    for (auto c = 0; c < i_RightVertexCount; ++c) {
      if (vi_forbiddenColors[c] != i_CurrentVertex) {
        m_vi_RightVertexColors[i_CurrentVertex] = c;
        break;
      }
    }
  }
  return m_vi_RightVertexColors;
}

/**
 * @brief Check the distance-2 colouring of the Bipartite graph
 * @tparam IntegerT
 * @param m_vi_LeftVertices
 * @param m_vi_RightVertices
 * @param m_vi_Edges
 * @param m_vi_RightVertexColors
 * @return
 */
template <typename IntegerT>
bool CheckPartialDistanceTwoColumnColoring(
    const std::vector<IntegerT> &m_vi_LeftVertices,
    const std::vector<IntegerT> &m_vi_RightVertices,
    const std::vector<IntegerT> &m_vi_Edges,
    const std::vector<IntegerT> &m_vi_RightVertexColors) {

  // for each of right vertices, find its D1 neighbour (left vertices)
  for (IntegerT i = 0; i < down_by_one(m_vi_RightVertices.size()); i++) {
    for (auto j = m_vi_RightVertices[i]; j < m_vi_RightVertices[i + 1]; j++) {
      // for each of the left vertices, find its D1 neighbour (right vertices
      // exclude the original right)
      for (auto k = m_vi_LeftVertices[m_vi_Edges[j]];
           k < m_vi_LeftVertices[m_vi_Edges[j] + 1]; k++) {
        if (m_vi_Edges[k] == i)
          continue;
        if (m_vi_RightVertexColors[m_vi_Edges[k]] ==
            m_vi_RightVertexColors[i]) {
          std::cout << "Right vertices " << i + 1 << " and "
                    << m_vi_Edges[k] + 1 << " (connected by left vertex "
                    << m_vi_Edges[j] + 1 << ") have the same color ("
                    << m_vi_RightVertexColors[i] << ")\n";
          return false;
        }
      }
    }
  }
  return true;
}

/**
 * @brief Convert a CSR graph to Bipartite graph
 * @tparam IntegerT
 * @param ip_RowIndex
 * @param i_RowCount
 * @param i_ColumnCount
 * @param ip_ColumnIndex
 * @param i_Start
 * @return
 */
template <typename IntegerT>
std::tuple<std::vector<IntegerT>, std::vector<IntegerT>, std::vector<IntegerT>>
BuildBPGraphFromCSRFormat(const IntegerT *const ip_RowIndex,
                          const IntegerT i_RowCount,
                          const IntegerT i_ColumnCount,
                          const IntegerT *const ip_ColumnIndex,
                          const IntegerT i_Start = 0) {

  std::map<IntegerT, std::vector<IntegerT>> colList;
  std::vector<IntegerT> m_vi_LeftVertices;
  m_vi_LeftVertices.reserve(i_RowCount + 1);
  std::vector<IntegerT> m_vi_RightVertices;
  m_vi_RightVertices.reserve(i_RowCount + 1);
  std::vector<IntegerT> m_vi_Edges;
  m_vi_Edges.reserve(2 * ip_RowIndex[i_RowCount]);

  m_vi_LeftVertices.push_back(0);
  for (IntegerT i = 0; i < i_RowCount; i++) {
    for (auto j = ip_RowIndex[i] - i_Start; j < ip_RowIndex[i + 1] - i_Start;
         j++) {
      m_vi_Edges.push_back(ip_ColumnIndex[j] - i_Start);
      colList[ip_ColumnIndex[j] - i_Start].push_back(i);
    }
    m_vi_LeftVertices.push_back(m_vi_Edges.size());
  }

  typename std::map<IntegerT, std::vector<IntegerT>>::iterator curr;
  m_vi_RightVertices.push_back(m_vi_Edges.size());
  for (IntegerT i = 0; i < i_ColumnCount; i++) {
    curr = colList.find(i);
    if (curr != colList.end()) {
      m_vi_Edges.insert(m_vi_Edges.end(), curr->second.begin(),
                        curr->second.end());
    }
    m_vi_RightVertices.push_back(m_vi_Edges.size());
  }
  return std::make_tuple(std::move(m_vi_LeftVertices),
                         std::move(m_vi_RightVertices), std::move(m_vi_Edges));
}

/**
 *
 * @tparam IntegerT
 * @param m_vi_LeftVertices
 * @param m_vi_RightVertices
 * @param m_vi_Edges
 * @return
 */
template <typename IntegerT>
std::tuple<IntegerT, IntegerT, IntegerT, IntegerT, IntegerT, IntegerT>
CalculateVertexDegrees(const std::vector<IntegerT> &m_vi_LeftVertices,
                       const std::vector<IntegerT> &m_vi_RightVertices,
                       const std::vector<IntegerT> &m_vi_Edges) {
  IntegerT i_LeftVertexCount(down_by_one(m_vi_LeftVertices.size()));
  IntegerT i_RightVertexCount(down_by_one(m_vi_RightVertices.size()));
  IntegerT i_TotalLeftVertexDegree(m_vi_Edges.size() / 2);
  IntegerT i_TotalRightVertexDegree(m_vi_Edges.size() / 2);

  IntegerT m_i_MinimumLeftVertexDegree =
      m_vi_LeftVertices[1] - m_vi_LeftVertices[0];
  IntegerT m_i_MaximumLeftVertexDegree =
      m_vi_LeftVertices[1] - m_vi_LeftVertices[0];
  for (IntegerT i = 0; i < i_LeftVertexCount; i++) {
    IntegerT i_VertexDegree = m_vi_LeftVertices[i + 1] - m_vi_LeftVertices[i];
    m_i_MaximumLeftVertexDegree =
        std::max(m_i_MaximumLeftVertexDegree, i_VertexDegree);
    m_i_MinimumLeftVertexDegree =
        std::min(m_i_MinimumLeftVertexDegree, i_VertexDegree);
  }

  IntegerT m_i_MinimumRightVertexDegree =
      m_vi_RightVertices[1] - m_vi_RightVertices[0];
  IntegerT m_i_MaximumRightVertexDegree =
      m_vi_RightVertices[1] - m_vi_RightVertices[0];
  for (IntegerT i = 0; i < i_RightVertexCount; i++) {
    IntegerT i_VertexDegree = m_vi_RightVertices[i + 1] - m_vi_RightVertices[i];
    m_i_MaximumRightVertexDegree =
        std::max(m_i_MaximumRightVertexDegree, i_VertexDegree);
    m_i_MinimumRightVertexDegree =
        std::min(m_i_MinimumRightVertexDegree, i_VertexDegree);
  }
  return std::make_tuple(
      m_i_MaximumLeftVertexDegree, m_i_MinimumLeftVertexDegree,
      i_TotalLeftVertexDegree, m_i_MaximumRightVertexDegree,
      m_i_MinimumRightVertexDegree, i_TotalRightVertexDegree);
}

#endif // TREEPART_GRAPH_COLOURING_HPP
