/*
 *  @file backend.cpp
 *  @brief The backend for setting up the data for the kernel
 *  @author Pavanakumar Mohanamuraly
 */

#include <cassert>
#include <hdf5.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "backend.h"
#include "colour.hpp"

int __nnodes = 0;
int __ncells = 0;
int __npedges = 0;
int __mcg = 10;
int __ncommon = 2;
int __max_cells_ingroup = 0;
int __ngcolours = 0;

std::vector<double> __x, __y, __z;
std::vector<int> __triangles;
std::vector<int> __pedges;
std::vector<int> __gofs, __gcolour;
std::vector<int> __cgofs, __cgroup;
std::vector<int> __epart;
std::vector<int> __eptr;
std::vector<int> __node_perm;
std::vector<int> __node_inv_perm;

/**
 * @brief
 *
 * @tparam T
 * @tparam uintT
 * @param data
 * @param perm
 * @param size
 */
template <typename T, typename uintT>
void InplacePermutationTriangles(size_t size, uintT *perm, T *data) {
  T temp[3];
  uintT j, k;
  for (uintT i = 0; i < size; ++i) {
    if (i != perm[i]) {
      temp[0] = data[3 * i];
      temp[1] = data[3 * i + 1];
      temp[2] = data[3 * i + 2];
      j = i;
      while (i != perm[j]) {
        k = perm[j];
        data[3 * j] = data[3 * k];
        data[3 * j + 1] = data[3 * k + 1];
        data[3 * j + 2] = data[3 * k + 2];
        perm[j] = j;
        j = k;
      }
      data[3 * j] = temp[0];
      data[3 * j + 1] = temp[1];
      data[3 * j + 2] = temp[2];
      perm[j] = j;
    }
  }
}

/**
 * @brief
 *
 * @tparam T
 * @tparam uintT
 * @param size
 * @param data1
 * @param data2
 * @param data3
 * @param perm
 */
template <typename T, typename uintT>
void InplacePermutation3(size_t size, uintT *perm, T *data1, T *data2,
                         T *data3) {
  T temp1, temp2, temp3;
  uintT j, k;
  for (uintT i = 0; i < size; ++i) {
    if (i != perm[i]) {
      temp1 = data1[i];
      temp2 = data2[i];
      temp3 = data3[i];
      j = i;
      while (i != perm[j]) {
        k = perm[j];
        data1[j] = data1[k];
        data2[j] = data2[k];
        data3[j] = data3[k];
        perm[j] = j;
        j = k;
      }
      data1[j] = temp1;
      data2[j] = temp2;
      data3[j] = temp3;
      perm[j] = j;
    }
  }
}

/**
 *
 * @param file_handle
 * @param path
 * @param file_dims
 * @return
 */
int GetDimensions(hid_t file_handle, const char *path, hsize_t *file_dims) {
  if (H5Lexists(file_handle, path, H5P_DEFAULT) <= 0)
    return -1;
  auto dset = H5Dopen(file_handle, path, H5P_DEFAULT);
  auto dspace = H5Dget_space(dset);
  auto rank = H5Sget_simple_extent_dims(dspace, file_dims, nullptr);
  H5Sclose(dspace);
  H5Dclose(dset);
  return rank;
}

/**
 * @brief Generic function to obtain HDF5 data type using template
 *
 * @tparam T
 * @return
 */
template <typename T> hid_t &GetDatatype();

/**
 *
 * @return
 */
template <> hid_t &GetDatatype<int>() { return H5T_NATIVE_INT; }

/**
 *
 * @return
 */
template <> hid_t &GetDatatype<double>() { return H5T_NATIVE_DOUBLE; }

/**
 *
 * @return
 */
template <> hid_t &GetDatatype<unsigned>() { return H5T_NATIVE_UINT; }

/**
 *
 * @tparam T
 * @param file
 * @param link
 * @param buf
 * @return
 */
template <typename T>
static int Read(hid_t &file, const char *link, std::vector<T> &buf) {
  if (H5Lexists(file, link, H5P_DEFAULT) <= 0) {
    std::stringstream cat;
    cat << "Error: Cannot find dataset " << link << " in hdf5 file";
    throw std::runtime_error(cat.str().c_str());
  }
  hsize_t dims[10];
  auto dset = H5Dopen(file, link, H5P_DEFAULT);
  auto dspace = H5Dget_space(dset);
  auto rank = H5Sget_simple_extent_dims(dspace, dims, nullptr);
  hsize_t bufsize = 1;
  for (int i = 0; i < rank; ++i)
    bufsize *= dims[i];
  buf.resize(bufsize);
  auto status = H5Dread(dset, GetDatatype<T>(), H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT, &buf[0]);
  assert(status >= 0);
  H5Sclose(dspace);
  H5Dclose(dset);
  return rank;
}

/**
 *
 * @tparam T_Real
 * @param x
 * @param y
 * @param s
 * @return
 */
template <typename T_Real>
T_Real area_triangle(const T_Real x[3], const T_Real y[3], T_Real *const s) {
  const auto fac = 1.0e0 / 3.0e0;
  auto val = 0.5 * (x[0] * y[1] - x[1] * y[0] + x[1] * y[2] - x[2] * y[1] +
                    x[2] * y[0] - x[0] * y[2]);
  s[0] = fac * val;
  s[1] = fac * val;
  s[2] = fac * val;
  return val;
}

/**
 *
 */
void ReorderNodesByGroup() {
  std::cout << "Reordering nodes by first touch of group policy\n";
  // std::map<unsigned, unsigned> map_reorder, inv_map_reorder;
  std::set<int> unique_node;
  __node_perm.resize(__nnodes);
  __node_inv_perm.resize(__nnodes);
  unsigned count = 0;
  for (unsigned kcolour = 0; kcolour < __ngcolours; ++kcolour) {
    for (unsigned k = __cgofs[kcolour]; k < __cgofs[kcolour + 1]; ++k) {
      unsigned kgroup = __cgroup[k];
      unsigned ibeg = __gofs[kgroup];
      unsigned iend = __gofs[kgroup + 1];
      for (unsigned i = ibeg; i < iend; ++i) {
        // Loop over the node ids of the triangle
        for (int id = 0; id < 3; ++id) {
          // Check if any of the ids are already present in the hash
          auto lid = __triangles[i * 3 + id];
          // If not present in hash add this to hash and continue
          if (unique_node.find(lid) == std::end(unique_node)) {
            unique_node.insert(lid);
            __node_perm[count] = lid;
            count++;
          }
          if (count - 1 == 0)
            std::cout << "First node in hash " << __node_perm[0] << "  count "
                      << count << "\n";
        }
      }
    }
  }
  std::cout << "Total nodes in mesh " << __nnodes << " total nodes in hash "
            << count << "\n";
  // Renumber the triangle node-ids
  for (unsigned i = 0; i < __nnodes; ++i)
    __node_inv_perm[__node_perm[i]] = i;
  for (auto &item : __triangles)
    item = __node_inv_perm[item];
  //  // Renumber periodicity
  //  for (auto &item : __pedges)
  //    item = map_reorder[item];
  // Reorder the x, y and z coordinates using the permutation
  std::vector<int> node_perm(__node_perm);
  InplacePermutation3(__nnodes, node_perm.data(), __x.data(), __y.data(),
                      __z.data());
}

/**
 *
 * @param mesh_filename
 * @param mcg
 */
void read_data(char *mesh_filename, const int mcg, const bool reorder) {
  std::vector<int> cell_perm;
  __mcg = mcg;
  std::string m_mesh_filename(mesh_filename);
  // Open mesh file and create the necessary data for kernel exec
  hsize_t file_dims[10];
  auto file = H5Fopen(m_mesh_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  GetDimensions(file, "/Coordinates/x", file_dims);
  Read(file, "/Coordinates/x", __x);
  Read(file, "/Coordinates/y", __y);
  __z.resize(__x.size());
  Read(file, "/Connectivity/tri->node", __triangles);
  __nnodes = __x.size();
  __ncells = __triangles.size() / 3;
  std::cout << "Nodes = " << __nnodes << "\nTriangles = " << __ncells << "\n";

  if (H5Lexists(file, "/Periodicity/periodic_node", H5P_DEFAULT) > 0) {
    Read(file, "/Periodicity/periodic_node", __pedges);
    __npedges = __pedges.size() / 2;
    std::cout << "Periodic = " << __npedges << "\n";
  }

  // Convert the Fortran 1 index to C style
  H5Fclose(file);

  __gofs.resize(__mcg + 1);
  __gcolour.resize(__mcg);
  __epart.resize(__ncells);
  cell_perm.resize(__ncells);
  __eptr.resize(__ncells + 1);

  int ierr = 1;
  for (size_t i = 0; i < __ncells; ++i) {
    __eptr[i] = ierr;
    ierr = ierr + 3;
  }
  __eptr[__ncells] = ierr;

  // Do the colouring of cells
  colour_api_call(
      &__nnodes,           /* number of nodes */
      &__ncells,           /* number of elements */
      &__mcg,              /* number of element groups */
      &__ncommon,          /* number of common nodes (dual-graph) */
      __eptr.data(),       /* element ptr */
      __triangles.data(),  /* element node index */
      __gofs.data(),       /* offset */
      __gcolour.data(),    /* colour */
      __epart.data(),      /* part array */
      cell_perm.data(),    /* permutation of the elements (size of ne) */
      &__max_cells_ingroup /* Max cells in a group */
  );

#if 0
  // Write output of colours for plotting
  write_colour(__nnodes, __ncells, __mcg, __x.data(), __y.data(), nullptr,
               __eptr.data(), __triangles.data(), __gofs.data(),
               __gcolour.data(), cell_perm.data());
#endif

#if 0
  for (const auto &item : __gcolour)
    std::cout << "colour " << item << "\n";
#endif

  // Form the colour groups and offsets
  __ngcolours = *(std::max_element(__gcolour.begin(), __gcolour.end())) + 1;
  __cgofs.resize(__ngcolours + 1);
  __cgroup.resize(__mcg);

  unsigned j = 0;
  for (unsigned kcolour = 0; kcolour < __ngcolours; ++kcolour) {
    for (unsigned kgroup = 0; kgroup < __mcg; ++kgroup) {
      if (__gcolour[kgroup] == kcolour) {
        __cgroup[j] = kgroup;
        ++__cgofs[kcolour + 1];
        ++j;
      }
    }
  }
  for (unsigned kcolour = 1; kcolour < __ngcolours + 1; ++kcolour)
    __cgofs[kcolour] += __cgofs[kcolour - 1];

#if 0
  for (unsigned kcolour = 0; kcolour < __ngcolours; ++kcolour) {
    std::cout << "colour " << kcolour << " begin " << __cgofs[kcolour]
              << " end " << __cgofs[kcolour + 1] << "\n";
    for (unsigned k = __cgofs[kcolour]; k < __cgofs[kcolour + 1]; ++k) {
      unsigned kgroup = __cgroup[k];
      std::cout << "group " << kgroup << "\n";
    }
  }
#endif

  for (auto &item : __triangles)
    --item;
  for (auto &item : __pedges)
    --item;
  for (auto &item : __gofs)
    --item;

  InplacePermutationTriangles(__ncells, cell_perm.data(), __triangles.data());

  if (reorder) {
    ReorderNodesByGroup();
  } else {
    __node_perm.resize(__nnodes);
    __node_inv_perm.resize(__nnodes);
    for (int i = 0; i < __nnodes; ++i) {
      __node_perm[i] = i;
      __node_inv_perm[i] = i;
    }
  }
}

/**
 *
 * @param nnodes
 * @param ncells
 * @param npedges
 * @param x
 * @param y
 * @param triangles
 * @param pedges
 */
void get_data_ptr(int *nnodes, int *ncells, int *npedges, int *ncolours,
                  int *maxgsize, double **x, double **y, int **group_offset,
                  int **cgroup_offset, int **cgroup, int **triangles,
                  int **pedges) {
  *nnodes = __nnodes;
  *ncells = __ncells;
  *npedges = __npedges;
  *ncolours = __ngcolours;
  *maxgsize = __max_cells_ingroup;
  *group_offset = __gofs.data();
  *cgroup_offset = __cgofs.data();
  *cgroup = __cgroup.data();
  *x = __x.data();
  *y = __y.data();
  *triangles = __triangles.data();
  *pedges = __pedges.data();
}

/**
 *
 */
void free_data() {
  __nnodes = 0;
  __ncells = 0;
  __npedges = 0;
  __x.clear();
  __x.shrink_to_fit();
  __y.clear();
  __y.shrink_to_fit();
  __triangles.clear();
  __triangles.shrink_to_fit();
  __pedges.clear();
  __pedges.shrink_to_fit();
}

/**
 *
 * @param i
 * @param value
 * @return
 */
double get_node_perm_value(int i, double *value) {
  return value[__node_inv_perm[i]];
}
