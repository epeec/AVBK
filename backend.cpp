/*
 * Gather-Scatter Mini-App
 * To compile : g++ -I.
 * -I/Users/mohanamuraly/NutsCFD_sandbox/NutsCFD/external/DIST/include
 *                  -L/Users/mohanamuraly/NutsCFD_sandbox/NutsCFD/external/DIST/lib
 * main.cpp -std=c++11 -lhdf5
 *
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
std::vector<double> __x, __y;
std::vector<int> __triangles;
std::vector<int> __pedges;
std::vector<int> __gofs, __gcolour;
std::vector<int> __epart;
std::vector<int> __perm;
std::vector<int> __eptr;

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
 * @param mesh_filename
 * @param mcg
 */
void read_data(char *mesh_filename, const int mcg) {
  __mcg = mcg;
  std::string m_mesh_filename(mesh_filename);
  // Open mesh file and create the necessary data for kernel exec
  hsize_t file_dims[10];
  auto file = H5Fopen(m_mesh_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  GetDimensions(file, "/Coordinates/x", file_dims);
  Read(file, "/Coordinates/x", __x);
  Read(file, "/Coordinates/y", __y);
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
  __perm.resize(__ncells);
  __eptr.resize(__ncells + 1);

  int ierr = 1;
  for (size_t i = 0; i < __ncells; ++i) {
    __eptr[i] = ierr;
    ierr = ierr + 3;
  }
  __eptr[__ncells] = ierr;

  // Do the colouring of cells
  colour_api_call(&__nnodes,          /* number of nodes */
                  &__ncells,          /* number of elements */
                  &__mcg,             /* number of element groups */
                  &__ncommon,         /* number of common nodes (dual-graph) */
                  __eptr.data(),      /* element ptr */
                  __triangles.data(), /* element node index */
                  __gofs.data(),      /* offset */
                  __gcolour.data(),   /* colour */
                  __epart.data(),     /* part array */
                  __perm.data(), /* permutation of the elements (size of ne) */
                  &ierr          /* Error code */
  );

#if 1
  // Write output of colours for plotting
  write_colour(__nnodes, __ncells, __mcg, __x.data(), __y.data(), nullptr,
               __eptr.data(), __triangles.data(), __gofs.data(),
               __gcolour.data(), __perm.data());
#endif

  for (auto &item : __triangles)
    --item;
  for (auto &item : __pedges)
    --item;
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
void get_data_ptr(int *nnodes, int *ncells, int *npedges, double **x,
                  double **y, int **triangles, int **pedges) {
  *nnodes = __nnodes;
  *ncells = __ncells;
  *npedges = __npedges;
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
