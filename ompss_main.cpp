/*
 *  @file main.cpp
 *  @brief Gather-Scatter Mini-App
 *  @author Pavanakumar Mohanamuraly
 */
#include "backend.h"
#include "kernel_ompss.h"
#include <iostream>

/**
 *
 * @param data
 * @param n
 */
void ZeroOutArray(const unsigned n, double *data) {
  for (unsigned i = 0; i < n; ++i)
    data[i] = 0.0;
}

/**
 *
 * @param nargs
 * @param args
 * @return
 */
int main(int nargs, char *args[]) {

  if (nargs != 2) {
    std::cout << "Need just one argument <input-hip-hdf5-mesh-file>\n";
    return -100;
  }

  // Read mesh data using C API call
  read_data(args[1], 1000, true);
  int nnodes, ncells, npedges, ncolours, maxgsize;
  int *cgroup_offset, *cgroup, *triangles, *pedges;
  int *group_offset;
  double *x, *y;
  get_data_ptr(&nnodes, &ncells, &npedges, &ncolours, &maxgsize, &x, &y,
               &group_offset, &cgroup_offset, &cgroup, &triangles, &pedges);

  // Allocate data for kernel
  double *vol_e = new double[ncells];
  double *vol_n = new double[nnodes];

  element_volumes_cgroup(nnodes, ncolours, maxgsize, group_offset, cgroup_offset,
                         cgroup, x, y, triangles, vol_e, vol_n);

#if 1
  for (size_t i = 0; i < nnodes; ++i)
    std::cout << "   (" << i << "): " << get_node_perm_value(i, vol_n) << ",\n";
#endif

//  // Clean up the memory
  free_data();
  delete[] vol_n;
  delete[] vol_e;
  return 0;
}
