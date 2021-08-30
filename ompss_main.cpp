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

  // Clean up the memory
  free_data();
  return 0;
}
