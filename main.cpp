/*
 *  @file main.cpp
 *  @brief Gather-Scatter Mini-App
 *  @author Pavanakumar Mohanamuraly
 */
#include "backend.h"
#include "kernel.h"
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
    abort();
  }

  // Read mesh data using C API call
  read_data(args[1], 100);
  int nnodes, ncells, npedges, ncolours, maxgsize;
  int *cgroup_offset, *cgroup, *triangles, *pedges;
  int *group_offset, *group;
  double *x, *y;
  get_data_ptr(&nnodes, &ncells, &npedges, &ncolours, &maxgsize, &x, &y,
               &group_offset, &group, &cgroup_offset, &cgroup, &triangles,
               &pedges);

  // Allocate data for kernel
  double *vol_e = new double[ncells];
  double *vol_n = new double[nnodes];
  double *work = new double[2 * npedges];

  // Zero out all the arrays
  ZeroOutArray(ncells, vol_e);
  ZeroOutArray(nnodes, vol_n);
  ZeroOutArray(2 * npedges, work);

#ifndef _OPENMP
  // Run the main kernel
  run_kernel(ncells, npedges, x, y, triangles, pedges, work, vol_e, vol_n);
#endif

#ifdef _OPENMP
  allocate_work_groups(maxgsize);
  element_volumes_cgroup(ncells, ncolours, group_offset, group, cgroup_offset,
                         cgroup, x, y, triangles, vol_e, vol_n);
  if (npedges > 0)
    apply_periodic_bc(npedges, pedges, work, vol_n);
#endif

#if 0
  for (size_t i = 0; i < ncells; ++i)
    std::cout << "   (" << i << "): " << vol_e[i] << ",\n";
#endif

  // Clean up the memory
  free_data();
  delete[] vol_n;
  delete[] vol_e;
  delete[] work;
  free_work_groups();
  return 0;
}
