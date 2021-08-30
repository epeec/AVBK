/*
 *  @file main.cpp
 *  @brief Gather-Scatter Mini-App
 *  @author Pavanakumar Mohanamuraly
 */
#include "backend.h"
#include "kernel.h"
#include <chrono>
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
  read_data(args[1], 100, true);
  int nnodes, ncells, npedges, ncolours, maxgsize;
  int *cgroup_offset, *cgroup, *triangles, *pedges;
  int *group_offset;
  double *x, *y;
  get_data_ptr(&nnodes, &ncells, &npedges, &ncolours, &maxgsize, &x, &y,
               &group_offset, &cgroup_offset, &cgroup, &triangles, &pedges);

  // Allocate data for kernel
  double *vol_e = new double[ncells];
  double *vol_n = new double[nnodes];
  double *work = new double[2 * npedges];
#ifdef _OPENMP
  allocate_work_groups(maxgsize);
#endif

#ifndef _OPENMP
  std::cout << "Running in serial mode\n";
  // Zero out all the arrays
  ZeroOutArray(ncells, vol_e);
  ZeroOutArray(nnodes, vol_n);
  ZeroOutArray(2 * npedges, work);
  // Run the main kernel
  run_kernel(ncells, npedges, x, y, triangles, pedges, work, vol_e, vol_n);
#endif

#ifdef _OPENMP
  std::cout << "Running in OpenMP mode\n";
  auto start = std::chrono::steady_clock::now();
  unsigned niter = 1;
  for (unsigned iter = 0; iter < niter; ++iter) {
    element_volumes_cgroup(nnodes, ncolours, group_offset, cgroup_offset,
                           cgroup, x, y, triangles, vol_e, vol_n);
    //    if (npedges > 0)
    //      apply_periodic_bc(npedges, pedges, work, vol_n);
  }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() / niter << "s\n";
#endif

#if 0
  for (size_t i = 0; i < nnodes; ++i)
    std::cout << "   (" << i << "): " << get_node_perm_value(i, vol_n) << ",\n";
#endif

  // Clean up the memory
  free_data();
  delete[] vol_n;
  delete[] vol_e;
  delete[] work;
  free_work_groups();
  return 0;
}
