#ifndef KERNEL_H

#define KERNEL_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#if defined(_OPENMP)
#include <omp.h>
#endif

double *work_x;
double *work_y;
double *work_scatter;
#pragma omp threadprivate(work_x, work_y, work_scatter)

/**
 *
 * @brief Calculate the area of the triangle given the three vertices
 * @param x
 * @param y
 * @param s
 * @return
 */
double area_triangle(const double x[3], const double y[3], double *const s) {
  const double fac = 1.0e0 / 3.0e0;
  double val = 0.5 * (x[0] * y[1] - x[1] * y[0] + x[1] * y[2] - x[2] * y[1] +
                      x[2] * y[0] - x[0] * y[2]);
  s[0] = fac * val;
  s[1] = fac * val;
  s[2] = fac * val;
  return val;
}

/**
 *
 * @brief Obtain the element and nodal volume values
 * @param nnodes
 * @param ncells
 * @param x
 * @param y
 * @param triangles
 * @param vol_e
 * @param vol_n
 *
 */
void element_volumes(const int ncells, const double *const x,
                     const double *const y, const int *const triangles,
                     double *const vol_e, double *const vol_n) {
  // Loop over all elements
  for (unsigned i = 0; i < ncells; ++i) {
    // Gather data to local nodes
    const int lid[3]{triangles[i * 3], triangles[i * 3 + 1],
                     triangles[i * 3 + 2]};
    const double lx[3]{x[*lid], x[*(lid + 1)], x[*(lid + 2)]};
    const double ly[3]{y[*lid], y[*(lid + 1)], y[*(lid + 2)]};
    // Compute local data (computationally demanding)
    double scatter_val[3];
    vol_e[i] = area_triangle(lx, ly, scatter_val);
    // Scatter to global nodes (accumulate)
    for (unsigned j = 0; j < 3; ++j)
      vol_n[lid[j]] += scatter_val[j];
  }
}

/**
 *
 * @param max_group_size
 *
 */
void allocate_work_groups(const unsigned max_group_size) {
#pragma omp parallel
  {
#if 0
    int nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    if (tid)
      std::cout << "Total threads = " << nt << "\n";
    std::cout << "Allocating arrays on thread-id " << tid << " of size "
              << max_group_size << "\n";
#endif
    work_x = (double *)calloc(3 * max_group_size, sizeof(double));
    work_y = (double *)calloc(3 * max_group_size, sizeof(double));
    work_scatter = (double *)calloc(3 * max_group_size, sizeof(double));
#if 0
    for (unsigned int i = 0; i < 3 * max_group_size; ++i) {
      work_x[i] = double(tid);
      work_y[i] = double(tid);
    }
#endif
  }
}

/**
 *
 */
void free_work_groups() {
#pragma omp parallel
  {
#if 0
    int tid = omp_get_thread_num();
    std::cout << "TID " << tid << " first data " << work_x[0] << ", "
              << work_y[0] << "\n";
#endif
    free(work_x);
    free(work_y);
    free(work_scatter);
  }
}

/**
 *
 * @param ncells
 * @param ncolours
 * @param group_offset
 * @param group
 * @param cgroup_offset
 * @param cgroup
 * @param x
 * @param y
 * @param triangles
 * @param vol_e
 * @param vol_n
 */
void element_volumes_cgroup(int ncells, int ncolours, const int *group_offset,
                            const int *group, const int *cgroup_offset,
                            const int *cgroup, const double *x, const double *y,
                            const int *triangles, double *vol_e,
                            double *vol_n) {
#pragma omp parallel shared(ncolours, group, group_offset, cgroup,             \
                            cgroup_offset, triangles, x, y, vol_e,             \
                            vol_n) default(none)
  {
#if 1
    int tid = omp_get_thread_num();
    std::ofstream fout("thread_" + std::to_string(tid) + ".dat");
    fout << "ncolours = " << ncolours << "\n";
#endif
    // Loop over all colours
    for (unsigned kcolour = 0; kcolour < ncolours; ++kcolour) {
      fout << "In colour " << kcolour << " cgroup_offset[kcolour] "
           << cgroup_offset[kcolour] << " cgroup_offset[kcolour + 1] "
           << cgroup_offset[kcolour + 1] << "\n";
#pragma omp for
      for (unsigned k = cgroup_offset[kcolour]; k < cgroup_offset[kcolour + 1];
           ++k) { // Loop over groups of the same colour (can be run on multiple
                  // threads)
        unsigned kgroup = cgroup[k];
        unsigned ibeg = group_offset[kgroup];
        unsigned iend = group_offset[kgroup + 1];
        unsigned isize = iend - ibeg;
        fout << "In group " << kgroup << " ibeg " << ibeg << " iend " << iend
             << " isize " << isize << "\n";
        // gather data at nodes
        for (unsigned i = ibeg; i < iend; ++i) {
#if 0
          fout << "TID " << tid << " colour " << kcolour << " group " << kgroup
               << " i " << i - ibeg << " size " << isize << "\n";
#endif
          const int lid[3]{triangles[i * 3], triangles[i * 3 + 1],
                           triangles[i * 3 + 2]};
          const double lx[3]{x[*lid], x[*(lid + 1)], x[*(lid + 2)]};
          const double ly[3]{y[*lid], y[*(lid + 1)], y[*(lid + 2)]};
          for (int inode = 0; inode < 3; ++inode) {
            work_x[3 * (i - ibeg) + inode] = lx[inode];
            work_y[3 * (i - ibeg) + inode] = ly[inode];
          }
        }

        // compute using gathered data
        for (unsigned i = 0; i < isize; ++i) {
          const double fac = 1.0e0 / 3.0e0;
          const double *const x_ptr = &work_x[3 * i];
          const double *const y_ptr = &work_y[3 * i];
          fout << *x_ptr << ", " << *y_ptr << "\n";
          double *const s = &work_scatter[3 * i];
          const double val = 0.5 * (x_ptr[0] * y_ptr[1] - x_ptr[1] * y_ptr[0] +
                                    x_ptr[1] * y_ptr[2] - x_ptr[2] * y_ptr[1] +
                                    x_ptr[2] * y_ptr[0] - x_ptr[0] * y_ptr[2]);
          s[0] = fac * val;
          s[1] = fac * val;
          s[2] = fac * val;
          vol_e[i + ibeg] = val;
        }

        // scatter and accumulation
        for (unsigned i = ibeg; i < iend; ++i) {
          const int lid[3]{triangles[i * 3], triangles[i * 3 + 1],
                           triangles[i * 3 + 2]};
          for (int inode = 0; inode < 3; ++inode)
            vol_n[lid[inode]] += work_scatter[3 * (i - ibeg) + inode];
        }
      }
    }
  }
}

/**
 * @brief Apply the periodic BC for nodal residue (using work array)
 * @param npedges
 * @param pedges
 * @param pwork
 * @param res
 */
void apply_periodic_bc(const int npedges, const int *const pedges,
                       double *const pwork, double *const res) {
  unsigned ibc_start = 0;
  unsigned ibc_end = 2 * npedges - 1;
  for (unsigned i = 0; i < npedges; ++i) {
    const unsigned i1 = pedges[ibc_start];
    const unsigned i2 = pedges[ibc_end];
    pwork[ibc_start] += res[i2];
    pwork[ibc_end] += res[i1];
    ibc_start++;
    ibc_end--;
  }

  ibc_start = 0;
  ibc_end = 2 * npedges - 1;
  for (unsigned i = 0; i < npedges; ++i) {
    const unsigned i1 = pedges[ibc_start];
    const unsigned i2 = pedges[ibc_end];
    res[i1] += pwork[ibc_start];
    res[i2] += pwork[ibc_end];
    ibc_start++;
    ibc_end--;
  }
}

/**
 * @brief Run the kernel
 * @param nnodes
 * @param ncells
 * @param npedges
 * @param x
 * @param y
 * @param triangles
 * @param pedges
 * @param work
 * @param vol_e
 * @param vol_n
 */
void run_kernel(const int ncells, const int npedges, const double *const x,
                const double *const y, const int *const triangles,
                const int *const pedges, double *const work,
                double *const vol_e, double *const vol_n) {

  element_volumes(ncells, x, y, triangles, vol_e, vol_n);
  if (npedges > 0)
    apply_periodic_bc(npedges, pedges, work, vol_n);
}

#endif
