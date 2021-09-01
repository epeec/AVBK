#include "kernel_ompss.h"
#include <cstdlib>

/**
 *
 * @param nnodes
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
void element_volumes_cgroup(int nnodes, int ncolours, int max_group_size,
                            const int *group_offset, const int *cgroup_offset,
                            const int *cgroup, const double *x, const double *y,
                            const int *triangles, double *vol_e,
                            double *vol_n) {
  //#pragma omp parallel shared(nnodes, ncolours, group_offset, cgroup,            \
//                            cgroup_offset, triangles, x, y, vol_e,             \
//                            vol_n) default(none)
  {
#if 0
    int tid = omp_get_thread_num();
    std::ofstream fout("thread_" + std::to_string(tid) + ".dat");
    fout << "ncolours = " << ncolours << "\n";
#endif
    // Zero-out the values
    //#pragma omp for
    for (unsigned i = 0; i < nnodes; ++i)
      vol_n[i] = 0.0;

    double *work_x = (double *)malloc(sizeof(double) * max_group_size);
    double *work_y = (double *)malloc(sizeof(double) * max_group_size);
    double *work_scatter = (double *)malloc(sizeof(double) * max_group_size);

    // Loop over all colours
    for (unsigned kcolour = 0; kcolour < ncolours; ++kcolour) {
#if 0
      fout << "In colour " << kcolour << " cgroup_offset[kcolour] "
      << cgroup_offset[kcolour] << " cgroup_offset[kcolour + 1] "
      << cgroup_offset[kcolour + 1] << "\n";
#endif
      //#pragma omp for
      for (unsigned k = cgroup_offset[kcolour]; k < cgroup_offset[kcolour + 1];
           ++k) { // Loop over groups of the same colour (can be run on multiple
        // threads)
        unsigned kgroup = cgroup[k];
        unsigned ibeg = group_offset[kgroup];
        unsigned iend = group_offset[kgroup + 1];
        unsigned isize = iend - ibeg;
#if 0
        fout << "In group " << kgroup << " ibeg " << ibeg << " iend " << iend
        << " isize " << isize << "\n";
#endif
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
#if 0
          fout << *x_ptr << ", " << *y_ptr << "\n";
#endif
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