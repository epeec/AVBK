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
                            const int *triangles, double *vol_e, double *vol_n) {
  {
    // Zero-out the values
    for (unsigned i = 0; i < nnodes; ++i)
      vol_n[i] = 0.0;

    // Loop over all colours
    for (unsigned kcolour = 0; kcolour < ncolours; ++kcolour) {
      unsigned offset_size = cgroup_offset[kcolour + 1] - cgroup_offset[kcolour];
      unsigned grp_size = 3 * max_group_size;

      // Loop over groups of the same colour (can be run on multiple threads)
      // (the taskloop is what we have to use here and not task)
#pragma oss taskloop
      for (unsigned k = cgroup_offset[kcolour]; k < cgroup_offset[kcolour + 1];
           ++k) {
         double *work_x = (double *)malloc(sizeof(double) * grp_size );
         double *work_y = (double *)malloc(sizeof(double) * grp_size );
         double *work_scatter = (double *)malloc(sizeof(double) * grp_size );

      	unsigned kgroup = cgroup[k];
        unsigned ibeg = group_offset[kgroup];
        unsigned iend = group_offset[kgroup + 1];
        unsigned isize = iend - ibeg;

        // Gather data at nodes
        for (unsigned i = ibeg; i < iend; ++i) {
          const int lid[3]{triangles[i * 3], triangles[i * 3 + 1],
                           triangles[i * 3 + 2]};
          const double lx[3]{x[*lid], x[*(lid + 1)], x[*(lid + 2)]};
          const double ly[3]{y[*lid], y[*(lid + 1)], y[*(lid + 2)]};
          for (int inode = 0; inode < 3; ++inode) {
            work_x[ ( 3 * (i - ibeg) + inode) ] = lx[inode];
            work_y[ ( 3 * (i - ibeg) + inode) ] = ly[inode];
          }
        }

        // Compute using gathered data
        for (unsigned i = 0; i < isize; ++i) {
          const double fac = 1.0e0 / 3.0e0;
          const double *const x_ptr = &work_x[ (3 * i) ];
          const double *const y_ptr = &work_y[ (3 * i) ];

          double *const s = &work_scatter[(3 * i) ];
          const double val = 0.5 * (x_ptr[0] * y_ptr[1] - x_ptr[1] * y_ptr[0] +
                                    x_ptr[1] * y_ptr[2] - x_ptr[2] * y_ptr[1] +
                                    x_ptr[2] * y_ptr[0] - x_ptr[0] * y_ptr[2]);
          s[0] = fac * val;
          s[1] = fac * val;
          s[2] = fac * val;
          vol_e[i + ibeg] = val;
        }

        // Scatter and accumulation
        for (unsigned i = ibeg; i < iend; ++i) {
          const int lid[3]{triangles[i * 3], triangles[i * 3 + 1],
                           triangles[i * 3 + 2]};
          for (int inode = 0; inode < 3; ++inode)
            vol_n[lid[inode]] += work_scatter[3 * (i - ibeg) + inode];
        }

        // Free the task local arrays that were allocated
        free(work_x);
        free(work_y);
        free(work_scatter);

      } // End of group of cells within colour loop
#pragma oss taskwait

    }  // End of colour loop
  }
}
}
