#ifndef KERNEL_H

#define KERNEL_H

void element_volumes_cgroup(int nnodes, int ncolours, int max_group_size,
                            const int *group_offset, const int *cgroup_offset,
                            const int *cgroup, const double *x, const double *y,
                            const int *triangles,  double *vol_e, double *vol_n);
#endif
