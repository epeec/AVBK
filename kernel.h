#ifndef KERNEL_H

#define KERNEL_H

/**
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
 * @brief Obtain the element and nodal volume values
 * @param nnodes
 * @param ncells
 * @param x
 * @param y
 * @param triangles
 * @param vol_e
 * @param vol_n
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
