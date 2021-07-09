
double area_triangle(const double x[3], const double y[3], double * __restrict__ const s) {
  const double fac = 1.0e0 / 3.0e0;
  double val = 0.5 * (x[0] * y[1] - x[1] * y[0] + x[1] * y[2] - x[2] * y[1] + x[2] * y[0] - x[0] * y[2]);
  s[0] = fac * val;
  s[1] = fac * val;
  s[2] = fac * val;
  return val;
}

// Obtain the element and nodal volume values
void element_volumes(const int nnodes,
                     const int ncells,
                     const double * __restrict__ const x,
                     const double * __restrict__ const y, 
                     const int * __restrict__ const triangles,
                     double * __restrict__ const vol_e,
                     double * __restrict__ const vol_n )
{
  // Loop over all elements
  for( size_t i=0; i<ncells; ++i) {
    // Gather data to local nodes
    const int lid[3]{ triangles[i * 3], triangles[i * 3 + 1], triangles[i * 3 + 2] };
    const double lx[3]{ x[*lid], x[*(lid + 1)], x[*(lid + 2)] };
    const double ly[3]{ y[*lid], y[*(lid + 1)], y[*(lid + 2)] };
    // Compute local data (computationally demanding)
    double scatter_val[3];
    vol_e[i] = area_triangle( lx, ly, scatter_val );
    // Scatter to global nodes (accumulate)
    for(size_t j=0; j<3; ++j)
      vol_n[lid[j]] += scatter_val[j];
  }
}

void apply_periodic_bc(const int npedges,
                       const int * __restrict__ const pedges,
                       double * __restrict__ const res) {
  size_t ibc_start = 0;
  size_t ibc_end   = npedges - 1;
  for( size_t i=0; i < npedges; ++i ) {
    const int i1 = pedges[ibc_start];
    const int i2 = pedges[ibc_end];
    const double new_res = res[i1] + res[i2];
    res[i1] = new_res;
    res[i2] = new_res;
  }
}

