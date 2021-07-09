#ifndef __BACKEND_HPP__

#define __BACKEND_HPP__

extern "C" {
  void read_data( char *mesh_filename );
  void get_data_ptr( int * nnodes,
                     int * ncells,
                     int * npedges,
                     double ** x,
                     double ** y,
                     int ** triangles,
                     int ** pedges);
  void free_data();
}

#endif

