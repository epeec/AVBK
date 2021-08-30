#ifndef __BACKEND_HPP__

#define __BACKEND_HPP__

extern "C" {
/**
 *
 * @param mesh_filename
 * @param mcg
 */
void read_data(char *mesh_filename, const int mcg);

/**
 *
 * @param nnodes
 * @param ncells
 * @param npedges
 * @param x
 * @param y
 * @param triangles
 * @param pedges
 */
void get_data_ptr(int *nnodes, int *ncells, int *npedges, int *ncolours,
                  int *maxgsize, double **x, double **y, int **group_offset,
                  int **group, int **cgroup_offset, int **cgroup,
                  int **triangles, int **pedges);
void free_data();
}

#endif
