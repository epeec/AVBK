/** @file colour.cpp
 *  @brief The main colouring API for AVBP OpenMP version. Taken
 *         from TreePart for quick and dirty colouring interface
 *         to AVBP OpenMP version
 *  @author Pavanakumar Mohanamuraly
*/

#ifndef COLOUR_HPP

#define COLOUR_HPP

#include "graph_colouring.hpp"
#include <metis.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <memory>
#include <set>
#include <fstream>

extern "C" {
void colour_api_call
(
    
  int32_t * nn,      /* number of nodes */
  int32_t * ne,      /* number of elements */
  int32_t * mcg,     /* number of element groups */
  int32_t * ncommon, /* number of common nodes (dual-graph) */
  int32_t * eptr,    /* element ptr */
  int32_t * eind,    /* element node index */
  int32_t * gofs,    /* offset */
  int32_t * gcolour, /* colour */
  int32_t * epart,   /* part array */
  int32_t * perm,    /* permutation of the elements (size of ne) */
  int32_t * ierr     /* Error code */
);
}


/*
 *  @brief Group cells into cache block sizes and colour them
 *
 *  perm (integer[ne])
 *  |-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-|
 *
 *  cell groups (integer[ng])
 *  |--------|---|------|----|----|----|-------|------|---------|
 *
 *  element type offset (obtained from eptr and eind itself integer[ntypes+1])
 *  |------TET---------------|------PYR--------|------HEX-------|
 *
 *  colour array (integer[ng]) stores the group ids that have the same colour
 *  colour offset (integer[ng + 1])
*/
int colour_api_call1
(
    
  const int32_t            nn, /* number of nodes */
  const int32_t            ne, /* number of elements */
  int32_t                 mcg, /* number of element groups */
  int32_t             ncommon, /* number of common nodes (dual-graph) */
  const int32_t * const  eptr, /* element ptr */
  const int32_t * const  eind, /* element node index */
  int32_t               *gofs, /* offset */
  int32_t            *gcolour, /* colour */
  int32_t              *epart, /* part array */
  int32_t * const        perm  /* permutation of the elements (size of ne) */
)
{
  // Use recursive bisection since max load imbalance is very low
  int32_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
  options[METIS_OPTION_NUMBERING] = 1; // Fortran style numbering
  options[METIS_OPTION_UFACTOR] = 1; // 1% load imbalance

  // Allocate the part arrays
  //std::cerr << "Number of cell groups " << mcg << "\n";
  // Call METIS to form the groups of cells
  {
    int32_t objval;
    std::vector<int32_t> npart(nn);
    //std::cerr << "ne = " << ne << "\n";
    //std::cerr << "nn = " << ne << "\n";
    //std::cerr << "eptr[0-1] = " << eptr[0] << ", " << eptr[1] << "\n";
    //std::cerr << "eind[0-2] = " << eind[0] << ", " << eind[1] << ", " << eind[2] << "\n";
    //std::cerr << "epart[0] = " << epart[0] << "\n";
    //std::cerr << "perm[0] = " << perm[0] << "\n";

    auto metis_ret = 
      METIS_PartMeshDual(const_cast<int32_t *>(&ne),
                         const_cast<int32_t *>(&nn),
                         const_cast<int32_t *>(eptr),
                         const_cast<int32_t *>(eind),
                         nullptr, nullptr,
                         &ncommon, &mcg, nullptr,
                         options, &objval,
                         &epart[0], &npart[0]);
    //std::cerr << "METIS returned = " << metis_ret << "\n";
    assert(metis_ret == METIS_OK);
  }

  // Node-Element graph
  std::vector<int32_t> ieptr(nn + 1), ieind;
  std::vector< std::set<int32_t> > iegraph(nn);
  for ( int32_t i = 0; i < ne; ++i )
    for ( int32_t j = eptr[i] - 1; j < eptr[i + 1] - 1; ++j )
      iegraph[ eind[j] - 1 ].insert(i + 1);
 
  ieptr[0] = 1;
  for ( int32_t i = 0; i < nn; ++i ) {
    ieptr[i + 1] = ieptr[i] + iegraph[i].size();
    std::copy( iegraph[i].begin(), iegraph[i].end(), std::back_inserter(ieind) );
  }
  iegraph.clear();

  // Form the connectivity between the
  // groups using ieptr and ieind
  std::vector<int32_t> cg_xadj, cg_adjncy;
  std::vector<int32_t> cg_colour;
  {
    std::vector<std::set<int32_t>> cg_graph(mcg);
    for ( int32_t i = 0; i < ne; ++i ) {
      auto my_grp_part = epart[i] - 1;
      for ( auto j = eptr[i] - 1; j < eptr[i + 1] - 1; ++j ) {
        auto inode = eind[j] - 1;
        for ( auto k = ieptr[inode] - 1; k < ieptr[inode + 1] - 1; ++k ) {
          auto ng_grp_part = epart[ ieind[k] - 1 ] - 1;
          cg_graph[my_grp_part].insert(ng_grp_part);
        }
      }
    }
    ieind.clear();
    ieptr.clear();

    // Convert graph to CSR format for colouring
    cg_xadj.resize(mcg + 1);
    cg_xadj[0] = 0;
    for ( int32_t i = 0; i < mcg; ++i ) {
      cg_xadj[i + 1] = cg_xadj[i] + cg_graph[i].size();
      std::copy( cg_graph[i].begin(), cg_graph[i].end(), std::back_inserter(cg_adjncy) );
    }
  }

  // Init the colouring interface and form the group colours
  {
    std::vector<int32_t> left, right, edges;
    std::tie(left, right, edges) = BuildBPGraphFromCSRFormat(&cg_xadj[0], mcg, mcg, &cg_adjncy[0]);
    cg_colour = std::move(PartialDistanceTwoColumnColoring(left, right, edges));
    CheckPartialDistanceTwoColumnColoring(left, right, edges, cg_colour);
  }

  // Colour the graph
  for ( int32_t i = 0; i < mcg; ++i )
    gcolour[i] = cg_colour[i];

  // Init the permutation (original order)
  for ( int32_t i = 0; i < ne; ++i )
    perm[i] = i;

  // Reorder the global ids (element-wise) as per the group id
  std::sort( perm, perm + ne,
             /* Sort lambda function */
            [&epart](const int32_t &a, const int32_t &b) {
              return epart[a] < epart[b];
            } );

  // Form the group offsets
  gofs[0] = 1;
  int32_t old_cg = epart[perm[0]], count = 1;
  for ( int32_t i = 0; i < ne; ++i ) {
    if ( old_cg != epart[perm[i]] ) {
      old_cg = epart[perm[i]];
      gofs[count++] = i + 1;
    }
  }
  gofs[mcg] = ne + 1;

  std::sort( epart, epart + ne );
  return 0;
}

/*
 *
*/
void write_colour
(
  int32_t nn,
  int32_t ne,
  int32_t mcg,
  double  *x,
  double  *y,
  double  *z,
  int32_t *eptr,
  int32_t *eind,
  int32_t *gofs,
  int32_t *gcolour,
  int32_t *perm
) {
  std::ofstream fout("tec.dat");

  // Title header
  fout << "VARIABLES =\"X\", \"Y\", ";
  if ( z != nullptr ) fout << "\"Z\", ";
  fout << "\"G_COLOUR\", \"CG_COLOUR\"\n";

  // Zone header
  fout << "ZONE DATAPACKING=BLOCK, NODES=" << nn;
  fout << ", ELEMENTS=" << ne << ", ZONETYPE=FETRIANGLE,";
  fout << " VARLOCATION=([3-4]=CELLCENTERED)\n";

  // Write coordinates
  for ( int32_t i = 0; i < nn; ++i ) fout << x[i] << "\n";
  for ( int32_t i = 0; i < nn; ++i ) fout << y[i] << "\n";
  if ( z != nullptr )
    for ( int32_t i = 0; i < nn; ++i )
      fout << z[i] << "\n";

  // Write group colour
  for ( int32_t i = 0; i < mcg; ++i )
    for ( int32_t j = gofs[i] - 1; j < gofs[i + 1] - 1; ++j )
      fout << i << "\n";
  for ( int32_t i = 0; i < mcg; ++i )
    for ( int32_t j = gofs[i] - 1; j < gofs[i + 1] - 1; ++j )
      fout << gcolour[i] << "\n";

  // Write the permutation to a tecplot file
  for ( int32_t ii = 0; ii < ne; ++ii ) {
    auto i = perm[ii];
    for ( auto j = eptr[i] - 1; j < eptr[i + 1] - 1; ++j ) {
      fout << eind[j] << " ";
    }
    fout << "\n";
  }
}


void colour_api_call
(
    
  int32_t * nn,      /* number of nodes */
  int32_t * ne,      /* number of elements */
  int32_t * mcg,     /* number of element groups */
  int32_t * ncommon, /* number of common nodes (dual-graph) */
  int32_t * eptr,    /* element ptr */
  int32_t * eind,    /* element node index */
  int32_t * gofs,    /* offset */
  int32_t * gcolour, /* colour */
  int32_t * epart,   /* part array */
  int32_t * perm,    /* permutation of the elements (size of ne) */
  int32_t * ierr     /* Error code */
)
{
  *ierr = colour_api_call1( *nn, *ne, *mcg, *ncommon, eptr, eind, gofs, gcolour, epart, perm );
}

#endif