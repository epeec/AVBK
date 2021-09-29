AVBK : A not-so Very Big Kernel
================================

This is a C++ mini-kernel that show cases the AVBP 
gather-computer-scatter paradigm using very simple
loop structure.

There are three types of loop structure

1) Plain old serial loop
2) OpenMP style loop
3) OmpSs style loop (under developmen)

We supply a very simple pre-processor that reads Hip
HDF5 meshes in 2D (triangles only) and creates the
group of cells datastrcuture. In addition we implement
a new colouring algorithms to remove data dependency
among groups and re-ordering of nodes/cells to have
better memory access after colouring.

Compile
=======

We use CMake for compiling the project and provide auto-install
script `external/install.sh` to install external libraries.

OmpSs and OpenMP must be detected automatically. In addtion,
we auto detect ASAN in debug mode to help you with debugging
those tough ones !


Test
====

We supply a simple hip script to generate triangle mesh.

Contact
=======

1) Pavanakumar Mohanamuraly
   Senior Researcher
   ALGO-COOP Team
   CERFACS

2) Gabriel Staffelbach
   Senior Researcher
   ALGO-COOP Team
   CERFACS

This project was supported by the EPEEC project under the EU Horizon 2020 project.

