YAFEMS 0.4 (Yet Another Finite Element Method Solver)
=====================================================

YAFEMS is a FEM solver programmed in Fortran program that reads a MED mesh file
produced by Salome (http://www.salome-platform.org/) with certain groups created
into the mesh, and with the help of an input text file (.yaf), performs a plane
stress or plane strain (yafems_2d) or 3D analysis (yafems_3d) and creates results
in plain text format and in MED format.

This way a complete analysis can be carried inside Salome, where a mesh file can
be created, exported and read by YAFEMS with the help of a .yaf input file and 
post-processed within Salome using ParaView.

For usage information, see the user manual included in the doc directory where 
the program usage is described in detail.

Version 0.4 new features:
- More performance boost due to the implementation of the Cuthill-McKee
algorithm. This algorithm reorders the stiffness matrix matrix and severely
reduces the bandwith. This leads to a huge decrease in time on the solving
stage. If the number of elements is not too big, the solving time can actually
increase a little, but those are the cases that solves in under a second, it
shouldn't be a problem. The piston case sees a reduction in solving time from
almost a minute to 9 seconds.
- Fixed some minor typos here and there.
- I've been trying to implement parallellization using OpenMP, but it seems that
in the loops I've tried it actually takes more time to create the threads,
execute them and then syncing then than the time it takes to do it with a
single thread. This has lead to postpone parallellization for now. Some
commented out code has been left in case I decide to give it a try in the future.
- Windows version released between 0.3 and 0.4.

================================================================================

YAFEMS is free software; you may redistribute it and/or modify it under 
the terms of the GNU General Public License (GPL) as published by the
Free Software Foundation. YAFEMS is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License (GPL) for details. 

================================================================================

You should have obtained this distribution from SourceForge.net:
https://sourceforge.net/projects/yafems

================================================================================
Included in this distribution:

LICENSE.txt  ... GNU GPL license
README.txt   ... this file
compile      ... small script to compile the source

doc/         ... documentation (in ODT and PDF format)
Examples/    ... examples of use
src/         ... Fortran source code 
bin/	     ... Executable files and DLL (only in the Windows version).

To compile and build under Linux, libmed lib is needed. .deb based distros
(Debian, Ubuntu and Linux Mint) have them in their repos, so it shouldn't be
a problem. Note that you'll need the development-header version too, apart from
libmed itself. In my system, there was a problem with the libmed package and a
bit of fiddling was necessary. More details in the manual.

A Fortran compiler is needed. I use gfortran.

For the Windows version, no additional software is needed for YAFEMS to work.
All DLLs and EXEs needed are provided, so no need to compile it from source.

A version of Salome is almost mandatory. It can be freely downloaded
at http://www.salome-platform.org/downloads/current-version , although a free
account is needed to be able to access the files. If you are interested in FEM,
Salome is a great tool for modeling and post-processing.

I also recommend the CAELinux distro (http://www.caelinux.com), that includes
Salome and a lot of FEM and scientific software. If you're a Windows person, a
virtual machine can be created to test CAELinux without having to install it on
your hard disk.


================================================================================

If you have any questions or comments please feel free to contact me.

YAFEMS author

(c) 2014-2015 Javier Marcelo Mora
Ingeniero de Caminos, Canales y Puertos
Universidad Politécnica de Madrid
mailto: javiermarcelomora@gmail.com

================================================================================
