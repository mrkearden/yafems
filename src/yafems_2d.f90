! This file is part of YAFEMS.
!
! COPYRIGHT (C) 2014-2015 Javier Marcelo Mora <javiermarcelomora@gmail.com>
! YAFEMS is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

 program yafems
! use omp_lib
 use rcm
 implicit none
 include 'med.hf'

! Declarations of variables for reading and parsing
! i, j, k, l, n: Used as loop indexes
! status_file: Receives error status from file opening and reading
! t_anaylis: 1 for plane stress, 2 for plane strain
! n_materials: Number of materials present in the analysis
! n_bc: Number of boundary conditions
! n_loads: Number of loads applied
! n_groups_mat: Number of faces in compound object
! n_eq: Number of equations of final system
! arg: Stores command line argument, ie: the analysis .yaf file
! mesh_file: Stores the name of the MED file
! temp: Temporary character string
! material_groups: Stores names of groups that has different materials in .yaf file
! bc_groups: Stores names of groups where boundary conditions are defined in .yaf file
! load_groups: Stores names of groups where boundary conditions are defined in .yaf file
! mat_def: Stores thickness, E, and Poisson for each material
! loads: Stores load values for each of the previously defined groups. X and Y values in global coordinates
! bc: Stores boundary condition restrictions in X and Y. 1 is fixed, 0 is free
! mat_of_group: Material number for every material_groups
! log_file: Log file name
! result_file: Name for the text result file
! result_med: Name for the MED result file
! huge_num: Large integer for testing purposes
! tima_a, time_b, start: Time variables for elapsed time
 integer :: i, j, k, l, n
 integer :: status_file, t_analysis, n_materials, n_bc, n_loads, n_groups_mat, n_eq
 character(len=200) :: arg, mesh_file, temp, log_file, result_file, result_med
 character(len=20), dimension(:), allocatable :: material_groups,bc_groups, load_groups
 real(kind=8), dimension(:,:), allocatable :: mat_def, loads
 real :: time_a, time_b, start
 integer, dimension(:,:), allocatable :: bc
 integer, dimension(:), allocatable :: mat_of_group
 integer, parameter :: huge_num=100000000

! Variables for MED loading and writing
! cret: Status of MED operation (fail, success)
! fid: file identification internal number
! mtype: type of mesh (structured or non structured)
! stype: Order of calculation sequence of mesh
! nstep: Number of calculation sequence
! atype: Type of coordinates of mesh
! chgt: Change indicator of previous calculation sequence
! tsf: Transformation indicator of previous calculation sequence
! num_groups: Number of groups of each family
! num_fam: Number of families
! num_groups_tot: Total number of groups counting all families
! num_nodes: Total number of nodes
! num_tria3: Total number of TRIA3 elements
! num_seg2: Total number of SEG2 elements
! sdim: Dimension of calculation space
! mdim: Dimension of mesh
! mesh_name: Name of mesh
! family_name: Name of family
! desc: Description of mesh
! dt_unit: Step unit associated with the calculation sequence
! nomcoo: Name of coordinate axis
! unicoo: Units of coordinates
! node_coord: Node coordinates, both x and y
! node_coord_x: X coordinate of nodes
! node_coord_y: Y coordinate of nodes
! s2: Connectivity table of SEG2 entities. Obtained from libmed
! t3: Connectivity table of TRIA3 entities. Obtained from libmed
! seg2: Connectivity table of SEG2 entities, processed for ease of use
! tria3: Connectivity table of TRIA3 entities, processed for ease of use
! group_name: Name of group in MED file
! group_number: Number of group in MED file
! num_fam_node: Number of family to which a node belongs
! num_fam_seg2: Number of family to which a SEG2 belongs
! num_fam_tria3: Number of family to which a TRIA3 belongs
! rounding: Value below which 0.0 is assumed
! field_disp: 'X' and 'Y'. Names of displacement field.
 integer :: cret, fid, mtype,stype, nstep, atype, chgt, tsf, num_groups, num_fam, num_groups_tot
 integer :: num_nodes, num_tria3, num_seg2, sdim, mdim
 character(len=64) :: mesh_name, family_name
 character(len=200) :: desc
 character(len=16) :: dt_unit
 character(len=16) :: nomcoo(2), unicoo(2)
 real(kind=8), allocatable, dimension (:) ::  node_coord, node_coord_x, node_coord_y
 integer, allocatable, dimension (:) :: s2, t3
 integer, allocatable, dimension(:,:) :: seg2, tria3
 character(len=80), allocatable, dimension (:) :: group_name
 integer, allocatable, dimension (:) :: group_number, num_fam_node, num_fam_seg2, num_fam_tria3
 real(kind=8),parameter :: rounding=1E-75_8
 character(len=20) :: field_name

! Variables for FEM calculations
! loop: Used in main loop to obtain K_global
! element_material: Material for each element
! bc_nodes: List of each node and its restrictions
! load_nodes: List of each node and its loads
! D: D material matrix, according to elasticity theory
! B: B matrix according to FEM theory
! LM: LM connectivity matrix.
! H_K: Intermediate matrix: H_K=0.5*Area*D*B
! K_elem: Stiffness matrix of element
! K_global: Global stiffness matrix for all nodes in skyline format
! K_graph: Graph matrix for K_global
! K_graph_degree: Graph node degree vector
! Q: Cuthill-McKee queue vector
! R: Cuthill-McKee result vector
! K_Visited: Vector that says if a node of K_graph has been visited
! no_nodes_Q: .TRUE. if no more nodes are in queue for RCM
! no_nodes_K: .TRUE. if no more nodes needs to be processed in K_graph
! maxa: Location of diagonal elements in K_global
! displacement: Final displacement vector
! length: Length of segment where a distributed load is applied
! area: Area of element
! x1, x2, x3, y1, y2, y3: Coordinates of the nodes of an element
! mk: half bandwith of global stiffness matrix
! mi: skyline of matrix
! column_height: (i-mi). Height of column i
! kn, kl, ku, kh, klt, ic, ki, nd, kk, c, b_g: Used in the solution using the Active Column Method
! stress: Stress result matrix for each element
! t1, t2: Temporary intermediate variables
! temp_field: Temporary array to write the values in the results MED, and other uses
! stress_nodes: Calculated sum stress in nodes from stress in element
! stress_nodes_mean: Calculated mean stress in nodes from stress in element
! total_elem: Used for calculation of stress in nodes
! von_mises: Von Mises stress in element
! von_mises_nodes: Von Mises stress in nodes
 integer :: loop, t_1, t_2, mk, kn, kl, ku, kh, klt, ic, ki, nd, kk
 integer, allocatable, dimension (:) :: element_material, mi, column_height, maxa, Q, R, K_graph_degree
 integer, allocatable, dimension (:,:) :: bc_nodes, LM, K_graph
 real(kind=8), allocatable, dimension (:,:) :: load_nodes, stress, stress_nodes, stress_nodes_mean
 real(kind=8), dimension (3,3) :: D
 real(kind=8), dimension (3,6) :: B
 real(kind=8), allocatable, dimension (:,:,:) :: H_K, K_elem
 real(kind=8) :: length, area, x1, x2, x3, y1, y2, y3, total_elem, c, b_g
 real(kind=8), allocatable, dimension (:) :: K_global, F_global, displacement, temp_field, von_mises, von_mises_nodes
 logical, allocatable, dimension (:) ::  K_Visited
 logical :: no_nodes_Q, no_nodes_K

! START OF PROGRAM
 call cpu_time(start)
 write(*,'(A)') '============================================================================='
 write(*,'(A)') '=                                   YAFEMS                                  ='
 write(*,'(A)') '=                                    v0.4                                   ='
 write(*,'(A)') '= GPL Copyright (C) 2015, Javier Marcelo Mora <javiermarcelomora@gmail.com> ='
 write(*,'(A)') '=          This is free software with absolutely no warranty.               ='
 write(*,'(A)') '=          For details, see the GPL license file, LICENSE.txt               ='
 write(*,'(A)') '============================================================================='
 write(*,*)
! write(*,'(a,i8)')'The number of processors available = ',omp_get_num_procs()
! write(*,'(a,i8)')'The number of threads available    = ',omp_get_max_threads()
! write(*,*)
 write(*,'(A)') 'Begin parsing of files and reading of initial data ...'
 call cpu_time(time_a)

! We read the command line arguments. We only need the name of the YAF file as argument.
 call get_command_argument(1, arg)
 if (LEN_TRIM(arg).eq.0) then
  write (*,*) 'Syntax: yafems input.yaf'
  call efexit(-1)        
 endif
      
! Next step is reading and parsing the input file. The file is opened first.
 open (10, iostat=status_file, file=TRIM(arg),action='read', status='old')
 if (status_file.ne.0) then
  write (*,*) "ERROR opening file '", TRIM(arg), "'"
  call efexit(-1)
 endif

! Parsing the file. Max lines of input file is 1000.
! Read mesh file name
 read(10,'(A)',iostat=status_file) mesh_file
 call test_file_error(status_file)

! Read type of analysis, number of materials, number of faces in compound,
! number of boundary conditions, number of loads
 read(10,*,iostat=status_file) t_analysis, n_materials, n_groups_mat, n_bc, n_loads
 call test_file_error(status_file)

! Allocate the physical properties variables
 allocate(mat_def(n_materials,3))
 allocate(material_groups(n_groups_mat))
 allocate(mat_of_group(n_groups_mat))
 allocate(bc_groups(n_bc))
 allocate(bc(n_bc,2))
 allocate(load_groups(n_loads))
 allocate(loads(n_loads,3))

! Populate the list of materials in mat_def
 do i=1,n_materials
  read(10,*,iostat=status_file) mat_def(i,1), mat_def(i,2), mat_def(i,3)
  call test_file_error(status_file)
 enddo

! Populate the list of groups with their assigned material
 do i=1,n_groups_mat
  read(10,'(A)',iostat=status_file) material_groups(i)
  call test_file_error(status_file)
  read(10,*,iostat=status_file) mat_of_group(i)
  call test_file_error(status_file)
 enddo

! Populate the boundary conditions list
 do i=1,n_bc
  read(10,'(A)',iostat=status_file) bc_groups(i)
  call test_file_error(status_file)
  read(10,*,iostat=status_file) bc(i,1), bc(i,2)
  call test_file_error(status_file)
 enddo

! Populate the load list
 do i=1,n_loads
  read(10,'(A)',iostat=status_file) load_groups(i)
  call test_file_error(status_file)
  read(10,*,iostat=status_file) loads(i,1), loads(i,2), loads(i,3)
  call test_file_error(status_file)
 enddo

! Read END stament
 read(10,'(A)',iostat=status_file) temp
 if (temp.ne.'END') then
  call test_file_error(status_file)
 endif
 close(10)

! The MED file is read using the MED library functions
! Opening of mesh file
 call mfiope(fid,TRIM(ADJUSTL(mesh_file)),MED_ACC_RDONLY, cret)
 call test_mesh_error(cret)

! Read info of mesh
 call mmhmii(fid,1,mesh_name,sdim,mdim,mtype,desc,dt_unit,stype,nstep,atype,nomcoo,unicoo,cret)
 call test_mesh_error(cret)

! Obtain number of nodes in mesh
 num_nodes=0
 call mmhnme(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,chgt,tsf,num_nodes,cret)
 call test_mesh_error(cret)

! Read node coordinates
 allocate(node_coord(num_nodes*sdim))
 call mmhcor(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_FULL_INTERLACE,node_coord,cret)
 call test_mesh_error(cret)

! Create two vectors for the x and y coordinates
 allocate(node_coord_x(num_nodes))
 allocate(node_coord_y(num_nodes))
 do i=1,num_nodes
  node_coord_x(i)=node_coord(2*i-1)
!  if (ABS(node_coord_x(i)).lt.rounding) node_coord_x(i)=0.0
  node_coord_y(i)=node_coord(2*i)
!  if (ABS(node_coord_y(i)).lt.rounding) node_coord_y(i)=0.0
 enddo
 deallocate(node_coord)

! Obtain the number of segments of type SEG2 in mesh
 num_seg2=0
 call mmhnme(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_SEG2,MED_CONNECTIVITY,MED_NODAL,chgt,tsf,num_seg2,cret)
 call test_mesh_error(cret)

! Read SEG2 connectivity data and put them in a matrix
 allocate(s2(num_seg2*2))
 allocate(seg2(num_seg2,2))
 call mmhcyr(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_SEG2,MED_NODAL,MED_FULL_INTERLACE,s2,cret)
 call test_mesh_error(cret)
 do i=1,num_seg2
  seg2(i,1)=s2(2*i-1)
  seg2(i,2)=s2(2*i)
 enddo
 deallocate(s2)

! Obtain the number of triangles of type TRIA3 in mesh
 num_tria3=0
 call mmhnme(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_TRIA3,MED_CONNECTIVITY,MED_NODAL,chgt,tsf,num_tria3,cret)
 call test_mesh_error(cret)

! Read TRIA3 connectivity data and put them in a matrix
 allocate(t3(num_tria3*3))
 allocate(tria3(num_tria3,3))
 call mmhcyr(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_TRIA3,MED_NODAL,MED_FULL_INTERLACE,t3,cret)
 call test_mesh_error(cret)
 do i=1,num_tria3
  tria3(i,1)=t3(3*i-2)
  tria3(i,2)=t3(3*i-1)
  tria3(i,3)=t3(3*i)      
 enddo
 deallocate(t3)

! Mesh group reading needed for loads, materials, and bc
! Number of families; each family can have several groups
 call mfanfa(fid,mesh_name,num_fam,cret)
 call test_mesh_error(cret)

! Read family and group names
! We read the total number first to allocate the group_name array
 num_groups_tot=0
 do i=1,num_fam
  call mfanfg(fid,mesh_name,i,num_groups,cret)
  call test_mesh_error(cret)
  num_groups_tot=num_groups_tot+num_groups
 enddo   
 allocate(group_name(num_groups_tot))
 allocate(group_number(num_groups_tot))

! We repeat for populating the group_name array
 k=0
 do i=1,num_fam
  call mfanfg(fid,mesh_name,i,num_groups,cret)
  call test_mesh_error(cret)
  k=k+num_groups
  if (num_groups.eq.0) cycle
  do j=1,num_groups
   call mfafai(fid,mesh_name,i,family_name,group_number(k-num_groups+j),group_name(k-num_groups+j),cret)
   call test_mesh_error(cret)
  enddo
 enddo
     
! Load number of family to which a node belongs
 allocate(num_fam_node(num_nodes))
 call mmhfnr(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_NODE,MED_NONE,num_fam_node,cret)
 call test_mesh_error(cret)

! Load number of family to which a seg2 belongs
 allocate(num_fam_seg2(num_seg2))
 call mmhfnr(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_SEG2,num_fam_seg2,cret)
 call test_mesh_error(cret)

! Load number of family to which a tria3 belongs
 allocate(num_fam_tria3(num_tria3))
 call mmhfnr(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_TRIA3,num_fam_tria3,cret)
 call test_mesh_error(cret)

! MED reading done. All needed data is now ready to use.
! Mesh file closing
 call mficlo(fid,cret)
 call test_mesh_error(cret)
      
! The following step is the creation of vectors with the material of
! each element, boundary conditions for each node and load for each node.
! Material for each element
 allocate(element_material(num_tria3))
 element_material(:)=0
 do i=1, n_groups_mat
  do j=1, num_groups_tot
   if (TRIM(material_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_tria3
     if (num_fam_tria3(k).eq.group_number(j)) then
      element_material(k)=mat_of_group(i)
     else
      cycle
     endif
    enddo
   else
    cycle
   endif
  enddo
 enddo
 deallocate(num_fam_tria3)

! Boundary conditions for each node. We can have restrictions in nodes or in segments (complete lines)
! We read the segment boundary conditions first
 allocate(bc_nodes(num_nodes,2))
 bc_nodes(:,:)=0
 do i=1, n_bc
  do j=1, num_groups_tot
   if (TRIM(bc_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_seg2
     if (num_fam_seg2(k).eq.group_number(j)) then
      if(bc_nodes(seg2(k,1),1).eq.0) bc_nodes(seg2(k,1),1)=bc(i,1)
      if(bc_nodes(seg2(k,1),2).eq.0) bc_nodes(seg2(k,1),2)=bc(i,2)
      if(bc_nodes(seg2(k,2),1).eq.0) bc_nodes(seg2(k,2),1)=bc(i,1)
      if(bc_nodes(seg2(k,2),2).eq.0) bc_nodes(seg2(k,2),2)=bc(i,2)
     else
      cycle
     endif
    enddo
   else
    cycle
   endif
  enddo
 enddo

! We read the node boundary conditions last
 do i=1, n_bc
  do j=1, num_groups_tot
   if (TRIM(bc_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_nodes
     if (num_fam_node(k).eq.group_number(j)) then
      if(bc_nodes(k,1).eq.0)bc_nodes(k,1)=bc(i,1)
      if(bc_nodes(k,2).eq.0)bc_nodes(k,2)=bc(i,2)
     else
      cycle
     endif
    enddo
   else
    cycle
   endif
  enddo
 enddo

! The bc_nodes matrix is transformed to the form explained in Finite Element Procedures
! by Bathe, section 12.2.2.1, matrix ID. Instead of 0 for fixed degrees of freedom, we use
! a large number so we can calculate mi properly, as mi depends on finding the minimum degree
! of freedom in LM and using 0 would invalidate the result.
 n_eq=0
 do i=1,num_nodes
  do j=1,2
   if (bc_nodes(i,j).eq.0) then
    n_eq=n_eq+1
    bc_nodes(i,j)=n_eq
   else
    bc_nodes(i,j)=huge_num
   endif 
  enddo
 enddo

! Loads for each node. We can have loads in nodes or in segments (complete lines)
! We read the segment loads first
 allocate(load_nodes(num_nodes,2))
 load_nodes(:,:)=0
 do i=1, n_loads
  do j=1, num_groups_tot
   if (TRIM(load_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_seg2
     if (num_fam_seg2(k).eq.group_number(j)) then
      length=sqrt((node_coord_x(seg2(k,1))-node_coord_x(seg2(k,2)))**2+(node_coord_y(seg2(k,1))-&
   &  node_coord_y(seg2(k,2)))**2)
      load_nodes(seg2(k,1),1)=load_nodes(seg2(k,1),1)+0.5*length*mat_def(nint(loads(i,3)),1)*loads(i,1)
      load_nodes(seg2(k,1),2)=load_nodes(seg2(k,1),2)+0.5*length*mat_def(nint(loads(i,3)),1)*loads(i,2)
      load_nodes(seg2(k,2),1)=load_nodes(seg2(k,2),1)+0.5*length*mat_def(nint(loads(i,3)),1)*loads(i,1)
      load_nodes(seg2(k,2),2)=load_nodes(seg2(k,2),2)+0.5*length*mat_def(nint(loads(i,3)),1)*loads(i,2)
     else
      cycle
     endif
    enddo
   else
    cycle
   endif
  enddo
 enddo
 deallocate(num_fam_seg2)

! We read the node loads last
 do i=1, n_loads
  do j=1, num_groups_tot
   if (TRIM(load_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_nodes
     if (num_fam_node(k).eq.group_number(j)) then
      load_nodes(k,1)=load_nodes(k,1)+loads(i,1)
      load_nodes(k,2)=load_nodes(k,2)+loads(i,2)
     else
      cycle
     endif
    enddo
   else
    cycle
   endif
  enddo
 enddo  
 deallocate(num_fam_node)

! F_global needs to be obtained from load_nodes
 allocate (F_global(n_eq))
 loop=1
 do i=1,num_nodes
  do j=1,2
   if (bc_nodes(i,j).eq.huge_num) then
    cycle
   else
    F_global(loop)=load_nodes(i,j)
    loop=loop+1
   endif
  enddo
 enddo
 
! In the following section we obtain matrices D, B, K_elem for each element and K_elem
! is added to the complete K_total matrix
 allocate(H_K(3,6,num_tria3))
 H_K(:,:,:)=0.0
 allocate(K_elem(6,6,num_tria3))
 K_elem(:,:,:)=0.0       
 allocate(LM(num_tria3,6))
 LM(:,:)=0

 call cpu_time(time_b)
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)
 write(*,'(A)') 'Creating the individual element matrices...'
 call cpu_time(time_a)

 do loop=1,num_tria3

! Material matrix D is obtained, depending on the analysis type defined.
  if(t_analysis.eq.1) then
   D(1,1)=(mat_def(element_material(loop),2)/(1-(mat_def(element_material(loop),3))**2))
   if (ABS(D(1,1)).lt.rounding) D(1,1)=0.0
   D(2,2)=D(1,1)
   D(1,2)=(mat_def(element_material(loop),2)/(1-(mat_def(element_material(loop),3))**2))*mat_def(element_material(loop),3)
!   if (ABS(D(1,2)).lt.rounding) D(1,2)=0.0
   D(2,1)=D(1,2)
   D(1,3)=0.0
   D(2,3)=0.0
   D(3,1)=0.0
   D(3,2)=0.0
   D(3,3)=(mat_def(element_material(loop),2)/(1-(mat_def(element_material(loop),3))**2))*&
   & ((1-mat_def(element_material(loop),3))/2.0)
!   if (ABS(D(3,3)).lt.rounding) D(3,3)=0.0
  else if (t_analysis.eq.2) then
   D(1,1)=(mat_def(element_material(loop),2)/((1+(mat_def(element_material(loop),3)))*&
   & (1-2*(mat_def(element_material(loop),3)))))*(1-mat_def(element_material(loop),3))
!   if (ABS(D(1,1)).lt.rounding) D(1,1)=0.0
   D(2,2)=D(1,1)
   D(1,2)=(mat_def(element_material(loop),2)/((1+(mat_def(element_material(loop),3)))*&
   & (1-2*(mat_def(element_material(loop),3)))))*mat_def(element_material(loop),3)
!   if (ABS(D(1,2)).lt.rounding) D(1,2)=0.0
   D(2,1)=D(1,2)
   D(1,3)=0.0
   D(2,3)=0.0
   D(3,1)=0.0
   D(3,2)=0.0
   D(3,3)=(mat_def(element_material(loop),2)/((1+(mat_def(element_material(loop),3)))*&
   & (1-2*(mat_def(element_material(loop),3)))))*(0.5-mat_def(element_material(loop),3))
!   if (ABS(D(3,3)).lt.rounding) D(3,3)=0.0
  else
   write(*,*) 'ERROR: type of analysis must be 1 or 2'
   call efexit(-1) 
  endif

! We obtain the B matrix. The area of the element is needed too, but we don't
! multiply B times (1/2A) in order to save calculation time later.
  x1=node_coord_x(tria3(loop,1))
  x2=node_coord_x(tria3(loop,2))
  x3=node_coord_x(tria3(loop,3))
  y1=node_coord_y(tria3(loop,1))
  y2=node_coord_y(tria3(loop,2))
  y3=node_coord_y(tria3(loop,3))
  area=ABS(0.5*(y1*x3-y3*x1+y3*x2-y2*x3+y2*x1-y1*x2))
  B(:,:)=0.0
  B(1,1)=y2-y3
  B(1,3)=y3-y1
  B(1,5)=y1-y2
  B(2,2)=x3-x2
  B(2,4)=x1-x3
  B(2,6)=x2-x1
  B(3,1)=x3-x2
  B(3,2)=y2-y3
  B(3,3)=x1-x3
  B(3,4)=y3-y1
  B(3,5)=x2-x1
  B(3,6)=y1-y2
!  do i=1,3
!   do j=1,6
!    if (ABS(B(i,j)).lt.rounding) B(i,j)=0.0
!   enddo
!  enddo

! We calculate a matrix H_K(3,6,num_tria3) that is going to serve us twice. Once for the K_elem matrix
! and for stress calculations later
  do i=1,3
   do j=1,6
    do k=1,3
     H_K(i,j,loop)=H_K(i,j,loop)+D(i,k)*(0.5/area)*B(k,j)
    enddo
   enddo
  enddo

! Next step is calculating the K matrix for each element. K_elem=0.5*thickness*B'*H_K
  do i=1,6
   do j=i,6
    do k=1,3
     K_elem(i,j,loop)=K_elem(i,j,loop)+0.5*(mat_def(element_material(loop),1))*B(k,i)*H_K(k,j,loop)
    enddo
    K_elem(j,i,loop)=K_elem(i,j,loop)
   enddo
  enddo

! Now we calculate LM vector for each element. See Finite Element Procedures by Bathe 12.2.2.1
  do i=1,6
   if ((i.eq.1).or.(i.eq.2)) LM(loop,i)=bc_nodes(tria3(loop,1),i)
   if ((i.eq.3).or.(i.eq.4)) LM(loop,i)=bc_nodes(tria3(loop,2),i-2)
   if ((i.eq.5).or.(i.eq.6)) LM(loop,i)=bc_nodes(tria3(loop,3),i-4)
  enddo
 enddo
! END OF ELEMENT LOOP
 call cpu_time(time_b)
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)
 deallocate(element_material)
 call cpu_time(time_a)
 write(*,'(A)') 'Creation of the global stiffness matrix...'

! Global stiffness matrix calculation in skyline format follows.

! mi defines the skyline of the matrix, and (i-mi) is the column height of the i-th column of K_global
 allocate(mi(n_eq))
 allocate(column_height(n_eq))
 do i=1,n_eq
  mi(i)=huge_num
  do j=1,num_tria3
   do k=1,6
    if (LM(j,k).ne.i) then
     cycle
    else
     mi(i)=MIN(mi(i), LM(j,1), LM(j,2), LM(j,3), LM(j,4), LM(j,5), LM(j,6))
    endif
    exit
   enddo
  enddo
 enddo

! Creation of column height matrix
 do i=1,n_eq
  column_height(i)=i-mi(i)
 enddo

! mk: Half bandwith of K_global
 mk=MAXVAL(column_height)

! maxa is created
 allocate(maxa(n_eq+1))
 maxa(1)=1
 do i=2,n_eq+1
  maxa(i)=SUM(column_height(1:i-1))+i
 enddo

! K_global population
 allocate(K_global(n_eq+SUM(column_height)))
 K_global(:)=0.0
 do loop=1,num_tria3
  do i=1,6
   if (LM(loop,i).eq.huge_num) then
    cycle
   endif
   do j=1,6
    if (LM(loop,j).eq.huge_num) then
     cycle
    else if (LM(loop,i).gt.LM(loop,j)) then
     cycle
    else
      K_global(maxa(LM(loop,j))+LM(loop,j)-LM(loop,i))=K_global(maxa(LM(loop,j))+LM(loop,j)-LM(loop,i))+K_elem(i,j,loop)
    endif
   enddo
  enddo
 enddo
 call cpu_time(time_b)
 write(*,'(A,i10)') 'K_global half-bandwidth = ', mk
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)

! Beginning of Cuthill-McKee algorithm.
! First step is creating the graph from the K_global matrix. K_global is what is called the adjacency matrix
! of the graph we are creating.
! A width of 200 is given to K_graph. This dimension should be enough, but in case a seg fault occurs,
! it should be increased
 call cpu_time(time_a)
 write(*,'(A)') 'Optimizing K_global...'
 allocate(K_graph(n_eq,200))
 K_graph(:,:)=0
 do i=1,n_eq
  k=1
  do j=MAX(i-mk, 1),MIN(i+mk, n_eq)
   if (i.le.j) then
    if ((i.ne.j).AND.(mi(j).le.i).AND.(K_global(maxa(j)+j-i).ne.0.0)) then
      K_graph(i,k)=j
      k=k+1
      cycle
    endif
   cycle
   endif
   if ((i.ne.j).AND.(mi(i).le.j).AND.(K_global(maxa(i)+i-j).ne.0.0)) then
     K_graph(i,k)=j
     k=k+1
   endif
  enddo
 enddo
 deallocate(K_global)

! We obtain the degree of each node of the graph and store them in K_graph_degree
 allocate(K_graph_degree(n_eq))
 do i=1,n_eq
  k=0
  do j=1,2*mk+1
   if (K_graph(i,j).ne.0)then
    k=k+1
   else
    K_graph_degree(i)=k
    exit
   endif
  enddo
 enddo

! Now that the graph for K_global is calculated, the Cuthill-McKee algorithm is performed
! BEWARE: The size of the queue can't be known, so a big enough size needs to be defined.
! If YAFEMS halts with a segmentation fault in this stage or runs for a long period of time,
! the dimension of Q needs to be changed and the rcm.f90 module needs to be modified to accomodate this
 allocate(R(n_eq))
 allocate(Q(2*n_eq))
 allocate(K_visited(n_eq))
 no_nodes_K=.FALSE.
 no_nodes_Q=.FALSE.
 R(:)=0
 K_Visited(:)=.FALSE.
 Q(:)=0
 do
  call next_unused_K_graph_node(i, n_eq, K_Visited, no_nodes_K)
  if (no_nodes_K.eqv..TRUE.) then
   exit
  else
   call add_node(i, R, K_Visited, Q, K_graph, K_graph_degree, n_eq, mk)
   do
    call process_queue(R, K_Visited, Q, K_graph, K_graph_degree, n_eq, mk, no_nodes_Q)
    if (no_nodes_Q.eqv..TRUE.) exit
   enddo
  endif
 enddo
 deallocate(Q)
 deallocate(K_visited)

! Now K_global needs to be rebuilt.
! First is rebuilding LM with the new node labeling
 do i=1,num_tria3
  do j=1,6
   if (LM(i,j).eq.0) exit
   call old_to_new_K(LM(i,j), LM(i,j), R, n_eq)
  enddo
 enddo

! New mi calculation
 do i=1,n_eq
  mi(i)=huge_num
  do j=1,num_tria3
   do k=1,6
    if (LM(j,k).ne.i) then
     cycle
    else
     mi(i)=MIN(mi(i), LM(j,1), LM(j,2), LM(j,3), LM(j,4), LM(j,5), LM(j,6))
    endif
    exit    
   enddo
  enddo
 enddo

! Creation of new column height matrix
 do i=1,n_eq
  column_height(i)=i-mi(i)
 enddo

! new mk: Half bandwith of K_global
 mk=MAXVAL(column_height)

! new maxa is created
 maxa(1)=1
 do i=2,n_eq+1
  maxa(i)=SUM(column_height(1:i-1))+i
 enddo
 deallocate(mi)

! Populate the new optimized K_global
 allocate(K_global(n_eq+SUM(column_height)))
 deallocate(column_height)
 K_global(:)=0.0
 do loop=1,num_tria3
  do i=1,6
   if (LM(loop,i).eq.huge_num) then
    cycle
   endif
   do j=1,6
    if (LM(loop,j).eq.huge_num) then
     cycle
    else if (LM(loop,i).gt.LM(loop,j)) then
     cycle
    else
      K_global(maxa(LM(loop,j))+LM(loop,j)-LM(loop,i))=K_global(maxa(LM(loop,j))+LM(loop,j)-LM(loop,i))+K_elem(i,j,loop)
    endif
   enddo
  enddo
 enddo

! Reordering of F_global as K_global has been reorganized.
 allocate(temp_field(n_eq))
 do i=1,n_eq
  call old_to_new_K(i, j, R, n_eq)
  temp_field(j)=F_global(i)
 enddo
 F_global(:)=temp_field(:)
 call cpu_time(time_b)
 write(*,'(A,i10,A)') 'Optimized K_global half-bandwidth = ', mk
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)

! The system K_global*U=F_global is solved using Gauss elimination (Active Column Solution).
! See Finite Element Procedures by Bathe, section 8.2.3
! L*D*L(T) factorization
 call cpu_time(time_a)
 write(*,'(A,I7)') 'Solving K*U=F system with number of equations = ',n_eq
 do n=1,n_eq
  kn=maxa(n) 
  kl=kn+1 
  ku=maxa(n+1)-1
  kh=ku-kl
  if (kh.lt.0) then
   if (k_global(kn).le.0) then
    write(*,*) "Error in system solving"
    call efexit(-1)
   else
    cycle
   endif
  endif
  if (kh.gt.0) then
   k=n-kh 
   ic=0
   klt=ku
   do j=1,kh 
    ic=ic+1 
    klt=klt-1 
    ki=maxa(k)
    nd=maxa(k+1)-ki-1
    if (nd.gt.0) then
     kk=min0(ic,nd)
     c=0.
     do l=1,kk
      c=c+k_global(ki+l)*k_global(klt+l)
     enddo
     k_global(klt)=k_global(klt)-c
    endif
    k=k+1
   enddo
   k=n
   b_g=0.
   do kk=kl,ku
    k=k-1
    ki=maxa(k) 
    c=k_global(kk)/k_global(ki)
    b_g=b_g + c*k_global(kk)
    k_global(kk)=c
   enddo
   k_global(kn)=k_global(kn)-b_g
   if (k_global(kn).le.0) then
    write(*,*) "Error in system solving"
    call efexit(-1)
   endif
  cycle
  endif
  if (kh.eq.0) then
   k=n
   b_g=0.
   do kk=kl,ku
    k=k-1
    ki=maxa(k) 
    c=k_global(kk)/k_global(ki)
    b_g=b_g + c*k_global(kk)
    k_global(kk)=c
   enddo
   k_global(kn)=k_global(kn)-b_g
   if (k_global(kn).le.0) then
    write(*,*) "Error in system solving"
    call efexit(-1)
   endif
  endif
 enddo 

! Right hand side load vector
 do n=1,n_eq 
  kl=maxa(n)+1
  ku=maxa(n+1)-1
  if ((ku-kl).ge.0) then
   k=n
   c=0.
   do kk=kl,ku 
    k=k-1
    c=c+k_global(kk)*f_global(k)
   enddo
    f_global(n)=f_global(n)-c 
  endif
 enddo

! Back substitute
 do n=1,n_eq 
  k=maxa(n)
  f_global(n)=f_global(n)/k_global(k)
 enddo
 if (n_eq.ne.1) then
  n=n_eq
  do l=2,n_eq
   kl=maxa(n)+1
   ku=maxa(n+1)-1
   if ((ku-kl).ge.0) then
    k=n
    do kk=kl,ku
     k=k-1
     f_global(k)=f_global(k)-k_global(kk)*f_global(n) 
    enddo
   endif
    n=n-1
  enddo
 endif
 deallocate(maxa)

! Reorder results to original ordering so the following calculations will be correct
 do i=1,n_eq
  call old_to_new_K(i, j, R, n_eq)
  temp_field(i)=F_global(j)
 enddo
 F_global(:)=temp_field(:)
 deallocate(temp_field)
 call cpu_time(time_b)
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)
 deallocate(R)

! Now F_global contains the displacement results, so we reorder them and put them in the displacement matrix
 allocate(displacement(2*num_nodes))
 loop=1
 do i=1,num_nodes
  do j=1,2
   if (bc_nodes(i,j).eq.huge_num) then
    displacement(2*i+j-2)=0.0
   else
    displacement(2*i+j-2)=F_global(loop)
    loop=loop+1
   endif
  enddo
 enddo
 deallocate(F_global)

! Reaction force calculation
! As K_global was changed during the solving process, we need to recalculate it again
! We also need the complete K matrix, including the fixed boundary condition nodes
 call cpu_time(time_a)
 write(*,'(A)') 'Reaction Force calculation. K reassembly and computation of forces...'
 deallocate(K_global)

! LM for total K matrix
 do loop=1,num_tria3
  LM(loop,1)=2*tria3(loop,1)-1
  LM(loop,2)=2*tria3(loop,1)
  LM(loop,3)=2*tria3(loop,2)-1
  LM(loop,4)=2*tria3(loop,2)
  LM(loop,5)=2*tria3(loop,3)-1
  LM(loop,6)=2*tria3(loop,3)
 enddo

! mi for total K matrix
 allocate(mi(2*num_nodes))
 do i=1,2*num_nodes
  mi(i)=huge_num
  do j=1,num_tria3
   do k=1,6
    if (LM(j,k).ne.i) then
     cycle
    else
     mi(i)=MIN(mi(i), LM(j,1), LM(j,2), LM(j,3), LM(j,4), LM(j,5), LM(j,6))
    endif
    exit
   enddo
  enddo
 enddo

! Creation of column height matrix
 allocate(column_height(2*num_nodes))
 do i=1,2*num_nodes
  column_height(i)=i-mi(i)
 enddo

! maxa is created
 allocate(maxa(2*num_nodes+1))
 maxa(1)=1
 do i=2,2*num_nodes+1
  maxa(i)=SUM(column_height(1:i-1))+i
 enddo

! New K_global of total system
 allocate(K_global(2*num_nodes+SUM(column_height)))
 deallocate(column_height)
 K_global(:)=0.0
 do loop=1,num_tria3
  do i=1,6
   do j=1,6
    if (LM(loop,i).gt.LM(loop,j)) then
     cycle
    else
      K_global(maxa(LM(loop,j))+LM(loop,j)-LM(loop,i))=K_global(maxa(LM(loop,j))+LM(loop,j)-LM(loop,i))+K_elem(i,j,loop)
    endif
   enddo
  enddo
 enddo
 deallocate(LM)
 deallocate(K_elem)

! Force calculation
 allocate(F_global(2*num_nodes))
 F_global(:)=0.0
 do i=1,2*num_nodes
  do j=1,2*num_nodes
   if (i.gt.j) then
    if (j.lt.mi(i)) then
     cycle
    else
     F_global(i)=F_global(i)+K_global(maxa(i)+i-j)*displacement(j)
     cycle
    endif
   endif
   if (i.lt.mi(j)) then
    cycle
   else 
    F_global(i)=F_global(i)+K_global(maxa(j)+j-i)*displacement(j)
   endif
  enddo
 enddo
 deallocate(mi)
 deallocate(maxa)
 deallocate(K_global)
 call cpu_time(time_b)
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)

! Stress calculation
 call cpu_time(time_a)
 write(*,'(A)') 'Stress calculation.'
 allocate(stress(num_tria3,3))
 allocate(stress_nodes(num_nodes,3))
 allocate(stress_nodes_mean(num_nodes,3))
 stress_nodes(:,:)=0.0
 stress_nodes_mean(:,:)=0.0
 do loop=1, num_tria3
  stress(loop,:)=0.0
  do i=1,3
   do j=1,6
    if ((j.eq.1).or.(j.eq.2)) then
     t_1=tria3(loop,1)
     t_2=j-2
    else if ((j.eq.3).or.(j.eq.4)) then
     t_1=tria3(loop,2)
     t_2=j-4
    else if ((j.eq.5).or.(j.eq.6)) then
     t_1=tria3(loop,3)
     t_2=j-6
    endif
    stress(loop,i)=stress(loop,i)+H_K(i,j,loop)*displacement(2*t_1+t_2)
   enddo
  enddo
 enddo
 deallocate(H_K)
 
! Stress in nodes are calculated using the mean of the elements on that node
 do i=1,num_nodes
  total_elem=0
  do j=1,num_tria3
   do k=1,3
    if (tria3(j,k).eq.i) then
     total_elem=total_elem+1
     do l=1,3
      stress_nodes(i,l)=stress_nodes(i,l)+stress(j,l)
     enddo
    endif
   enddo
  enddo
  do l=1,3
   stress_nodes_mean(i,l)=stress_nodes(i,l)/total_elem
  enddo
 enddo
 call cpu_time(time_b)
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)

! A log file is written with several useful data for reference
 call cpu_time(time_a)
 write(*,'(A)') 'Writing of result files...'
 do i=1,200
  if (arg(i:i).eq.'.') then
   temp=arg(1:(i-1))
   exit
  endif
 enddo
 log_file=TRIM(TRIM(temp)//'.log')
 open (11, iostat=status_file, file=log_file, action='write', status='replace')
 if (status_file.ne.0) then
  write (*,*) "ERROR creating log file '", log_file, "'."
  call efexit(-1)
 endif
 write(11,'(A)') '============================================================================='
 write(11,'(A)') '=                                   YAFEMS                                  ='
 write(11,'(A)') '=                                    v0.4                                   ='
 write(11,'(A)') '= GPL Copyright (C) 2015, Javier Marcelo Mora <javiermarcelomora@gmail.com> ='
 write(11,'(A)') '=          This is free software with absolutely no warranty.               ='
 write(11,'(A)') '=          For details, see the GPL license file, LICENSE.txt               ='
 write(11,'(A)') '============================================================================='
 write(11,*)
 write(11,'(A)') '======================================='
 write(11,'(A)') '= GENERAL INFORMATION FROM INPUT FILE ='
 write(11,'(A)') '======================================='
 write(11,'(A,A)') 'Name of mesh file to read: ',TRIM(ADJUSTL(mesh_file))
 write(11,'(A,I1)') 'Analysis type=',t_analysis
 write(11,'(A,I2)') 'Number of materials=',n_materials
 write(11,'(A,I2)') 'Number of boundary conditions=',n_bc
 write(11,'(A,I2)') 'Number of loads=',n_loads
 write(11,'(A,I2)') 'Number of faces in compound=',n_groups_mat
 write(11,*)
 write(11,'(A)')'Materials'
 do i=1,n_materials
  write(11,'(A,I2)') 'Material #',i
  write(11,'(A,F10.3,5X,A,F15.0,5X,A,F10.3)') 'Thickness=',mat_def(i,1),'E=',mat_def(i,2),'Poisson=', mat_def(i,3)
 enddo
 write(11,*)
 write(11,'(A)') 'Groups and materials'
 do i=1,n_groups_mat
  write(11,'(A,I3)') TRIM(ADJUSTL(material_groups(i))),mat_of_group(i)
 enddo
 write(11,*)
 write(11,'(A)') 'Boundary conditions'
 do i=1,n_bc
  write(11,'(A,2X,I2,2X,I2)') TRIM(ADJUSTL(bc_groups(i))), bc(i,1), bc(i,2)
 enddo
 write(11,*)
 write(11,'(A)') 'Loads'
 do i=1,n_loads
  write(11,'(A,F20.3,F20.3)') TRIM(ADJUSTL(load_groups(i))),loads(i,1), loads(i,2)      
 enddo
 deallocate(material_groups,bc_groups,load_groups,bc,mat_of_group)
 write(11,*)
 write(11,'(A)') '======================================='
 write(11,'(A)') '=  GENERAL INFORMATION FROM MED FILE  ='
 write(11,'(A)') '======================================='
 write(11,*) 
 write(11,'(A,A)') 'Name of mesh: ', TRIM(mesh_name)
 write(11,'(A,I6)') 'Number of nodes: ', num_nodes
 write(11,'(A,I6)') 'Number of segments (seg2 entities):', num_seg2
 write(11,'(A,I6)') 'Number of elements (tria3 entities):', num_tria3
 write(11,'(A,I3)') 'Number of families:', num_fam
 write(11,'(A,I3)') 'Total number of groups:', num_groups_tot
 do i=1,num_groups_tot
  write(11,'(A,A,A,I3)') 'Name of group: ', TRIM(group_name(i)),'. Number assigned to group: ', group_number(i)
 enddo
 deallocate(group_name, group_number)
 write(11,'(A)')
 write(11,'(A)') '========================================'
 write(11,'(A)') '=              NODE LIST               ='
 write(11,'(A)') '========================================'
 write(11,'(A)') '| NODE |       X       |       Y       |'
 write(11,'(A)') '========================================'
 do i=1,num_nodes
  write(11,'(A,I6,A,F15.7,A,F15.7,A)') '|',i,'|',node_coord_x(i),'|',node_coord_y(i),'|'
 enddo
 deallocate(node_coord_x, node_coord_y)
 write(11,'(A)') '========================================'
 write(11,'(A)')
 write(11,'(A)') '========================================'
 write(11,'(A)') '=              SEG2 LIST               ='
 write(11,'(A)') '========================================'
 write(11,'(A)') '| NUMB |     NODE 1    |    NODE 2     |'
 write(11,'(A)') '========================================'
 do i=1,num_seg2
  write(11,'(A,I6,A,I15,A,I15,A)') '|',i,'|',seg2(i,1),'|',seg2(i,2),'|'
 enddo
 deallocate(seg2)
 write(11,'(A)') '========================================'
 write(11,'(A)')
 write(11,'(A)') '========================================'
 write(11,'(A)') '=             TRIA3 LIST               ='
 write(11,'(A)') '========================================'
 write(11,'(A)') '| ELEM |  NODE 1  |  NODE 2  |  NODE 3 |'
 write(11,'(A)') '========================================'
 do i=1,num_tria3
  write(11,'(A,I6,A,I10,A,I10,A,I9,A)') '|',i,'|',tria3(i,1),'|',tria3(i,2),'|',tria3(i,3),'|'
 enddo
 deallocate(tria3)
 write(11,'(A)') '========================================'
 write(11,'(A)')
 write(11,'(A)') '========================================'
 write(11,'(A)') '=       BOUNDARY CONDITION LIST        ='
 write(11,'(A)') '========================================'
 write(11,'(A)') '| NODE |     COND X    |    COND Y     |'
 write(11,'(A)') '========================================'
 do i=1,num_nodes
  write(11,'(A,I6,A,I15,A,I15,A)') '|',i,'|',bc_nodes(i,1),'|',bc_nodes(i,2),'|'
 enddo
 deallocate(bc_nodes)
 write(11,'(A)') '========================================'
 write(11,'(A)')
 write(11,'(A)') '========================================'
 write(11,'(A)') '=              LOAD LIST               ='
 write(11,'(A)') '========================================'
 write(11,'(A)') '| NODE |       X       |       Y       |'
 write(11,'(A)') '========================================'
 do i=1,num_nodes
  write(11,'(A,I6,A,F15.2,A,F15.2,A)') '|',i,'|',load_nodes(i,1),'|',load_nodes(i,2),'|'
 enddo
 write(11,'(A)') '========================================'
 close(11)

! Two result files are written. One is a text file, readable with the .result.txt extension and one .rmed MED result
! file. The MED file can be loaded in a graphical post-processing program such as Paraview.
 do i=1,200
  if (arg(i:i).eq.'.') then
   temp=arg(1:(i-1))
   exit
  endif
 enddo
 result_file=TRIM(TRIM(temp)//'.result.txt')
 open (12, iostat=status_file, file=result_file, action='write', status='replace')
 if (status_file.ne.0) then
  write (*,*) "ERROR creating result file '", result_file, "'."
  call efexit(-1)
 endif
 write(12,'(A)') '============================================================================='
 write(12,'(A)') '=                                   YAFEMS                                  ='
 write(12,'(A)') '=                                    v0.4                                   ='
 write(12,'(A)') '= GPL Copyright (C) 2015, Javier Marcelo Mora <javiermarcelomora@gmail.com> ='
 write(12,'(A)') '=          This is free software with absolutely no warranty.               ='
 write(12,'(A)') '=          For details, see the GPL license file, LICENSE.txt               ='
 write(12,'(A)') '============================================================================='
 write(12,*)
 write(12,'(A)') '================================================'
 write(12,'(A)') '=                 DISPLACEMENTS                ='
 write(12,'(A)') '================================================'
 write(12,'(A)') '| NODE |         X         |         Y         |'
 write(12,'(A)') '================================================'
 do i=1,num_nodes
  write(12,'(A,I6,A,E19.11,A,E19.11,A)') '|',i,'|',displacement(2*i-1),'|',displacement(2*i),'|'
 enddo
 write(12,'(A)') '================================================'
 write(12,*)
 write(12,'(A)') '================================================================'
 write(12,'(A)') '=                           STRESSES                           ='
 write(12,'(A)') '================================================================'
 write(12,'(A)') '| ELEM |      SX     |     SY      |     TXY     |  VON MISES  |'
 write(12,'(A)') '================================================================'
 allocate(von_mises(num_tria3))
 do i=1,num_tria3
  von_mises(i)=DSQRT((stress(i,1))**2-(stress(i,1))*(stress(i,2))+(stress(i,2))**2+3*(stress(i,3))**2)
 enddo
 do i=1,num_tria3
  write(12,'(A,I6,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A)') '|',i,'|',stress(i,1),'|',stress(i,2),'|',stress(i,3),'|',&
  & von_mises(i),'|'
 enddo
 write(12,'(A)') '================================================================'
 write(12,*)
 write(12,'(A)') '================================================================'
 write(12,'(A)') '=                      STRESSES ON NODES                       ='
 write(12,'(A)') '================================================================'
 write(12,'(A)') '| NODE |      SX     |     SY      |     TXY     |  VON MISES  |'
 write(12,'(A)') '================================================================'
 allocate(von_mises_nodes(num_nodes))
 do i=1,num_nodes
  von_mises_nodes(i)=DSQRT((stress_nodes_mean(i,1))**2-(stress_nodes_mean(i,1))*(stress_nodes_mean(i,2))&
  &+(stress_nodes_mean(i,2))**2+3*(stress_nodes_mean(i,3))**2)
 enddo
 do i=1,num_nodes
  write(12,'(A,I6,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A)') '|',i,'|',stress_nodes_mean(i,1),'|',stress_nodes_mean(i,2),&
  & '|',stress_nodes_mean(i,3),'|',von_mises_nodes(i),'|'
 enddo
 write(12,'(A)') '================================================================'
 write(12,*)
 write(12,'(A)') '================================================'
 write(12,'(A)') '=               FORCES ON NODES                ='
 write(12,'(A)') '================================================'
 write(12,'(A)') '| NODE |         X         |         Y         |'
 write(12,'(A)') '================================================'
 do i=1,num_nodes
  write(12,'(A,I6,A,E19.11,A,E19.11,A)') '|',i,'|',F_global(2*i-1),'|',F_global(2*i),'|'
 enddo
 write(12,'(A)') '================================================'
 close(12)

! A MED file is created and results are written to it.
! Creation of mesh result file
 do i=1,200
  if (arg(i:i).eq.'.') then
   temp=arg(1:(i-1))
   exit
  endif
 enddo 
 result_med=TRIM(TRIM(temp)//'.rmed')
 
! If result_med file already exists, delete it
  open(unit=15, iostat=status_file, file=TRIM(result_med), status='old')
  if (status_file.eq.0) close(15, status='delete')
! endif

! The original MED file is copied with the rmed extension
 call system("cp " // trim(mesh_file) // " " // trim(result_med)) 

! The following routine copies the files using POSIX libs, which I can't test
! as I'm developing in Linux with gfortran

! call PXFLINK (mesh_file, len_trim (mesh_file), result_med, len_trim (result_med), status_file)
! if (status_file.NE.0) then
!  write(*,*) "ERROR: Can't copy med file to rmed"
!  call efexit(-1) 
! endif 
! call PXFUNLINK (mesh_file, len_trim (mesh_file), status_file)
! if (status_file.NE.0) then
!  write(*,*) "ERROR: Can't copy med file to rmed"
!  call efexit(-1) 
! endif 
 
 call mfiope(fid,result_med,MED_ACC_RDWR, cret)
! This line doesn't apply. See the comment 5 lines after this one
! call mfiope(fid,result_med,MED_ACC_CREAT, cret)
 call test_rmesh_error(cret)

! Link result file mesh with original mesh
! This doesn't work for whatever reason, so the solution is copying the original mesh file
! and writing the results there. This has the advantage that mesh and results are on the same
! file.
!
! call mlnliw (fid,mesh_name,mesh_file, cret)
! call test_rmesh_error(cret)

! Create total displacement result field in file
 field_name='TOTAL_DISPLACEMENT'
 call mfdcre(fid,'TOTAL_DISPLACEMENT',MED_FLOAT64,2,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write total displacement into field
 call mfdrvw(fid, 'TOTAL_DISPLACEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,displacement,cret)
 call test_rmesh_error(cret)

! Create X-displacement result field in file
 field_name='UX'
 call mfdcre(fid,'X_DISPLACEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write X-displacement into field
 allocate(temp_field(num_nodes))
 do i=1,num_nodes
  temp_field(i)=displacement(2*i-1)
 enddo
 call mfdrvw(fid, 'X_DISPLACEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create Y-displacement result field in file
 field_name='UY'
 call mfdcre(fid,'Y_DISPLACEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Y-displacement into field
 do i=1,num_nodes
  temp_field(i)=displacement(2*i)
 enddo
 call mfdrvw(fid, 'Y_DISPLACEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create force result field in file
 field_name='TOTAL_FORCE'
 call mfdcre(fid,'TOTAL_FORCE',MED_FLOAT64,2,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write total forces into field
 call mfdrvw(fid, 'TOTAL_FORCE', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT, num_nodes,F_global,cret)
 call test_rmesh_error(cret)

! Create Force-X result field in file
 field_name='FX'
 call mfdcre(fid,'FORCE_X',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Force-X into field
 do i=1,num_nodes
  temp_field(i)=F_global(2*i-1)
 enddo
 call mfdrvw(fid, 'FORCE_X', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create Force-Y result field in file
 field_name='FY'
 call mfdcre(fid,'FORCE_Y',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Force-Y into field
 do i=1,num_nodes
  temp_field(i)=F_global(2*i)
 enddo
 call mfdrvw(fid, 'FORCE_Y', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)
 deallocate(temp_field)

! Create SX stress result field in file
 field_name='SX_ELEMENT'
 call mfdcre(fid,'SX_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SX stress into field
 allocate(temp_field(num_tria3))
 do i=1,num_tria3
  temp_field(i)=stress(i,1)
 enddo
 call mfdrvw(fid,'SX_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TRIA3, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tria3,temp_field,cret)
 call test_rmesh_error(cret)

! Create SY stress result field in file
 field_name='SY_ELEMENT'
 call mfdcre(fid,'SY_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SY stress into field
 do i=1,num_tria3
  temp_field(i)=stress(i,2)
 enddo
 call mfdrvw(fid,'SY_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TRIA3, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tria3,temp_field,cret)
 call test_rmesh_error(cret)

! Create TXY stress result field in file
 field_name='TXY_ELEMENT'
 call mfdcre(fid,'TXY_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TXY stress into field
 do i=1,num_tria3
  temp_field(i)=stress(i,3)
 enddo
 call mfdrvw(fid,'TXY_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TRIA3, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tria3,temp_field,cret)
 call test_rmesh_error(cret)

! Create Von Mises stress result field in file
 field_name='Von_Mises_ELEMENT'
 call mfdcre(fid,'Von_Mises_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Von_Mises stress into field
 call mfdrvw(fid,'Von_Mises_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TRIA3, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tria3,von_mises,cret)
 call test_rmesh_error(cret)

! Create SX stress result field in file
 deallocate(temp_field)
 allocate(temp_field(num_nodes))
 field_name='SX_NODES'
 call mfdcre(fid,'SX_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SX stress into field
 do i=1,num_nodes
  temp_field(i)=stress_nodes_mean(i,1)
 enddo
 call mfdrvw(fid,'SX_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create SY stress result field in file
 field_name='SY_NODES'
 call mfdcre(fid,'SY_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SY stress into field
 do i=1,num_nodes
  temp_field(i)=stress_nodes_mean(i,2)
 enddo
 call mfdrvw(fid,'SY_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create TXY stress result field in file
 field_name='TXY_NODES'
 call mfdcre(fid,'TXY_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TXY stress into field
 do i=1,num_nodes
  temp_field(i)=stress_nodes_mean(i,3)
 enddo
 call mfdrvw(fid,'TXY_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create Von Mises stress result field in file
 field_name='Von_Mises_NODES'
 call mfdcre(fid,'Von_Mises_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Von_Mises stress into field
 call mfdrvw(fid,'Von_Mises_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,von_mises_nodes,cret)
 call test_rmesh_error(cret)

! Mesh file closing
 call mficlo(fid,cret)
 call test_rmesh_error(cret)
 deallocate(stress,temp_field,displacement,F_global)
 call cpu_time(time_b)
 write(*,'(A,f10.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)
 call cpu_time(time_b)
 write(*,'(A,f10.3,A)') 'TOTAL COMPUTATION TIME= ',time_b-start,' seconds'
 write(*,*)

! Exit program
 call efexit(cret)
 end program yafems

! Some error testing subroutines used
 subroutine test_file_error(status_file)
  integer status_file
  if (status_file.ne.0) then
   write(*,*) 'ERROR: Bad syntax of file. See template.yaf'
   call efexit(-1) 
  endif
 end

 subroutine test_mesh_error(cret)
  integer cret
   if (cret.ne.0) then
    write(*,*) 'ERROR reading the mesh file.'
    call efexit(-1) 
   endif
 end

 subroutine test_rmesh_error(cret)
  integer cret
   if (cret.ne.0) then
    write(*,*) 'ERROR writing result mesh file.'
    call efexit(-1) 
   endif
 end

