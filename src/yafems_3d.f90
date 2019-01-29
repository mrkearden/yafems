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
! i, j, k, l, m, n: Used as loop indexes
! status_file: Receives error status from file opening and reading
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
! mat_def: E, and Poisson for each material
! loads: Stores load values for each of the previously defined groups. X and Y values in global coordinates
! bc: Stores boundary condition restrictions in X and Y. 1 is fixed, 0 is free
! mat_of_group: Material number for every material_groups
! log_file: Log file name
! result_file: Name for the text result file
! result_med: Name for the MED result file
! huge_num: Large integer for testing purposes
! tima_a, time_b, start: Time variables for elapsed time
 integer :: i, j, k, l, m, n
 integer :: status_file, n_materials, n_bc, n_loads, n_groups_mat, n_eq
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
! num_tetra4: Total number of TETRA4 elements
! sdim: Dimension of calculation space
! mdim: Dimension of mesh
! mesh_name: Name of mesh
! family_name: Name of family
! desc: Description of mesh
! dt_unit: Step unit associated with the calculation sequence
! nomcoo: Name of coordinate axis
! unicoo: Units of coordinates
! node_coord: Node coordinates, all x, y and z
! node_coord_x: X coordinate of nodes
! node_coord_y: Y coordinate of nodes
! node_coord_z: Z coordinate of nodes
! s2: Connectivity table of SEG2 entities. Obtained from libmed
! t3: Connectivity table of TRIA3 entities. Obtained from libmed
! tt4: Connectivity table of TETRA4 entities. Obtained from libmed
! seg2: Connectivity table of SEG2 entities, processed for ease of use
! tria3: Connectivity table of TRIA3 entities, processed for ease of use
! tetra4: Connectivity table of TETRA4 entities, processed for ease of use
! group_name: Name of group in MED file
! group_number: Number of group in MED file
! num_fam_node: Number of family to which a node belongs
! num_fam_seg2: Number of family to which a SEG2 belongs
! num_fam_tria3: Number of family to which a TRIA3 belongs
! num_fam_tetra4: Number of family to which a TRIA3 belongs
! rounding: Value below which 0.0 is assumed
! field_disp: 'X' and 'Y'. Names of displacement field.
 integer :: cret, fid, mtype,stype, nstep, atype, chgt, tsf, num_groups, num_fam, num_groups_tot
 integer :: num_nodes, num_tria3, num_seg2, num_tetra4, sdim, mdim
 character(len=64) :: mesh_name, family_name
 character(len=200) :: desc
 character(len=16) :: dt_unit
 character(len=16) :: nomcoo(2), unicoo(2)
 real(kind=8), allocatable, dimension (:) :: node_coord, node_coord_x, node_coord_y, node_coord_z
 integer, allocatable, dimension (:) :: s2, t3, tt4
 integer, allocatable, dimension(:,:) :: seg2, tria3, tetra4
 character(len=80), allocatable, dimension (:) :: group_name
 integer, allocatable, dimension (:) :: group_number, num_fam_node, num_fam_seg2, num_fam_tria3, num_fam_tetra4
 real(kind=8),parameter :: rounding=1E-75_8
 character(len=20) :: field_name

! Variables for FEM calculations
! loop: Used in main loop to obtain K_global
! element_material: Material for each element
! bc_nodes: List of each node and its restrictions
! load_nodes: List of each node and its loads
! D: D material matrix, according to elasticity theory
! D_fact: Temporary variable containing the factor that multiplies D
! B: B matrix according to FEM theory
! LM: LM connectivity matrix.
! H_K: Intermediate matrix for calculating K_elem
! K_elem: Stiffness matrix of element
! K_global: Global stiffness matrix for all nodes
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
! area: Area of face
! volume: Volume of element
! mk: half bandwith of global stiffness matrix
! mi: skyline of matrix
! column_height: (i-mi). Height of column i
! kn, kl, ku, kh, klt, ic, ki, nd, kk, c, b_g: Real numbers used in the solution using the Active Column Method
! x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4: Coordinates of the nodes of an element
! x_1, x_2, x_3, y_1, y_2, y_3: Temporary variables to calculate area of TRIA3 face
! b_1, b_2, b_3, b_4, c_1, c_2, c_3, c_4, d_1, d_2, d_3, d_4: Elements of B matrix
! stress: Stress result matrix for each element
! t1, t2: Temporary intermediate variables
! temp_field: Temporary array to write the values in the results MED and other uses
! stress_nodes: Calculated sum stress in nodes from stress in element
! stress_nodes_mean: Calculated mean stress in nodes from stress in element
! total_elem: Used for calculation of stress in nodes
! von_mises: Von Mises stress in element
! von_mises_nodes: Von Mises stress in nodes
 integer :: loop, t_1, t_2, mk, kn, kl, ku, kh, klt, ic, ki, nd, kk
 integer, allocatable, dimension (:) :: element_material, column_height, mi, maxa, Q, R, K_graph_degree
 integer, allocatable, dimension (:,:) :: bc_nodes, LM, K_graph
 real(kind=8), allocatable, dimension (:,:) :: load_nodes, stress, stress_nodes, stress_nodes_mean
 real(kind=8), dimension (6,6) :: D
 real(kind=8), dimension (6,12) :: B
 real(kind=8), allocatable, dimension (:,:,:) :: H_K, K_elem
 real(kind=8) :: length, area, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, x_1, x_2, x_3, y_1, y_2, y_3
 real(kind=8) :: b_1, b_2, b_3, b_4, c_1, c_2, c_3, c_4, d_1, d_2, d_3, d_4, total_elem, volume, D_fact, c, b_g
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

! Read number of materials, number of faces in compound,
! number of boundary conditions, number of loads
 read(10,*,iostat=status_file) n_materials, n_groups_mat, n_bc, n_loads
 call test_file_error(status_file)

! Allocate the physical properties variables
 allocate(mat_def(n_materials,2))
 allocate(material_groups(n_groups_mat))
 allocate(mat_of_group(n_groups_mat))
 allocate(bc_groups(n_bc))
 allocate(bc(n_bc,3))
 allocate(load_groups(n_loads))
 allocate(loads(n_loads,4))

! Populate the list of materials in mat_def
 do i=1,n_materials
  read(10,*,iostat=status_file) mat_def(i,1), mat_def(i,2)
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
  read(10,*,iostat=status_file) bc(i,1), bc(i,2), bc(i,3)
  call test_file_error(status_file)
 enddo

! Populate the load list
 do i=1,n_loads
  read(10,'(A)',iostat=status_file) load_groups(i)
  call test_file_error(status_file)
  read(10,*,iostat=status_file) loads(i,1), loads(i,2), loads(i,3), loads(i,4)
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

! Create three vectors for the x, y and z coordinates
 allocate(node_coord_x(num_nodes))
 allocate(node_coord_y(num_nodes))
 allocate(node_coord_z(num_nodes))
 do i=1,num_nodes
  node_coord_x(i)=node_coord(3*i-2)
!  if (ABS(node_coord_x(i)).lt.rounding) node_coord_x(i)=0.0
  node_coord_y(i)=node_coord(3*i-1)
!  if (ABS(node_coord_y(i)).lt.rounding) node_coord_y(i)=0.0
  node_coord_z(i)=node_coord(3*i)
!  if (ABS(node_coord_z(i)).lt.rounding) node_coord_z(i)=0.0
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

! Obtain the number of tetrahedrons of type TETRA4 in mesh
 num_tetra4=0
 call mmhnme(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_TETRA4,MED_CONNECTIVITY,MED_NODAL,chgt,tsf,num_tetra4,cret)
 call test_mesh_error(cret)

! Read TETRA4 connectivity data and put them in a matrix
 allocate(tt4(num_tetra4*4))
 allocate(tetra4(num_tetra4,4))
 call mmhcyr(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_TETRA4,MED_NODAL,MED_FULL_INTERLACE,tt4,cret)
 call test_mesh_error(cret)
 do i=1,num_tetra4
  tetra4(i,1)=tt4(4*i-3)
  tetra4(i,2)=tt4(4*i-2)
  tetra4(i,3)=tt4(4*i-1)
  tetra4(i,4)=tt4(4*i)
 enddo
 deallocate(tt4)

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

! Load number of family to which a tetra4 belongs
 allocate(num_fam_tetra4(num_tetra4))
 call mmhfnr(fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,MED_TETRA4,num_fam_tetra4,cret)
 call test_mesh_error(cret)

! MED reading done. All needed data is now ready to use.
! Mesh file closing
 call mficlo(fid,cret)
 call test_mesh_error(cret)
      
! The following step is the creation of vectors with the material of
! each element, boundary conditions for each node and load for each node.
! Material for each element
 allocate(element_material(num_tetra4))
 element_material(:)=0
 do i=1, n_groups_mat
  do j=1, num_groups_tot
   if (TRIM(material_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_tetra4
     if (num_fam_tetra4(k).eq.group_number(j)) then
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
 deallocate(num_fam_tetra4)

! Boundary conditions for each node. We can have restrictions in nodes, segments (complete lines) and faces
! We read the face boundary conditions first
 allocate(bc_nodes(num_nodes,3))
 bc_nodes(:,:)=0
 do i=1, n_bc
  do j=1, num_groups_tot
   if (TRIM(bc_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_tria3
     if (num_fam_tria3(k).eq.group_number(j)) then
      do l=1,3
       do m=1,3
        if(bc_nodes(tria3(k,l),m).eq.0) bc_nodes(tria3(k,l),m)=bc(i,m)
       enddo
      enddo
     else
      cycle
     endif
    enddo
   else
    cycle
   endif
  enddo
 enddo

! We read the segment boundary conditions second
 do i=1, n_bc
  do j=1, num_groups_tot
   if (TRIM(bc_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_seg2
     if (num_fam_seg2(k).eq.group_number(j)) then
      do l=1,2
       do m=1,3
        if(bc_nodes(seg2(k,l),m).eq.0) bc_nodes(seg2(k,l),m)=bc(i,m)
       enddo
      enddo
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
      do l=1,3
       if(bc_nodes(k,l).eq.0)bc_nodes(k,l)=bc(i,l)
      enddo
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
  do j=1,3
   if (bc_nodes(i,j).eq.0) then
    n_eq=n_eq+1
    bc_nodes(i,j)=n_eq
   else
    bc_nodes(i,j)=huge_num
   endif 
  enddo
 enddo

! Loads for each node. We can have loads in nodes, in segments (complete lines), or faces
! We read the face loads first
 allocate(load_nodes(num_nodes,3))
 load_nodes(:,:)=0
 do i=1, n_loads
  do j=1, num_groups_tot
   if (TRIM(load_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_tria3
     if (num_fam_tria3(k).eq.group_number(j)) then
      x_1=node_coord_x(tria3(k,2))-node_coord_x(tria3(k,1))
      x_2=node_coord_y(tria3(k,2))-node_coord_y(tria3(k,1))
      x_3=node_coord_z(tria3(k,2))-node_coord_z(tria3(k,1))
      y_1=node_coord_x(tria3(k,3))-node_coord_x(tria3(k,1))
      y_2=node_coord_y(tria3(k,3))-node_coord_y(tria3(k,1))
      y_3=node_coord_z(tria3(k,3))-node_coord_z(tria3(k,1))
      area=0.5_8*sqrt((x_2*y_3-x_3*y_2)**2+(x_3*y_1-x_1*y_3)**2+(x_1*y_2-x_2*y_1)**2)
      do l=1,3
       do m=1,3
        load_nodes(tria3(k,l),m)=load_nodes(tria3(k,l),m)+(area/3.0_8)*loads(i,m)
       enddo
      enddo
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

! We read the segment loads second
 do i=1, n_loads
  do j=1, num_groups_tot
   if (TRIM(load_groups(i)).eq.TRIM(group_name(j))) then
    do k=1,num_seg2
     if (num_fam_seg2(k).eq.group_number(j)) then
      length=sqrt((node_coord_x(seg2(k,1))-node_coord_x(seg2(k,2)))**2+(node_coord_y(seg2(k,1))-&
   &  node_coord_y(seg2(k,2)))**2+(node_coord_z(seg2(k,1))-node_coord_z(seg2(k,2)))**2)
      do l=1,2
       do m=1,3
        load_nodes(seg2(k,l),m)=load_nodes(seg2(k,l),m)+load_nodes(seg2(k,l),m)+0.5_8*length*loads(i,m)
       enddo
      enddo
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
      do l=1,3
       load_nodes(k,l)=load_nodes(k,l)+loads(i,l)
      enddo
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
  do j=1,3
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
 allocate(H_K(6,12,num_tetra4))
 H_K(:,:,:)=0.0
 allocate(K_elem(12,12,num_tetra4))
 K_elem(:,:,:)=0.0       
 allocate(LM(num_tetra4,12))
 LM(:,:)=0
 
 call cpu_time(time_b)
 write(*,'(A,f6.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)
 write(*,'(A)') 'Creating the individual element matrices...'
 call cpu_time(time_a)

 do loop=1,num_tetra4

! Material matrix D is obtained
  D_fact=mat_def(element_material(loop),1)/(1+mat_def(element_material(loop),2))
  D(:,:)=0.0
  D(1,1)=D_fact*(1-mat_def(element_material(loop),2))/(1-2*mat_def(element_material(loop),2))
  D(2,2)=D(1,1)
  D(3,3)=D(1,1)
  D(4,4)=0.5_8*D_fact
  D(5,5)=0.5_8*D_fact
  D(6,6)=0.5_8*D_fact
  D(1,2)=D_fact*(mat_def(element_material(loop),2))/(1-2*mat_def(element_material(loop),2))
  D(2,1)=D(1,2)
  D(1,3)=D(1,2)
  D(3,1)=D(1,2)
  D(2,3)=D(1,2)
  D(3,2)=D(1,2)

! We obtain the B matrix. The volume of the element is needed too.
  x1=node_coord_x(tetra4(loop,1))
  x2=node_coord_x(tetra4(loop,2))
  x3=node_coord_x(tetra4(loop,3))
  x4=node_coord_x(tetra4(loop,4))
  y1=node_coord_y(tetra4(loop,1))
  y2=node_coord_y(tetra4(loop,2))
  y3=node_coord_y(tetra4(loop,3))
  y4=node_coord_y(tetra4(loop,4))
  z1=node_coord_z(tetra4(loop,1))
  z2=node_coord_z(tetra4(loop,2))
  z3=node_coord_z(tetra4(loop,3))
  z4=node_coord_z(tetra4(loop,4))
  b_1=(y4-y2)*(z3-z2)-(y3-y2)*(z4-z2)
  b_2=(y3-y1)*(z4-z3)-(y3-y4)*(z1-z3)
  b_3=(y2-y4)*(z1-z4)-(y1-y4)*(z2-z4)
  b_4=(y1-y3)*(z2-z1)-(y1-y2)*(z3-z1)
  c_1=(z4-z2)*(x3-x2)-(z3-z2)*(x4-x2)
  c_2=(z3-z1)*(x4-x3)-(z3-z4)*(x1-x3)
  c_3=(z2-z4)*(x1-x4)-(z1-z4)*(x2-x4)
  c_4=(z1-z3)*(x2-x1)-(z1-z2)*(x3-x1)
  d_1=(x4-x2)*(y3-y2)-(x3-x2)*(y4-y2)
  d_2=(x3-x1)*(y4-y3)-(x3-x4)*(y1-y3)
  d_3=(x2-x4)*(y1-y4)-(x1-x4)*(y2-y4)
  d_4=(x1-x3)*(y2-y1)-(x1-x2)*(y3-y1)
  volume=(1.0_8/6.0_8)*ABS(-(((z3-z4)*y2-((z2-z4)*y3-(z2-z3)*y4))*x1-(((z3-z4)*y1-((z1-z4)*y3-(z1-z3)*y4))*&
  &x2-(((z2-z4)*y1-((z1-z4)*y2-(z1-z2)*y4))*x3-((z2-z3)*y1-((z1-z3)*y2-(z1-z2)*y3))*x4))))
  B(:,:)=0.0
  B(1,1)=b_1
  B(1,4)=b_2
  B(1,7)=b_3
  B(1,10)=b_4
  B(2,2)=c_1
  B(2,5)=c_2
  B(2,8)=c_3
  B(2,11)=c_4
  B(3,3)=d_1
  B(3,6)=d_2
  B(3,9)=d_3
  B(3,12)=d_4
  B(4,1)=c_1
  B(4,2)=b_1
  B(4,4)=c_2
  B(4,5)=b_2
  B(4,7)=c_3
  B(4,8)=b_3
  B(4,10)=c_4
  B(4,11)=b_4
  B(5,2)=d_1
  B(5,3)=c_1
  B(5,5)=d_2
  B(5,6)=c_2
  B(5,8)=d_3
  B(5,9)=c_3
  B(5,11)=d_4
  B(5,12)=c_4
  B(6,1)=d_1
  B(6,3)=b_1
  B(6,4)=d_2
  B(6,6)=b_2
  B(6,7)=d_3
  B(6,9)=b_3
  B(6,10)=d_4
  B(6,12)=b_4

!  do i=1,3
!   do j=1,6
!    if (ABS(B(i,j)).lt.rounding) B(i,j)=0.0
!   enddo
!  enddo

! We calculate a matrix H_K(6,12,num_tetra4) that is going to serve us twice. Once for the K_elem matrix
! and for stress calculations later
  do i=1,6
   do j=1,12
    do k=1,6
     H_K(i,j,loop)=H_K(i,j,loop)+D(i,k)*(1.0_8/(6.0_8*volume))*B(k,j)
    enddo
   enddo
  enddo

! Next step is calculating the K matrix for each element.
  do i=1,12
   do j=i,12
    do k=1,6
     K_elem(i,j,loop)=K_elem(i,j,loop)+(1.0_8/6.0_8)*B(k,i)*H_K(k,j,loop)
    enddo
    K_elem(j,i,loop)=K_elem(i,j,loop)
   enddo
  enddo

! Now we calculate LM vector for each element. See Finite Element Procedures by Bathe 12.2.2.1
  do i=1,12
   if ((i.eq.1).or.(i.eq.2).or.(i.eq.3)) LM(loop,i)=bc_nodes(tetra4(loop,1),i)
   if ((i.eq.4).or.(i.eq.5).or.(i.eq.6)) LM(loop,i)=bc_nodes(tetra4(loop,2),i-3)
   if ((i.eq.7).or.(i.eq.8).or.(i.eq.9)) LM(loop,i)=bc_nodes(tetra4(loop,3),i-6)
   if ((i.eq.10).or.(i.eq.11).or.(i.eq.12)) LM(loop,i)=bc_nodes(tetra4(loop,4),i-9)
  enddo
 enddo
! END OF ELEMENT LOOP
 call cpu_time(time_b)
 write(*,'(A,f6.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
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
  do j=1,num_tetra4
   do k=1,12
    if (LM(j,k).ne.i) then
     cycle
    else
     mi(i)=MIN(mi(i), LM(j,1), LM(j,2), LM(j,3), LM(j,4), LM(j,5), LM(j,6),&
     & LM(j,7), LM(j,8), LM(j,9), LM(j,10), LM(j,11), LM(j,12))
    endif
    exit    
   enddo
  enddo
 enddo

! Creation of column height matrix
 do i=1,n_eq
  column_height(i)=i-mi(i)
 enddo

! Half-bandwith of global stiffness matrix
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
 do loop=1,num_tetra4
  do i=1,12
   if (LM(loop,i).eq.huge_num) then
    cycle
   endif
   do j=1,12
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
 write(*,'(A,f6.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
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

! Now that the graph for K_global is calculated, the Cuhill-McKee algorithm is performed
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
 do i=1,num_tetra4
  do j=1,12
   if (LM(i,j).eq.0) exit
   call old_to_new_K(LM(i,j), LM(i,j), R, n_eq)
  enddo
 enddo

! New mi calculation
 do i=1,n_eq
  mi(i)=huge_num
  do j=1,num_tetra4
   do k=1,12
    if (LM(j,k).ne.i) then
     cycle
    else
     mi(i)=MIN(mi(i), LM(j,1), LM(j,2), LM(j,3), LM(j,4), LM(j,5), LM(j,6),&
     & LM(j,7), LM(j,8), LM(j,9), LM(j,10), LM(j,11), LM(j,12))
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
 do loop=1,num_tetra4
  do i=1,12
   if (LM(loop,i).eq.huge_num) then
    cycle
   endif
   do j=1,12
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
 write(*,'(A,I6)') 'Solving K*U=F system with number of equations = ',n_eq
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
 allocate(displacement(3*num_nodes))
 loop=1
 do i=1,num_nodes
  do j=1,3
   if (bc_nodes(i,j).eq.huge_num) then
    displacement(3*i+j-3)=0.0
   else
    displacement(3*i+j-3)=F_global(loop)
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
 do loop=1,num_tetra4
  LM(loop,1)=3*tetra4(loop,1)-2
  LM(loop,2)=3*tetra4(loop,1)-1
  LM(loop,3)=3*tetra4(loop,1)
  LM(loop,4)=3*tetra4(loop,2)-2
  LM(loop,5)=3*tetra4(loop,2)-1
  LM(loop,6)=3*tetra4(loop,2)
  LM(loop,7)=3*tetra4(loop,3)-2
  LM(loop,8)=3*tetra4(loop,3)-1
  LM(loop,9)=3*tetra4(loop,3)
  LM(loop,10)=3*tetra4(loop,4)-2
  LM(loop,11)=3*tetra4(loop,4)-1
  LM(loop,12)=3*tetra4(loop,4)
 enddo

! mi for total K matrix
 allocate(mi(3*num_nodes))
 do i=1,3*num_nodes
  mi(i)=huge_num
  do j=1,num_tetra4
   do k=1,12
    if (LM(j,k).ne.i) then
     cycle
    else
     mi(i)=MIN(mi(i), LM(j,1), LM(j,2), LM(j,3), LM(j,4), LM(j,5), LM(j,6),&
     & LM(j,7), LM(j,8), LM(j,9), LM(j,10), LM(j,11), LM(j,12))
    endif
    exit    
   enddo
  enddo
 enddo

! Creation of column height matrix
 allocate(column_height(3*num_nodes))
 do i=1,3*num_nodes
  column_height(i)=i-mi(i)
 enddo

! maxa is created
 allocate(maxa(3*num_nodes+1))
 maxa(1)=1
 do i=2,3*num_nodes+1
  maxa(i)=SUM(column_height(1:i-1))+i
 enddo

! New K_global of total system
 allocate(K_global(3*num_nodes+SUM(column_height)))
 deallocate(column_height)
 K_global(:)=0.0
 do loop=1,num_tetra4
  do i=1,12
   do j=1,12
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
 allocate(F_global(3*num_nodes))
 F_global(:)=0.0
 do i=1,3*num_nodes
  do j=1,3*num_nodes
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
 write(*,'(A,f6.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
 write(*,*)


! Stress calculation
 call cpu_time(time_a)
 write(*,'(A)') 'Stress calculation.'
 allocate(stress(num_tetra4,6))
 allocate(stress_nodes(num_nodes,6))
 allocate(stress_nodes_mean(num_nodes,6))
 stress_nodes(:,:)=0.0
 stress_nodes_mean(:,:)=0.0
 do loop=1, num_tetra4
  stress(loop,:)=0.0
  do i=1,6
   do j=1,12
    if ((j.eq.1).or.(j.eq.2).or.(j.eq.3)) then
     t_1=tetra4(loop,1)
     t_2=j-3
    else if ((j.eq.4).or.(j.eq.5).or.(j.eq.6)) then
     t_1=tetra4(loop,2)
     t_2=j-6
    else if ((j.eq.7).or.(j.eq.8).or.(j.eq.9)) then
     t_1=tetra4(loop,3)
     t_2=j-9
    else if ((j.eq.10).or.(j.eq.11).or.(j.eq.12)) then
     t_1=tetra4(loop,4)
     t_2=j-12
    endif
    stress(loop,i)=stress(loop,i)+H_K(i,j,loop)*displacement(3*t_1+t_2)
   enddo
  enddo
 enddo
 deallocate(H_K)

! Stress in nodes are calculated using the mean of the elements on that node
 do i=1,num_nodes
  total_elem=0
  do j=1,num_tetra4
   do k=1,4
    if (tetra4(j,k).eq.i) then
     total_elem=total_elem+1
     do l=1,6
      stress_nodes(i,l)=stress_nodes(i,l)+stress(j,l)
     enddo
    endif
   enddo
  enddo
  do l=1,6
   stress_nodes_mean(i,l)=stress_nodes(i,l)/total_elem
  enddo
 enddo
 call cpu_time(time_b)
 write(*,'(A,f6.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
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
 write(11,'(A,I2)') 'Number of materials=',n_materials
 write(11,'(A,I2)') 'Number of boundary conditions=',n_bc
 write(11,'(A,I2)') 'Number of loads=',n_loads
 write(11,'(A,I2)') 'Number of faces in compound=',n_groups_mat
 write(11,*)
 write(11,'(A)')'Materials'
 do i=1,n_materials
  write(11,'(A,I2)') 'Material #',i
  write(11,'(A,F15.0,5X,A,F10.3)') 'E=',mat_def(i,1),'Poisson=', mat_def(i,2)
 enddo
 write(11,*)
 write(11,'(A)') 'Groups and materials'
 do i=1,n_groups_mat
  write(11,'(A,I3)') TRIM(ADJUSTL(material_groups(i))),mat_of_group(i)
 enddo
 write(11,*)
 write(11,'(A)') 'Boundary conditions'
 do i=1,n_bc
  write(11,'(A,2X,I2,2X,I2,2X,I2)') TRIM(ADJUSTL(bc_groups(i))), bc(i,1), bc(i,2), bc(i,3)
 enddo
 write(11,*)
 write(11,'(A)') 'Loads'
 do i=1,n_loads
  write(11,'(A,F20.3,F20.3,F20.3)') TRIM(ADJUSTL(load_groups(i))),loads(i,1), loads(i,2), loads(i,3)
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
 write(11,'(A,I6)') 'Number of faces (tria3 entities):', num_tria3
 write(11,'(A,I6)') 'Number of elements (tetra4 entities):', num_tetra4
 write(11,'(A,I3)') 'Number of families:', num_fam
 write(11,'(A,I3)') 'Total number of groups:', num_groups_tot
 do i=1,num_groups_tot
  write(11,'(A,A,A,I3)') 'Name of group: ', TRIM(group_name(i)),'. Number assigned to group: ', group_number(i)
 enddo
 deallocate(group_name, group_number)
 write(11,'(A)')
 write(11,'(A)') '========================================================'
 write(11,'(A)') '=                      NODE LIST                       ='
 write(11,'(A)') '========================================================'
 write(11,'(A)') '| NODE |       X       |       Y       |       Z       |'
 write(11,'(A)') '========================================================'
 do i=1,num_nodes
  write(11,'(A,I6,A,F15.7,A,F15.7,A,F15.7,A)') '|',i,'|',node_coord_x(i),'|',node_coord_y(i),'|', node_coord_z(i),'|'
 enddo
 deallocate(node_coord_x, node_coord_y, node_coord_z)
 write(11,'(A)') '========================================================'
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
 write(11,'(A)') '| ELEM |  NODE 1  |  NODE 2  |  NODE3  |'
 write(11,'(A)') '========================================'
 do i=1,num_tria3
  write(11,'(A,I6,A,I10,A,I10,A,I9,A)') '|',i,'|',tria3(i,1),'|',tria3(i,2),'|',tria3(i,3),'|'
 enddo
 deallocate(tria3)
 write(11,'(A)') '========================================'
 write(11,'(A)')
 write(11,'(A)') '=================================================='
 write(11,'(A)') '=                   TETRA4 LIST                  ='
 write(11,'(A)') '=================================================='
 write(11,'(A)') '| ELEM |  NODE 1  |  NODE 2  |  NODE 3 |  NODE 4 |'
 write(11,'(A)') '=================================================='
 do i=1,num_tetra4
  write(11,'(A,I6,A,I10,A,I10,A,I9,A,I9,A)') '|',i,'|',tetra4(i,1),'|',tetra4(i,2),'|',tetra4(i,3),'|',tetra4(i,4),'|'
 enddo
 deallocate(tetra4)
 write(11,'(A)') '=================================================='
 write(11,'(A)')
 write(11,'(A)') '========================================================'
 write(11,'(A)') '=               BOUNDARY CONDITION LIST                ='
 write(11,'(A)') '========================================================'
 write(11,'(A)') '| NODE |     COND X    |    COND Y     |    COND Z     |'
 write(11,'(A)') '========================================================'
 do i=1,num_nodes
  write(11,'(A,I6,A,I15,A,I15,A,I15,A)') '|',i,'|',bc_nodes(i,1),'|',bc_nodes(i,2),'|',bc_nodes(i,3),'|'
 enddo
 deallocate(bc_nodes)
 write(11,'(A)') '========================================================'
 write(11,'(A)')
 write(11,'(A)') '========================================================'
 write(11,'(A)') '=                        LOAD LIST                     ='
 write(11,'(A)') '========================================================'
 write(11,'(A)') '| NODE |       X       |       Y       |       Z       |'
 write(11,'(A)') '========================================================'
 do i=1,num_nodes
  write(11,'(A,I6,A,F15.2,A,F15.2,A,F15.2,A)') '|',i,'|',load_nodes(i,1),'|',load_nodes(i,2),'|',load_nodes(i,3),'|'
 enddo
 write(11,'(A)') '========================================================'
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
 write(12,'(A)') '===================================================================='
 write(12,'(A)') '=                           DISPLACEMENTS                          ='
 write(12,'(A)') '===================================================================='
 write(12,'(A)') '| NODE |         X         |         Y         |         Z         |'
 write(12,'(A)') '===================================================================='
 do i=1,num_nodes
  write(12,'(A,I6,A,E19.11,A,E19.11,A,E19.11,A)') '|',i,'|',displacement(3*i-2),'|',displacement(3*i-1),'|',displacement(3*i),'|'
 enddo
 write(12,'(A)') '===================================================================='
 write(12,*)
 write(12,'(A)') '=========================================================================================================='
 write(12,'(A)') '=                                                STRESSES                                                ='
 write(12,'(A)') '=========================================================================================================='
 write(12,'(A)') '| ELEM |      SX     |     SY      |     SZ      |     TXY     |     TYZ     |     TXZ     |  VON MISES  |'
 write(12,'(A)') '=========================================================================================================='
 allocate(von_mises(num_tetra4))
 do i=1,num_tetra4
  von_mises(i)=DSQRT(0.5_8*((stress(i,1)-stress(i,2))**2+(stress(i,2)-stress(i,3))**2+(stress(i,3)-stress(i,1))**2+&
  & 6.*(stress(i,4)**2+stress(i,5)**2+stress(i,6)**2)))
 enddo
 do i=1,num_tetra4
  write(12,'(A,I6,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A)') '|',i,'|',stress(i,1),'|',stress(i,2),'|',&
   & stress(i,3),'|', stress(i,4),'|',stress(i,5),'|',stress(i,6),'|',von_mises(i),'|'
 enddo
 write(12,'(A)') '=========================================================================================================='
 write(12,*)
 write(12,'(A)') '=========================================================================================================='
 write(12,'(A)') '=                                            STRESSES ON NODES                                           ='
 write(12,'(A)') '=========================================================================================================='
 write(12,'(A)') '| ELEM |      SX     |     SY      |     SZ      |     TXY     |     TYZ     |     TXZ     |  VON MISES  |'
 write(12,'(A)') '=========================================================================================================='
 allocate(von_mises_nodes(num_nodes))
 do i=1,num_nodes
  von_mises_nodes(i)=DSQRT(0.5_8*((stress_nodes_mean(i,1)-stress_nodes_mean(i,2))**2+(stress_nodes_mean(i,2)-&
  & stress_nodes_mean(i,3))**2+(stress_nodes_mean(i,3)-stress_nodes_mean(i,1))**2+6.*(stress_nodes_mean(i,4)**2+&
  & stress_nodes_mean(i,5)**2+stress_nodes_mean(i,6)**2)))
 enddo
 do i=1,num_nodes
  write(12,'(A,I6,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A)') '|',i,'|',stress_nodes_mean(i,1),'|',&
  & stress_nodes_mean(i,2),'|', stress_nodes_mean(i,3),'|',stress_nodes_mean(i,4),'|',stress_nodes_mean(i,5),'|',&
  & stress_nodes_mean(i,6),'|',von_mises_nodes(i),'|'
 enddo
 write(12,'(A)') '=========================================================================================================='
 write(12,*)
 write(12,'(A)') '===================================================================='
 write(12,'(A)') '=                          FORCES ON NODES                         ='
 write(12,'(A)') '===================================================================='
 write(12,'(A)') '| NODE |         X         |         Y         |         Z         |'
 write(12,'(A)') '===================================================================='
 do i=1,num_nodes
  write(12,'(A,I6,A,E19.11,A,E19.11,A,E19.11,A)') '|',i,'|',F_global(3*i-2),'|',F_global(3*i-1),'|',F_global(3*i),'|'
 enddo
 write(12,'(A)') '===================================================================='
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
 call mfdcre(fid,'TOTAL_DISPLACEMENT',MED_FLOAT64,3,field_name,'',dt_unit,mesh_name,cret)
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
  temp_field(i)=displacement(3*i-2)
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
  temp_field(i)=displacement(3*i-1)
 enddo
 call mfdrvw(fid, 'Y_DISPLACEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create Z-displacement result field in file
 field_name='UZ'
 call mfdcre(fid,'Z_DISPLACEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Z-displacement into field
 do i=1,num_nodes
  temp_field(i)=displacement(3*i)
 enddo
 call mfdrvw(fid, 'Z_DISPLACEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create force result field in file
 field_name='TOTAL_FORCE'
 call mfdcre(fid,'TOTAL_FORCE',MED_FLOAT64,3,field_name,'',dt_unit,mesh_name,cret)
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
  temp_field(i)=F_global(3*i-2)
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
  temp_field(i)=F_global(3*i-1)
 enddo
 call mfdrvw(fid, 'FORCE_Y', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create Force-Z result field in file
 field_name='FZ'
 call mfdcre(fid,'FORCE_Z',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Force-Z into field
 do i=1,num_nodes
  temp_field(i)=F_global(3*i)
 enddo
 call mfdrvw(fid, 'FORCE_Z', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)
 deallocate(temp_field)

! Create SX stress result field in file
 field_name='SX_ELEMENT'
 call mfdcre(fid,'SX_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SX stress into field
 allocate(temp_field(num_tetra4))
 do i=1,num_tetra4
  temp_field(i)=stress(i,1)
 enddo
 call mfdrvw(fid,'SX_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TETRA4, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tetra4,temp_field,cret)
 call test_rmesh_error(cret)

! Create SY stress result field in file
 field_name='SY_ELEMENT'
 call mfdcre(fid,'SY_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SY stress into field
 do i=1,num_tetra4
  temp_field(i)=stress(i,2)
 enddo
 call mfdrvw(fid,'SY_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TETRA4, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tetra4,temp_field,cret)
 call test_rmesh_error(cret)

! Create SZ stress result field in file
 field_name='SZ_ELEMENT'
 call mfdcre(fid,'SZ_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SZ stress into field
 do i=1,num_tetra4
  temp_field(i)=stress(i,3)
 enddo
 call mfdrvw(fid,'SZ_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TETRA4, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tetra4,temp_field,cret)
 call test_rmesh_error(cret)

! Create TXY stress result field in file
 field_name='TXY_ELEMENT'
 call mfdcre(fid,'TXY_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TXY stress into field
 do i=1,num_tetra4
  temp_field(i)=stress(i,4)
 enddo
 call mfdrvw(fid,'TXY_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TETRA4, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tetra4,temp_field,cret)
 call test_rmesh_error(cret)

! Create TYZ stress result field in file
 field_name='TYZ_ELEMENT'
 call mfdcre(fid,'TYZ_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TYZ stress into field
 do i=1,num_tetra4
  temp_field(i)=stress(i,5)
 enddo
 call mfdrvw(fid,'TYZ_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TETRA4, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tetra4,temp_field,cret)
 call test_rmesh_error(cret)

! Create TXZ stress result field in file
 field_name='TXZ_ELEMENT'
 call mfdcre(fid,'TXZ_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TXZ stress into field
 do i=1,num_tetra4
  temp_field(i)=stress(i,6)
 enddo
 call mfdrvw(fid,'TXZ_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TETRA4, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tetra4,temp_field,cret)
 call test_rmesh_error(cret)

! Create Von Mises stress result field in file
 field_name='Von_Mises_ELEMENT'
 call mfdcre(fid,'Von_Mises_ELEMENT',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write Von_Mises stress into field
 call mfdrvw(fid,'Von_Mises_ELEMENT', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_CELL,MED_TETRA4, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_tetra4,von_mises,cret)
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

! Create SZ stress result field in file
 field_name='SZ_NODES'
 call mfdcre(fid,'SZ_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write SZ stress into field
 do i=1,num_nodes
  temp_field(i)=stress_nodes_mean(i,3)
 enddo
 call mfdrvw(fid,'SZ_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create TXY stress result field in file
 field_name='TXY_NODES'
 call mfdcre(fid,'TXY_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TXY stress into field
 do i=1,num_nodes
  temp_field(i)=stress_nodes_mean(i,4)
 enddo
 call mfdrvw(fid,'TXY_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create TYZ stress result field in file
 field_name='TYZ_NODES'
 call mfdcre(fid,'TYZ_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TYZ stress into field
 do i=1,num_nodes
  temp_field(i)=stress_nodes_mean(i,5)
 enddo
 call mfdrvw(fid,'TYZ_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
 & MED_ALL_CONSTITUENT,num_nodes,temp_field,cret)
 call test_rmesh_error(cret)

! Create TXZ stress result field in file
 field_name='TXZ_NODES'
 call mfdcre(fid,'TXZ_NODES',MED_FLOAT64,1,field_name,'',dt_unit,mesh_name,cret)
 call test_rmesh_error(cret)

! Write TXZ stress into field
 do i=1,num_nodes
  temp_field(i)=stress_nodes_mean(i,6)
 enddo
 call mfdrvw(fid,'TXZ_NODES', MED_NO_DT, MED_NO_IT,MED_NO_DT, MED_NODE,MED_NONE, MED_FULL_INTERLACE, &
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
 write(*,'(A,f6.3,A)') 'Done. Time used= ',time_b-time_a,' seconds'
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

