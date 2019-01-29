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
 module rcm
 contains

! Finds the next unused node in K_graph and returns it in "node"
 subroutine next_unused_K_graph_node(node, dim_K, K_Visited, no_nodes_K)
 implicit none
 integer, intent(out) :: node
 integer, intent(in) :: dim_K
 logical, dimension(dim_K), intent(in) :: K_visited
 logical, intent(out) :: no_nodes_K
 integer :: i
 do i=1, dim_K
  if (K_visited(i).eqv..FALSE.) then
   node=i
   return
  endif
 enddo
 no_nodes_K=.TRUE.
 end subroutine next_unused_K_graph_node

! Adds the node "node" to R and puts the nodes conected to "node" in Q (queue)
 subroutine add_node(node, R, K_visited, Q, K_graph, K_graph_degree, dimen, mk)
 implicit none
 integer, intent(in) :: node, dimen, mk
 integer, dimension(dimen), intent(inout) :: R
 integer, dimension(2*dimen), intent(inout) :: Q
 logical, dimension(dimen) , intent(inout):: K_visited
 integer, dimension(dimen,2*mk+1),intent(in) :: K_graph
 integer, dimension(dimen), intent(in) :: K_graph_degree
 integer :: i, j, k
 if (K_visited(node).eqv..TRUE.) return
 do i=dimen,1,-1
  if (R(i).eq.0) then
   R(i)=node
   K_visited(node)=.TRUE.
   exit
  endif
 enddo
 do j=1,2*dimen
  if (Q(j).eq.0) then
    Q(j:j+K_graph_degree(node)-1)=K_graph(node,1:K_graph_degree(node))
    do k=j,j+K_graph_degree(node)-1
     if (K_visited(Q(k)).eqv..TRUE.) Q(k:k+K_graph_degree(node))=Q(k+1:k+K_graph_degree(node)+1)
    enddo
    exit
  endif
 enddo
 end subroutine add_node

! Processes an element of the Q (queue) and add it to R if it is not already present
 subroutine process_queue(R, K_Visited, Q, K_graph, K_graph_degree, dimen, mk, no_nodes_Q)
 implicit none
 integer, intent(in) :: dimen, mk
 integer, dimension(dimen), intent(inout) :: R
 integer, dimension(2*dimen), intent(inout) :: Q
 logical, dimension(dimen), intent(inout) :: K_visited
 integer, dimension(dimen,2*mk+1), intent(in) :: K_graph
 integer, dimension(dimen), intent(in) :: K_graph_degree
 logical, intent(out) :: no_nodes_Q
 integer :: i, j
 do i=dimen,1,-1
  if (R(i).eq.Q(1)) then
   Q(:2*dimen-1)=Q(2:2*dimen)
   return
  endif
 enddo
 call add_node (Q(1), R, K_visited, Q, K_graph, K_graph_degree, dimen, mk) 
 Q(:2*dimen-1)=Q(2:2*dimen)
 Q(2*dimen)=0
 if (Q(1).eq.0) then
  no_nodes_Q=.TRUE.
 endif
 end subroutine process_queue

! Receives a node position in the old K matrix and return the position it will have in the new matrix
 subroutine old_to_new_K(i, j, R, dimen)
 integer, intent(in) :: i, dimen
 integer, intent(out) :: j
 integer, dimension(dimen), intent(in) :: R
 integer loop
 do loop=1,dimen
  if (i.eq.R(loop)) then
   j=loop
   exit
  endif
 enddo
 end subroutine old_to_new_K
   
 end module rcm
 