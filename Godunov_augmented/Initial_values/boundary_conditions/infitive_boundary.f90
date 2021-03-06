module infinitive_boundary
! Infinitive boundary constion.

use solution
! 'solution.f90'

implicit none

public infinitive_boundary_cell, infinitive_boundary_HS


contains

!******************************************************************************************
! Treatment of boundary
subroutine infinitive_boundary_HS

implicit none

integer :: i,k
   
do i=-ighost+1,0   
 d_rho(i)=0.0d0 
 d_v(i)=0.0d0
 d_p(i)=0.0d0 
enddo
   
do i=1,ighost
 d_rho(cells_number+i)=0.0d0  
 d_v(cells_number+i)=0.0d0
 d_p(cells_number+i)=0.0d0  
enddo

end subroutine infinitive_boundary_HS
!******************************************************************************************


!******************************************************************************************
!the treatment of boundary
subroutine infinitive_boundary_cell

implicit none

integer :: i,k

do i=-ighost+1,0
 do k=1,3
  u(i,k)=u(1,k)   
 enddo 
 uu(i)=uu(0)  
 v_aver(i)=v_aver(1)
 p_aver(i)=p_aver(1)
enddo

do i=1,ighost
 do k=1,3
  u(cells_number+i,k)=u(cells_number,k) 
 enddo 
 uu(cells_number+i)=uu(cells_number)
 v_aver(cells_number+i)=v_aver(cells_number) 
 p_aver(cells_number+i)=p_aver(cells_number)    
enddo

end subroutine infinitive_boundary_cell
!******************************************************************************************


end module infinitive_boundary