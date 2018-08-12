module initialcondition

!use consant
use function_definitions

use typedefinition

implicit none

! real*8,parameter :: xl=-10.0d0,xm1=0.0d0,xm2=0.0d0,xr=10.0d0,time=5.0d0
real*8,parameter :: xl=-1.0d0,xm1=0.0d0,xm2=0.0d0,xr=1.0d0,time=0.4d0
! real*8,parameter :: xl=-2.0d0,xm1=0.0d0,xm2=0.0d0,xr=2.0d0,time=0.8d0
! real*8,parameter :: xl=-1.0d0,xm1=0.0d0,xm2=0.0d0,xr=0.0d0,time=0.8d0
! real*8,parameter :: xl=-0.9d0,xm1=0.0d0,xm2=0.0d0,xr=0.1d0,time=0.6d0
real*8,parameter :: rhol0=1.0d0,vl0=0.0d0,pl0=1.0d0
real*8,parameter :: rhom0=1.0d0,vm0=0.0d0,pm0=1.0d0
real*8,parameter :: rhor0=0.4263194230d0,vr0=0.9274526117d0,pr0=0.3031301790d0
! real*8,parameter :: rhol0=0.4263194230d0,vl0=-0.9274526117d0,pl0=0.3031301790d0
! real*8,parameter :: rhom0=0.4263194230d0,vm0=-0.9274526117d0,pm0=0.3031301790d0
! real*8,parameter :: rhor0=1.0d0,vr0=0.0d0,pr0=1.0d0
contains

!******************************************************************************************
!intialize u and U
subroutine initialize

implicit none

integer :: i,k
real*8 :: rho,v,p


dx=(xr-xl)/nx
do i=-ighost,nx+ighost 
   x(i)=xl+dx*i+0.5d0*dx
enddo
do i=-ighost,nx+ighost 
   if(x(i)>xm2) then
      rho=rhor0    
      v=vr0          
      p=pr0  
   elseif(x(i)>xm1) then
      rho=rhom0    
      v=vm0          
      p=pm0
   else
      rho=rhol0     
      v=vl0      
      p=pl0
   endif
   phys%u(i,1)=rho
   phys%u(i,2)=f_m(rho,v,p)
   phys%u(i,3)=f_e(rho,v,p)
   phys%uu(i)=uu(rho,v,p)
   phys%v_aver(i)=v
   phys%p_aver(i)=p
enddo
do i=0,nx
   phys%rho_boundary(i)=0.5d0*(phys%u(i,1)+phys%u(i+1,1))
   phys%v_boundary(i)=0.5d0*(f_v(phys%u(i,1),phys%u(i,2),phys%u(i,3))+f_v(phys%u(i+1,1),phys%u(i+1,2),phys%u(i+1,3)))
   phys%p_boundary(i)=0.5d0*(f_p(phys%u(i,1),phys%u(i,2),phys%u(i,3))+f_p(phys%u(i+1,1),phys%u(i+1,2),phys%u(i+1,3)))
enddo

do i=-ighost,-1   
   phys%rho_boundary(i)=phys%rho_boundary(0) 
   phys%v_boundary(i)=phys%v_boundary(0) 
   phys%p_boundary(i)=phys%p_boundary(0) 
enddo

do i=1,ighost     
   phys%rho_boundary(nx+i)=phys%rho_boundary(nx) 
   phys%v_boundary(nx+i)=phys%v_boundary(nx) 
   phys%p_boundary(nx+i)=phys%p_boundary(nx) 
enddo

end subroutine initialize
!******************************************************************************************
!******************************************************************************************
!the treatment of boundary
subroutine boundary_d

implicit none

integer :: i,k

do i=-ighost,-1   
   phys%d_rho(i)=phys%d_rho(0)  
   phys%d_v(i)=phys%d_v(0) 
   phys%d_p(i)=phys%d_p(0) 
enddo

do i=1,ighost
   phys%d_rho(nx+i)=phys%d_rho(nx)  
   phys%d_v(nx+i)=phys%d_v(nx)
   phys%d_p(nx+i)=phys%d_p(nx)  
enddo

end subroutine boundary_d
!******************************************************************************************





!******************************************************************************************
!the treatment of boundary
subroutine boundary

implicit none

integer :: i,k

do i=-ighost,-1
   do k=1,3
      phys%u(i,k)=phys%u(0,k)   
   enddo 
   phys%uu(i)=phys%uu(0)  
   phys%v_aver(i)=phys%v_aver(0)
   phys%p_aver(i)=phys%p_aver(0)
enddo

do i=1,ighost
   do k=1,3
      phys%u(nx+i,k)=phys%u(nx,k) 
   enddo 
   phys%uu(nx+i)=phys%uu(nx)
   phys%v_aver(nx+i)=phys%v_aver(nx) 
   phys%p_aver(nx+i)=phys%p_aver(nx)    
enddo

end subroutine boundary
!******************************************************************************************

!******************************************************************************************
!get dt using CFL condition number
subroutine cfldt

implicit none

integer :: i
real*8 :: rho,v,p,temp1,temp2,temp,rlamdamax

rlamdamax=0.0d0
do i=0,nx
   rho=phys%u(i,1)   
   v=f_v(phys%u(i,1),phys%u(i,2),phys%u(i,3))
   p=f_p(phys%u(i,1),phys%u(i,2),phys%u(i,3))
   temp=abs(v)+dsqrt(gamma*p/rho)
   if(temp>rlamdamax) rlamdamax=temp
enddo
dt=dx/rlamdamax*cfl
end subroutine cfldt 
!******************************************************************************************
end module initialcondition






