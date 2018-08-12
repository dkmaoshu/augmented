program Riemann_initialize

use function_definitions
! 'function_definitions.f90'

use type_definitions
! 'type_definitions.f90'

use out_put
! 'out_put.f90'

implicit none

! Geometrical values
!real(8) :: xl, xr

! Physical values
real(8) :: rhol, vl, pl, rhor, vr, pr

! real*8,parameter :: xl=-10.0d0,xm1=0.0d0,xm2=0.0d0,xr=10.0d0,time=5.0d0
!real*8,parameter :: xl=-1.0d0,xm1=0.0d0,xm2=0.0d0,xr=1.0d0,time=0.4d0
! real*8,parameter :: xl=-2.0d0,xm1=0.0d0,xm2=0.0d0,xr=2.0d0,time=0.8d0
! real*8,parameter :: xl=-1.0d0,xm1=0.0d0,xm2=0.0d0,xr=0.0d0,time=0.8d0
! real*8,parameter :: xl=-0.9d0,xm1=0.0d0,xm2=0.0d0,xr=0.1d0,time=0.6d0
!real*8,parameter :: rhol0=1.0d0,vl0=0.0d0,pl0=1.0d0
!real*8,parameter :: rhom0=1.0d0,vm0=0.0d0,pm0=1.0d0
!real*8,parameter :: rhor0=0.4263194230d0,vr0=0.9274526117d0,pr0=0.3031301790d0
! real*8,parameter :: rhol0=0.4263194230d0,vl0=-0.9274526117d0,pl0=0.3031301790d0
! real*8,parameter :: rhom0=0.4263194230d0,vm0=-0.9274526117d0,pm0=0.3031301790d0
! real*8,parameter :: rhor0=1.0d0,vr0=0.0d0,pr0=1.0d0

integer :: i, k, case_number
real*8 :: rho, v, p

!xl=-1.0d0; xr=1.0d0

print*, 'INPUT CASE NUMBER'
read(*,'(i5)') case_number

select case(case_number)

 case(1)
  rhol=1.0d0; vl=0.0d0; pl=1.0d0
  rhor=0.125d0; vr=0.0d0; pr=0.1d0
   
 case(2) 
  rhol=1.0d0; vl=0.0d0; pl=1.0d0
  rhor=0.4263194230d0; vr=0.9274526117d0; pr=0.3031301790d0
   
 case(3)
  rhol=0.4263194230d0; vl=-0.9274526117d0; pl=0.3031301790d0
  rhor=1.0d0; vr=0.0d0; pr=1.0d0

 case(4) 
  rhol=0.2655737113d0; vl=0.9274526117d0; pl=0.3031301790d0
  rhor=0.125d0; vr=0.0d0; pr=0.1d0
    
 case default
    
end select

! Fix grid.
dx=(xr-xl)/dfloat(cells_number)
do i=-ighost, cells_number+ighost 
 x(i)=xl+dfloat(i)*dx
enddo

do i=-ighost, cells_number+ighost 
 if(x(i).gt.xm1) then
  rho=rhor    
  v=vr          
  p=pr  
 else
  rho=rhol
  v=vl      
  p=pl
 endif
 phys%u(i,1)=rho
 phys%u(i,2)=f_m(rho,v,p)
 phys%u(i,3)=f_e(rho,v,p)
 phys%uu(i)=entropy_per_volume(rho,v,p)
 phys%v_aver(i)=v
 phys%p_aver(i)=p
end do
do i=0, cells_number
 phys%rho_boundary(i)=0.5d0*(phys%u(i,1)+phys%u(i+1,1))
 phys%v_boundary(i)=0.5d0*(f_v(phys%u(i,1),phys%u(i,2),phys%u(i,3))+f_v(phys%u(i+1,1),phys%u(i+1,2),phys%u(i+1,3)))
 phys%p_boundary(i)=0.5d0*(f_p(phys%u(i,1),phys%u(i,2),phys%u(i,3))+f_p(phys%u(i+1,1),phys%u(i+1,2),phys%u(i+1,3)))
enddo

do i=-ighost, -1   
 phys%rho_boundary(i)=phys%rho_boundary(0) 
 phys%v_boundary(i)=phys%v_boundary(0) 
 phys%p_boundary(i)=phys%p_boundary(0) 
enddo

do i=1, ighost     
 phys%rho_boundary(cells_number+i)=phys%rho_boundary(cells_number) 
 phys%v_boundary(cells_number+i)=phys%v_boundary(cells_number) 
 phys%p_boundary(cells_number+i)=phys%p_boundary(cells_number) 
enddo

call output_solution


!******************************************************************************************
end program Riemann_initialize