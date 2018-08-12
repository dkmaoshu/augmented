module half_step_reconstruction
! Compute the half steps for the reconstruction.

use solution
! 'solution.f90'

use physics
! 'physics.f90'

public  compute_half_steps
private ABS_mass_HS, ABS_velocity_HS, Vaver_via_Paver_drho, Paver_fix_Vaver_dv, fix_Vaver_drho, conservation_check


contains


subroutine compute_half_steps

implicit none

integer :: i
real*8, dimension(-ighost:cells_number+ighost) :: d_rho_ub,d_rho_e,d_v_e,d_p_e,temp_v
real*8 :: rho,m,e,error
real*8 :: rho_minus,rho_plus
real*8 :: v_minus,v,v_plus
real*8 :: p_minus,p,p_plus
real*8 :: drhol,drhor,drho,dpl,dpr,dp,dv
real*8 :: temp,temp1,temp1_v,temp1_vv,temp_p,temp_pp,temp1_p,temp2_p,temp1_pp,temp_rho
real*8 :: f_entropy_left_max,dp_max,dp_max_modify
real*8 :: d_rho_tvd,d_v_tvd,d_p_tvd,dp1,dp2
real*8 :: d_rho_temp, d_v_temp,v2
real*8 :: paver
real*8 :: drho0,dp0,alpha,uu_max,dp_star,f_max,f_d_rho_max

 real*8 :: dv2

 real(8) :: xxx
 integer :: counter
   
flag_tvd_rho=0
flag_tvd_v=0
   
do i=0, cells_number 
  
 rho=u(i,1)
 m=u(i,2)
 e=u(i,3) 
  
 rho_minus=u(i-1,1)
 rho_plus=u(i+1,1)
  
 v_minus=f_v(u(i-1,1),u(i-1,2),u(i-1,3)) 
 v=f_v(u(i,1),u(i,2),u(i,3)) 
 v_plus=f_v(u(i+1,1),u(i+1,2),u(i+1,3)) 
  
 p_minus=f_p(u(i-1,1),u(i-1,2),u(i-1,3)) 
 p=f_p(u(i,1),u(i,2),u(i,3)) 
 p_plus=f_p(u(i+1,1),u(i+1,2),u(i+1,3)) 
  
 if(p_aver(i).gt.p) then
  if(p_aver(i).gt.p+tol)  call error_message 
  p_aver(i)=p
 endif
  
 d_rho_tvd=rminmod(rho-rho_minus,rho_plus-rho)
 call ABS_mass_HS(rho,v,p,v_aver(i),p_aver(i),d_rho_tvd,d_rho_temp) 
! d_rho_e(i)=sgn(d_rho_tvd)*dsqrt(drho2) 
 if(dabs(d_rho_tvd)<d_rho_temp) then
  d_rho(i)=d_rho_tvd
  flag_tvd_rho(i)=1
  call Vaver_via_Paver_drho(rho,d_rho(i),v,p,p_aver(i),v_aver(i))

!****************************************************       
! 	   print*,	'Q和密度台阶修改', i,x(i)
!****************************************************
 else
  d_rho(i)=sgn(d_rho_tvd)*d_rho_temp
 endif
  
 d_v_tvd=rminmod(rminmod(v_aver(i)-v_minus,v_plus-v_aver(i)),rminmod(v-v_minus,v_plus-v))
 call ABS_velocity_HS(rho,v,p,v_aver(i),p_aver(i),d_v_temp) 
! d_v_e(i)=sgn(v_plus-v_minus)*dsqrt(dv2)
 if(dabs(d_v_tvd)<d_v_temp) then
!  temp1_p=p-0.5d0*(gamma-1.0d0)*rho*(d_v_tvd**2.0d0-(v_aver(i)-v)**2.0d0)  
  call Paver_fix_Vaver_dv(rho,v,v_aver(i),d_v_tvd,p,temp1_p)
!  temp2_p=p-0.5d0*(gamma-1.0d0)*rho*(v_aver(i)-v)**2.0d0*(rho**2.0d0-d_rho_tvd**2.0d0)/(d_rho_tvd**2.0d0+tol)
  call Paver_fix_Vaver_drho(rho,d_rho_tvd,v,v_aver(i),p,temp2_p)
  if(temp1_p<temp2_p) then
   p_aver(i)=temp1_p
   d_v(i)=d_v_tvd
   call ABS_mass_HS(rho,v,p,v_aver(i),p_aver(i),d_rho_tvd,d_rho_temp) 
   d_rho(i)=sgn(d_rho_tvd)*d_rho_temp
    
!****************************************************       
!   print*,	'P和速度台阶修改', i,x(i)
!****************************************************
  else
   d_v(i)=d_v_tvd
!   temp_pp=0.0d0
   call compute_p_drho(rho,d_v_tvd,d_rho(i),p,temp_pp) 
   p_aver(i)=temp_pp
   call Vaver_via_Paver_drho(rho,d_rho(i),v,p,p_aver(i),v_aver(i))
   
!****************************************************       
!    print*,'P，Q和速度台阶修改', i,x(i)
!****************************************************
  endif 
  flag_tvd_v(i)=1
 else
  d_v(i)=sgn(d_v_tvd)*d_v_temp
  if(dabs(d_v(i)*d_rho(i)-rho*(v-v_aver(i))).gt.1.0d-11) call error_message 
 endif
   
!  xxx=0.5d0*(f_m(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+f_m(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
              
 d_p_tvd=rminmod(rminmod(p_aver(i)-p_minus,p_plus-p_aver(i)),rminmod(p-p_minus,p_plus-p))
 f_max=entropy_per_volume(rho,v,p_aver(i))
 dp_max=p_aver(i)*d_rho(i)/rho
 f_d_rho_max=f_entropy_left(rho,v_aver(i),p_aver(i),uu(i),d_rho(i),dp_max)
  
 ! 3.1
 if(f_max>uu(i)) then
 ! 3.1.1
  if(f_d_rho_max>uu(i)) then
 ! 3.1.1.1  
   temp=f_entropy_left(rho,v_aver(i),p_aver(i),uu(i),d_rho(i),d_p_tvd)
   if(temp<uu(i)) then
    dpl=dp_max
    dpr=d_p_tvd
    call iterative_method_dp(rho,v_aver(i),p_aver(i),uu(i),d_rho(i),dpl,dpr,dp_star) 
! 3.1.1.1.1
    if(((dp_star*d_p_tvd)>0.0d0).and.(dabs(dp_star)<dabs(d_p_tvd))) then
     d_p(i)=dp_star			
! 3.1.1.1.2和3.1.1.1.3			
    else
     if((dp_star*d_p_tvd)<0.0d0) then
      drho0=d_rho(i)
      dp0=0.0d0
      if(dabs(drho0)<tol) then     
       alpha=tol
      else
       alpha=dp0*1.0d8/(drho0*1.0d8)
      endif
     else
      drho0=d_rho(i)
      dp0=d_p_tvd
      if(dabs(drho0)<tol) then     
       alpha=tol
      else
       alpha=dp0*1.0d8/(drho0*1.0d8)
      endif
     endif
     drhol=0.0d0
     drhor=drho0
     call iterative_method_drho_cs(rho,m,e,drhol,drhor,p_aver(i),alpha,uu(i),drho) 
     d_rho(i)=drho
     d_p(i)=alpha*drho
     call compute_dv2_p(rho,d_rho(i),p,p_aver(i),dv2) 
     d_v(i)=sgn(v_plus-v_minus)*dsqrt(dv2)
     call Vaver_via_Paver_drho(rho,d_rho(i),v,p,p_aver(i),v_aver(i))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   print*,	' Q修改 3.1.1.1.2-3', i,x(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
! 3.1.1.2  
   else
    error=10000.d0
! 迭代       		   
      
     counter=0
      
     xxx=d_rho(i)*d_v(i)-u(i,1)*(f_v(u(i,1),u(i,2),u(i,3))-v_aver(i))
      
    do while(dabs(error)>=tol.and.f_dp(rho,p_aver(i),uu(i),d_rho(i),d_p_tvd)>=0.0d0)
     temp_rho=d_rho(i)
     temp1=f_entropy_left(rho,v_aver(i),p_aver(i),uu(i),d_rho_tvd,d_p_tvd)  
! 3.1.1.2.1 
     if(temp1<uu(i)) then
      drho=0.0d0
      drhol=d_rho(i)
      drhor=d_rho_tvd
      call iterative_method_drho_crho(rho,drhol,drhor,p_aver(i),uu(i),d_p_tvd,drho)   
      d_rho(i)=drho
      d_p(i)=d_p_tvd
! 3.1.1.2.2
     else
      d_rho(i)=d_rho_tvd
      d_p(i)=d_p_tvd
      uu(i)=temp1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  print*,	' Entropy和压力台阶修改3.1.1.2.2', i,x(i)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     endif 
      
      xxx=d_rho(i)*d_v(i)-u(i,1)*(f_v(u(i,1),u(i,2),u(i,3))-v_aver(i))
      
     if(dabs(d_rho(i))<tol)then
      if(d_rho(i)>=0.0d0)then
       d_rho(i)=tol
      else
       d_rho(i)=-tol
      endif
     endif
!     call compute_p_q(rho,v,p,d_rho(i),v_aver(i),paver) 
     call Paver_fix_Vaver_drho(rho,d_rho(i),v,v_aver(i),p,paver)
     p_aver(i)=paver
     
      xxx=d_rho(i)*d_v(i)-u(i,1)*(f_v(u(i,1),u(i,2),u(i,3))-v_aver(i))
     
     call ABS_velocity_HS(rho,v,p,v_aver(i),p_aver(i),dv2) 
     d_v(i)=sgn(d_v_tvd)*dsqrt(dv2)
     
      xxx=d_rho(i)*d_v(i)-u(i,1)*(f_v(u(i,1),u(i,2),u(i,3))-v_aver(i)) !!! The problem happens in the above instruction!!!!!!!!!
     
     error=dabs(temp_rho-d_rho(i))
     if(p_aver(i)>p) then
      p_aver(i)=p
     endif
      
      counter=counter+1	  
       
!      xxx=0.5d0*(f_m(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+f_m(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
              
      continue
      
    enddo
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  print*,	' P和Q修改 3.1.1.2.1-2', i,x(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!需要处理 
   endif
! 3.1.2
  else
! 3.1.2.1
   if(((dp_max*d_p_tvd)>0.0d0).and.(dabs(dp_max)<dabs(d_p_tvd))) then
    drho0=d_rho(i)
    dp0=dp_max
    if(dabs(drho0)<tol) then     
     alpha=tol
    else
     alpha=dp0*1.0d8/(drho0*1.0d8)
    endif		
! 3.1.2.2和3.1.2.3			
   else
    if((dp_max*d_p_tvd)<0.0d0) then
     drho0=d_rho(i)
     dp0=0.0d0
     if(dabs(drho0)<tol) then     
      alpha=tol
     else
      alpha=dp0*1.0d8/(drho0*1.0d8)
     endif
    else
     drho0=d_rho(i)
     dp0=d_p_tvd
     if(dabs(drho0)<tol) then     
      alpha=tol
     else
      alpha=dp0*1.0d8/(drho0*1.0d8)
     endif
    endif
   endif
    
   drhol=0.0d0
   drhor=drho0
   call iterative_method_drho_cs(rho,m,e,drhol,drhor,p_aver(i),alpha,uu(i),drho) 
   d_rho(i)=drho
   d_p(i)=alpha*drho
   call compute_dv2_p(rho,d_rho(i),p,p_aver(i),dv2) 
   d_v(i)=sgn(v_plus-v_minus)*dsqrt(dv2)
     
   call Vaver_via_Paver_drho(rho,d_rho(i),v,p,p_aver(i),v_aver(i))
   		 
!   if(v_aver(i)>=v) then
!    v_aver(i)=v+dsqrt(v2)
!   else
!    v_aver(i)=v-dsqrt(v2)
!   endif 
  		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print*,' Q修改3.1.2.1-3', i,x(i)      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		      
  endif
! 3.2
 else
  d_rho(i)=tol
  d_p(i)=tol
   
  p_aver(i)=rho**gamma*dexp(uu(i)/rho)
  v_aver(i)=v
    
  if(p_aver(i)>p) then
   p_aver(i)=p
  endif
  call compute_dv2_p(rho,d_rho(i),p,p_aver(i),dv2)
  d_v(i)=sgn(v_plus-v_minus)*dsqrt(dv2)
  call Vaver_via_Paver_drho(rho,d_rho(i),v, p,p_aver(i),v_aver(i))
     		 
!  if(v_aver(i)>=v) then
!   v_aver(i)=v+dsqrt(v2)
!  else
!   v_aver(i)=v-dsqrt(v2)
!  endif    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       print*,' P和Q修改3.2 ', i,x(i)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 endif

enddo

! continue

call conservation_check

end subroutine compute_half_steps


!***************************************************************************************
subroutine ABS_mass_HS(rho,v,p,vaver,paver,d_rho_tvd,d_rho_temp) 
! Compute the absolute value of mass HS from momentum and energy conservation.

implicit none

real*8, intent(in) :: rho, v, p, vaver, paver, d_rho_tvd
real*8, intent(out) :: d_rho_temp

real*8 :: temp_p, temp_v,temp

temp_p=p-paver
temp_v=(vaver-v)**2.0

d_rho_temp=0.5d0*rho**3.0d0*temp_v+tol*d_rho_tvd**2.0d0
d_rho_temp=d_rho_temp/(temp_p/(gamma-1.0d0)+0.5d0*rho*temp_v+tol)

d_rho_temp=dsqrt(d_rho_temp)

end subroutine ABS_mass_HS


!***************************************************************************************
subroutine ABS_velocity_HS(rho,v,p,vaver,paver,d_v_temp) 
! Compute the absolute value of velocity HS from momentum and energy conservation.

implicit none

real*8, intent(in) :: rho, v, p, vaver, paver
real*8, intent(out) :: d_v_temp

d_v_temp=2.0d0*(p-paver)/((gamma-1.0d0)*rho)+(vaver-v)**2.0d0
d_v_temp=dsqrt(d_v_temp)


end subroutine ABS_velocity_HS
!***************************************************************************************


subroutine Vaver_via_Paver_drho(rho,drho,v,p,paver,vaver) 
! Modify velocity average with fixed pressure average and mass HS.

implicit none

real*8, intent(in) :: rho, drho, v, p, paver
real*8, intent(out) :: vaver

real*8 :: difference

difference=2.0d0*drho**2.0d0/((gamma-1.0d0)*rho*(rho**2.0d0-drho**2.0d0))*(p-paver)

if(vaver.gt.v) then
 vaver=v+dsqrt(difference)
else
 vaver=v-dsqrt(difference)
endif

end subroutine Vaver_via_Paver_drho


subroutine Paver_fix_Vaver_dv(rho,v,vaver,dv,p,paver)
! Modify pressure average with fixed velocity average and velocity HS.

implicit none

real(8), intent(in) :: rho, v, vaver, dv, p
real(8), intent(out) :: paver

paver=p-0.5d0*(gamma-1.0d0)*rho*(dv**2.0d0-(vaver-v)**2.0d0)  

end subroutine Paver_fix_Vaver_dv


subroutine Paver_fix_Vaver_drho(rho,drho,v,vaver,p,paver)
! Modify pressure average with fixed velocity average and velocity HS.

implicit none

real(8), intent(in) :: rho, drho, v, vaver, p
real(8), intent(out) :: paver

paver=p-0.5d0*(gamma-1.0d0)*rho*(vaver-v)**2.0d0*(rho**2.0d0-drho**2.0d0)/(drho**2.0d0+tol*tol)

end subroutine Paver_fix_Vaver_drho


!***************************************************************************************
subroutine compute_dv2_p(rho,drho,p,paver,dv2) 

implicit none

real*8 :: rho,drho,p,paver,dv2

dv2=2.0d0*rho/((gamma-1.0d0)*(rho**2.0d0-drho**2.0d0))*(p-paver)

end subroutine compute_dv2_p


!***************************************************************************************
subroutine compute_p_drho(rho,dv,drho,p,paver) 

implicit none

real*8 :: rho,dv,drho,p,paver

paver=p-dv**2.0d0*(gamma-1.0d0)*(rho**2.0d0-drho**2.0d0)/(2.0d0*rho)

end subroutine compute_p_drho


!!***************************************************************************************
!subroutine compute_p_q(rho,v,p,drho,vaver,paver) 
!
!implicit none
!
!real*8 :: rho,v,p,drho,vaver,paver
!
!paver=p-0.5d0*rho*(gamma-1.0d0)*(rho**2.0d0-drho**2.0d0)*(vaver-v)**2.0d0/drho**2.0d0
!
!end subroutine compute_p_q


!***************************************************************************************
subroutine iterative_method_drho_cs(rho,m,e,drhol,drhor,paver,alpha,uue,drho) 

implicit none

integer :: iterate_step
real*8 :: rho,m,e,uue,drho
real*8 :: drhol,f_valuel
real*8 :: f_valuem
real*8 :: drhor,f_valuer
real*8 :: xl,xm,xr,xm0,error,paver,alpha

xl=drhol
xr=drhor
error=1.0d3
iterate_step=0

f_valuel=f_drho(rho,m,e,paver,alpha,uue,xl)
f_valuer=f_drho(rho,m,e,paver,alpha,uue,xr)
xm=0.5d0*(xl+xr)

if(f_valuer<=0.0d0) then
 if(f_valuel*f_valuer<=0.0d0) then
  do while(error>tol)
   xm0=xm
   f_valuel=f_drho(rho,m,e,paver,alpha,uue,xl)
   f_valuer=f_drho(rho,m,e,paver,alpha,uue,xr)
   f_valuem=f_drho(rho,m,e,paver,alpha,uue,xm)
   if((f_valuel*f_valuem)>0.0d0.and.(f_valuer*f_valuem)>0.0d0) then
    print*,'It is wrong!'
    stop		  
   endif
   if((f_valuel*f_valuem)<0.0d0.and.(f_valuer*f_valuem)<0.0d0) then
    print*,'It is wrong!'
    stop		  
   endif
   if((f_valuel*f_valuem)<=0.0d0) then
    xr=xm
   else
    xl=xm
   endif
   xm=0.5d0*(xl+xr)
   error=dabs(xm-xm0)
   iterate_step=iterate_step+1
	     !print*,error,iterate_step
  enddo
  drho=xm
 else
  drho=0.0d0
!       drho=drhor
 endif
else
 drho=drhor
endif

end subroutine iterative_method_drho_cs
!***************************************************************************************


!***************************************************************************************
subroutine iterative_method_drho_crho(rho,drhol,drhor,paver,uue,dp,drho) 

implicit none

integer :: iterate_step
real*8 :: rho,drhol,drhor,paver,uue,dp,drho
real*8 :: f_valuel,f_valuem,f_valuer
real*8 :: xl,xm,xr,xm0,error

xl=drhol
xr=drhor
error=1.0d3
iterate_step=0

f_valuel=f_dp(rho,paver,uue,xl,dp)
f_valuer=f_dp(rho,paver,uue,xr,dp)
xm=0.5d0*(xl+xr)

if(f_valuer<=0.0d0) then
 if(f_valuel*f_valuer<=0.0d0) then
  do while(error>tol)
   xm0=xm
   f_valuel=f_dp(rho,paver,uue,xl,dp)
   f_valuer=f_dp(rho,paver,uue,xr,dp)
    
   f_valuem=f_dp(rho,paver,uue,xm,dp)
   if((f_valuel*f_valuem)>0.0d0.and.(f_valuer*f_valuem)>0.0d0) then
    print*,'It is wrong!'
    stop		  
   endif
   if((f_valuel*f_valuem)<0.0d0.and.(f_valuer*f_valuem)<0.0d0) then
    print*,'It is wrong!'
    stop		  
   endif
   if((f_valuel*f_valuem)<=0.0d0) then
    xr=xm
   else
    xl=xm
   endif
   xm=0.5d0*(xl+xr)
   error=dabs(xm-xm0)
   iterate_step=iterate_step+1
	     !print*,error,iterate_step
  enddo
  drho=xm
 else
!       print*,"cuowu"
  drho=0.0d0
!       drho=drhor
 endif
else
   drho=drhor
endif

end subroutine iterative_method_drho_crho
!***************************************************************************************


!***************************************************************************************
subroutine iterative_method_dp(rho,vaver,paver,uue,drho,dpl,dpr,dp) 

implicit none

integer :: iterate_step
real*8 :: rho,vaver,paver,uue,drho,dp
real*8 :: dpl,f_valuel
real*8 :: dpm,f_valuem
real*8 :: dpr,f_valuer
real*8 :: xl,xm,xr,xm0,error

xl=dpl
xr=dpr

error=1.0d3
iterate_step=0

f_valuel=f_dp(rho,paver,uue,drho,xl)
f_valuer=f_dp(rho,paver,uue,drho,xr)
xm=0.5d0*(xl+xr)


! if(f_valuer<=0.0d0) then
if(f_valuel*f_valuer<=0.0d0) then
 do while(error>tol)
  xm0=xm
  f_valuel=f_dp(rho,paver,uue,drho,xl)
  f_valuer=f_dp(rho,paver,uue,drho,xr)
   
  f_valuem=f_dp(rho,paver,uue,drho,xm)
  if((f_valuel*f_valuem)>0.0d0.and.(f_valuer*f_valuem)>0.0d0) then
   print*,'It is wrong!'
   stop		  
  end if
  if((f_valuel*f_valuem)<0.0d0.and.(f_valuer*f_valuem)<0.0d0) then
   print*,'It is wrong!'
   stop		  
  endif
  if((f_valuel*f_valuem)<=0.0d0) then
   xr=xm
  else
   xl=xm
  endif
  xm=0.5d0*(xl+xr)
  error=dabs(xm-xm0)
  iterate_step=iterate_step+1
     !print*,error,iterate_step
 enddo
 dp=xm
else
 dp=0.0d0
endif

! continue

end subroutine iterative_method_dp


subroutine conservation_check
! Check the conservation and entropy non-decrease after the reconstruction.

implicit none

integer :: i
real(8) :: xxx, yyy, zzz

real(8) :: toll

! real(8) :: ttt

toll=10.0d0*tol

do i=1, cells_number-1
   
!  ttt=d_rho(i)*d_v(i)-u(i,1)*(f_v(u(i,1),u(i,2),u(i,3))-v_aver(i))
   
 xxx=0.5d0*(f_m(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+f_m(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
 yyy=0.5d0*(f_e(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+f_e(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
 zzz=0.5d0*(entropy_per_volume(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+entropy_per_volume(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
 if(dabs(xxx-u(i,2)).gt.toll) then
  print*, i, '   momentum'; call error_message
 end if
 if(dabs(yyy-u(i,3)).gt.toll) then
  print*, i, '   energy'; call error_message
 end if
 if(dabs(zzz-uu(i)).gt.toll) then
  print*, i, '   entropy'; call error_message
 end if
end do

end subroutine conservation_check


!***************************************************************************************
end module half_step_reconstruction