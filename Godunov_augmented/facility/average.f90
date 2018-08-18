module average

use solution
! 'solution.f90'

use physics
!' physics.f90'

use RiemannSolver
! 'RiemannSolver.f90'

implicit none


contains
!***************************************************************************************


!***************************************************************************************
subroutine average_vandp

implicit none

integer :: i,IW
real*8 :: RS,US,PS
real*8 :: vl_aver,vr_aver,pl_aver,pr_aver,entropyl,entropyr,uul_aver,uur_aver
real*8 :: rarefaction_rho_l,rarefaction_rho_r,shock_rho_l,shock_rho_r,contact_rho_l,contact_rho_r
real*8,dimension(-ighost+1:cells_number+ighost) :: temp_v_aver,temp_p_aver
real*8,dimension(-ighost+1:cells_number+ighost) :: R1,U1,P1,R2,U2,P2

real(8) :: u_star, p_star

temp_v_aver=v_aver;temp_p_aver=p_aver

do i=-ighost+1, cells_number+ighost 

 R1(i)=u(i,1)-d_rho(i)
 U1(i)=v_aver(i)-d_v(i)
 P1(i)=p_aver(i)-d_p(i)
   
 R2(i)=u(i,1)+d_rho(i)
 U2(i)=v_aver(i)+d_v(i)
 P2(i)=p_aver(i)+d_p(i)

!   endif
enddo

do i=1, cells_number

 rarefaction_rho(i)=0.0d0;shock_rho(i)=0.0d0;contact_rho(i)=0.0d0
 v_aver(i)=0.0d0;p_aver(i)=0.0d0;entropy(i)=0.0d0
   
!***************************************************************************************
 rarefaction_rho_l=0.0d0;rarefaction_rho_r=0.0d0;shock_rho_l=0.0d0;shock_rho_r=0.0d0;contact_rho_l=0.0d0;contact_rho_r=0.0d0
 vl_aver=0.0d0;vr_aver=0.0d0;pl_aver=0.0d0;pr_aver=0.0d0;uul_aver=0.0d0;uur_aver=0.0d0;entropyl=0.0d0;entropyr=0.0d0
   
 CALL RIEMANN0(R2(i-1),U2(i-1),P2(i-1),R1(i),U1(i),P1(i),u_star,p_star,IW)
 call aver_velocityandpressure(vl_aver,vr_aver,pl_aver,pr_aver,uul_aver,uur_aver,entropyl,entropyr,rarefaction_rho_l,rarefaction_rho_r,shock_rho_l,shock_rho_r,contact_rho_l,contact_rho_r)
 v_aver(i)=v_aver(i)+vr_aver
 p_aver(i)=p_aver(i)+pr_aver
 entropy(i)=entropy(i)+entropyr
!!! uu(i)=uu(i)+uur_aver

! print*,i,"×ó",iw
!***************************************************************************************
 rarefaction_rho_l=0.0d0;rarefaction_rho_r=0.0d0;shock_rho_l=0.0d0;shock_rho_r=0.0d0;contact_rho_l=0.0d0;contact_rho_r=0.0d0
 vl_aver=0.0d0;vr_aver=0.0d0;pl_aver=0.0d0;pr_aver=0.0d0;uul_aver=0.0d0;uur_aver=0.0d0;entropyl=0.0d0;entropyr=0.0d0
   
 CALL RIEMANN0(R1(i),U1(i),P1(i),R2(i),U2(i),P2(i),u_star,p_star,IW)
 call aver_velocityandpressure(vl_aver,vr_aver,pl_aver,pr_aver,uul_aver,uur_aver,entropyl,entropyr,rarefaction_rho_l,rarefaction_rho_r,shock_rho_l,shock_rho_r,contact_rho_l,contact_rho_r)
 v_aver(i)=v_aver(i)+vl_aver+vr_aver
 p_aver(i)=p_aver(i)+pl_aver+pr_aver
 entropy(i)=entropy(i)+entropyl+entropyr 
!!! uu(i)=uu(i)+uul_aver+uur_aver

! print*,i,"ÖÐ",iw
!***************************************************************************************
 rarefaction_rho_l=0.0d0;rarefaction_rho_r=0.0d0;shock_rho_l=0.0d0;shock_rho_r=0.0d0;contact_rho_l=0.0d0;contact_rho_r=0.0d0
 vl_aver=0.0d0;vr_aver=0.0d0;pl_aver=0.0d0;pr_aver=0.0d0;uul_aver=0.0d0;uur_aver=0.0d0;entropyl=0.0d0;entropyr=0.0d0
    
 CALL RIEMANN0(R2(i),U2(i),P2(i),R1(i+1),U1(i+1),P1(i+1),u_star,p_star,IW)
 call aver_velocityandpressure(vl_aver,vr_aver,pl_aver,pr_aver,uul_aver,uur_aver,entropyl,entropyr,rarefaction_rho_l,rarefaction_rho_r,shock_rho_l,shock_rho_r,contact_rho_l,contact_rho_r)

 v_aver(i)=v_aver(i)+vl_aver
 p_aver(i)=p_aver(i)+pl_aver
 entropy(i)=entropy(i)+entropyl
!!!! uu(i)=uu(i)+uul_aver

! print*,i,"ÓÒ",iw

enddo

end subroutine average_vandp
!***************************************************************************************


!***************************************************************************************                                                  
subroutine aver_velocityandpressure(vl_aver,vr_aver,pl_aver,pr_aver,uul_aver,uur_aver,entropyl,entropyr,rarefaction_rho_l,rarefaction_rho_r,shock_rho_l,shock_rho_r,contact_rho_l,contact_rho_r)

implicit none

integer :: i
real*8 :: G1, G2, G3, G4, G5, G6, G7, G8, G9
real*8 :: DL, UL, PL, CL, DR, UR, PR, CR, PM, UM
real*8 :: SML,SMR
real*8 :: C,SHL,CML,STL,PML,SL,PMR,SR,SHR,CMR,STR
real*8 :: a,b
real*8 :: vl_aver,vr_aver,pl_aver,pr_aver,uul_aver,uur_aver
real*8 :: intl_rho,intl_v,intl_p,intll_rho_sonic,intll_v_sonic,intll_p_sonic,intlr_rho_sonic,intlr_v_sonic,intlr_p_sonic
real*8 :: intr_rho,intr_v,intr_p,intrl_rho_sonic,intrl_v_sonic,intrl_p_sonic,intrr_rho_sonic,intrr_v_sonic,intrr_p_sonic
real*8 :: rhoml_s,rhomr_s,rhoml_r,rhomr_r
real*8 :: entropy_l,entropy_r,entropym_l,entropym_r
real*8 :: rhoml,rhomr,entropyl,entropyr
real*8 :: rho_l_sonic,rho_r_sonic
real*8 :: rarefaction_rho_l,rarefaction_rho_r,shock_rho_l,shock_rho_r,contact_rho_l,contact_rho_r
real*8, dimension(-1:1) :: intp_v,intp_p !,gp

COMMON /GAMMAS/ G1, G2, G3, G4, G5, G6, G7, G8, G9
COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR, PM, UM

!gp(-1)=-dsqrt(0.6d0)
!gp(0)=0.0d0
!gp(1)=dsqrt(0.6d0)

! if(pm>pl) then
!    print*,pm,pr
! endif

entropy_l=entropy_per_mass(dl,ul,pl) ! Left entropy per mass
entropy_r=entropy_per_mass(dr,ur,pr) ! Right entropy per mass

!shock

sl=ul-cl*dsqrt(g2*pm/pl+g1) ! Left shock speed
sr=ur+cr*dsqrt(g2*pm/pr+g1) ! Right shock speed

!rarefaction

cml=cl*(pm/pl)**g1 ! Left middle sound speed
cmr=cr*(pm/pr)**g1 ! Right middle sound speed

shl=ul-cl  ! Left characteristic speed
stl=um-cml ! left middle characteristic speed

shr=ur+cr  ! Right characteristic speed
str=um+cmr ! Right middle characteristic speed

rhoml_r=dl*(pm/pl)**g8 ! Left middle mass density
rhomr_r=dr*(pm/pr)**g8 ! Right middle mass density
 
!the density between the shock and contact discontinuity
rhoml_s=dl*((pm/pl+g6)/(pm/pl*g6+1.0d0))
rhomr_s=dr*((pm/pr+g6)/(pm/pr*g6+1.0d0))
entropym_l=entropy_per_mass(rhoml_s,um,pm)
entropym_r=entropy_per_mass(rhomr_s,um,pm)
!the density between the shock and contact discontinuity
rhoml=dl*((pm/pl+g6)/(pm/pl*g6+1.0d0))
rhomr=dr*((pm/pr+g6)/(pm/pr*g6+1.0d0))
!!the density of the sonic point
!rho_l_sonic=dl*(g5+g6*ul/cl)**g4
!rho_r_sonic=dr*(g5-g6*ur/cr)**g4

call compute_rarefaction_left(shl,stl,intl_rho,intl_v,intl_p)
call compute_rarefaction_left(shl,0.0d0,intll_rho_sonic,intll_v_sonic,intll_p_sonic)
call compute_rarefaction_left(0.0d0,stl,intlr_rho_sonic,intlr_v_sonic,intlr_p_sonic)
call compute_rarefaction_right(str,shr,intr_rho,intr_v,intr_p)
call compute_rarefaction_right(str,0.0d0,intrl_rho_sonic,intrl_v_sonic,intrl_p_sonic)
call compute_rarefaction_right(0.0d0,shr,intrr_rho_sonic,intrr_v_sonic,intrr_p_sonic)

if(um>0) then
! The contact discontinuity goes to the right.

 if(pm>pl) then
! The left wave is a shock.   

  if(sl>0.0d0) then
! The left shock goes to the right; therefore, the let half of the Riemann solution is occupied by the initial left state.

   uul_aver=0.25d0*dl*entropy_l
   vl_aver=0.25d0*ul                                      
   pl_aver=0.25d0*pl
   entropyr=entropyr+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
   if(pm>pr) then
    entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
    uur_aver=sl*dt/dx*dl*entropy_l+(um-sl)*dt/dx*rhoml_s*entropym_l+(sr-um)*dt/dx*rhomr_s*entropym_r+(0.25d0-sr*dt/dx)*dr*entropy_r
    vr_aver=sl*dt/dx*ul+(sr-sl)*dt/dx*um+(0.25d0-sr*dt/dx)*ur
    pr_aver=sl*dt/dx*pl+(sr-sl)*dt/dx*pm+(0.25d0-sr*dt/dx)*pr
   else
    uur_aver=sl*dt/dx*dl*entropy_l+(um-sl)*dt/dx*rhoml_s*entropym_l+((um+cmr-um)*dt/dx*rhomr_r+(0.25d0-(ur+cr)*dt/dx)*dr+intr_rho)*entropy_r
    vr_aver=sl*dt/dx*ul+(um+cmr-sl)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v
    pr_aver=sl*dt/dx*pl+(um+cmr-sl)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p
   endif
  else
   entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
   uul_aver=(0.25d0-dabs(sl)*dt/dx)*dl*entropy_l+dabs(sl)*dt/dx*rhoml_s*entropym_l
   vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+dabs(sl)*dt/dx*um
   pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+dabs(sl)*dt/dx*pm
!    shock_rho_l=shock_rho_l+dabs(dl-rhoml_s)
   if(pm>pr) then
    entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
    uur_aver=um*dt/dx*rhoml_s*entropym_l+(sr-um)*dt/dx*rhomr_s*entropym_r+(0.25d0-sr*dt/dx)*dr*entropy_r
    vr_aver=sr*dt/dx*um+(0.25d0-sr*dt/dx)*ur
    pr_aver=sr*dt/dx*pm+(0.25d0-sr*dt/dx)*pr  
!     shock_rho_r=shock_rho_r+dabs(rhomr_s-dr)	
!     contact_rho_r=contact_rho_r+dabs(rhomr_s-rhoml_s) 
   else
    uur_aver=um*dt/dx*rhoml_s*entropym_l+(((um+cmr)-um)*dt/dx*rhomr_r+(0.25d0-(ur+cr)*dt/dx)*dr+intr_rho)*entropy_r
    vr_aver=(um+cmr)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v
    pr_aver=(um+cmr)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p
!     rarefaction_rho_r=rarefaction_rho_r+dabs(rhomr_r-dr)		
!     contact_rho_r=contact_rho_r+dabs(rhomr_r-rhoml_s)			
   endif
  endif
 else

! The left wave is a rarefaction wave.
  if((ul-cl)>0.0d0) then
! The left edge of the left rarefaction wave goes to the right; therefore, the left half of the Riemann solution is occupied by the initial left state.

   uul_aver=0.25d0*dl*entropy_l
   vl_aver=0.25d0*ul
   pl_aver=0.25d0*pl
   if(pm>pr) then
    entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
    uur_aver=((ul-cl)*dt/dx*dl+(um-(um-cml))*dt/dx*rhoml_r+intl_rho)*entropy_l+(sr-um)*dt/dx*rhomr_s*entropym_r+(0.25d0-sr*dt/dx)*dr*entropy_r
    vr_aver=(ul-cl)*dt/dx*ul+(sr-(um-cml))*dt/dx*um+(0.25d0-sr*dt/dx)*ur+intl_v
    pr_aver=(ul-cl)*dt/dx*pl+(sr-(um-cml))*dt/dx*pm+(0.25d0-sr*dt/dx)*pr+intl_p         
   else
    uur_aver=((ul-cl)*dt/dx*dl+(um-(um-cml))*dt/dx*rhoml_r+intl_rho)*entropy_l+((um+cmr-um)*dt/dx*rhomr_r+intr_rho+(0.25d0-(ur+cr)*dt/dx)*dr)*entropy_r
    vr_aver=(ul-cl)*dt/dx*ul+((um+cmr)-(um-cml))*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intl_v+intr_v
    pr_aver=(ul-cl)*dt/dx*pl+((um+cmr)-(um-cml))*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intl_p+intr_p       		    
   endif
  else
   if((um-cml)<0.0d0) then
! The right edge of the left rarefaction wave goes to the left; therefore, the whole left rarefaction wave goes to the left.
! Because the contact discontinuity now goes to the right, the left half of the Riemann solution consists of the left state, the left rarefaction wave and part of the left middle state. 

!    uul_aver=((0.25d0-dabs(ul-cl)*dt/dx)*dl+intl_rho+dabs(um-cml)*dt/dx*rhoml_r)*entropy_l 
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+dabs(um-cml)*dt/dx*um+intl_v 
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+dabs(um-cml)*dt/dx*pm+intl_p
      		 
    if(pm>pr) then
! The right wave is a shock. Since the contact discontinuity now goes to the right, the shock also goes to the right.

     entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
!     uur_aver=um*dt/dx*rhoml_r*entropy_l+(sr-um)*dt/dx*rhomr_s*entropym_r+(0.25d0-sr*dt/dx)*dr*entropy_r
     vr_aver=sr*dt/dx*um+(0.25d0-sr*dt/dx)*ur
     pr_aver=sr*dt/dx*pm+(0.25d0-sr*dt/dx)*pr

    else
! The right wave is a rarefaction wave. Since the contact discontinuity now goes to the right, the rarefaction wave also goes to the right.

!     uur_aver=um*dt/dx*rhoml_r*entropy_l+(((um+cmr)-um)*dt/dx*rhomr_r+intr_rho+(0.25d0-(ur+cr)*dt/dx)*dr)*entropy_r			
     vr_aver=(um+cmr)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v
     pr_aver=(um+cmr)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p

    endif 

   else
    uul_aver=((0.25d0-dabs(ul-cl)*dt/dx)*dl+intll_rho_sonic)*entropy_l
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+intll_v_sonic
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+intll_p_sonic
!     rarefaction_rho_l=rarefaction_rho_l+dabs(rho_l_sonic-dl)
!     rarefaction_rho_r=rarefaction_rho_r+dabs(rho_l_sonic-rhoml_r)
    if(pm>pr) then
     entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
     uur_aver=(intlr_rho_sonic+(um-dabs(um-cml))*dt/dx*rhoml_r)*entropy_l+(sr-um)*dt/dx*rhomr_s*entropym_r+(0.25d0-sr*dt/dx)*dr*entropy_r
     vr_aver=(sr-dabs(um-cml))*dt/dx*um+(0.25d0-sr*dt/dx)*ur+intlr_v_sonic
     pr_aver=(sr-dabs(um-cml))*dt/dx*pm+(0.25d0-sr*dt/dx)*pr+intlr_p_sonic
!      shock_rho_r=shock_rho_r+dabs(rhomr_s-dr)
!      contact_rho_r=contact_rho_r+dabs(rhoml_r-rhomr_s)		    
    else
     uur_aver=(intlr_rho_sonic+(um-dabs(um-cml))*dt/dx*rhoml_r)*entropy_l+(((um+cmr)-um)*dt/dx*rhomr_r+intr_rho+(0.25d0-(ur+cr)*dt/dx)*dr)*entropy_r
     vr_aver=((um+cmr)-(um-cml))*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intlr_v_sonic+intr_v
     pr_aver=((um+cmr)-(um-cml))*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intlr_p_sonic+intr_p	
!      rarefaction_rho_r=rarefaction_rho_r+dabs(rhomr_r-dr)
!      contact_rho_r=contact_rho_r+dabs(rhoml_r-rhomr_r)		 
    endif
   endif
  endif
 end if
else
 if(pm>pr) then
  if(sr<0.0d0) then
   uur_aver=0.25d0*dr*entropy_r 
   vr_aver=0.25d0*ur    
   pr_aver=0.25d0*pr 
   entropyl=entropyl+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
!    shock_rho_l=shock_rho_l+dabs(dr-rhomr_s)
   if(pm>pl) then
    entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
    uul_aver=(0.25d0-dabs(sl)*dt/dx)*dl*entropy_l+(dabs(sl)-dabs(um))*dt/dx*rhoml_s*entropym_l+(dabs(um)-dabs(sr))*dt/dx*rhomr_s*entropym_r+dabs(sr)*dt/dx*dr*entropy_r
    vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+(dabs(sl)-dabs(sr))*dt/dx*um+dabs(sr)*dt/dx*ur 
    pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+(dabs(sl)-dabs(sr))*dt/dx*pm+dabs(sr)*dt/dx*pr 
!     shock_rho_l=shock_rho_l+dabs(dl-rhoml_s)	
!     contact_rho_l=contact_rho_l+dabs(rhomr_s-rhoml_s) 	 
   else
    uul_aver=((0.25d0-dabs(ul-cl)*dt/dx)*dl+intl_rho+(dabs(um-cml)-dabs(um))*dt/dx*rhoml_r)*entropy_l+(dabs(um)-dabs(sr))*dt/dx*rhomr_s*entropym_r+dabs(sr)*dt/dx*dr*entropy_r
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+(dabs(um-cml)-dabs(sr))*dt/dx*um+dabs(sr)*dt/dx*ur+intl_v
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+(dabs(um-cml)-dabs(sr))*dt/dx*pm+dabs(sr)*dt/dx*pr+intl_p 
!     rarefaction_rho_l=rarefaction_rho_l+dabs(rhoml_r-dl)
!     contact_rho_l=contact_rho_l+dabs(rhoml_r-rhomr_s) 		 
   endif
  else
   entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
   uur_aver=sr*dt/dx*rhomr_s*entropym_r+(0.25d0-sr*dt/dx)*dr*entropy_r 
   vr_aver=sr*dt/dx*um+(0.25d0-sr*dt/dx)*ur    
   pr_aver=sr*dt/dx*pm+(0.25d0-sr*dt/dx)*pr 
!    shock_rho_r=shock_rho_r+dabs(dr-rhomr_s)
   if(pm>pl) then
    entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
    uul_aver=(0.25d0-dabs(sl)*dt/dx)*dl*entropy_l+(dabs(sl)-dabs(um))*dt/dx*rhoml_s*entropym_l+dabs(um)*dt/dx*rhomr_s*entropym_r
    vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+dabs(sl)*dt/dx*um 
    pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+dabs(sl)*dt/dx*pm
!     shock_rho_l=shock_rho_l+dabs(dl-rhoml_s)
!     contact_rho_l=contact_rho_l+dabs(rhomr_s-rhoml_s) 	
   else
    uul_aver=((0.25d0-dabs(ul-cl)*dt/dx)*dl+intl_rho+(dabs(um-cml)-dabs(um))*dt/dx*rhoml_r)*entropy_l+dabs(um)*dt/dx*rhomr_s*entropym_r
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+dabs(um-cml)*dt/dx*um+intl_v
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+dabs(um-cml)*dt/dx*pm+intl_p
!     rarefaction_rho_l=rarefaction_rho_l+dabs(rhoml_r-dl)
!     contact_rho_l=contact_rho_l+dabs(rhoml_r-rhomr_s) 
   endif
  endif
 else
  if((ur+cr)<0.0d0) then
   uur_aver=0.25d0*dr*entropy_r
   vr_aver=0.25d0*ur    
   pr_aver=0.25d0*pr 
!    rarefaction_rho_l=rarefaction_rho_l+dabs(rhomr_r-dr)
   if(pm>pl) then
    entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
    uul_aver=(0.25d0-dabs(sl)*dt/dx)*dl*entropy_l+(dabs(sl)-dabs(um))*dt/dx*rhoml_s*entropym_l+((dabs(um)-dabs(um+cmr))*dt/dx*rhomr_r+intr_rho+dabs(ur+cr)*dt/dx*dr)*entropy_r
    vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+(dabs(sl)-dabs(um+cmr))*dt/dx*um+dabs(ur+cr)*dt/dx*ur+intr_v
    pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+(dabs(sl)-dabs(um+cmr))*dt/dx*pm+dabs(ur+cr)*dt/dx*pr+intr_p	    
!     shock_rho_l=shock_rho_l+dabs(rhoml_s-dl)
!     contact_rho_l=contact_rho_l+dabs(rhoml_s-rhomr_r) 	 
   else
    uul_aver=((0.25d0-dabs(ul-cl)*dt/dx)*dl+intl_rho+(dabs(um-cml)-dabs(um))*dt/dx*rhoml_r)*entropy_l+((dabs(um)-dabs(um+cmr))*dt/dx*rhomr_r+intr_rho+dabs(ur+cr)*dt/dx*dr)*entropy_r
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+(dabs(um-cml)-dabs(um+cmr))*dt/dx*um+dabs(ur+cr)*dt/dx*ur+intl_v+intr_v
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+(dabs(um-cml)-dabs(um+cmr))*dt/dx*pm+dabs(ur+cr)*dt/dx*pr+intl_p+intr_p
!     rarefaction_rho_l=rarefaction_rho_l+dabs(rhoml_r-dl)
!     contact_rho_l=contact_rho_l+dabs(rhoml_r-rhomr_r) 	 		  
   endif	     
  else
   if((um+cmr)>0.0d0) then
    uur_aver=((um+cmr)*dt/dx*rhomr_r+intr_rho+(0.25d0-(ur+cr)*dt/dx)*dr)*entropy_r 
    vr_aver=(um+cmr)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v  
    pr_aver=(um+cmr)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p
!     rarefaction_rho_r=rarefaction_rho_r+dabs(rhomr_r-dr)
    if(pm>pl) then
     entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
     uul_aver=(0.25d0-dabs(sl)*dt/dx)*dl*entropy_l+(dabs(sl)-dabs(um))*dt/dx*rhoml_s*entropym_l+dabs(um)*dt/dx*rhomr_r*entropy_r
     vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+dabs(sl)*dt/dx*um
     pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+dabs(sl)*dt/dx*pm
!      shock_rho_l=shock_rho_l+dabs(dl-rhoml_s)
!      contact_rho_l=contact_rho_l+dabs(rhoml_s-rhomr_r) 
    else
     uul_aver=((0.25d0-dabs(ul-cl)*dt/dx)*dl+intl_rho+(dabs(um-cml)-dabs(um))*dt/dx*rhoml_r)*entropy_l+dabs(um)*dt/dx*rhomr_r*entropy_r
     vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+dabs(um-cml)*dt/dx*um+intl_v
     pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+dabs(um-cml)*dt/dx*pm+intl_p
!      rarefaction_rho_l=rarefaction_rho_l+dabs(rhoml_r-dl)
!      contact_rho_l=contact_rho_l+dabs(rhoml_r-rhomr_r) 
    endif	
   else 
    uur_aver=(intrr_rho_sonic+(0.25d0-(ur+cr)*dt/dx)*dr)*entropy_r
    vr_aver=(0.25d0-(ur+cr)*dt/dx)*ur+intrr_v_sonic    
    pr_aver=(0.25d0-(ur+cr)*dt/dx)*pr+intrr_p_sonic
!     rarefaction_rho_r=rarefaction_rho_r+dabs(rho_r_sonic-dr)
!     rarefaction_rho_l=rarefaction_rho_l+dabs(rho_r_sonic-rhomr_r)
    if(pm>pl) then
     entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
     uul_aver=(0.25d0-dabs(sl)*dt/dx)*dl*entropy_l+(dabs(sl)-dabs(um))*dt/dx*rhoml_s*entropym_l+((dabs(um)-dabs(um+cmr))*dt/dx*rhomr_r+intrl_rho_sonic)*entropy_r   
     vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+(dabs(sl)-dabs(um+cmr))*dt/dx*um+intrl_v_sonic   
     pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+(dabs(sl)-dabs(um+cmr))*dt/dx*pm+intrl_p_sonic   		    
!      shock_rho_l=shock_rho_l+dabs(dl-rhoml_s)
!      contact_rho_l=contact_rho_l+dabs(rhoml_s-rhomr_r) 
    else
     uul_aver=((0.25d0-dabs(ul-cl)*dt/dx)*dl+intl_rho+(dabs(um-cml)-dabs(um))*dt/dx*rhoml_r)*entropy_l+((dabs(um)-dabs(um+cmr))*dt/dx*rhomr_r+intrl_rho_sonic)*entropy_r   
     vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+(dabs(um-cml)-dabs(um+cmr))*dt/dx*um+intrl_v_sonic+intl_v   
     pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+(dabs(um-cml)-dabs(um+cmr))*dt/dx*pm+intrl_p_sonic+intl_p    	
!      rarefaction_rho_l=rarefaction_rho_l+dabs(rhoml_r-dl)
!      contact_rho_l=contact_rho_l+dabs(rhoml_r-rhomr_r) 
    endif  	      
   endif
  endif
 endif
end if

end  subroutine aver_velocityandpressure
!***************************************************************************************


!***************************************************************************************      
subroutine compute_rarefaction_left(a,b,value_intl_rho,value_intl_v,value_intl_p) 

integer :: i
real*8 :: G1, G2, G3, G4, G5, G6, G7, G8, G9
real*8 :: DL, UL, PL, CL, DR, UR, PR, CR, PM, UM
real*8 :: a,b,value_intl_rho,value_intl_v,value_intl_p
!real*8,dimension(-1:1) :: gp,intp_rho,intp_v,intp_p

COMMON /GAMMAS/ G1, G2, G3, G4, G5, G6, G7, G8, G9
COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR, PM, UM

value_intl_rho=dt/dx*dl*cl*((g5+g6/cl*(ul-a))**(1.0d0/g6)-(g5+g6/cl*(ul-b))**(1.0d0/g6))
value_intl_v=dt/dx*g5*0.5d0*((cl+g7*ul+b)**2.0d0-(cl+g7*ul+a)**2.0d0)
value_intl_p=dt/dx*(gamma+1.0d0)/(3.0d0*gamma-1.0d0)*pl*cl*((g5+g6/cl*(ul-a))**(g3+1.0d0)-(g5+g6/cl*(ul-b))**(g3+1.0d0))

end subroutine compute_rarefaction_left
!***************************************************************************************


!***************************************************************************************
subroutine compute_rarefaction_right(a,b,value_intr_rho,value_intr_v,value_intr_p) 

integer :: i
real*8 :: G1, G2, G3, G4, G5, G6, G7, G8, G9
real*8 :: DL, UL, PL, CL, DR, UR, PR, CR, PM, UM
real*8 :: a,b,value_intr_rho,value_intr_v,value_intr_p
real*8,dimension(-1:1) :: gp,intp_rho,intp_v,intp_p

COMMON /GAMMAS/ G1, G2, G3, G4, G5, G6, G7, G8, G9
COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR, PM, UM

value_intr_rho=dt/dx*dr*cr*((g5-g6/cr*(ur-b))**(1.0d0/g6)-(g5-g6/cr*(ur-a))**(1.0d0/g6))
value_intr_v=dt/dx*g5*0.5d0*((-cr+g7*ur+b)**2.0d0-(-cr+g7*ur+a)**2.0d0)
value_intr_p=dt/dx*(gamma+1.0d0)/(3.0d0*gamma-1.0d0)*pr*cr*((g5-g6/cr*(ur-b))**(g3+1.0d0)-(g5-g6/cr*(ur-a))**(g3+1.0d0))

end subroutine compute_rarefaction_right
!***************************************************************************************


!***************************************************************************************  
end module average