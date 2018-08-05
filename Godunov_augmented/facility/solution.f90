module solution
! The definition of numerical solution. It includes the original quantities, augmented quantities, half-steps. 
! The numerical fluxes are also defined.

use grid_and_parameters

real*8, dimension(-ighost+1:cells_number+ighost,3) :: u
real*8, dimension(-ighost+1:cells_number+ighost) :: uu, v_aver, p_aver
real*8, dimension(-ighost+1:cells_number+ighost) :: d_rho, d_v, d_p, entropy_increase
real*8, dimension(-ighost+1:cells_number+ighost) :: d_rho_tvd, d_v_tvd, d_p_tvd !, d_rho_contact, s_rho, s_entropy, s_p,uu_aver

real*8,dimension(-ighost+1:cells_number+ighost-2,3) :: flux_u 
real*8,dimension(-ighost+1:cells_number+ighost-2) :: flux_uu   

end module solution