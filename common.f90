module mod_common

  use mod_param

  implicit none

  real :: lx,ly,lz ! domain size
  real :: vol_frac ! volume fraction
  real :: vol_box, dt
  real :: t_w0,t_w1, t_nb0,t_nb1,t_nb_tot, t_comp0,t_comp1,t_comp_tot  ! timing

  integer :: istep

  real :: F_elst  ! electrostatic repulsion force scale

  ! particle ID, radius
  integer, dimension(np) :: p_id
  real,    dimension(np) :: radii
  
  ! positions, velocities, accelerations (both translational and rotational)
  real, dimension(3*np) :: pp,po, pv,pw, pa,pb   

  ! matrices containing pairwise info of interacting neighbors
  integer, dimension(np)      :: num_nb           ! number of neighbors for each particle
  integer, dimension(  30,np) :: nb_ind           ! neighbor-j indice (assuming num_nb < 30 #NeedVerification)
  real,    dimension(  np,np) :: d_ij_mx          ! distance matrix
  real,    dimension(3,np,np) :: r_ij_mx, v_ij_mx ! position/velocity vector matrix
  real,    dimension(3,   np) :: F_stk_mx                            ! Stokes drag matrix
  real,    dimension(3,np,np) :: F_rep_mx, F_lub_mx,F_col_mx,F_f_mx  ! rep/lub/col/f force matrix
  integer, dimension(  np,np) :: fr_prev_mx, fr_curr_mx              ! frictional contact flags
  real,    dimension(3,np,np) :: xi_prev_mx, xi_curr_mx              ! tangential stretch vectors
  real                        :: F_crit                              ! critical-load friction
  
  ! bulk stress tensor
  real, dimension(3,3) :: stress

  ! Stokes coefficients, mass, moment of inertia
  real, dimension(2) :: F_stk_factor, T_stk_factor, masses, mois

  ! lubrication coefficients
  real, dimension(-1:1) :: R_xiia1,R_xiia2, R_yiia,R_yiib,R_yjib,R_yiic,R_yijc,R_yjjc
    
  ! periodic image shift
  real :: y_shift, v_shift

  ! time in the simulation
  real :: t9

  ! output intervals
  integer :: iout0, iout1, iout2, iout3
  
end module mod_common



module mod_common_output

  use mod_param

  implicit none

  integer :: l_col,l_bad  ! number of collision/bad pairs
  integer :: club,ccol    ! lub/col counter of particle 1
  
  real :: h_avr, h_max  ! particle overlap

  real, dimension(3) :: F_stk_1,F_rep_1,F_lub_1,F_col_1,F_f_1  ! forces of particle 1
  
  real, dimension(3,3) :: s_rep, s_lub, s_ctt  ! repulsion-/lubrication-/contact-force stress tensor

end module mod_common_output
