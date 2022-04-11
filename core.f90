module mod_core

  use OMP_LIB
  use mod_param
  use mod_common
  use mod_output
  
  implicit none
  
  private
  public time_marching

contains
  
  subroutine time_marching(istart,iend)

    integer, intent(in) :: istart,iend

    call init_marching
    
    do istep = istart+1,iend
       t9 = t9 + dt
       call progress(istart,iend)
       call integrator
       call data_output(t9)
    enddo
    
    return
  end subroutine time_marching
  

  !========================================================================

  
  subroutine init_marching

    ! Initialize neighbor-list, force/torque, and output.

    call gather_nb_info(   pp,pv)
    call comp_force_torque(pp,pv,po,pw, pa,pb)
    if (.not. restart_from_last) then
       call data_output(0.)
    endif
    
    return
  end subroutine init_marching

  !---------------------------------------------------------

  subroutine integrator
    
    ! The main stepper that handles temporal integration.
    
    real, dimension(3*np) :: pv_mid,pw_mid, pv_new,pw_new
    real, dimension(3*np) :: pp_new,po_new, pa_new,pb_new

    select case (integrate_method)
       
    case(1) ! Forward Euler (not tested)

       pp_new = pp + dt*pv + (0.5*dt**2)*pa
       po_new = po + dt*pw + (0.5*dt**2)*pb

       pv_new = pv + dt*pa
       pw_new = pw + dt*pb

       call gather_nb_info(   pp_new,pv_new)
       call comp_force_torque(pp_new,pv_new,po_new,pw_new, pa_new,pb_new)

    case(2) ! (modified) Velocity-Verlet

       pp_new = pp + dt*pv + (0.5*dt**2)*pa
       po_new = po + dt*pw + (0.5*dt**2)*pb

       pv_mid = pv + (0.5*dt)*pa
       pw_mid = pw + (0.5*dt)*pb
       
       call gather_nb_info(   pp_new,pv_mid)
       call comp_force_torque(pp_new,pv_mid,po_new,pw_mid, pa_new,pb_new)
    
       pv_new = pv + dt*(pa + pa_new)/2.
       pw_new = pw + dt*(pb + pb_new)/2.
       
    end select

    pp = pp_new
    po = po_new
    pv = pv_new
    pw = pw_new
    pa = pa_new
    pb = pb_new       

    fr_prev_mx = fr_curr_mx
    xi_prev_mx = xi_curr_mx

    call periodicity
    
    return
  end subroutine integrator

  !------------------------------------------------------
  
  subroutine gather_nb_info(pos,vel)
  
    use mod_common_output
    
    ! For each particle (i), find its interacting neighbors (j>i), and
    !
    ! (1) store the neighbor indices in the matrix nb_ind(k,i),
    !     where k is the counting index;
    ! (2) store the distance between i and j in the matrix d_ij_mx(j,i);
    ! (3) store the relative position vector in the matrix r_ij_mx(1:3,j,i);
    ! (4) store the relative velocity vector in the matrix v_ij_mx(1:3,j,i).
  
    real, dimension(3*np), intent(in)  :: pos,vel
  
    integer :: i,j,k
    real    :: a_i,a_j, d_ij,h_ij,hh
  
    real, dimension(3) :: pos_i,pos_j, vel_i,vel_j

    !call cpu_time(t_nb0)
    t_nb0 = omp_get_wtime() 
    
    l_col = 0; l_bad = 0     ! pair counts
    h_max = 0.; h_avr = 0.   ! maximal and average surface gaps
    
    !$OMP parallel default(private) &
    !$OMP shared(radii,pos,vel,nb_ind,num_nb,d_ij_mx,r_ij_mx,v_ij_mx) &
    !$OMP reduction(+:l_col,l_bad,h_avr) &
    !$OMP reduction(max:h_max) 
    !$OMP do schedule(dynamic)
    do i = 1,np
       a_i        = radii(i)
       pos_i(1:3) = pos(3*i-2 : 3*i)
       vel_i(1:3) = vel(3*i-2 : 3*i)
       
       k = 0
       do j = i+1,np
          a_j        = radii(j)
          pos_j(1:3) = pos(3*j-2 : 3*j)
          vel_j(1:3) = vel(3*j-2 : 3*j)
          
          call comp_least_dist(a_i,a_j,pos_i,pos_j,vel_j,d_ij)

          h_ij = d_ij-(a_i+a_j)  ! surface gap
          
          interacting_neighbor: if (h_ij < h_lub_out) then             
             k = k+1
             nb_ind(     k,i) = j                        ! index
             d_ij_mx(    j,i) = d_ij                     ! distance
             r_ij_mx(1:3,j,i) = pos_j(1:3) - pos_i(1:3)  ! position vector (from i to j)
             v_ij_mx(1:3,j,i) = vel_i(1:3) - vel_j(1:3)  ! velocity vector (i relative to j)

             contact_pair: if (h_ij < 0.) then  
                l_col = l_col+1
                hh = abs(h_ij)/a_1
                h_avr = h_avr +hh
                if (hh > h_max) h_max = hh                
                if (hh > 0.05) then
                   l_bad = l_bad +1
                   !#call output_diagnosis(a_i,a_j,hh,pos_i,pos_i+r_ij) 
                endif                
             endif contact_pair
          endif interacting_neighbor
          
       enddo
       num_nb(i) = k  ! number of neighbors
       
    enddo
    !$OMP end parallel
    
    !call cpu_time(t_nb1)
    t_nb1 = omp_get_wtime() 
    t_nb_tot = t_nb_tot + (t_nb1-t_nb0)
    
    return
  end subroutine gather_nb_info
  
  !------------------------------------------------------
  
  subroutine comp_force_torque(pos,vel,ori,omg, acc,bta)

    use mod_common_output
  
    ! Input:  pos    particle center position
    !         vel    particle translational velocity
    !         ori    particle orientation
    !         omg    particle angular velocity
    !
    ! Output: acc    particle translational acceleration
    !         bta    particle angular acceleration
    !
    ! (Implementing OpenMP)
  
    real, dimension(3*np), intent(in)  :: pos,vel,ori,omg
    real, dimension(3*np), intent(out) :: acc,bta
  
    integer :: i,j,k
    integer :: ix,iy,iz, jx,jy,jz
  
    real :: a_i,a_j  
    
    real, dimension(3) :: pos_i, vel_i,  omg_i
    real, dimension(3) :: r_ij,  vel_ij, omg_j   
    
    !call cpu_time(t_comp0)
    t_comp0 = omp_get_wtime() 
    
    !---- Initialization ----!
    
    acc(:) = 0.; bta(:) = 0. ! translational/angular accelerations
    fr_curr_mx(:,:) = 0      ! friction flag (assume no frictional contact until detacted)
    xi_curr_mx(:,:,:) = 0.   ! tangential stretch

    !---- Looping ----!
    
    !$OMP parallel default(private) &
    !$OMP shared(radii,pos,vel,omg,nb_ind,num_nb,r_ij_mx,v_ij_mx) &
    !$OMP reduction(+:acc,bta)
    !$OMP do schedule(static) 
    center_particle_i: do i = 1,np
       
       ix = 3*i-2
       iy = ix+1
       iz = iy+1

       a_i = radii(i)

       pos_i(1:3) = pos(ix:iz)
       vel_i(1:3) = vel(ix:iz)
       omg_i(1:3) = omg(ix:iz)
              
       call comp_Stokes(i,pos_i,vel_i,omg_i, acc,bta) ! single-body force/torque
  
       interacting_neighbor_j: do k = 1,num_nb(i)  ! higher-index neighbors (j>i)
  
          j = nb_ind(k,i)

          jx = 3*j-2
          jy = jx+1
          jz = jy+1

          a_j = radii(j)

          r_ij(  1:3) = r_ij_mx(1:3,j,i) ! position vector (from i to j)
          vel_ij(1:3) = v_ij_mx(1:3,j,i) ! relative velocities (i relative to j)
          omg_j( 1:3) = omg(jx:jz)

          call comp_pair_interactions(i,a_i,omg_i, j,a_j,omg_j, r_ij,vel_ij, acc,bta)
          
       enddo interacting_neighbor_j
    enddo center_particle_i
    !$OMP end parallel

    !call cpu_time(t_comp1)
    t_comp1 = omp_get_wtime() 
    t_comp_tot = t_comp_tot + (t_comp1-t_comp0)
    
    return
  end subroutine comp_force_torque

  !----
  
  subroutine comp_Stokes(i,pos_i,vel_i,omg_i, acc,bta)

    ! Compute the Stokes drag and torque

    integer,               intent(in)     :: i
    real, dimension(3),    intent(in)     :: pos_i, vel_i, omg_i
    real, dimension(3*np), intent(inout)  :: acc, bta

    integer :: ix,iy,iz
    
    real, dimension(3) :: vel_infty, omg_infty, F_stk, T_stk 

    ix = 3*i-2; iy = ix+1; iz = iy+1
       
    ! undisturbed translational/angular velocities
    vel_infty(1) = 0.
    if (oscillatory_shear) then
       vel_infty(2) = (pos_i(3)-lz/2.)*(amp0*lz*freq/lz)*cos(freq*t9)
    else
       vel_infty(2) = (pos_i(3)-lz/2.)*shear_rate
    end if
    vel_infty(3) = 0.
    if (oscillatory_shear) then
       omg_infty(1) = -(amp0*lz*freq/lz)*cos(freq*t9)/2.
    else
       omg_infty(1) = -shear_rate/2.
    end if
    omg_infty(2) = 0.
    omg_infty(3) = 0.

    F_stk(1:3) = F_stk_factor(p_id(i))*(vel_i(1:3)-vel_infty(1:3))
    T_stk(1:3) = T_stk_factor(p_id(i))*(omg_i(1:3)-omg_infty(1:3))
  
    ! add Stokes drag/torque (one-body) to the total force/torque
    
    acc(ix:iz) = acc(ix:iz)-F_stk(1:3)/masses(p_id(i))
    bta(ix:iz) = bta(ix:iz)-T_stk(1:3)/mois(  p_id(i))

    F_stk_mx(:,i) = F_stk(:)  ! statistics in a big matrix

    return
  end subroutine comp_Stokes

  !----
  
  subroutine comp_pair_interactions(i,a_i,omg_i, j,a_j,omg_j, r_ij,vel_ij, acc,bta)

    ! Compute the lubrication, collision, friction, and physiochemical interactions

    integer,               intent(in)     :: i,j
    real,                  intent(in)     :: a_i,a_j
    real, dimension(3),    intent(in)     :: omg_i,omg_j, r_ij, vel_ij
    real, dimension(3*np), intent(inout)  :: acc, bta

    integer :: ix,iy,iz, jx,jy,jz, ind
    
    real :: d_ij, h_ij, n_x, n_y, n_z
    real :: hh, delta,delta_inv,log_delta_inv
    real :: X_iiA,Y_iiA,Y_iiB,Y_jiB,Y_iiC,Y_ijC,Y_jjC  
    real :: F_col_mag, F_f_test, F_f_slid, buf0
    
    real, dimension(3) :: n_ij
    real, dimension(3) :: P_n_x, P_n_y, P_n_z
    real, dimension(3) :: P_t_x, P_t_y, P_t_z 
    real, dimension(3) :: vel_ij_n, vel_ij_t 
    real, dimension(3) :: omgxn_i, omgxn_j, buf_vec   
    real, dimension(3) :: F_lub,F_col,F_f, T_lub,T_col    

    ix = 3*i-2; iy = ix+1; iz = iy+1
    jx = 3*j-2; jy = jx+1; jz = jy+1
    
    !-- Geometrical/kinematic terms --!
    
    d_ij = d_ij_mx(j,i)     ! particle center distance
    h_ij = d_ij -(a_i+a_j)  ! particle surface gap (can be negative)

    ! unit normal vectors (from i to j)
    n_ij(1:3) = r_ij(1:3)/d_ij  
    n_x=n_ij(1);n_y=n_ij(2);n_z=n_ij(3)
  
    ! normal projection matrix (nn)
    P_n_x(1) = n_x*n_x; P_n_x(2) = n_x*n_y; P_n_x(3) = n_x*n_z
    P_n_y(1) = n_x*n_y; P_n_y(2) = n_y*n_y; P_n_y(3) = n_y*n_z
    P_n_z(1) = n_x*n_z; P_n_z(2) = n_y*n_z; P_n_z(3) = n_z*n_z
  
    ! tangential projection matrix (I-nn)
    P_t_x = -P_n_x; P_t_y = -P_n_y; P_t_z = -P_n_z
    P_t_x(1) = 1. + P_t_x(1)
    P_t_y(2) = 1. + P_t_y(2)
    P_t_z(3) = 1. + P_t_z(3)

    ! normal velocity (x,y,z components)
    vel_ij_n(1) = dot_product(P_n_x,vel_ij) 
    vel_ij_n(2) = dot_product(P_n_y,vel_ij) 
    vel_ij_n(3) = dot_product(P_n_z,vel_ij) 

    ! tangential velocity (x,y,z components) [The rotation is NOT added (intentional).]
    vel_ij_t(1) = dot_product(P_t_x,vel_ij) 
    vel_ij_t(2) = dot_product(P_t_y,vel_ij) 
    vel_ij_t(3) = dot_product(P_t_z,vel_ij) 
  
    omgxn_i = cross(omg_i,n_ij); omgxn_j = cross(omg_j,n_ij) 

    
    !-- Electrostatic and/or van der Waals --!

    if (elst_repl) then
       call comp_elst_vdw(i,j,a_i,a_j,h_ij,n_ij, acc)
    endif
    

    !-- Lubrication (always on) --!
                
    delta = 2.*max(h_ij,h_lub_inn)/(a_i+a_j)  ! regularized nondimensional gap
    delta_inv = 1./delta
    log_delta_inv = log(delta_inv)
  
    ! scalar resistances
    ind = p_id(j) - p_id(i) ! the index of the R arrays (= 1 or 0 or -1)

    X_iiA = (a_i )*(R_xiia2(ind)*log_delta_inv + R_xiia1(ind)*delta_inv)
    Y_iiA = (a_i   )*R_yiia(ind)*log_delta_inv
    Y_iiB = (a_i**2)*R_yiib(ind)*log_delta_inv
    Y_jiB = (a_j**2)*R_yjib(ind)*log_delta_inv
    Y_iiC = (a_i**3)*R_yiic(ind)*log_delta_inv
    Y_ijC = (a_i**3)*R_yijc(ind)*log_delta_inv
    Y_jjC = (a_j**3)*R_yjjc(ind)*log_delta_inv
                
    ! lubrication force (all three components)
    F_lub(1:3) = X_iiA*(-vel_ij_n(1:3)) + Y_iiA*(-vel_ij_t(1:3)) &
         +       Y_iiB*omgxn_i(1:3)     + Y_jiB*omgxn_j(1:3)

    ! add lubrication force to the total force (i,j) (equal and opposite)
    acc(ix:iz) = acc(ix:iz) + F_lub(1:3)/masses(p_id(i))
    acc(jx:jz) = acc(jx:jz) - F_lub(1:3)/masses(p_id(j))

    F_lub_mx(1:3,j,i) = F_lub(1:3)
    
    ! lubrication torque (all three components)
    buf_vec = cross(n_ij,vel_ij)

    ! acting on i by j
    T_lub(1) = Y_iiB*buf_vec(1) -Y_iiC*dot_product(P_t_x,omg_i)-Y_ijC*dot_product(P_t_x,omg_j)
    T_lub(2) = Y_iiB*buf_vec(2) -Y_iiC*dot_product(P_t_y,omg_i)-Y_ijC*dot_product(P_t_y,omg_j)
    T_lub(3) = Y_iiB*buf_vec(3) -Y_iiC*dot_product(P_t_z,omg_i)-Y_ijC*dot_product(P_t_z,omg_j)

    bta(ix:iz) = bta(ix:iz) + T_lub(1:3)/mois(p_id(i))

    ! acting on j by i
    T_lub(1) = Y_jiB*buf_vec(1) -Y_ijC*dot_product(P_t_x,omg_i)-Y_jjC*dot_product(P_t_x,omg_j)
    T_lub(2) = Y_jiB*buf_vec(2) -Y_ijC*dot_product(P_t_y,omg_i)-Y_jjC*dot_product(P_t_y,omg_j)
    T_lub(3) = Y_jiB*buf_vec(3) -Y_ijC*dot_product(P_t_z,omg_i)-Y_jjC*dot_product(P_t_z,omg_j)

    bta(jx:jz) = bta(jx:jz) + T_lub(1:3)/mois(p_id(j))

    
    !-- Collision/Friction --!

    if (h_ij < 0.) then

       ! normal collision froce: spring and dashpot (all three components, acting on i by j)
       F_col(1:3) = k_n*h_ij*n_ij(1:3) - dp_n*dot_product(n_ij,vel_ij)*n_ij(1:3)

       friction: if (frictional) then

          F_col_mag = sqrt_3d(F_col)

          if (F_col_mag > F_crit) then
             
             fr_curr_mx(j,i) = 1 ! activate friction (pair i,j)

             if (fr_prev_mx(j,i) == 0) then
                xi_curr_mx(:,j,i) = 0. ! reset accumulated tangential stretch
             else
                ! tangential velocity vector (the complete expression, i relative to j)
                buf_vec = a_i*omgxn_i + a_j*omgxn_j ! relative tang. vel. due to rotation
                vel_ij_t(1) = vel_ij_t(1) + dot_product(P_t_x, buf_vec) 
                vel_ij_t(2) = vel_ij_t(2) + dot_product(P_t_y, buf_vec) 
                vel_ij_t(3) = vel_ij_t(3) + dot_product(P_t_z, buf_vec)
  
                ! tangential stretch (i relative to j)
                buf_vec(1:3) = xi_prev_mx(1:3,j,i) + vel_ij_t(1:3)*dt
                xi_curr_mx(1,j,i) = dot_product(P_t_x, buf_vec)
                xi_curr_mx(2,j,i) = dot_product(P_t_y, buf_vec)
                xi_curr_mx(3,j,i) = dot_product(P_t_z, buf_vec)
  
                ! friction force (all three components, acting on i by j)
                F_f(1:3) = - k_t*xi_curr_mx(1:3,j,i)
                F_f_test = sqrt_3d(F_f)
                F_f_slid = mu_c*(F_col_mag - F_crit)
                      
                if (F_f_test > F_f_slid) then             ! sliding
                   F_f(1:3) = F_f(1:3)/F_f_test*F_f_slid  ! apply Coulomb's law of friction
                   xi_curr_mx(1:3,j,i) = - F_f(1:3)/k_t   ! saturate tangential stretch
                endif

                F_f_mx(1:3,j,i) = F_f(1:3)

                ! final collision force (all three components, acting on i by j)
                F_col(1:3) = F_col(1:3) + F_f(1:3)

                ! collision torque (acting on i by j)
                T_col = cross(n_ij,F_f)  ! [now I only calculate it when fric, not true in general] 
                bta(ix:iz) = bta(ix:iz) + T_col(1:3)*a_i/mois(p_id(i))
                bta(jx:jz) = bta(jx:jz) + T_col(1:3)*a_j/mois(p_id(j))

             endif
          endif

       endif friction
       
       ! add contact force to the total force (i,j)
       acc(ix:iz) = acc(ix:iz) + F_col(1:3)/masses(p_id(i))
       acc(jx:jz) = acc(jx:jz) - F_col(1:3)/masses(p_id(j))

       F_col_mx(1:3,j,i) = F_col(1:3)
       
    endif
       

    return
  end subroutine comp_pair_interactions

  !----

  subroutine comp_elst_vdw(i,j,a_i,a_j,h_ij,n_ij, acc)

    ! Compute the electrostatic repulsion and/or van der Waals attraction

    integer,            intent(in)    :: i,j
    real,               intent(in)    :: h_ij, a_i,a_j
    real, dimension(3), intent(in)    :: n_ij
    real, dimension(3), intent(inout) :: acc

    integer :: ix,iy,iz, jx,jy,jz
    
    real               :: buf0, F_rep_mag, eps
    real, dimension(3) :: F_rep

    eps = (3.16e-3)*a_1  ! regularisation term for small separations
    
    if (h_ij < h_lub_out) then  ! activate when gap is small
       
       ix = 3*i-2; iy = ix+1; iz = iy+1
       jx = 3*j-2; jy = jx+1; jz = jy+1
    
       buf0 = (2.*a_i*a_j)/(a_i+a_j)  !! harmonic mean
       F_rep_mag = 0.!!-F_elst*(buf0/a_1)*exp(-max(h_ij,0.)/Debye_len) >>> temporarily disabled

       if (van_der_W) then
          !F_rep_mag = F_rep_mag + C_Hamaker*buf0/(12.*(h_ij**2+0.01*buf0**2))
          F_rep_mag = F_rep_mag + C_Hamaker*buf0/(12.*(h_ij**2+ eps**2))
       endif

       ! acting on i by j (saturates below contact)
       F_rep(1:3) = F_rep_mag*n_ij(1:3)

       ! add to the total force (i,j)
       acc(ix:iz) = acc(ix:iz) + F_rep(1:3)/masses(p_id(i))
       acc(jx:jz) = acc(jx:jz) - F_rep(1:3)/masses(p_id(j))

       F_rep_mx(1:3,j,i) = F_rep(1:3)

    endif

    return
  end subroutine comp_elst_vdw
  
  !------------------------------------------------------------------

  subroutine return_image(xi,yi, xj,yj,zj, xj_im,yj_im,zj_im)

    ! Return three images of j based on quadrant location of i
    ! (applies to doubly periodic BC, regardless of z-coord).

    real, intent(in) :: xi,yi, xj,yj,zj
    real, dimension(3), intent(out) :: xj_im,yj_im,zj_im

    xj_im(1) = xj+lx*sign(1., xi-lx/2.)
    yj_im(1) = yj
    zj_im(1) = zj
    xj_im(2) = xj
    yj_im(2) = yj+ly*sign(1., yi-ly/2.)
    zj_im(2) = zj
    xj_im(3) = xj+lx*sign(1., xi-lx/2.)
    yj_im(3) = yj+ly*sign(1., yi-ly/2.)
    zj_im(3) = zj

    return
  end subroutine return_image

  !------------------------------------------------------------------
  
  subroutine comp_least_dist(a_i,a_j,pos_i,pos_j,vel_j,d_ij)
    
    ! Compute the least center-to-center distance between particles i and j (including 
    ! possible images of j), and output the possibly-modified position and velocity of j.
    ! 
    ! Note that i cannot interact with more than one image of j (including j itself),
    ! because the interaction is short-range.

    real,               intent(in)    :: a_i, a_j
    real, dimension(3), intent(in)    :: pos_i
    real, dimension(3), intent(inout) :: pos_j,vel_j
    real,               intent(out)   :: d_ij

    integer :: k,jj
    real :: drim
    real :: xi,yi,zi, xj,yj,zj
    real, dimension(7) :: xj_im,yj_im,zj_im, vj_im_shift, d_im

    xi = pos_i(1); yi = pos_i(2); zi = pos_i(3)
    xj = pos_j(1); yj = pos_j(2); zj = pos_j(3)

    ! distance between i and the original j
    d_ij = sqrt( (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2 )

    ! check images of j
    xj_im(:) = xj
    yj_im(:) = yj
    zj_im(:) = zj
    vj_im_shift(:) = 0.
    jj = 0
    drim = a_i + a_j + h_lub_out
    
    if (d_ij > a_i + a_j + h_lub_out) then

       ! x- and y-dir (doubly-periodic)
       if ((xi < drim) .or. (xi > lx-drim) .or. (yi < drim) .or. (yi > ly-drim)) then
          call return_image(xi,yi, xj,yj,zj, xj_im(1:3),yj_im(1:3),zj_im(1:3))
       endif

       ! z-dir (shear-periodic)
       if (zi < drim) then
          
          xj_im(4) = xj
          yj_im(4) = mod(ly+yj-y_shift, ly)
          zj_im(4) = zj-lz

          vj_im_shift(4:7) = -v_shift
          
          call return_image(xi,yi, xj_im(4),yj_im(4),zj_im(4), xj_im(5:7),yj_im(5:7),zj_im(5:7))
          
       elseif (zi > lz-drim) then
          
          xj_im(4) = xj
          yj_im(4) = mod(ly+yj+y_shift, ly)
          zj_im(4) = zj+lz

          vj_im_shift(4:7) = v_shift
          
          call return_image(xi,yi, xj_im(4),yj_im(4),zj_im(4), xj_im(5:7),yj_im(5:7),zj_im(5:7))
          
       endif

       ! compute distances between i and the two images of j
       do k = 1,7
          d_im(k) = sqrt( (xi-xj_im(k))**2+(yi-yj_im(k))**2+(zi-zj_im(k))**2 )
          if (d_im(k) < a_i + a_j + h_lub_out) jj = k
       enddo
       
    endif

    ! update position/distance
    if (jj /= 0) then
       pos_j(1) = xj_im(jj)
       pos_j(2) = yj_im(jj)
       pos_j(3) = zj_im(jj)
       vel_j(2) = vel_j(2) + vj_im_shift(jj) ! shift the streamwise velocity component
       d_ij = d_im(jj)
    endif
    
    return
  end subroutine comp_least_dist

  !------------------------------------------------------------------

  subroutine periodicity

    ! Lees-Edwards boundary condition
    
    integer :: i, ix,iy,iz

    if (oscillatory_shear) then
       y_shift = mod(amp0*lz*sin(freq*t9), ly) ! bounded in [0,ly)
       v_shift = (amp0*lz*freq)*cos(freq*t9)
    else
       y_shift = mod(shear_rate*lz*istep*dt, ly) ! bounded in [0,ly)
    end if
    
    do i = 1,np
       
       ix = 3*i-2; iy = ix+1; iz = iy+1

       ! x-dir
       pp(ix) = mod(lx+pp(ix), lx) 

       ! y-dir
       if (pp(iz) > lz) then
          pp(iy) = mod(ly+pp(iy)-y_shift, ly)
          pv(iy) = pv(iy) - v_shift
       elseif (pp(iz) < 0.) then
          pp(iy) = mod(ly+pp(iy)+y_shift, ly)
          pv(iy) = pv(iy) + v_shift
       else
          pp(iy) = mod(ly+pp(iy), ly) 
       endif

       ! z-dir
       pp(iz) = mod(lz+pp(iz), lz)
       
    enddo

    return
  end subroutine periodicity

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  Function cross(a, b)
    
    real, dimension(3), intent(in) :: a, b
    real, dimension(3) :: cross
    
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

  End Function cross
  
  
  Function sqrt_3d(vec) result(mag_vec)
    
    real, dimension(3), intent(in) :: vec
    real :: mag_vec
    
    mag_vec = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

  End Function sqrt_3d

  
end module mod_core
