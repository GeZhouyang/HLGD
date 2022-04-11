module mod_output
  
  use mod_param
  use mod_common
  use mod_common_output
  
  implicit none
  
  private
  public progress, data_output, output_vtk, output_diagnosis, output_timing

contains

  subroutine progress(istart,iend)

    integer, intent(in) :: istart,iend
    integer :: i,j 

    open(17,file=trim(datadir)//'log.txt',position='append')
    if (istep .eq. istart+1) write(17,'(A)') 'Begin time marching...'
    
    i = iend/10
    do j = 1,10
       if (istep == i*j) write(17,'(A,I3,A)') 'progress: ',10*j,'%'
    enddo
    close(17)

    return
  end subroutine progress

  !-------------------------------------------------------
  
  subroutine data_output(t)

    real, intent(in) :: t
    integer :: N_period
    real :: T_os, t_p

    if (mod(istep,iout0) == 0) then
       !!call output_F1
       !call output_all_F(istep)
       call output_gap(istep)
       call output_fabric(istep)
    endif
       
    if (mod(istep,iout1) == 0) then
       !call output_num_nb
       call output_overlap
       call output_stress(t)
    endif
          
    !if (mod(istep,iout2) == 0) then
    !   call output_vtk(istep)
    !endif
    T_os = 2.*pi/freq
    N_period = int(t/T_os)+1 
    t_p = N_period*T_os
    if (t <= t_p .and. t+dt > t_p) then
       call output_vtk(N_period)
    endif

    
    if (mod(istep,iout3) == 0) then
       call output_restart(istep)
    endif
    
    return
  end subroutine data_output

  
  !=====================================================================


  subroutine output_num_nb

    integer :: i
    character(len=7) :: iterchar
    character(len=25) :: filename

    write(iterchar,'(i7.7)') istep/100
    filename = trim(datadir)//'num_nb'//iterchar//'h.vtk'

    open( 21, file=filename, status='replace')
    do i = 1,np
       write(21,'(2I5)') i, num_nb(i)
    enddo
    close(21)

    return
  end subroutine output_num_nb
  
  !-------------------------------------------------------
  
  subroutine output_overlap
        
    open( 21,file=trim(datadir)//'overlap.txt',position='append')
    write(21,'(I10,2ES10.2,2I6)') istep, h_avr/max(l_col,1),h_max,l_bad,l_col
    close(21)

    return
  end subroutine output_overlap
  
  !-------------------------------------------------------
  
  subroutine output_diagnosis(a_i,a_j,hh,pos_i,pos_j)

    real,               intent(in) :: a_i,a_j,hh
    real, dimension(3), intent(in) :: pos_i,pos_j
    
    open( 22,file=trim(datadir)//'diagnosis.txt',position='append')
    write(22,'(I9,9F8.3)') istep, a_i,a_j, hh, &
         pos_i(1),pos_i(2),pos_i(3), pos_j(1),pos_j(2),pos_j(3) 
    close(22)

    return
  end subroutine output_diagnosis

  !-------------------------------------------------------
  
  subroutine output_max_v(vel,velinf)

    real, dimension(3), intent(in) :: vel
    real, intent(in) :: velinf
    real :: c1,c2

    ! non-dimensionalization
    c1 = shear_rate*a_1
    !c2 = F_stk_factor(1)*shear_rate*a_1
    
    open( 22,file=trim(datadir)//'max_vel.txt',position='append')
    write(22,'(I9,4ES12.3)') istep, vel(1)/c1,vel(2)/c1,vel(3)/c1, velinf/c1 
    close(22)

    return
  end subroutine output_max_v

  !-------------------------------------------------------
  
  subroutine output_F1

    ! Compute and output the forces acting on particle 1.
    ! Note, the loops have to be EXACTLY the same as in comp_force_torque.

    integer :: i,j,k
    real :: a_i,a_j, d_ij,h_ij, c2

    club = 0; ccol = 0  ! force-pair counter
    F_lub_1(:) = 0.; F_col_1(:) = 0.; F_f_1(1:3) = 0.

    i = 1; a_i = radii(i)
       
    F_stk_1(1:3) = F_stk_mx(1:3,i)

    do k = 1,num_nb(i)
       j = nb_ind(k,i)

       a_j = radii(j)
       d_ij = d_ij_mx(j,i); h_ij = d_ij -(a_i+a_j)  

       lubrication: if (.true.) then !if (h_ij > 0.) then !
          club = club +1
          F_lub_1(1:3) = F_lub_1(1:3) + F_lub_mx(1:3,j,i)
       endif lubrication

       collision: if (h_ij < 0.) then
          ccol = ccol +1
          F_col_1(1:3) = F_col_1(1:3) + F_col_mx(1:3,j,i)
          friction: if (k_t > 0.) then
             if (fr_curr_mx(j,i) == 1) then
                if (fr_prev_mx(j,i) == 1) then
                   F_f_1(1:3) = F_f_1(1:3) + F_f_mx(1:3,j,i)
                endif
             endif
          endif friction
       endif collision
    enddo
       
    ! non-dimensionalization
    c2 = F_stk_factor(1)*shear_rate*a_1
    
    open( 22,file=trim(datadir)//'max_F.txt',position='append')
    write(22,'(I9,3ES12.3,I3,3ES12.3,I3,6ES12.3)') istep, &
         F_stk_1(1)/c2,F_stk_1(2)/c2,F_stk_1(3)/c2, club, &
         F_lub_1(1)/c2,F_lub_1(2)/c2,F_lub_1(3)/c2, ccol, &
         F_col_1(1)/c2,F_col_1(2)/c2,F_col_1(3)/c2, &
         F_f_1(  1)/c2,F_f_1(  2)/c2,F_f_1(  3)/c2
    close(22)

    return
  end subroutine output_F1

  !-------------------------------------------------------
  
  subroutine output_all_F(istep)

    ! Output the magnitudes of Stokes, lubrication, collision, and friction forces 
    ! acting on all particles.

    integer, intent(in) :: istep
    
    integer :: i,j,k
    real :: a_i,a_j, d_ij,h_ij, c2, buf1,buf2,buf3,buf4,buf5
    character(len=7) :: iterchar
    character(len=40) :: filename

    
    write(iterchar,'(i7.7)') istep/1000
    filename = trim(datadir)//'force'//trim(dsub_num)//'/all_F'//iterchar//'k.txt'

    !club = 0; ccol = 0  ! force-pair counter

    do i = 1,np

       a_i = radii(i)
       
       F_stk_1(1:3) = F_stk_mx(1:3,i)
       buf1 = sqrt_3d(F_stk_1)

       F_rep_1(:) = 0.; F_lub_1(:) = 0.; F_col_1(:) = 0.; F_f_1(1:3) = 0.
       
       do k = 1,num_nb(i)
          j = nb_ind(k,i)

          a_j = radii(j)
          d_ij = d_ij_mx(j,i); h_ij = d_ij -(a_i+a_j)  

          electrostatic_repulsion: if (elst_repl) then
             F_rep_1(1:3) = F_rep_1(1:3) + F_rep_mx(1:3,j,i)
          endif electrostatic_repulsion

          lubrication: if (.true.) then !if (h_ij > 0.) then !
             !club = club +1
             F_lub_1(1:3) = F_lub_1(1:3) + F_lub_mx(1:3,j,i)
          endif lubrication

          collision: if (h_ij < 0.) then
             !ccol = ccol +1
             F_col_1(1:3) = F_col_1(1:3) + F_col_mx(1:3,j,i)
             friction: if (k_t > 0.) then
                if (fr_curr_mx(j,i) == 1) then
                   if (fr_prev_mx(j,i) == 1) then
                      F_f_1(1:3) = F_f_1(1:3) + F_f_mx(1:3,j,i)
                   endif
                endif
             endif friction
          endif collision
       enddo

       buf2 = sqrt_3d(F_lub_1)
       buf3 = sqrt_3d(F_col_1)
       buf4 = sqrt_3d(F_f_1)
       buf5 = sqrt_3d(F_rep_1)

       open( 22,file=filename,position='append')
       write(22,'(I5,F5.1,5ES14.5)') i, a_i, buf1,buf2,buf3,buf4,buf5
       close(22)

    enddo

    return
  end subroutine output_all_F

  !-------------------------------------------------------
  
  subroutine output_gap(istep)

    ! Output the surface gap between interaction particles.

    integer, intent(in) :: istep
    
    integer :: i,j,k
    real :: a_i,a_j, d_ij,h_ij
    character(len=7) :: iterchar
    character(len=40) :: filename

    
    write(iterchar,'(i7.7)') istep/1000
    filename = trim(datadir)//'gap'//trim(dsub_num)//'/gap'//iterchar//'k.txt'

    !club = 0; ccol = 0  ! force-pair counter

    do i = 1,np

       a_i = radii(i)
       
       do k = 1,num_nb(i)
          j = nb_ind(k,i)

          a_j = radii(j)
          d_ij = d_ij_mx(j,i); h_ij = d_ij -(a_i+a_j)  

          open( 29,file=filename,position='append')
          write(29,'(2I5,ES14.5)') i,j, h_ij
          close(29)
          
       enddo
    enddo

    return
  end subroutine output_gap

  !-------------------------------------------------------
  
  subroutine output_fabric(istep)

    ! Output the suspension fabric of a given microstructure.

    integer, intent(in) :: istep
    
    integer :: i,j,k, club
    real    :: rx,ry,rz
    real, dimension(6) :: A_h
    character(len=40) :: filename


    club = 0     ! force-pair counter
    A_h(:) = 0.  ! suspension fabric matrix

    do i = 1,np       
       do k = 1,num_nb(i) ! loop thru all points since they are all within h_lub_out
          j = nb_ind(k,i)

          club = club +1

          rx = r_ij_mx(1,j,i)
          ry = r_ij_mx(2,j,i)
          rz = r_ij_mx(3,j,i)

          A_h(1) = A_h(1) + rx*rx
          A_h(2) = A_h(2) + rx*ry
          A_h(3) = A_h(3) + rx*rz
          A_h(4) = A_h(4) + ry*ry
          A_h(5) = A_h(5) + ry*rz
          A_h(6) = A_h(6) + rz*rz
          
       enddo
    enddo

    if (club > 0) then
       A_h(1) = A_h(1)/club
       A_h(2) = A_h(2)/club
       A_h(3) = A_h(3)/club
       A_h(4) = A_h(4)/club
       A_h(5) = A_h(5)/club
       A_h(6) = A_h(6)/club
    endif

    filename = trim(datadir)//'fabric.txt'

    open( 29,file=filename,position='append')
    write(29,'(I12,6ES14.5)') istep, A_h(1),A_h(2),A_h(3),A_h(4),A_h(5),A_h(6)
    close(29)

    return
  end subroutine output_fabric
  
  !-------------------------------------------------------
  
  subroutine output_stress(t)

    real, intent(in) :: t

    ! Compute and output the bulk stress tensor.
    ! Note, the loops have to be EXACTLY the same as in comp_force_torque.

    integer :: i,j,k, m,n
    real    :: a_i,a_j, d_ij,h_ij 
    real, dimension(3) :: r_ij, F_rep,F_lub,F_col
    real, dimension(3,3) :: ss

    s_rep(:,:)  = 0.
    s_lub(:,:)  = 0.
    s_ctt(:,:)  = 0.
    stress(:,:) = 0.  
    
    do i = 1,np
       a_i = radii(i)
       
       do k = 1,num_nb(i)
          j = nb_ind(k,i)

          a_j = radii(j)
          r_ij(1:3) = r_ij_mx(1:3,j,i)
          d_ij = d_ij_mx(j,i); h_ij = d_ij -(a_i+a_j)  

          electrostatic_repulsion: if (elst_repl) then
             ! stresses due to electrostatic repulsion (symmetric)
             F_rep(1:3) = F_rep_mx(1:3,j,i)
             do m = 1,3
                s_rep(m,1:3) = s_rep(m,1:3) + F_rep(m)*r_ij(1:3)
             enddo
          endif electrostatic_repulsion

          lubrication: if (.true.) then !if (h_ij > 0.) then !
             ! stresses due to pair-wise lubrication (symmetric)
             F_lub(1:3) = F_lub_mx(1:3,j,i)
             do n = 1,3
                do m = 1,3
                   s_lub(m,n) = s_lub(m,n) +(F_lub(m)*r_ij(n) + F_lub(n)*r_ij(m))*0.5
                enddo
             enddo
          endif lubrication

          collision: if (h_ij < 0.) then
             ! stresses due to contact force (symmetric)
             F_col(1:3) = F_col_mx(1:3,j,i)
             do m = 1,3
                s_ctt(m,1:3) = s_ctt(m,1:3) + F_col(m)*r_ij(1:3)
             enddo
          endif collision
       enddo
    enddo

    ! bulk stress tensor (adding the per-particle hydro. stresslet)
    stress(:,:) = (s_rep(:,:) + s_lub(:,:) + s_ctt(:,:))/vol_box
    if (oscillatory_shear) then
       stress(2,3) = stress(2,3) + mu_f*(amp0*freq*cos(freq*t))*(1.+2.5*vol_frac) 
       stress(3,2) = stress(3,2) + mu_f*(amp0*freq*cos(freq*t))*(1.+2.5*vol_frac)
    else
       stress(2,3) = stress(2,3) + mu_f*shear_rate*(1.+2.5*vol_frac) 
       stress(3,2) = stress(3,2) + mu_f*shear_rate*(1.+2.5*vol_frac) 
    endif
    
    ss = stress
    open( 23,file=trim(datadir)//'bulk_stress_all.txt',position='append')
    write(23,'(I10,9ES12.3)') istep, &
         ss(1,1),ss(1,2),ss(1,3), &
         ss(2,1),ss(2,2),ss(2,3), &
         ss(3,1),ss(3,2),ss(3,3)
    close(23)

    ss = s_lub/vol_box
    open( 23,file=trim(datadir)//'bulk_stress_lub.txt',position='append')
    write(23,'(I10,9ES12.3)') istep, &
         ss(1,1),ss(1,2),ss(1,3), &
         ss(2,1),ss(2,2),ss(2,3), &
         ss(3,1),ss(3,2),ss(3,3)
    close(23)

    if (elst_repl) then
       ss = s_rep/vol_box
       open( 23,file=trim(datadir)//'bulk_stress_rep.txt',position='append')
       write(23,'(I10,9ES12.3)') istep, &
            ss(1,1),ss(1,2),ss(1,3), &
            ss(2,1),ss(2,2),ss(2,3), &
            ss(3,1),ss(3,2),ss(3,3)
       close(23)
    endif

    return
  end subroutine output_stress
  
  !-------------------------------------------------------
  
  !subroutine output_stress(abc, ss)
  !
  !  character(len=3), intent(in) :: abc
  !  
  !  real, dimension(3,3), intent(in) :: ss
  !
  !  open( 23,file=trim(datadir)//'bulk_stress_'//abc//'.txt',position='append')
  !  write(23,'(I10,9ES12.3)') istep, &
  !       ss(1,1),ss(1,2),ss(1,3), &
  !       ss(2,1),ss(2,2),ss(2,3), &
  !       ss(3,1),ss(3,2),ss(3,3)
  !  close(23)
  !
  !  return
  !end subroutine output_stress

  !-------------------------------------------------------

  subroutine output_vtk(iter)

    ! An unelegant vtk writer.
    
    integer, intent(in) :: iter
    
    integer :: i,j, fh
    character(len=2) :: digit_count
    character(len=7) :: iterchar
    character(len=40) :: filename

    
    fh = 11 ! file handle
    
    !write(iterchar,'(i7.7)') iter/1000
    !filename = trim(datadir)//'para'//trim(dsub_num)//'/para'//iterchar//'k.vtk'
    write(iterchar,'(i7.7)') iter
    filename = trim(datadir)//'para'//trim(dsub_num)//'/para'//iterchar//'.vtk'
    
    ! get the digit length of np
    write(digit_count,'(I1)') floor(log10(1.*np))+1 
    
    open(fh, file=filename, status='replace')

    write(fh,'(A)') '# vtk DataFile Version 3.0'
    write(fh,'(A)') 'a stupid file'
    write(fh,'(A)') 'ASCII'
    write(fh,'(A)') 'DATASET UNSTRUCTURED_GRID' 

    !-- particle positions --!

    write(fh,'(A,I'//trim(digit_count)//',A)') 'POINTS ', np, ' float'
    
    j=1
    do i = 1,np
       write(fh,'(3F12.6)') pp(j),pp(j+1),pp(j+2) 
       j=j+3
    enddo

    !-- vtk cell type (1) for points --!
    
    write(fh,*)
    write(fh,'(A,I'//trim(digit_count)//')') 'CELL_TYPES ', np 

    do i = 1,np
       write(fh,'(I1)') 1
    enddo

    !-- some scalar value for each point --!
    
    write(fh,*)
    write(fh,'(A,I'//trim(digit_count)//')') 'POINT_DATA ', np
    write(fh,'(A,I1)') 'SCALARS index float ', 1
    write(fh,'(A)') 'LOOKUP_TABLE default'
    
    do i = 1,np
       write(fh,'(F5.1)') radii(i) 
    enddo

    close(fh)
    
    return
  end subroutine output_vtk

  !-------------------------------------------------------

  subroutine output_timing

    real :: twt,tun ! [total wall time] and [t unit]
    character(len=5) :: cun  ! time unit
    
    twt = t_w1 - t_w0

    if (twt < 60.) then
       tun = 1.
       cun = ' sec.'
    elseif (twt < 3600.) then
       tun = 60.
       cun = ' min.'
    else
       tun = 3600.
       cun = ' hr.'
    endif

    open( 17,file=trim(datadir)//'log.txt',position='append')
    write(17,'(A)') ''
    write(17,'(A,ES9.3,A)') 'Total computing time = ',twt/tun,cun
    write(17,'(A,F5.1,A,F5.1,A)') '(Neighbor',t_nb_tot/twt*100., &
         '%, F&T compute', t_comp_tot/twt*100.,'%.)'
    close(17)

    return
  end subroutine output_timing

  !-------------------------------------------------------

  subroutine output_restart(iter)

    ! write variables necessary for restarting a simulation

    integer, intent(in) :: iter
    
    integer :: i,j, fh
    character(len=7) :: iterchar
    character(len=25) :: filename

    
    fh = 11 ! file handle
    
    write(iterchar,'(i7.7)') iter/1000
    filename = trim(datadir)//'restart'//iterchar//'k.txt'
    
    open(fh, file=filename, status='replace')

    write(fh,'(A,I10)') 'This is a restarting file from istep = ', istep

    write(fh,'(1ES16.8)') t9
    write(fh,'(4ES16.8)') vol_frac, lx,ly,lz
    write(fh,'(2ES16.8)') y_shift, v_shift    
    
    do i = 1,np
       write(fh,'(I6,ES16.8,I3)') i, radii(i), p_id(i)
    enddo

    j=1
    do i = 1,np
       write(fh,'(I6,3ES16.8)') i, pp(j),pp(j+1),pp(j+2) 
       j=j+3
    enddo

    j=1
    do i = 1,np
       write(fh,'(I6,3ES16.8)') i, pv(j),pv(j+1),pv(j+2) 
       j=j+3
    enddo

    j=1
    do i = 1,np
       write(fh,'(I6,3ES16.8)') i, po(j),po(j+1),po(j+2) 
       j=j+3
    enddo

    j=1
    do i = 1,np
       write(fh,'(I6,3ES16.8)') i, pw(j),pw(j+1),pw(j+2) 
       j=j+3
    enddo

    do i = 1,np
       do j = 1,np
          write(fh,'(2I6,I3,3ES16.8)') i,j, fr_prev_mx(j,i), xi_prev_mx(1,j,i),xi_prev_mx(2,j,i),xi_prev_mx(3,j,i)   
       enddo
    enddo

    write(fh,'(A)') 'End.'
    close(fh)
    

    return
  end subroutine output_restart




  !-----------------------------------------------------

  Function sqrt_3d(vec) result(mag_vec)
    
    real, dimension(3), intent(in) :: vec
    real :: mag_vec
    
    mag_vec = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

  End Function sqrt_3d

  

end module mod_output
