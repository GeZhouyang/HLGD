program SLD

  use OMP_LIB

  use mod_param
  use mod_common
  use mod_zero
  use mod_core
  use mod_output

  implicit none

  integer :: istart,iend
  real :: t_simu, t_nondi
  
  istart = 1500*1000*1000
  t_nondi = 50.  ! uneffective now, search for iend in zero.f90
  
  if (oscillatory_shear) then
     t_simu = t_nondi*(2.*pi/freq)
  else
     t_simu = t_nondi/shear_rate
  endif
  
  call cpu_time(t_w0)
  t_w0 = omp_get_wtime() 

  call launcher(istart,t_simu,iend)  ! initialization
  call time_marching(istart,iend)    ! the main loop

  call cpu_time(t_w1)
  t_w1 = omp_get_wtime()

  call output_timing
  
  stop
end program SLD
