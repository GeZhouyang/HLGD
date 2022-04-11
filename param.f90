module mod_param

  implicit none

  !-- Constants --!
  
  real, parameter :: pi = acos(-1.)

  !-- Miscellaneous --!
    
  character(len=2),  parameter :: case_num = "1/"
  character(len=6),  parameter :: datadir = 'data'//trim(case_num)
  character(len=50), parameter :: init_fn = 'D3N500VF0.5Bidi1.4_0.5Cubic_1_.dat'
  logical,           parameter :: restart_from_last = .true.
  logical,           parameter :: start_oscl_from_stdy = .false.
  !character(len=50), parameter :: restart_fname = 'restart_500/restart0025000k.txt'
  character(len=50), parameter :: restart_fname = 'data1/restart1500000k.txt'
  character(len=1),  parameter :: dsub_num = "z"  ! index of sub folders in /data?, eg. /gap?
  
  !-- Physical paramters --!
  
  integer, parameter ::   np = 500      ! number of particles
  
  real, parameter ::  a_1 = 1.          ! radius of group 1 sphere [length] (the smaller one)
  real, parameter ::  a_2 = 1.4         ! radius of group 2 sphere [length]
  real, parameter ::  rho = 1.          ! particle density [mass/length^3]
  real, parameter ::  k_n = 1e4         ! normal spring constant (collision) [mass/time^2]
  real, parameter :: dp_n = k_n*0.0     ! normal damping constant (collision) [mass/time]
  real, parameter ::  k_t = k_n*2./7.   ! tangential spring constant (collision) [mass/time^2]
  real, parameter :: mu_c = 0.5         ! sliding friction coefficient [1]
                                        
  real, parameter :: mu_f = 1.          ! fluid dynamic viscosity [mass/(time*length)]
  real, parameter :: shear_rate = 1e-2  ! [1/time] (ineffective if oscillatory_shear)

  logical, parameter :: oscillatory_shear = .true.
  real,    parameter :: amp0 = 5e-2      ! oscillation amplitude (at z=Lz) normalized by Lz  [1]
  real,    parameter :: freq = 1e-1      ! oscillation frequency [1/time]

  logical, parameter :: elst_repl = .true.     ! electrostatic repulsion >>> temporarily disabled
  real,    parameter :: elst_sr   = 1e2        ! electrostatic repulsion scaled shear rate [1]
  real,    parameter :: Debye_len = 0.05*a_1   ! Debye length (appearing in the electrostatic repulsion model) [length]
  logical, parameter :: van_der_W = .true.     ! van der Waals (VdW) attraction
  real,    parameter :: C_Hamaker = 7e-8       ! Hamaker constant [energy]

  real, parameter :: fric_sr = 1e-2                 ! friction-scaled shear rate [1], redundant if fully_frictional

  logical, parameter :: frictional = .true.         ! frictional contact
  logical, parameter :: fully_frictional = .true.  ! always activate friction (i.e. zero critical loading)
  
  !-- Numerical parameters --!
  
  real, parameter :: h_lub_out = min(a_1,a_2)*0.2    ! lubrication cutoff (outer range) !!!!!!!!!!!!!!############
  real, parameter :: h_lub_inn = min(a_1,a_2)*0.001   ! lubrication cutoff (inner range)

  integer, parameter :: integrate_method = 2
  ! 1: Forward Euler
  ! 2: Velocity-Verlet
  
end module mod_param
