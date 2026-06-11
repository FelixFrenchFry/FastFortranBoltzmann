module settings
    ! imports
    use iso_fortran_env, only: int32, int64, real32, real64
    implicit none

    ! floating point precision
#ifdef FFB_FP64
    integer(int32), parameter :: FP = real64
    character(len=*), parameter :: FP_DTYPE = "real64"
#else
    integer(int32), parameter :: FP = real32
    character(len=*), parameter :: FP_DTYPE = "real32"
#endif

    ! sim size and duration
#ifdef FFB_USE_CMAKE_SETTINGS
    integer(int32), parameter :: N_X = FFB_N_X
    integer(int32), parameter :: N_Y = FFB_N_Y
    integer(int32), parameter :: N_STEPS = FFB_N_STEPS
#else
    integer(int32), parameter :: N_X = 14400
    integer(int32), parameter :: N_Y = 14400
    integer(int32), parameter :: N_STEPS = 300
#endif
    integer(int64), parameter :: N_CELLS = int(N_X, int64) * int(N_Y, int64)
    integer(int32), parameter :: N_DIRS = 9

    ! sim parameters
#ifdef FFB_USE_CMAKE_SETTINGS
    real(FP), parameter :: RHO_0 = FFB_RHO_0
    real(FP), parameter :: OMEGA = FFB_OMEGA
    real(FP), parameter :: U_MAX = FFB_U_MAX
    real(FP), parameter :: N_SIN = FFB_N_SIN
    real(FP), parameter :: U_WALL = FFB_U_WALL
    real(FP), parameter :: U_LID = FFB_U_LID
    real(FP), parameter :: RHO_IN = FFB_RHO_IN
    real(FP), parameter :: RHO_OUT = FFB_RHO_OUT
#else
    real(FP), parameter :: RHO_0 = 1.0_FP
    real(FP), parameter :: OMEGA = 1.5_FP
    real(FP), parameter :: U_MAX = 0.1_FP
    real(FP), parameter :: N_SIN = 2.0_FP
    real(FP), parameter :: U_WALL = 0.1_FP
    real(FP), parameter :: U_LID = 0.1_FP
    real(FP), parameter :: RHO_IN = 1.001_FP
    real(FP), parameter :: RHO_OUT = 0.999_FP
#endif

    ! constants for sim modes
    integer(int32), parameter :: SIM_SHEAR_WAVE = 1
    integer(int32), parameter :: SIM_COUETTE_FLOW = 2
    integer(int32), parameter :: SIM_POISEUILLE_FLOW = 3
    integer(int32), parameter :: SIM_SLIDING_LID = 4
#ifdef FFB_USE_CMAKE_SETTINGS
    integer(int32), parameter :: SIM_MODE = FFB_SIM_MODE ! selected sim mode
#else
    integer(int32), parameter :: SIM_MODE = 4 ! selected sim mode
#endif

    ! algorithm selection
#ifdef FFB_USE_CMAKE_SETTINGS
    logical, parameter :: USE_UNROLLED_KERNELS = FFB_USE_UNROLLED_KERNELS
#else
    logical, parameter :: USE_UNROLLED_KERNELS = .true.
#endif

    ! export settings
#ifdef FFB_USE_CMAKE_SETTINGS
    logical, parameter :: EXPORT_MACROS = FFB_EXPORT_MACROS
    logical, parameter :: EXPORT_ENDPOINT_STATES = FFB_EXPORT_ENDPOINT_STATES
    integer(int32), parameter :: EXPORT_INTERVAL = FFB_EXPORT_INTERVAL
    character(len=*), parameter :: EXPORT_NUM = FFB_EXPORT_NUM
#else
    logical, parameter :: EXPORT_MACROS = .false.
    logical, parameter :: EXPORT_ENDPOINT_STATES = .true.
    integer(int32), parameter :: EXPORT_INTERVAL = 10000
    character(len=*), parameter :: EXPORT_NUM = "run_000"
#endif

    ! progress display settings
#ifdef FFB_USE_CMAKE_SETTINGS
    logical, parameter :: INTERACTIVE_PROGRESS = FFB_INTERACTIVE_PROGRESS
    integer(int32), parameter :: PROGRESS_INTERVAL = FFB_PROGRESS_INTERVAL
#else
    logical, parameter :: INTERACTIVE_PROGRESS = .true.
    integer(int32), parameter :: PROGRESS_INTERVAL = 1
#endif

    ! D2Q9 lattice velocities and weights
    integer(int32), parameter :: C_X(N_DIRS) = [ 0,  1,  0, -1,  0,  1, -1, -1,  1 ]
    integer(int32), parameter :: C_Y(N_DIRS) = [ 0,  0,  1,  0, -1,  1,  1, -1, -1 ]
    real(FP), parameter :: C_X_FP(N_DIRS) = real(C_X, FP) ! fp-version for compute
    real(FP), parameter :: C_Y_FP(N_DIRS) = real(C_Y, FP) ! fp-version for compute
    ! ---------
    ! | 7 3 6 |
    ! | 4 1 2 |
    ! | 8 5 9 |
    ! ---------
    ! 1: ( 0,  0) = rest
    ! 2: ( 1,  0) = east
    ! 3: ( 0,  1) = north
    ! 4: (-1,  0) = west
    ! 5: ( 0, -1) = south
    ! 6: ( 1,  1) = north-east
    ! 7: (-1,  1) = north-west
    ! 8: (-1, -1) = south-west
    ! 9: ( 1, -1) = south-east
    real(FP), parameter :: W(N_DIRS) = [ &
        4.0_FP/9.0_FP, &
        1.0_FP/9.0_FP, &
        1.0_FP/9.0_FP, &
        1.0_FP/9.0_FP, &
        1.0_FP/9.0_FP, &
        1.0_FP/36.0_FP, &
        1.0_FP/36.0_FP, &
        1.0_FP/36.0_FP, &
        1.0_FP/36.0_FP]

    ! misc
    real(FP), parameter :: PI = 3.141592653589793238462643383279502884197_FP

contains

    pure function sim_mode_to_string( &
        sim_mode &
        ) result(name)
        ! inputs
        integer(int32), intent(in) :: sim_mode

        ! output
        character(len=:), allocatable :: name

        select case (sim_mode)
        case (SIM_SHEAR_WAVE)
            name = "shear_wave"
        case (SIM_COUETTE_FLOW)
            name = "couette_flow"
        case (SIM_POISEUILLE_FLOW)
            name = "poiseuille_flow"
        case (SIM_SLIDING_LID)
            name = "sliding_lid"
        case default
            name = "unknown_sim_mode"
        end select
    end function sim_mode_to_string

end module settings
