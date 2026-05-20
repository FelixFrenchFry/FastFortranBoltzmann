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
    integer(int32), parameter :: N_X = 14400
    integer(int32), parameter :: N_Y = 14400
    integer(int32), parameter :: N_STEPS = 300
    integer(int64), parameter :: N_CELLS = int(N_X, int64) * int(N_Y, int64)
    integer(int32), parameter :: N_DIRS = 9

    ! constants for sim modes
    integer(int32), parameter :: SIM_SHEAR_WAVE = 1
    integer(int32), parameter :: SIM_COUETTE_FLOW = 2
    integer(int32), parameter :: SIM_POISEUILLE_FLOW = 3
    integer(int32), parameter :: SIM_SLIDING_LID = 4
    integer(int32), parameter :: SIM_MODE = 4 ! selected sim mode

    ! kernel selection
    logical, parameter :: USE_UNROLLED_KERNELS = .true.
    logical, parameter :: USE_UNIVERSAL_KERNELS = .false.
    logical, parameter :: USE_PULL_SHIFT_KERNELS = .false.

    ! export settings
    logical, parameter :: EXPORT_RHO = .true.
    logical, parameter :: EXPORT_U_X = .true.
    logical, parameter :: EXPORT_U_Y = .true.
    logical, parameter :: EXPORT_U_MAG = .true.
    integer(int32), parameter :: EXPORT_INTERVAL = 100000
    logical, parameter :: EXPORT_INITIAL_STATE = .true.
    logical, parameter :: EXPORT_FINAL_STATE = .true.
    character(len=*), parameter :: OUTPUT_DIR_NAME = "output"
    character(len=*), parameter :: EXPORT_NUM = "run_000"

    ! progress display settings
    logical, parameter :: INTERACTIVE_PROGRESS = .true.
    integer(int32), parameter :: PROGRESS_INTERVAL = 1

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

    ! sim parameter sets for each sim mode
    type :: shear_wave_params_t
        real(FP) :: rho_0 ! rest density
        real(FP) :: omega ! relaxation factor
        real(FP) :: u_max ! initial velocity
        real(FP) :: n_sin ! num sin periods
    end type shear_wave_params_t

    type :: couette_flow_params_t
        real(FP) :: rho_0 ! rest density
        real(FP) :: omega ! relaxation factor
        real(FP) :: u_wall ! top wall velocity
    end type couette_flow_params_t

    type :: poiseuille_flow_params_t
        real(FP) :: rho_0 ! rest density
        real(FP) :: omega ! relaxation factor
        real(FP) :: rho_in ! inlet density
        real(FP) :: rho_out ! outlet density
    end type poiseuille_flow_params_t

    type :: sliding_lid_params_t
        real(FP) :: rho_0 ! rest density
        real(FP) :: omega ! relaxation factor
        real(FP) :: u_wall ! top wall velocity
    end type sliding_lid_params_t

    ! parameter set for shear wave
    type(shear_wave_params_t), parameter :: SW_PARAMS = shear_wave_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        u_max = 0.1_FP, &
        n_sin = 2.0_FP &
    )

    ! parameter set for couette flow
    type(couette_flow_params_t), parameter :: CF_PARAMS = couette_flow_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        u_wall = 0.1_FP &
    )

    ! parameter set for poiseuille flow
    type(poiseuille_flow_params_t), parameter :: PF_PARAMS = poiseuille_flow_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        rho_in = 1.001_FP, &
        rho_out = 0.999_FP &
    )

    ! parameter set for sliding lid
    type(sliding_lid_params_t), parameter :: SL_PARAMS = sliding_lid_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        u_wall = 0.1_FP &
    )

contains

    pure function sim_mode_to_string( &
        sim_mode &
        ) result(name)
        ! read-only inputs
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
