module settings
    ! imports
    use iso_fortran_env, only: int32, int64, real32, real64
    implicit none

    ! sim size and duration
    integer(int32), parameter :: N_X = 500
    integer(int32), parameter :: N_Y = 500
    integer(int32), parameter :: N_STEPS = 10000
    integer(int64), parameter :: N_CELLS = int(N_X, int64) * int(N_Y, int64)
    integer(int32), parameter :: N_DIRS = 9

    ! constants for sim modes
    integer(int32), parameter :: SIM_SHEAR_WAVE = 1
    integer(int32), parameter :: SIM_COUETTE_FLOW = 2
    integer(int32), parameter :: SIM_POISEUILLE_FLOW = 3
    integer(int32), parameter :: SIM_SLIDING_LID = 4
    integer(int32), parameter :: SIM_MODE = 1 ! selected sim mode

    ! floating point precision
#ifdef FFB_FP64
    integer(int32), parameter :: FP = real64
    character(len=*), parameter :: FP_DTYPE = "real64"
#else
    integer(int32), parameter :: FP = real32
    character(len=*), parameter :: FP_DTYPE = "real32"
#endif

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
