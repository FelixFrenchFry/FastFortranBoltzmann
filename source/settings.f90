module settings
    ! imports
    use iso_fortran_env, only: int32, int64, real32
    implicit none

    ! sim size and duration
    integer(int32), parameter :: N_X = 500
    integer(int32), parameter :: N_Y = 500
    integer(int32), parameter :: N_STEPS = 120000
    integer(int64), parameter :: N_CELLS = int(N_X, int64) * int(N_Y, int64)
    integer(int32), parameter :: N_DIRS = 9

    ! constants for sim modes
    integer(int32), parameter :: SIM_SHEAR_WAVE = 1
    integer(int32), parameter :: SIM_COUETTE_FLOW = 2
    integer(int32), parameter :: SIM_POISEUILLE_FLOW = 3
    integer(int32), parameter :: SIM_SLIDING_LID = 4
    integer(int32), parameter :: sim_mode = 2

    ! misc
    real(real32), parameter :: PI = 3.1415927410125732421875_real32

    ! sim parameter sets for each sim mode
    type :: shear_wave_params_t
        real(real32) :: rho_0 ! rest density
        real(real32) :: omega ! relaxation factor
        real(real32) :: u_max ! initial velocity
        real(real32) :: n_sin ! num sin periods
    end type shear_wave_params_t

    type :: couette_flow_params_t
        real(real32) :: rho_0 ! rest density
        real(real32) :: omega ! relaxation factor
        real(real32) :: u_wall ! top wall velocity
    end type couette_flow_params_t

    type :: poiseuille_flow_params_t
        real(real32) :: rho_0 ! rest density
        real(real32) :: omega ! relaxation factor
        real(real32) :: rho_in ! inlet density
        real(real32) :: rho_out ! outlet density
    end type poiseuille_flow_params_t

    type :: sliding_lid_params_t
        real(real32) :: rho_0 ! rest density
        real(real32) :: omega ! relaxation factor
        real(real32) :: u_wall ! top wall velocity
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
