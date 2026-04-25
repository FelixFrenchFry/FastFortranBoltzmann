module settings
    ! imports
    use iso_fortran_env, only: int32
    implicit none

    integer(int32), parameter :: SIM_SHEAR_WAVE = 1
    integer(int32), parameter :: SIM_COUETTE_FLOW = 2
    integer(int32), parameter :: SIM_POISEUILLE_FLOW = 3
    integer(int32), parameter :: SIM_SLIDING_LID = 4
    integer(int32), parameter :: sim_mode = 1

contains

    pure function sim_mode_to_string() result(name)
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
