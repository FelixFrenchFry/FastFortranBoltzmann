module export
    ! imports
    use iso_fortran_env, only: int32, real32
    implicit none

    private

    ! export modes
    integer(int32), parameter, public :: EXPORT_NONE = 0
    integer(int32), parameter, public :: EXPORT_DENSITY = 1
    integer(int32), parameter, public :: EXPORT_VELOCITY_X = 2
    integer(int32), parameter, public :: EXPORT_VELOCITY_Y = 3
    integer(int32), parameter, public :: EXPORT_VELOCITY_MAG = 4

    public :: export_mode_to_string
    public :: should_export_step
    public :: export_selected_data

contains

    pure function export_mode_to_string( &
        export_mode &
        ) result(name)
        ! read-only inputs
        integer(int32), intent(in) :: export_mode

        ! output
        character(len=:), allocatable :: name

        select case (export_mode)
        case (EXPORT_NONE)
            name = "none"
        case (EXPORT_DENSITY)
            name = "density"
        case (EXPORT_VELOCITY_X)
            name = "velocity_x"
        case (EXPORT_VELOCITY_Y)
            name = "velocity_y"
        case (EXPORT_VELOCITY_MAG)
            name = "velocity_mag"
        case default
            name = "unknown"
        end select
    end function export_mode_to_string


    pure function should_export_step( &
        N_STEPS, step, export_mode, export_interval, export_initial_state, export_final_state &
        ) result(write_step)
        ! read-only inputs
        integer(int32), intent(in) :: N_STEPS
        integer(int32), intent(in) :: step
        integer(int32), intent(in) :: export_mode
        integer(int32), intent(in) :: export_interval
        logical, intent(in) :: export_initial_state
        logical, intent(in) :: export_final_state

        ! output
        logical :: write_step

        ! no output if export is disabled
        if (export_mode == EXPORT_NONE .or. export_interval <= 0) then
            write_step = .false.

        ! optional output at initial state
        else if (step == 0) then
            write_step = export_initial_state

        ! regular interval output plus optional output at final state
        else
            write_step = mod(step, export_interval) == 0 .or. &
                (export_final_state .and. step == N_STEPS)
        end if
    end function should_export_step


    subroutine export_selected_data( &
        export_mode, output_dir_name, export_num, suffix_num, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: export_mode
        character(len=*), intent(in) :: output_dir_name
        character(len=*), intent(in) :: export_num
        integer(int32), intent(in) :: suffix_num
        real(real32), intent(in) :: rho(:,:)
        real(real32), intent(in) :: u_x(:,:)
        real(real32), intent(in) :: u_y(:,:)

        ! temp
        real(real32), allocatable :: velocity_mag(:,:)

        ! choose and export selected scalar field
        select case (export_mode)
        case (EXPORT_NONE)
            return
        case (EXPORT_DENSITY)
            call export_scalar_field(rho, "density", output_dir_name, export_num, suffix_num)
        case (EXPORT_VELOCITY_X)
            call export_scalar_field(u_x, "velocity_x", output_dir_name, export_num, suffix_num)
        case (EXPORT_VELOCITY_Y)
            call export_scalar_field(u_y, "velocity_y", output_dir_name, export_num, suffix_num)
        case (EXPORT_VELOCITY_MAG)
            allocate(velocity_mag(size(u_x, 1), size(u_x, 2)))
            velocity_mag = sqrt(u_x * u_x + u_y * u_y) ! element-wise sqrt of velocity magnitude
            call export_scalar_field(velocity_mag, "velocity_mag", output_dir_name, export_num, suffix_num)
        case default
            error stop "unknown export mode in export_selected_data()"
        end select
    end subroutine export_selected_data


    subroutine export_scalar_field( &
        field, field_name, output_dir_name, export_num, suffix_num &
        )
        ! read-only inputs
        real(real32), intent(in) :: field(:,:)
        character(len=*), intent(in) :: field_name
        character(len=*), intent(in) :: output_dir_name
        character(len=*), intent(in) :: export_num
        integer(int32), intent(in) :: suffix_num

        ! temp
        character(len=:), allocatable :: output_path
        character(len=:), allocatable :: file_path

        ! assemble output path and filename with its data type and step suffix
        output_path = trim(output_dir_name) // "/" // trim(export_num)
        file_path = output_path // "/" // trim(field_name) // format_step_suffix(suffix_num) // ".bin"

        call ensure_output_directory(output_path)
        call write_binary_field(field, file_path)
    end subroutine export_scalar_field


    pure function format_step_suffix( &
        suffix_num &
        ) result(suffix)
        ! read-only inputs
        integer(int32), intent(in) :: suffix_num

        ! output
        character(len=10) :: suffix

        ! temp
        character(len=9) :: suffix_digits

        write(suffix_digits, '(I9.9)') suffix_num
        suffix = "_" // suffix_digits
    end function format_step_suffix


    subroutine ensure_output_directory( &
        output_path &
        )
        ! read-only inputs
        character(len=*), intent(in) :: output_path

        ! temp
        character(len=:), allocatable :: command
        integer :: cmdstat
        integer :: exitstat

        ! create output directory if needed
        command = 'mkdir -p "' // trim(output_path) // '"'
        call execute_command_line(command, exitstat=exitstat, cmdstat=cmdstat)

        if (cmdstat /= 0 .or. exitstat /= 0) then
            error stop "failed to create output directory"
        end if
    end subroutine ensure_output_directory


    subroutine write_binary_field( &
        field, file_path &
        )
        ! read-only inputs
        real(real32), intent(in) :: field(:,:)
        character(len=*), intent(in) :: file_path

        ! temp
        integer :: unit
        integer :: io_stat

        ! open raw binary stream file
        open(newunit=unit, file=trim(file_path), access="stream", form="unformatted", &
            status="replace", action="write", iostat=io_stat)

        if (io_stat /= 0) then
            error stop "could not open binary output file"
        end if

        ! write raw real32 field data to file
        write(unit, iostat=io_stat) field

        if (io_stat /= 0) then
            close(unit)
            error stop "could not write binary output file"
        end if

        close(unit)
    end subroutine write_binary_field

end module export
