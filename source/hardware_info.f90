module hardware_info
    ! imports
    use iso_fortran_env, only: int64
    implicit none

#ifndef FFB_COMPILER_ID
#define FFB_COMPILER_ID "unknown"
#endif
#ifndef FFB_COMPILER_VERSION
#define FFB_COMPILER_VERSION "unknown"
#endif
#ifndef FFB_FORTRAN_FLAGS
#define FFB_FORTRAN_FLAGS "unknown"
#endif

    private

    public :: hardware_info_t
    public :: collect_hardware_info
    public :: print_hardware_summary
    public :: write_hardware_metadata

    type :: hardware_info_t

        ! machine info
        character(len=:), allocatable :: cpu_model
        character(len=:), allocatable :: logical_threads
        character(len=:), allocatable :: slurm_nodes

        ! build info
        character(len=:), allocatable :: compiler
        character(len=:), allocatable :: compiler_flags

    end type hardware_info_t

    character(len=*), parameter :: unknown_value = "unknown"
    character(len=*), parameter :: compiler_id = FFB_COMPILER_ID
    character(len=*), parameter :: compiler_version = FFB_COMPILER_VERSION
    character(len=*), parameter :: fortran_flags = FFB_FORTRAN_FLAGS

contains

    subroutine collect_hardware_info( &
        info &
        )
        ! output
        type(hardware_info_t), intent(out) :: info

        ! locals
        character(len=256) :: line
        integer :: exitstat

        call set_unknown_hardware_info(info)

        call read_command_output( &
            "sh -c ""sed -n 's/^model name[[:space:]]*:[[:space:]]*//p' /proc/cpuinfo | head -n 1""", &
            info%cpu_model)
        call read_command_output( &
            "sh -c ""grep -c '^processor' /proc/cpuinfo""", &
            info%logical_threads)

        call get_environment_variable("SLURM_NNODES", line, status=exitstat)
        if (exitstat == 0 .and. len_trim(line) > 0) then
            info%slurm_nodes = trim(adjustl(line))
        else
            info%slurm_nodes = "1"
        end if

        info%compiler = trim(compiler_id) // " " // trim(compiler_version)

        if (len_trim(fortran_flags) > 0) then
            info%compiler_flags = trim(fortran_flags)
        else
            info%compiler_flags = "none"
        end if
    end subroutine collect_hardware_info


    subroutine print_hardware_summary( &
        info &
        )
        ! inputs
        type(hardware_info_t), intent(in) :: info

        print '(A)', "--- [ hardware ] ---------------------------------------------------------------"
        print '(A,T27,A,A)', "cpu model", "= ", trim(info%cpu_model)
        print '(A,T27,A,A)', "logical threads/node", "= ", trim(info%logical_threads)
        print '(A,T27,A,A)', "compute nodes", "= ", trim(info%slurm_nodes)
        print '(A,T27,A,A)', "compiler", "= ", trim(info%compiler)
        print '(A,T27,A,A)', "flags", "= ", trim(info%compiler_flags)
    end subroutine print_hardware_summary


    subroutine write_hardware_metadata( &
        unit, info &
        )
        ! inputs
        integer, intent(in) :: unit
        type(hardware_info_t), intent(in) :: info

        write(unit, '(A)') '  "hardware": {'
        write(unit, '(A,A,A)') '    "cpu_model": "', json_escape(info%cpu_model), '",'
        write(unit, '(A,A,A)') '    "logical_threads_per_node": ', trim(integer_text_to_json(info%logical_threads)), ','
        write(unit, '(A,A,A)') '    "compute_nodes": ', trim(integer_text_to_json(info%slurm_nodes)), ','
        write(unit, '(A,A,A)') '    "compiler": "', json_escape(info%compiler), '",'
        write(unit, '(A,A,A)') '    "compiler_flags": "', json_escape(info%compiler_flags), '"'
        write(unit, '(A)') '  },'
        write(unit, '(A)') ""
    end subroutine write_hardware_metadata


    subroutine set_unknown_hardware_info( &
        info &
        )
        ! output
        type(hardware_info_t), intent(out) :: info

        info%cpu_model = unknown_value
        info%logical_threads = unknown_value
        info%slurm_nodes = unknown_value
        info%compiler = unknown_value
        info%compiler_flags = unknown_value
    end subroutine set_unknown_hardware_info


    subroutine read_command_output( &
        command, command_output &
        )
        ! inputs
        character(len=*), intent(in) :: command

        ! output
        character(len=:), allocatable, intent(out) :: command_output

        ! locals
        character(len=*), parameter :: tmp_file = ".ffb_hardware_info.tmp"
        character(len=512) :: line
        integer :: unit
        integer :: io_stat
        integer :: cmdstat
        integer :: exitstat

        command_output = unknown_value

        call execute_command_line(trim(command) // " > " // tmp_file // " 2>/dev/null", &
            exitstat=exitstat, cmdstat=cmdstat)

        if (cmdstat == 0) then
            open(newunit=unit, file=tmp_file, form="formatted", status="old", action="read", iostat=io_stat)

            if (io_stat == 0) then
                read(unit, '(A)', iostat=io_stat) line
                if (io_stat == 0 .and. len_trim(line) > 0) then
                    command_output = trim(adjustl(line))
                end if
                close(unit)
            end if
        end if

        call execute_command_line("rm -f " // tmp_file)
    end subroutine read_command_output


    function integer_text_to_json( &
        text &
        ) result(value_json)
        ! inputs
        character(len=*), intent(in) :: text

        ! output
        character(len=32) :: value_json

        ! locals
        integer(int64) :: value
        integer :: io_stat

        read(text, *, iostat=io_stat) value

        if (io_stat == 0) then
            write(value_json, '(I0)') value
        else
            value_json = "null"
        end if

        value_json = adjustl(value_json)
    end function integer_text_to_json


    function json_escape( &
        text &
        ) result(escaped)
        ! inputs
        character(len=*), intent(in) :: text

        ! output
        character(len=:), allocatable :: escaped

        ! locals
        integer :: i

        escaped = ""

        do i = 1, len_trim(text)
            select case (text(i:i))
            case ('"')
                escaped = escaped // '\"'
            case ('\')
                escaped = escaped // '\\'
            case default
                escaped = escaped // text(i:i)
            end select
        end do
    end function json_escape

end module hardware_info
