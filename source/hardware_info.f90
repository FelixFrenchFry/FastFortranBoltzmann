module hardware_info
    ! imports
    use iso_fortran_env, only: int64, real64
    implicit none

#ifndef FFB_COMPILER_ID
#define FFB_COMPILER_ID "unknown"
#endif
#ifndef FFB_COMPILER_VERSION
#define FFB_COMPILER_VERSION "unknown"
#endif
#ifndef FFB_BUILD_PRESET
#define FFB_BUILD_PRESET "unknown"
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
        character(len=:), allocatable :: hostname
        character(len=:), allocatable :: cpu_model
        character(len=:), allocatable :: logical_threads
        character(len=:), allocatable :: memory_total

        ! build info
        character(len=:), allocatable :: compiler
        character(len=:), allocatable :: build_preset
        character(len=:), allocatable :: compiler_flags

    end type hardware_info_t

    character(len=*), parameter :: unknown_value = "unknown"
    character(len=*), parameter :: compiler_id = FFB_COMPILER_ID
    character(len=*), parameter :: compiler_version = FFB_COMPILER_VERSION
    character(len=*), parameter :: cmake_build_preset = FFB_BUILD_PRESET
    character(len=*), parameter :: fortran_flags = FFB_FORTRAN_FLAGS

contains

    subroutine collect_hardware_info( &
        info &
        )
        ! output
        type(hardware_info_t), intent(out) :: info

        ! locals
        character(len=:), allocatable :: mem_total_kb_text
        integer(int64) :: mem_total_kb
        integer :: io_stat

        call set_unknown_hardware_info(info)

        call read_command_output("hostname", info%hostname)
        call read_command_output( &
            "sh -c ""sed -n 's/^model name[[:space:]]*:[[:space:]]*//p' /proc/cpuinfo | head -n 1""", &
            info%cpu_model)
        call read_command_output( &
            "sh -c ""grep -c '^processor' /proc/cpuinfo""", &
            info%logical_threads)
        call read_command_output( &
            "sh -c ""grep '^MemTotal:' /proc/meminfo | tr -s ' ' | cut -d' ' -f2""", &
            mem_total_kb_text)

        if (len_trim(mem_total_kb_text) > 0 .and. trim(mem_total_kb_text) /= unknown_value) then
            read(mem_total_kb_text, *, iostat=io_stat) mem_total_kb
            if (io_stat == 0) then
                info%memory_total = format_memory_gb(mem_total_kb)
            end if
        end if

        info%compiler = trim(compiler_id) // " " // trim(compiler_version)
        info%build_preset = trim(cmake_build_preset)

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

        print '(A)', "--- [ hardware ] ----------------------------------------------------------"
        print '(A,T27,A,A)', "hostname", "= ", trim(info%hostname)
        print '(A,T27,A,A)', "cpu model", "= ", trim(info%cpu_model)
        print '(A,T27,A,A)', "logical threads", "= ", trim(info%logical_threads)
        print '(A,T27,A,A)', "memory total", "= ", trim(info%memory_total)
        print '(A,T27,A,A)', "compiler", "= ", trim(info%compiler)
        print '(A,T27,A,A)', "build preset", "= ", trim(info%build_preset)
        print '(A,T27,A,A)', "compiler flags", "= ", trim(info%compiler_flags)
    end subroutine print_hardware_summary


    subroutine write_hardware_metadata( &
        unit, info &
        )
        ! inputs
        integer, intent(in) :: unit
        type(hardware_info_t), intent(in) :: info

        write(unit, '(A)') '  "hardware": {'
        write(unit, '(A,A,A)') '    "hostname": "', json_escape(info%hostname), '",'
        write(unit, '(A,A,A)') '    "cpu_model": "', json_escape(info%cpu_model), '",'
        write(unit, '(A,A,A)') '    "logical_threads": ', trim(integer_text_to_json(info%logical_threads)), ','
        write(unit, '(A,A,A)') '    "memory_total": "', json_escape(info%memory_total), '",'
        write(unit, '(A,A,A)') '    "compiler": "', json_escape(info%compiler), '",'
        write(unit, '(A,A,A)') '    "build_preset": "', json_escape(info%build_preset), '",'
        write(unit, '(A,A,A)') '    "compiler_flags": "', json_escape(info%compiler_flags), '"'
        write(unit, '(A)') '  },'
        write(unit, '(A)') ""
    end subroutine write_hardware_metadata


    subroutine set_unknown_hardware_info( &
        info &
        )
        ! output
        type(hardware_info_t), intent(out) :: info

        info%hostname = unknown_value
        info%cpu_model = unknown_value
        info%logical_threads = unknown_value
        info%memory_total = unknown_value
        info%compiler = unknown_value
        info%build_preset = unknown_value
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


    function format_memory_gb( &
        mem_total_kb &
        ) result(text)
        ! inputs
        integer(int64), intent(in) :: mem_total_kb

        ! output
        character(len=32) :: text

        write(text, '(F0.1,A)') real(mem_total_kb, real64) / 1024.0_real64 / 1024.0_real64, " GB"
        text = adjustl(text)
    end function format_memory_gb


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
