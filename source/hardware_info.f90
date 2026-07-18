module hardware_info
    ! imports
    use iso_fortran_env, only: int32, int64
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
#ifndef FFB_GIT_COMMIT
#define FFB_GIT_COMMIT "unknown"
#endif

    private

    public :: hardware_info_t
    public :: collect_hardware_info
    public :: collect_image_host_names
    public :: print_hardware_summary
    public :: print_image_host_table
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
    character(len=*), parameter :: latest_commit = FFB_GIT_COMMIT
    character(len=256) :: image_host_names[*]

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


    subroutine collect_image_host_names()
        ! locals
        character(len=256) :: host_name
        character(len=256) :: slurm_job_id
        character(len=256) :: tmp_file
        character(len=32) :: image_id_text
        character(len=:), allocatable :: command_host_name
        integer :: env_status

        host_name = ""

        write(image_id_text, '(I0)') this_image()
        slurm_job_id = ""
        call get_environment_variable("SLURM_JOB_ID", slurm_job_id, status=env_status)

        if (len_trim(slurm_job_id) > 0) then
            tmp_file = "/tmp/.ffb_hostname_" // trim(slurm_job_id) // "_" // &
                trim(image_id_text) // ".tmp"
        else
            tmp_file = "/tmp/.ffb_hostname_" // trim(image_id_text) // ".tmp"
        end if

        command_host_name = unknown_value
        call read_command_output("hostname", command_host_name, tmp_file)

        if (len_trim(command_host_name) > 0 .and. trim(command_host_name) /= unknown_value) then
            host_name = command_host_name
        else
            call get_environment_variable("HOSTNAME", host_name, status=env_status)

            if (len_trim(host_name) == 0) then
                host_name = ""
                call get_environment_variable("SLURMD_NODENAME", host_name, status=env_status)
            end if
        end if

        if (len_trim(host_name) == 0) then
            image_host_names = unknown_value
        else
            image_host_names = trim(adjustl(host_name))
        end if

        sync all
    end subroutine collect_image_host_names


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
        print '(A,T27,A,A)', "latest commit", "= ", trim(latest_commit)
    end subroutine print_hardware_summary


    subroutine print_image_host_table( &
        n_images, n_images_x &
        )
        ! inputs
        integer(int32), intent(in) :: n_images
        integer(int32), intent(in) :: n_images_x

        ! locals
        integer(int32) :: image_id
        integer(int32) :: image_x
        integer(int32) :: image_y

        if (this_image() /= 1) then
            return
        end if

        print '(A)', "image | [    X /    Y ] | host"
        print '(A)', "--------------------------------------------------------------------------------"

        do image_id = 1, n_images
            image_x = modulo(image_id - 1, n_images_x) + 1
            image_y = (image_id - 1) / n_images_x + 1
            print '(I5,A,I4,A,I4,A,A)', image_id, " | [ ", image_x, " / ", image_y, " ] | ", &
                trim(image_host_names[image_id])
        end do

        call print_image_host_buckets(n_images)
    end subroutine print_image_host_table


    subroutine print_image_host_buckets( &
        n_images &
        )
        ! inputs
        integer(int32), intent(in) :: n_images

        ! locals
        character(len=*), parameter :: table_separator = "--------------------------------------------------------------------------------"
        integer(int32), parameter :: table_width = len(table_separator)
        integer(int32) :: image_id
        integer(int32) :: bucket_image_id
        character(len=32) :: image_id_text
        character(len=:), allocatable :: line
        character(len=:), allocatable :: continuation_prefix
        character(len=:), allocatable :: image_separator
        logical :: is_first_bucket
        logical :: line_has_image

        print '(A)', ""
        print '(A)', "host | images"
        print '(A)', table_separator

        do image_id = 1, n_images
            is_first_bucket = .true.

            do bucket_image_id = 1, image_id - 1
                if (trim(image_host_names[bucket_image_id]) == &
                    trim(image_host_names[image_id])) then
                    is_first_bucket = .false.
                    exit
                end if
            end do

            if (.not. is_first_bucket) then
                cycle
            end if

            line = trim(image_host_names[image_id]) // " | [ "
            continuation_prefix = repeat(" ", len_trim(image_host_names[image_id])) // " | "
            line_has_image = .false.

            do bucket_image_id = image_id, n_images
                if (trim(image_host_names[bucket_image_id]) == &
                    trim(image_host_names[image_id])) then
                    write(image_id_text, '(I0)') bucket_image_id

                    if (line_has_image) then
                        image_separator = ", "
                    else
                        image_separator = ""
                    end if

                    if (len(line) + len(image_separator) + len_trim(image_id_text) + &
                        len(" ]") + 1 > table_width) then
                        if (line_has_image) then
                            print '(A)', trim(line) // ","
                            line = continuation_prefix
                        else
                            print '(A)', trim(image_host_names[image_id]) // " |"
                            line = " | [ "
                            continuation_prefix = " | "
                        end if

                        image_separator = ""
                        line_has_image = .false.
                    end if

                    line = line // image_separator // trim(image_id_text)
                    line_has_image = .true.
                end if
            end do

            print '(A)', trim(line) // " ]"
        end do

        print '(A)', ""
    end subroutine print_image_host_buckets


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
        command, command_output, tmp_file &
        )
        ! inputs
        character(len=*), intent(in) :: command
        character(len=*), intent(in), optional :: tmp_file

        ! output
        character(len=:), allocatable, intent(out) :: command_output

        ! locals
        character(len=*), parameter :: default_tmp_file = ".ffb_hardware_info.tmp"
        character(len=512) :: line
        integer :: unit
        integer :: io_stat
        integer :: cmdstat
        integer :: exitstat
        character(len=:), allocatable :: tmp_file_resolved

        if (present(tmp_file)) then
            tmp_file_resolved = trim(tmp_file)
        else
            tmp_file_resolved = default_tmp_file
        end if

        command_output = unknown_value

        call execute_command_line(trim(command) // " > " // trim(tmp_file_resolved) // " 2>/dev/null", &
            exitstat=exitstat, cmdstat=cmdstat)

        if (cmdstat == 0) then
            open(newunit=unit, file=trim(tmp_file_resolved), form="formatted", status="old", action="read", iostat=io_stat)

            if (io_stat == 0) then
                read(unit, '(A)', iostat=io_stat) line
                if (io_stat == 0 .and. len_trim(line) > 0) then
                    command_output = trim(adjustl(line))
                end if
                close(unit)
            end if
        end if

        call execute_command_line("rm -f " // trim(tmp_file_resolved))
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
