module export
    ! imports
    use iso_fortran_env, only: int32, int64, real64
    use domain, only: domain_t
    use hardware_info, only: hardware_info_t, write_hardware_metadata
    use settings, only: N_X, N_Y, N_STEPS, N_CELLS, N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, FP, FP_DTYPE, &
        DIST_FUNC_LAYOUT, USE_UNROLLED_KERNELS, &
        RHO_0, OMEGA, U_MAX, N_SIN, U_WALL, U_LID, RHO_IN, RHO_OUT, sim_mode_to_string, &
        calculate_reynolds
    implicit none

    private

    public :: should_export_step
    public :: export_selected_data_distributed
    public :: export_metadata

contains

    pure function should_export_step( &
        step, export_endpoint_states, export_interval &
        ) result(write_step)
        ! read-only inputs
        integer(int32), intent(in) :: step
        logical, intent(in) :: export_endpoint_states
        integer(int32), intent(in) :: export_interval

        ! output
        logical :: write_step

        ! no output if export is disabled
        if (export_interval <= 0) then
            write_step = .false.

        ! optional output at initial endpoint state
        else if (step == 0) then
            write_step = export_endpoint_states

        ! regular interval output plus optional output at final endpoint state
        else
            write_step = mod(step, export_interval) == 0 .or. &
                (export_endpoint_states .and. step == N_STEPS)
        end if
    end function should_export_step


    subroutine export_selected_data_distributed( &
        domain_info, export_num, suffix_num, rho, u_x, u_y &
        )
        ! read-only inputs
        type(domain_t), intent(in) :: domain_info
        character(len=*), intent(in) :: export_num
        integer(int32), intent(in) :: suffix_num
        real(FP), intent(in) :: rho(:,:)
        real(FP), intent(in) :: u_x(:,:)
        real(FP), intent(in) :: u_y(:,:)

        ! temp
        real(FP), allocatable :: velocity_mag(:,:)

        ! export macro scalar fields
        call export_scalar_field_distributed( &
            domain_info, rho, "density", export_num, suffix_num)
        call export_scalar_field_distributed( &
            domain_info, u_x, "velocity_x", export_num, suffix_num)
        call export_scalar_field_distributed( &
            domain_info, u_y, "velocity_y", export_num, suffix_num)

        allocate(velocity_mag(size(u_x, 1), size(u_x, 2)))
        velocity_mag = sqrt(u_x * u_x + u_y * u_y) ! element-wise sqrt of velocity magnitude
        call export_scalar_field_distributed( &
            domain_info, velocity_mag, "velocity_mag", export_num, suffix_num)
    end subroutine export_selected_data_distributed


    subroutine export_metadata( &
        machine_info, domain_info, sim_mode, &
        export_macros, export_endpoint_states, export_interval, export_num, &
        dist_function_buffers_bytes, macro_field_buffers_bytes, &
        total_buffer_bytes, total_bytes_per_cell &
        )
        ! read-only inputs
        type(hardware_info_t), intent(in) :: machine_info
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: sim_mode
        logical, intent(in) :: export_macros
        logical, intent(in) :: export_endpoint_states
        integer(int32), intent(in) :: export_interval
        character(len=*), intent(in) :: export_num
        integer(int64), intent(in) :: dist_function_buffers_bytes
        integer(int64), intent(in) :: macro_field_buffers_bytes
        integer(int64), intent(in) :: total_buffer_bytes
        real(real64), intent(in) :: total_bytes_per_cell

        ! temp
        character(len=:), allocatable :: output_path
        character(len=:), allocatable :: file_path
        integer :: unit
        integer :: io_stat
        real(real64) :: gb_per_byte
        real(real64) :: halo_cell_percent

        ! assemble output path and metadata filename
        output_path = "output/" // trim(export_num)
        file_path = output_path // "/config.json"
        gb_per_byte = 1.0e-9_real64
        halo_cell_percent = 100.0_real64 * &
            real((domain_info%n_x + 2) * (domain_info%n_y + 2) - &
            domain_info%n_x * domain_info%n_y, real64) / &
            real(domain_info%n_x * domain_info%n_y, real64)

        call ensure_output_directory(output_path)

        ! document run configuration as .json
        open(newunit=unit, file=trim(file_path), form="formatted", status="replace", &
            action="write", iostat=io_stat)

        if (io_stat /= 0) then
            error stop "error: could not open metadata output file"
        end if

        write(unit, '(A)') "{"
        call write_hardware_metadata(unit, machine_info)
        write(unit, '(A,A,A)') '  "SIM_MODE": "', trim(sim_mode_to_string(sim_mode)), '",'

        select case (sim_mode)
        case (SIM_SHEAR_WAVE)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real_to_json(RHO_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real_to_json(OMEGA)), ','
            write(unit, '(A,A,A)') '  "u_max": ', trim(real_to_json(U_MAX)), ','
            write(unit, '(A,A,A)') '  "n_sin": ', trim(real_to_json(N_SIN)), ','

        case (SIM_COUETTE_FLOW)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real_to_json(RHO_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real_to_json(OMEGA)), ','
            write(unit, '(A,A,A)') '  "u_wall": ', trim(real_to_json(U_WALL)), ','
            write(unit, '(A,A,A)') '  "reynolds": ', trim(real_to_json(calculate_reynolds(sim_mode))), ','

        case (SIM_POISEUILLE_FLOW)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real_to_json(RHO_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real_to_json(OMEGA)), ','
            write(unit, '(A,A,A)') '  "rho_in": ', trim(real_to_json(RHO_IN)), ','
            write(unit, '(A,A,A)') '  "rho_out": ', trim(real_to_json(RHO_OUT)), ','
            continue

        case (SIM_SLIDING_LID)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real_to_json(RHO_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real_to_json(OMEGA)), ','
            write(unit, '(A,A,A)') '  "u_lid": ', trim(real_to_json(U_LID)), ','
            write(unit, '(A,A,A)') '  "reynolds": ', trim(real_to_json(calculate_reynolds(sim_mode))), ','

        case default
            error stop "error: unknown sim mode in export_metadata()"
        end select

        write(unit, '(A)') ""
        write(unit, '(A,I0,A)') '  "N_X": ', N_X, ','
        write(unit, '(A,I0,A)') '  "N_Y": ', N_Y, ','
        write(unit, '(A,I0,A)') '  "N_STEPS": ', N_STEPS, ','
        write(unit, '(A,I0,A)') '  "N_CELLS": ', N_CELLS, ','
        write(unit, '(A,I0,A)') '  "N_DIRS": ', N_DIRS, ','
        write(unit, '(A)') ""
        write(unit, '(A,A,A)') '  "dist_func_layout": "', trim(DIST_FUNC_LAYOUT), '",'
        write(unit, '(A,A,A)') '  "use_unrolled_kernels": ', trim(logical_to_json(USE_UNROLLED_KERNELS)), ','
        write(unit, '(A)') ""
    #ifdef FFB_FP64
        write(unit, '(A)') '  "FFB_FP64": true,'
    #else
        write(unit, '(A)') '  "FFB_FP64": false,'
    #endif
    #ifdef FFB_DENSITY_CHECKS
        write(unit, '(A)') '  "FFB_DENSITY_CHECKS": true,'
    #else
        write(unit, '(A)') '  "FFB_DENSITY_CHECKS": false,'
    #endif
        write(unit, '(A)') ""
        write(unit, '(A,A,A)') '  "export_macros": ', trim(logical_to_json(export_macros)), ','
        write(unit, '(A,A,A)') '  "export_endpoint_states": ', trim(logical_to_json(export_endpoint_states)), ','
        write(unit, '(A,I0,A)') '  "export_interval": ', export_interval, ','
        write(unit, '(A)') ""
        write(unit, '(A)') '  "domain_decomposition": {'
        write(unit, '(A,I0,A)') '    "coarray_images": ', domain_info%n_images, ','
        write(unit, '(A,I0,A)') '    "image_grid_x": ', domain_info%n_images_x, ','
        write(unit, '(A,I0,A)') '    "image_grid_y": ', domain_info%n_images_y, ','
        write(unit, '(A,I0,A)') '    "local_n_x": ', domain_info%n_x, ','
        write(unit, '(A,I0,A)') '    "local_n_y": ', domain_info%n_y, ','
        write(unit, '(A,A)') '    "halo_cells_percent": ', trim(real64_to_json(halo_cell_percent))
        write(unit, '(A)') '  },'
        write(unit, '(A)') ""
        write(unit, '(A)') '  "memory_usage": {'
        write(unit, '(A,I0,A)') '    "dist_function_buffers_per_cell_B": ', &
            nint(real(dist_function_buffers_bytes, real64) / real(N_CELLS, real64), int64), ','
        write(unit, '(A,A,A)') '    "dist_function_buffers_total_GB": ', &
            trim(real64_to_json(real(dist_function_buffers_bytes, real64) * gb_per_byte)), ','
        write(unit, '(A,I0,A)') '    "macro_field_buffers_per_cell_B": ', &
            nint(real(macro_field_buffers_bytes, real64) / real(N_CELLS, real64), int64), ','
        write(unit, '(A,A,A)') '    "macro_field_buffers_total_GB": ', &
            trim(real64_to_json(real(macro_field_buffers_bytes, real64) * gb_per_byte)), ','
        write(unit, '(A,I0,A)') '    "total_per_cell_B": ', nint(total_bytes_per_cell, int64), ','
        write(unit, '(A,A)') '    "total_GB": ', trim(real64_to_json(real(total_buffer_bytes, real64) * gb_per_byte))
        write(unit, '(A)') '  },'
        write(unit, '(A)') ""
        write(unit, '(A,A,A)') '  "output_dir": "', trim(output_path), '",'
        write(unit, '(A,A,A)') '  "file_dtype": "', FP_DTYPE, '"'
        write(unit, '(A)') "}"

        close(unit)
    end subroutine export_metadata


    subroutine export_scalar_field_distributed( &
        domain_info, local_field, field_name, export_num, suffix_num &
        )
        ! read-only inputs
        type(domain_t), intent(in) :: domain_info
        real(FP), intent(in) :: local_field(:,:)
        character(len=*), intent(in) :: field_name
        character(len=*), intent(in) :: export_num
        integer(int32), intent(in) :: suffix_num

        ! temp
        integer(int32) :: image_id
        integer(int32) :: image_x
        integer(int32) :: image_y
        integer(int32) :: x_global_start
        integer(int32) :: x_global_end
        integer(int32) :: y_global_start
        integer(int32) :: y_global_end
        character(len=:), allocatable :: output_path
        character(len=:), allocatable :: file_path
        real(FP), allocatable :: global_field(:,:)
        real(FP), allocatable :: export_buffer(:,:)[:]

        if (size(local_field, 1) /= domain_info%n_x .or. &
            size(local_field, 2) /= domain_info%n_y) then
            error stop "error: local distributed export field has wrong shape"
        end if

        allocate(export_buffer(domain_info%n_x, domain_info%n_y)[*])
        export_buffer(:, :) = local_field(:, :)

        sync all

        if (this_image() == 1) then
            allocate(global_field(N_X, N_Y))

            do image_id = 1, domain_info%n_images
                image_x = modulo(image_id - 1, domain_info%n_images_x) + 1
                image_y = (image_id - 1) / domain_info%n_images_x + 1

                x_global_start = (image_x - 1) * domain_info%n_x + 1
                x_global_end = image_x * domain_info%n_x
                y_global_start = (image_y - 1) * domain_info%n_y + 1
                y_global_end = image_y * domain_info%n_y

                global_field(x_global_start:x_global_end, y_global_start:y_global_end) = &
                    export_buffer(:, :)[image_id]
            end do

            output_path = "output/" // trim(export_num)
            file_path = output_path // "/" // trim(field_name) // format_step_suffix(suffix_num) // ".bin"

            call ensure_output_directory(output_path)
            call write_binary_field(global_field, file_path)
        end if

        sync all
    end subroutine export_scalar_field_distributed


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


    pure function logical_to_json( &
        value &
        ) result(text)
        ! read-only inputs
        logical, intent(in) :: value

        ! output
        character(len=5) :: text

        if (value) then
            text = "true"
        else
            text = "false"
        end if
    end function logical_to_json


    function real_to_json( &
        value &
        ) result(text)
        ! read-only inputs
        real(FP), intent(in) :: value

        ! output
        character(len=32) :: text

    #ifdef FFB_FP64
        write(text, '(ES24.16)') value
    #else
        write(text, '(ES16.8)') value
    #endif
        text = adjustl(text)
    end function real_to_json


    function real64_to_json( &
        value &
        ) result(text)
        ! read-only inputs
        real(real64), intent(in) :: value

        ! output
        character(len=32) :: text

        write(text, '(F32.6)') value
        text = adjustl(text)
    end function real64_to_json


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
            error stop "error: failed to create output directory"
        end if
    end subroutine ensure_output_directory


    subroutine write_binary_field( &
        field, file_path &
        )
        ! read-only inputs
        real(FP), intent(in) :: field(:,:)
        character(len=*), intent(in) :: file_path

        ! temp
        integer :: unit
        integer :: io_stat

        ! open raw binary stream file
        open(newunit=unit, file=trim(file_path), access="stream", form="unformatted", &
            status="replace", action="write", iostat=io_stat)

        if (io_stat /= 0) then
            error stop "error: could not open binary output file"
        end if

        ! write raw floating point field data to file
        write(unit, iostat=io_stat) field

        if (io_stat /= 0) then
            close(unit)
            error stop "error: could not write binary output file"
        end if

        close(unit)
    end subroutine write_binary_field


end module export
