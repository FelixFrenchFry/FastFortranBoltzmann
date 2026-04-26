module export
    ! imports
    use iso_fortran_env, only: int32, int64, real32
    use settings, only: SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, &
        shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t, sim_mode_to_string
    implicit none

    private

    public :: should_export_step
    public :: export_selected_data
    public :: export_metadata

contains

    pure function should_export_step( &
        N_STEPS, step, export_interval, export_initial_state, export_final_state &
        ) result(write_step)
        ! read-only inputs
        integer(int32), intent(in) :: N_STEPS
        integer(int32), intent(in) :: step
        integer(int32), intent(in) :: export_interval
        logical, intent(in) :: export_initial_state
        logical, intent(in) :: export_final_state

        ! output
        logical :: write_step

        ! no output if export is disabled
        if (export_interval <= 0) then
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
        export_rho, export_u_x, export_u_y, export_u_mag, &
        output_dir_name, export_num, suffix_num, rho, u_x, u_y &
        )
        ! read-only inputs
        logical, intent(in) :: export_rho
        logical, intent(in) :: export_u_x
        logical, intent(in) :: export_u_y
        logical, intent(in) :: export_u_mag
        character(len=*), intent(in) :: output_dir_name
        character(len=*), intent(in) :: export_num
        integer(int32), intent(in) :: suffix_num
        real(real32), intent(in) :: rho(:,:)
        real(real32), intent(in) :: u_x(:,:)
        real(real32), intent(in) :: u_y(:,:)

        ! temp
        real(real32), allocatable :: velocity_mag(:,:)

        ! export selected scalar fields
        if (export_rho) then
            call export_scalar_field(rho, "density", output_dir_name, export_num, suffix_num)
        end if

        if (export_u_x) then
            call export_scalar_field(u_x, "velocity_x", output_dir_name, export_num, suffix_num)
        end if

        if (export_u_y) then
            call export_scalar_field(u_y, "velocity_y", output_dir_name, export_num, suffix_num)
        end if

        if (export_u_mag) then
            allocate(velocity_mag(size(u_x, 1), size(u_x, 2)))
            velocity_mag = sqrt(u_x * u_x + u_y * u_y) ! element-wise sqrt of velocity magnitude
            call export_scalar_field(velocity_mag, "velocity_mag", output_dir_name, export_num, suffix_num)
        end if
    end subroutine export_selected_data


    subroutine export_metadata( &
        sim_mode, shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
        N_X, N_Y, N_STEPS, N_CELLS, N_DIRS, pi, &
        export_rho, export_u_x, export_u_y, export_u_mag, export_interval, output_dir_name, export_num, &
        export_initial_state, export_final_state &
        )
        ! read-only inputs
        character(len=*), intent(in) :: output_dir_name
        character(len=*), intent(in) :: export_num
        integer(int32), intent(in) :: sim_mode
        type(shear_wave_params_t), intent(in) :: shear_wave_params
        type(couette_flow_params_t), intent(in) :: couette_flow_params
        type(poiseuille_flow_params_t), intent(in) :: poiseuille_flow_params
        type(sliding_lid_params_t), intent(in) :: sliding_lid_params
        integer(int32), intent(in) :: N_X
        integer(int32), intent(in) :: N_Y
        integer(int32), intent(in) :: N_STEPS
        integer(int64), intent(in) :: N_CELLS
        integer(int32), intent(in) :: N_DIRS
        real(real32), intent(in) :: pi
        logical, intent(in) :: export_rho
        logical, intent(in) :: export_u_x
        logical, intent(in) :: export_u_y
        logical, intent(in) :: export_u_mag
        integer(int32), intent(in) :: export_interval
        logical, intent(in) :: export_initial_state
        logical, intent(in) :: export_final_state

        ! temp
        character(len=:), allocatable :: output_path
        character(len=:), allocatable :: file_path
        integer :: unit
        integer :: io_stat
        real(real32) :: k

        ! assemble output path and metadata filename
        output_path = trim(output_dir_name) // "/" // trim(export_num)
        file_path = output_path // "/config.json"

        call ensure_output_directory(output_path)

        ! document run configuration as .json
        open(newunit=unit, file=trim(file_path), form="formatted", status="replace", &
            action="write", iostat=io_stat)

        if (io_stat /= 0) then
            error stop "error: could not open metadata output file"
        end if

        write(unit, '(A)') "{"
        write(unit, '(A,A,A)') '  "sim_mode": "', trim(sim_mode_to_string(sim_mode)), '",'

        select case (sim_mode)
        case (SIM_SHEAR_WAVE)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real32_to_json(shear_wave_params%rho_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real32_to_json(shear_wave_params%omega)), ','
            write(unit, '(A,A,A)') '  "u_max": ', trim(real32_to_json(shear_wave_params%u_max)), ','
            write(unit, '(A,A,A)') '  "n_sin": ', trim(real32_to_json(shear_wave_params%n_sin)), ','

        case (SIM_COUETTE_FLOW)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real32_to_json(couette_flow_params%rho_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real32_to_json(couette_flow_params%omega)), ','
            write(unit, '(A,A,A)') '  "u_wall": ', trim(real32_to_json(couette_flow_params%u_wall)), ','

        case (SIM_POISEUILLE_FLOW)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real32_to_json(poiseuille_flow_params%rho_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real32_to_json(poiseuille_flow_params%omega)), ','
            write(unit, '(A,A,A)') '  "rho_in": ', trim(real32_to_json(poiseuille_flow_params%rho_in)), ','
            write(unit, '(A,A,A)') '  "rho_out": ', trim(real32_to_json(poiseuille_flow_params%rho_out)), ','
            continue

        case (SIM_SLIDING_LID)
            write(unit, '(A,A,A)') '  "rho_0": ', trim(real32_to_json(sliding_lid_params%rho_0)), ','
            write(unit, '(A,A,A)') '  "omega": ', trim(real32_to_json(sliding_lid_params%omega)), ','
            write(unit, '(A,A,A)') '  "u_wall": ', trim(real32_to_json(sliding_lid_params%u_wall)), ','

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
        write(unit, '(A,A,A)') '  "export_rho": ', trim(logical_to_json(export_rho)), ','
        write(unit, '(A,A,A)') '  "export_u_x": ', trim(logical_to_json(export_u_x)), ','
        write(unit, '(A,A,A)') '  "export_u_y": ', trim(logical_to_json(export_u_y)), ','
        write(unit, '(A,A,A)') '  "export_u_mag": ', trim(logical_to_json(export_u_mag)), ','
        write(unit, '(A,I0,A)') '  "export_interval": ', export_interval, ','
        write(unit, '(A,A,A)') '  "export_initial_state": ', trim(logical_to_json(export_initial_state)), ','
        write(unit, '(A,A,A)') '  "export_final_state": ', trim(logical_to_json(export_final_state)), ','
        write(unit, '(A)') ""
        write(unit, '(A,A,A)') '  "output_dir_name": "', trim(output_dir_name), '",'
        write(unit, '(A,A,A)') '  "export_num": "', trim(export_num), '",'
        write(unit, '(A)') '  "file_dtype": "real32"'
        write(unit, '(A)') "}"

        close(unit)
    end subroutine export_metadata


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


    function real32_to_json( &
        value &
        ) result(text)
        ! read-only inputs
        real(real32), intent(in) :: value

        ! output
        character(len=24) :: text

        write(text, '(ES16.8)') value
        text = adjustl(text)
    end function real32_to_json


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
        real(real32), intent(in) :: field(:,:)
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

        ! write raw real32 field data to file
        write(unit, iostat=io_stat) field

        if (io_stat /= 0) then
            close(unit)
            error stop "error: could not write binary output file"
        end if

        close(unit)
    end subroutine write_binary_field

end module export
