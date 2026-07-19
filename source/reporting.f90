module reporting
    ! imports
    use iso_fortran_env, only: int32, int64, real64, output_unit
    use domain, only: domain_t, print_domain_summary
    use hardware_info, only: hardware_info_t, print_hardware_summary, print_image_host_table
    use settings, only: N_X, N_Y, N_STEPS, N_CELLS, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, &
        DIST_FUNC_LAYOUT, USE_UNROLLED_KERNELS, &
        RHO_0, OMEGA, U_MAX, N_SIN, U_WALL, U_LID, RHO_IN, RHO_OUT, sim_mode_to_string, &
        calculate_reynolds
    implicit none
    private

    public :: print_run_summary
    public :: print_launch_timestamp
    public :: print_progress_status
    public :: print_finish_timestamp
    public :: print_execution_summary

contains

    subroutine print_run_summary( &
        machine_info, domain_info, sim_mode, &
        export_macros, export_endpoint_states, export_interval, &
        export_num, dist_function_buffers_bytes, macro_field_buffers_bytes, &
        total_buffer_bytes, total_bytes_per_cell &
        )
        ! inputs
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
        real(real64) :: gb_per_byte

        gb_per_byte = 1.0e-9_real64

        print '(A)', ""
        call print_hardware_summary(machine_info)

        print '(A)', ""
        call print_image_host_table(domain_info%n_images, domain_info%n_images_x)

        print '(A)', "--- [ simulation parameters ] --------------------------------------------------"
        print '(A,T27,A,A)',     "SIM_MODE", "= ", trim(sim_mode_to_string(sim_mode))

        select case (sim_mode)
        case (SIM_SHEAR_WAVE)
            print '(A,T27,A,F8.6)', "rho_0", "= ", RHO_0
            print '(A,T27,A,F8.6)', "omega", "= ", OMEGA
            print '(A,T27,A,F8.6)', "u_max", "= ", U_MAX
            print '(A,T27,A,F8.6)', "n_sin", "= ", N_SIN
        case (SIM_COUETTE_FLOW)
            print '(A,T27,A,F8.6)', "rho_0", "= ", RHO_0
            print '(A,T27,A,F8.6)', "omega", "= ", OMEGA
            print '(A,T27,A,F8.6)', "u_wall", "= ", U_WALL
            print '(A,T27,A,F0.6)', "reynolds", "= ", calculate_reynolds(sim_mode)
        case (SIM_POISEUILLE_FLOW)
            print '(A,T27,A,F8.6)', "rho_0", "= ", RHO_0
            print '(A,T27,A,F8.6)', "omega", "= ", OMEGA
            print '(A,T27,A,F8.6)', "rho_in", "= ", RHO_IN
            print '(A,T27,A,F8.6)', "rho_out", "= ", RHO_OUT
        case (SIM_SLIDING_LID)
            print '(A,T27,A,F8.6)', "rho_0", "= ", RHO_0
            print '(A,T27,A,F8.6)', "omega", "= ", OMEGA
            print '(A,T27,A,F8.6)', "u_lid", "= ", U_LID
            print '(A,T27,A,F0.6)', "reynolds", "= ", calculate_reynolds(sim_mode)
        case default
            error stop "error: unknown sim mode in print_run_summary()"
        end select

        ! parameter info
        print '(A)', ""
        print '(A)', "--- [ other parameters ] -------------------------------------------------------"
        print '(A,T27,A,I0)',    "N_X_TOTAL", "= ", N_X
        print '(A,T27,A,I0)',    "N_Y_TOTAL", "= ", N_Y
        print '(A,T27,A,I0)',    "N_STEPS", "= ", N_STEPS
        print '(A,T27,A,A)',     "dist_func_layout", "= ", DIST_FUNC_LAYOUT
        print '(A,T27,A,L1)',    "use_unrolled_kernels", "= ", USE_UNROLLED_KERNELS
        print '(A,T27,A,L1)',    "export_macros", "= ", export_macros
        print '(A,T27,A,L1)',    "export_endpoint_states", "= ", export_endpoint_states
        print '(A,T27,A,I0)',    "export_interval", "= ", export_interval
        print '(A,T27,A,A)',     "output_dir", "= ", "output/" // trim(export_num)

        call print_domain_summary(domain_info)
        print '(A)', ""

        ! memory info
        print '(A,T47,A,T50,A,T64,A,T67,A)', "memory usage", "|", "per cell [B]", "|", "all cells [GB]"
        print '(A)', "--------------------------------------------------------------------------------"
        print '(A,T47,A,T50,I12,T64,A,T67,F14.3)', "dist function buffers", "|", &
            nint(real(dist_function_buffers_bytes, real64) / real(N_CELLS, real64), int64), "|", &
            real(dist_function_buffers_bytes, real64) * gb_per_byte
        print '(A,T47,A,T50,I12,T64,A,T67,F14.3)', "macro field buffers", "|", &
            nint(real(macro_field_buffers_bytes, real64) / real(N_CELLS, real64), int64), "|", &
            real(macro_field_buffers_bytes, real64) * gb_per_byte
        print '(A,T47,A,T50,I12,T64,A,T67,F14.3)', "total", "|", &
            nint(total_bytes_per_cell, int64), "|", real(total_buffer_bytes, real64) * gb_per_byte
        print *
    end subroutine print_run_summary


    subroutine print_launch_timestamp()
        ! temp
        character(len=10) :: current_time

        call date_and_time(time=current_time)
        print '(A)', "[" // current_time(1:2) // ":" // current_time(3:4) // ":" // current_time(5:6) // "] &
            launched ------------------------------------------------------------"
    end subroutine print_launch_timestamp


    subroutine print_progress_status( &
        step, clock_start, clock_rate, interactive_progress &
        )
        ! inputs
        integer(int32), intent(in) :: step
        integer(int64), intent(in) :: clock_start
        integer(int64), intent(in) :: clock_rate
        logical, intent(in) :: interactive_progress

        ! temp
        integer(int64) :: clock_now
        integer(int64) :: eta_total_seconds
        integer(int64) :: eta_hours
        integer(int64) :: eta_minutes
        integer(int64) :: eta_seconds_int
        real(real64) :: elapsed_now
        real(real64) :: eta_seconds
        real(real64) :: sim_percent
        real(real64) :: avg_millisec_per_step
        character(len=10) :: current_time

        call system_clock(clock_now)
        call date_and_time(time=current_time)

        elapsed_now = real(clock_now - clock_start, real64) / real(clock_rate, real64)
        avg_millisec_per_step = 1000.0_real64 * elapsed_now / real(step, real64)
        eta_seconds = elapsed_now * real(N_STEPS - step, real64) / real(step, real64)
        sim_percent = 100.0_real64 * real(step, real64) / real(N_STEPS, real64)

        eta_total_seconds = int(eta_seconds + 0.999_real64, int64)
        eta_hours = eta_total_seconds / 3600_int64
        eta_minutes = modulo(eta_total_seconds, 3600_int64) / 60_int64
        eta_seconds_int = modulo(eta_total_seconds, 60_int64)

        if (interactive_progress) then ! progress output with carriage return to overwrite previous line
            write(output_unit,'(A,3A,F6.2,A,I2.2,A,I2.2,A,I2.2,A,I0,A,I0,A,F0.3,A)', advance='no') &
                achar(13), "[", current_time(1:2)//":"//current_time(3:4)//":"//current_time(5:6), "] ", &
                sim_percent, " %  (T-", &
                eta_hours, ":", eta_minutes, ":", eta_seconds_int, ")  ", &
                step, "/", N_STEPS, " steps  |  avg step: ", avg_millisec_per_step, " ms   "

            flush(output_unit)
        else ! fallback to non-interactive progress output
            print '(3A,F6.2,A,I2.2,A,I2.2,A,I2.2,A,I0,A,I0,A,F0.3,A)', &
                "[", current_time(1:2)//":"//current_time(3:4)//":"//current_time(5:6), "] ", &
                sim_percent, " %  (T-", &
                eta_hours, ":", eta_minutes, ":", eta_seconds_int, ")  ", &
                step, "/", N_STEPS, " steps  |  avg step: ", avg_millisec_per_step, " ms   "
        end if
    end subroutine print_progress_status


    subroutine print_finish_timestamp( &
        interactive_progress &
        )
        ! inputs
        logical, intent(in) :: interactive_progress

        ! temp
        character(len=10) :: current_time

        if (interactive_progress) then
            print *
        end if

        call date_and_time(time=current_time)
        print '(A)', "[" // current_time(1:2) // ":" // current_time(3:4) // ":" // current_time(5:6) // "] &
            finished ------------------------------------------------------------"
    end subroutine print_finish_timestamp


    subroutine print_execution_summary( &
        best_seconds, worst_seconds, best_image_ids, worst_image_ids, seconds_per_step, mlups &
        )
        ! inputs
        real(real64), intent(in) :: best_seconds(:)
        real(real64), intent(in) :: worst_seconds(:)
        integer(int32), intent(in) :: best_image_ids(:)
        integer(int32), intent(in) :: worst_image_ids(:)
        real(real64), intent(in) :: seconds_per_step
        real(real64), intent(in) :: mlups
        character(len=32) :: total_time_text
        character(len=32) :: step_time_text

        if (size(best_seconds) < 5 .or. size(worst_seconds) < 5 .or. &
            size(best_image_ids) < 5 .or. size(worst_image_ids) < 5) then
            error stop "error: execution timing arrays are too small"
        end if

        call format_compact_real(worst_seconds(5), total_time_text)
        call format_compact_real(seconds_per_step * 1000.0_real64, step_time_text)

        print '(A)', ""
        print '(A)', "image execution time spread |       best [sec] (image) |     worst [sec] (image)"
        print '(A)', "--------------------------------------------------------------------------------"

        call print_timing_spread_row( &
            "kernel compute", best_seconds(1), best_image_ids(1), worst_seconds(1), worst_image_ids(1))
        call print_timing_spread_row( &
            "halo sync", best_seconds(2), best_image_ids(2), worst_seconds(2), worst_image_ids(2))
        call print_timing_spread_row( &
            "halo transfer", best_seconds(3), best_image_ids(3), worst_seconds(3), worst_image_ids(3))
        call print_timing_spread_row( &
            "other", best_seconds(4), best_image_ids(4), worst_seconds(4), worst_image_ids(4))
        call print_timing_spread_row( &
            "total", best_seconds(5), best_image_ids(5), worst_seconds(5), worst_image_ids(5))

        print '(A)', ""
        print '(A)', "--- [ perf metrics ] -----------------------------------------------------------"
        print '(A,T24,A,I0,A,I0,A,I0,A)', "sim size [X/Y/N]", "= [ ", N_X, " / ", N_Y, " / ", N_STEPS, " ]"
        print '(A,T24,A,A,A)',            "total time", "= ", trim(total_time_text), " sec"
        print '(A,T24,A,A,A)',            "step time", "= ", trim(step_time_text), " ms"
        print '(A,T24,A,I0)',             "MLUPS", "= ", int(mlups, int64)
    end subroutine print_execution_summary


    subroutine format_compact_real(value, value_text)
        ! inputs
        real(real64), intent(in) :: value

        ! outputs
        character(len=*), intent(out) :: value_text

        ! temp
        character(len=32) :: raw_text

        write(raw_text, '(F0.3)') value
        raw_text = adjustl(raw_text)

        if (raw_text(1:1) == ".") then
            value_text = "0" // trim(raw_text)
        else if (len_trim(raw_text) >= 2 .and. raw_text(1:2) == "-.") then
            value_text = "-0" // trim(raw_text(2:))
        else
            value_text = trim(raw_text)
        end if
    end subroutine format_compact_real


    subroutine format_timing_measurement(seconds, image_id, value_text)
        ! inputs
        real(real64), intent(in) :: seconds
        integer(int32), intent(in) :: image_id

        ! output
        character(len=*), intent(out) :: value_text

        ! temp
        character(len=64) :: raw_text
        integer(int32) :: value_length

        write(raw_text, '(F0.3," (",I0,")")') seconds, image_id
        value_length = len_trim(raw_text)
        if (value_length > len(value_text)) then
            error stop "error: timing measurement text is too wide"
        end if

        value_text = repeat(" ", len(value_text))
        value_text(len(value_text) - value_length + 1:) = trim(raw_text)
    end subroutine format_timing_measurement


    subroutine print_timing_spread_row( &
        row_name, best_seconds, best_image_id, worst_seconds, worst_image_id &
        )
        ! inputs
        character(len=*), intent(in) :: row_name
        real(real64), intent(in) :: best_seconds
        integer(int32), intent(in) :: best_image_id
        real(real64), intent(in) :: worst_seconds
        integer(int32), intent(in) :: worst_image_id

        ! temp
        character(len=27) :: row_text
        character(len=24) :: best_text
        character(len=23) :: worst_text

        if (len_trim(row_name) > len(row_text)) then
            error stop "error: timing row name is too wide"
        end if

        row_text = " "
        row_text(1:len_trim(row_name)) = trim(row_name)
        call format_timing_measurement(best_seconds, best_image_id, best_text)
        call format_timing_measurement(worst_seconds, worst_image_id, worst_text)

        print '(A)', row_text // " | " // best_text // " | " // worst_text
    end subroutine print_timing_spread_row


end module reporting
