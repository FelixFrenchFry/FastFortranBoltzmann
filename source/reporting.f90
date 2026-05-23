module reporting
    ! imports
    use iso_fortran_env, only: int32, int64, real64, output_unit
    use domain, only: domain_t, print_domain_summary
    use hardware_info, only: hardware_info_t, print_hardware_summary
    use settings, only: N_X, N_Y, N_STEPS, N_CELLS, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, &
        USE_UNROLLED_KERNELS, USE_UNIVERSAL_KERNELS, &
        RHO_0, OMEGA, U_MAX, N_SIN, U_WALL, U_LID, RHO_IN, RHO_OUT, sim_mode_to_string
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
        export_rho, export_u_x, export_u_y, export_u_mag, export_interval, export_initial_state, &
        export_final_state, export_num, dist_function_buffers_bytes, macro_field_buffers_bytes, &
        total_buffer_bytes, total_bytes_per_cell &
        )
        ! inputs
        type(hardware_info_t), intent(in) :: machine_info
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: sim_mode
        logical, intent(in) :: export_rho
        logical, intent(in) :: export_u_x
        logical, intent(in) :: export_u_y
        logical, intent(in) :: export_u_mag
        integer(int32), intent(in) :: export_interval
        logical, intent(in) :: export_initial_state
        logical, intent(in) :: export_final_state
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
        print '(A)', "--- [ simulation parameters ] ---------------------------------------------"
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
        case (SIM_POISEUILLE_FLOW)
            print '(A,T27,A,F8.6)', "rho_0", "= ", RHO_0
            print '(A,T27,A,F8.6)', "omega", "= ", OMEGA
            print '(A,T27,A,F8.6)', "rho_in", "= ", RHO_IN
            print '(A,T27,A,F8.6)', "rho_out", "= ", RHO_OUT
        case (SIM_SLIDING_LID)
            print '(A,T27,A,F8.6)', "rho_0", "= ", RHO_0
            print '(A,T27,A,F8.6)', "omega", "= ", OMEGA
            print '(A,T27,A,F8.6)', "u_lid", "= ", U_LID
        case default
            error stop "error: unknown sim mode in print_run_summary()"
        end select

        ! parameter info
        print '(A)', ""
        print '(A)', "--- [ other parameters ] --------------------------------------------------"
        print '(A,T27,A,I0)',    "N_X_TOTAL", "= ", N_X
        print '(A,T27,A,I0)',    "N_Y_TOTAL", "= ", N_Y
        print '(A,T27,A,I0)',    "N_STEPS", "= ", N_STEPS
        print '(A,T27,A,L1)',    "use_unrolled_kernels", "= ", USE_UNROLLED_KERNELS
        print '(A,T27,A,L1)',    "use_universal_kernels", "= ", USE_UNIVERSAL_KERNELS
        print '(A,T27,A,L1)',    "distributed_coarrays", "= ", .true.
        print '(A,T27,A,L1)',    "export_rho", "= ", export_rho
        print '(A,T27,A,L1)',    "export_u_x", "= ", export_u_x
        print '(A,T27,A,L1)',    "export_u_y", "= ", export_u_y
        print '(A,T27,A,L1)',    "export_u_mag", "= ", export_u_mag
        print '(A,T27,A,I0)',    "export_interval", "= ", export_interval
        print '(A,T27,A,L1)',    "export_initial_state", "= ", export_initial_state
        print '(A,T27,A,L1)',    "export_final_state", "= ", export_final_state
        print '(A,T27,A,A)',     "output_dir_name", "= ", "output"
        print '(A,T27,A,A)',     "export_num", "= ", export_num

        call print_domain_summary(domain_info)
        print '(A)', ""

        ! memory info
        print '(A,T42,A,T45,A,T59,A,T62,A)', "memory usage", "|", "per cell [B]", "|", "all cells [GB]"
        print '(A)', "---------------------------------------------------------------------------"
        print '(A,T42,A,T45,I12,T59,A,T62,F14.3)', "dist function buffers", "|", &
            nint(real(dist_function_buffers_bytes, real64) / real(N_CELLS, real64), int64), "|", &
            real(dist_function_buffers_bytes, real64) * gb_per_byte
        print '(A,T42,A,T45,I12,T59,A,T62,F14.3)', "macro field buffers", "|", &
            nint(real(macro_field_buffers_bytes, real64) / real(N_CELLS, real64), int64), "|", &
            real(macro_field_buffers_bytes, real64) * gb_per_byte
        print '(A,T42,A,T45,I12,T59,A,T62,F14.3)', "total", "|", &
            nint(total_bytes_per_cell, int64), "|", real(total_buffer_bytes, real64) * gb_per_byte
        print *
    end subroutine print_run_summary


    subroutine print_launch_timestamp()
        ! temp
        character(len=10) :: current_time

        call date_and_time(time=current_time)
        print '(A)', "[" // current_time(1:2) // ":" // current_time(3:4) // ":" // current_time(5:6) // "] &
            launched -------------------------------------------------------"
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
            finished -------------------------------------------------------"
    end subroutine print_finish_timestamp


    subroutine print_execution_summary( &
        kernel_compute_seconds, halo_exchange_seconds, other_seconds, &
        elapsed_seconds, seconds_per_step, mlups &
        )
        ! inputs
        real(real64), intent(in) :: kernel_compute_seconds
        real(real64), intent(in) :: halo_exchange_seconds
        real(real64), intent(in) :: other_seconds
        real(real64), intent(in) :: elapsed_seconds
        real(real64), intent(in) :: seconds_per_step
        real(real64), intent(in) :: mlups

        print '(A)', ""
        print '(A,T42,A,T46,A,T59,A,T67,A)', "execution time", "|", "total [sec]", "|", "share [%]"
        print '(A)', "---------------------------------------------------------------------------"

        call print_execution_time_row("kernel compute", kernel_compute_seconds, elapsed_seconds)
        call print_execution_time_row("halo exchange", halo_exchange_seconds, elapsed_seconds)
        call print_execution_time_row("other", other_seconds, elapsed_seconds)
        call print_execution_time_row("total", elapsed_seconds, elapsed_seconds)

        print '(A)', ""
        print '(A)', "--- [ perf metrics ] ------------------------------------------------------"
        print '(A,I0,A,I0,A,I0,A)', "sim size [X/Y/N]:      [ ", N_X, " / ", N_Y, " / ", N_STEPS, " ]"
        print '(A,F12.3,A)',        "total time:     ", elapsed_seconds, " sec"
        print '(A,F12.3,A)',        "step time:      ", seconds_per_step * 1000.0_real64, " ms"
        print '(A,F12.3)',          "MLUPS:          ", mlups
    end subroutine print_execution_summary


    subroutine print_execution_time_row( &
        row_name, total_seconds, total_loop_seconds &
        )
        ! inputs
        character(len=*), intent(in) :: row_name
        real(real64), intent(in) :: total_seconds
        real(real64), intent(in) :: total_loop_seconds

        ! temp
        real(real64) :: time_share

        if (total_loop_seconds > 0.0_real64) then
            time_share = 100.0_real64 * total_seconds / total_loop_seconds
        else
            time_share = 0.0_real64
        end if

        print '(A,T42,A,T45,F12.3,T59,A,T64,F10.3,A)', &
            row_name, "|", total_seconds, "|", time_share, " %"
    end subroutine print_execution_time_row


end module reporting
