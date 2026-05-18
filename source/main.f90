program main
    ! imports
    use iso_fortran_env, only: int32, int64, real64, output_unit
    use domain, only: domain_t, initialize_domain, print_domain_summary
    use exchange, only: halo_buffers_t, allocate_halo_buffers, exchange_halos, &
        exchange_halos_direct_put, finish_halos_direct_put
    use export, only: should_export_step, export_selected_data, export_selected_data_distributed, export_metadata
    use hardware_info, only: hardware_info_t, collect_hardware_info, print_hardware_summary
    use initialization, only: initialize_sim_condition, apply_condition_shear_wave_local
    use settings, only: N_X, N_Y, N_STEPS, N_CELLS, N_DIRS, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, &
        USE_UNROLLED_KERNELS, USE_PULL_SHIFT_KERNELS, USE_DIRECT_COARRAY_HALOS, &
        shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t, sim_mode_to_string
    use shear_wave, only: fuzed_unrolled_pull_streaming_collision_local_SW, &
        fuzed_unrolled_pull_streaming_collision_range_SW
    use simulation, only: execute_full_sim_step, swap_distribution_function_buffers
    implicit none

    ! misc
    integer(int32) :: step
    integer(int32) :: active_buffer
    integer(int32) :: next_buffer
    logical :: write_macro_fields
    logical :: use_distributed_shear_wave
    logical :: direct_coarray_halos_enabled
    type(domain_t) :: domain_info
    type(halo_buffers_t) :: halo_buffers
    type(hardware_info_t) :: machine_info

    ! parameter set for shear wave
    type(shear_wave_params_t), parameter :: shear_wave_params = shear_wave_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        u_max = 0.1_FP, &
        n_sin = 3.0_FP &
    )

    ! parameter set for couette flow
    type(couette_flow_params_t), parameter :: couette_flow_params = couette_flow_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        u_wall = 0.1_FP &
    )

    ! parameter set for poiseuille flow
    type(poiseuille_flow_params_t), parameter :: poiseuille_flow_params = poiseuille_flow_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        rho_in = 1.001_FP, &
        rho_out = 0.999_FP &
    )

    ! parameter set for sliding lid
    type(sliding_lid_params_t), parameter :: sliding_lid_params = sliding_lid_params_t( &
        rho_0 = 1.0_FP, &
        omega = 1.5_FP, &
        u_wall = 0.1_FP &
    )

    ! export settings
    logical, parameter :: export_rho = .false.
    logical, parameter :: export_u_x = .false.
    logical, parameter :: export_u_y = .false.
    logical, parameter :: export_u_mag = .false.
    integer(int32), parameter :: export_interval = 10000
    logical, parameter :: export_initial_state = .true.
    logical, parameter :: export_final_state = .true.
    character(len=*), parameter :: output_dir_name = "output"
    character(len=*), parameter :: export_num = "run_000"

    ! progress display settings
    logical, parameter :: interactive_progress = .true.
    integer(int32), parameter :: progress_interval = 1

    ! metrics
    integer(int64) :: clock_start
    integer(int64) :: clock_end
    integer(int64) :: clock_rate
    integer(int64) :: clock_now
    integer(int64) :: clock_section_start
    integer(int64) :: clock_section_end
    integer(int64) :: eta_total_seconds
    integer(int64) :: eta_hours
    integer(int64) :: eta_minutes
    integer(int64) :: eta_seconds_int
    integer(int64) :: bytes_fp
    integer(int64) :: dist_function_buffers_bytes
    integer(int64) :: macro_field_buffers_bytes
    integer(int64) :: total_buffer_bytes
    real(real64) :: elapsed_seconds
    real(real64) :: seconds_per_step
    real(real64) :: avg_millisec_per_step
    real(real64) :: elapsed_now
    real(real64) :: eta_seconds
    real(real64) :: sim_percent
    real(real64) :: mlups
    real(real64) :: gb_per_byte
    real(real64) :: total_bytes_per_cell
    real(real64) :: kernel_compute_seconds
    real(real64) :: halo_exchange_seconds
    real(real64) :: buffer_swap_seconds
    real(real64) :: export_seconds
    real(real64) :: progress_seconds
    real(real64) :: measured_seconds
    real(real64) :: other_seconds
    real(real64) :: execution_time_values(7)[*]
    character(len=10) :: current_time
    integer(int32) :: image_id
    integer(int32) :: timing_image_id

    ! allocate sim data structures (double-buffered distribution functions)
    real(FP), allocatable :: f(:, :, :) ! read-version of distribution functions f(x, y, dir)
    real(FP), allocatable :: f_next(:, :, :) ! write-version version of f(x, y, dir)
    real(FP), allocatable :: f_direct(:, :, :, :)[:] ! double-buffered coarray version of f(x, y, dir)
    real(FP), allocatable :: rho(:,:)
    real(FP), allocatable :: u_x(:,:)
    real(FP), allocatable :: u_y(:,:)

    ! setup domain decomposition
    call initialize_domain(domain_info)

    use_distributed_shear_wave = SIM_MODE == SIM_SHEAR_WAVE
    direct_coarray_halos_enabled = use_distributed_shear_wave .and. USE_DIRECT_COARRAY_HALOS
    active_buffer = 1_int32
    next_buffer = 2_int32

    if (domain_info%n_images > 1 .and. SIM_MODE /= SIM_SHEAR_WAVE) then
        error stop "error: distributed coarray execution is only implemented for shear wave yet"
    end if

    if (use_distributed_shear_wave .and. USE_PULL_SHIFT_KERNELS) then
        error stop "error: distributed pull-shift shear wave is not implemented yet"
    end if

    if (use_distributed_shear_wave .and. .not. USE_UNROLLED_KERNELS) then
        error stop "error: distributed regular shear wave is not implemented yet"
    end if

    if (this_image() == 1) then
        call collect_hardware_info(machine_info)
    end if

    if (direct_coarray_halos_enabled) then
        allocate(f_direct(0:domain_info%n_x_local+1, 0:domain_info%n_y_local+1, N_DIRS, 2)[*])
        allocate(rho(domain_info%n_x_local, domain_info%n_y_local))
        allocate(u_x(domain_info%n_x_local, domain_info%n_y_local))
        allocate(u_y(domain_info%n_x_local, domain_info%n_y_local))
    else if (use_distributed_shear_wave) then
        allocate(f(0:domain_info%n_x_local+1, 0:domain_info%n_y_local+1, N_DIRS))
        allocate(f_next(0:domain_info%n_x_local+1, 0:domain_info%n_y_local+1, N_DIRS))
        allocate(rho(domain_info%n_x_local, domain_info%n_y_local))
        allocate(u_x(domain_info%n_x_local, domain_info%n_y_local))
        allocate(u_y(domain_info%n_x_local, domain_info%n_y_local))
        call allocate_halo_buffers(domain_info, halo_buffers)
    else
        allocate(f(N_X, N_Y, N_DIRS))
        allocate(f_next(N_X, N_Y, N_DIRS))
        allocate(rho(N_X, N_Y))
        allocate(u_x(N_X, N_Y))
        allocate(u_y(N_X, N_Y))
    end if

    ! compute memory metrics for persistent main sim buffers
    bytes_fp = int(storage_size(0.0_FP), int64) / 8_int64
    if (direct_coarray_halos_enabled) then
        dist_function_buffers_bytes = size(f_direct, kind=int64) * bytes_fp * &
            int(domain_info%n_images, int64)
    else
        dist_function_buffers_bytes = (size(f, kind=int64) + size(f_next, kind=int64)) * bytes_fp * &
            int(domain_info%n_images, int64)
    end if
    macro_field_buffers_bytes = (size(rho, kind=int64) + size(u_x, kind=int64) + size(u_y, kind=int64)) * bytes_fp * &
        int(domain_info%n_images, int64)
    total_buffer_bytes = dist_function_buffers_bytes + macro_field_buffers_bytes
    total_bytes_per_cell = real(total_buffer_bytes, real64) / real(N_CELLS, real64)
    gb_per_byte = 1.0e-9_real64

    ! inital condition
    if (direct_coarray_halos_enabled) then
        call apply_condition_shear_wave_local( &
            shear_wave_params%rho_0, shear_wave_params%u_max, shear_wave_params%n_sin, &
            domain_info%n_x_local, domain_info%n_y_local, domain_info%y_global_start, &
            f_direct(:, :, :, active_buffer), rho, u_x, u_y)
    else if (use_distributed_shear_wave) then
        call apply_condition_shear_wave_local( &
            shear_wave_params%rho_0, shear_wave_params%u_max, shear_wave_params%n_sin, &
            domain_info%n_x_local, domain_info%n_y_local, domain_info%y_global_start, f, rho, u_x, u_y)
    else
        call initialize_sim_condition(shear_wave_params, couette_flow_params, poiseuille_flow_params, &
            sliding_lid_params, f, rho, u_x, u_y)
    end if

    ! print sim info
    if (this_image() == 1) then
        print '(A)', ""
        call print_hardware_summary(machine_info)

        print '(A)', ""
        print '(A)', "--- [ simulation parameters ] ---------------------------------------------"
        print '(A,T27,A,A)',     "SIM_MODE", "= ", trim(sim_mode_to_string(SIM_MODE))

        select case (SIM_MODE)
        case (SIM_SHEAR_WAVE)
            print '(A,T27,A,F8.6)', "rho_0", "= ", shear_wave_params%rho_0
            print '(A,T27,A,F8.6)', "omega", "= ", shear_wave_params%omega
            print '(A,T27,A,F8.6)', "u_max", "= ", shear_wave_params%u_max
            print '(A,T27,A,F8.6)', "n_sin", "= ", shear_wave_params%n_sin
        case (SIM_COUETTE_FLOW)
            print '(A,T27,A,F8.6)', "rho_0", "= ", couette_flow_params%rho_0
            print '(A,T27,A,F8.6)', "omega", "= ", couette_flow_params%omega
            print '(A,T27,A,F8.6)', "u_wall", "= ", couette_flow_params%u_wall
        case (SIM_POISEUILLE_FLOW)
            print '(A,T27,A,F8.6)', "rho_0", "= ", poiseuille_flow_params%rho_0
            print '(A,T27,A,F8.6)', "omega", "= ", poiseuille_flow_params%omega
            print '(A,T27,A,F8.6)', "rho_in", "= ", poiseuille_flow_params%rho_in
            print '(A,T27,A,F8.6)', "rho_out", "= ", poiseuille_flow_params%rho_out
        case (SIM_SLIDING_LID)
            print '(A,T27,A,F8.6)', "rho_0", "= ", sliding_lid_params%rho_0
            print '(A,T27,A,F8.6)', "omega", "= ", sliding_lid_params%omega
            print '(A,T27,A,F8.6)', "u_wall", "= ", sliding_lid_params%u_wall
        case default
            error stop "error: unknown sim mode in main print block"
        end select

        ! parameter info
        print '(A)', ""
        print '(A)', "--- [ other parameters ] --------------------------------------------------"
        print '(A,T27,A,I0)',    "N_X_TOTAL", "= ", N_X
        print '(A,T27,A,I0)',    "N_Y_TOTAL", "= ", N_Y
        print '(A,T27,A,I0)',    "N_STEPS", "= ", N_STEPS
        print '(A,T27,A,L1)',    "use_unrolled_kernels", "= ", USE_UNROLLED_KERNELS
        print '(A,T27,A,L1)',    "use_pull_shift_kernels", "= ", USE_PULL_SHIFT_KERNELS
        print '(A,T27,A,L1)',    "use_direct_coarray_halos", "= ", direct_coarray_halos_enabled
        print '(A,T27,A,L1)',    "distributed_coarrays", "= ", .true.
        print '(A,T27,A,L1)',    "export_rho", "= ", export_rho
        print '(A,T27,A,L1)',    "export_u_x", "= ", export_u_x
        print '(A,T27,A,L1)',    "export_u_y", "= ", export_u_y
        print '(A,T27,A,L1)',    "export_u_mag", "= ", export_u_mag
        print '(A,T27,A,I0)',    "export_interval", "= ", export_interval
        print '(A,T27,A,L1)',    "export_initial_state", "= ", export_initial_state
        print '(A,T27,A,L1)',    "export_final_state", "= ", export_final_state
        print '(A,T27,A,A)',     "output_dir_name", "= ", output_dir_name
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
    end if

    ! export metadata
    if (this_image() == 1) then
        call export_metadata(machine_info, shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
            export_rho, export_u_x, export_u_y, export_u_mag, export_interval, &
            output_dir_name, export_num, export_initial_state, export_final_state)
    end if

    ! export initial condition
    if (should_export_step(0_int32, export_interval, &
        export_initial_state, export_final_state)) then
        if (use_distributed_shear_wave) then
            call export_selected_data_distributed(domain_info, export_rho, export_u_x, export_u_y, export_u_mag, &
                output_dir_name, export_num, 0_int32, rho, u_x, u_y)
        else if (this_image() == 1) then
            call export_selected_data(export_rho, export_u_x, export_u_y, export_u_mag, &
                output_dir_name, export_num, 0_int32, rho, u_x, u_y)
        end if
    end if

    ! print sim launch timestamp
    if (this_image() == 1) then
        call date_and_time(time=current_time)
        print '(A)', "[" // current_time(1:2) // ":" // current_time(3:4) // ":" // current_time(5:6) // "] &
            launched -------------------------------------------------------"
    end if

    if (use_distributed_shear_wave) then
        sync all
    end if

    kernel_compute_seconds = 0.0_real64
    halo_exchange_seconds = 0.0_real64
    buffer_swap_seconds = 0.0_real64
    export_seconds = 0.0_real64
    progress_seconds = 0.0_real64

    call system_clock(clock_start, clock_rate)

    ! simulation loop
    do step = 1, N_STEPS

        ! decide if density and velocity fields need to be stored in this step
        write_macro_fields = SIM_MODE == SIM_POISEUILLE_FLOW .or. &
            should_export_step(step, export_interval, export_initial_state, export_final_state) .and. &
            (export_rho .or. export_u_x .or. export_u_y .or. export_u_mag)

        if (direct_coarray_halos_enabled) then
            call system_clock(clock_section_start)
            call exchange_halos_direct_put( &
                domain_info, domain_info%n_x_local, domain_info%n_y_local, active_buffer, f_direct)
            call system_clock(clock_section_end)
            halo_exchange_seconds = halo_exchange_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

            call system_clock(clock_section_start)
            call fuzed_unrolled_pull_streaming_collision_range_SW( &
                2_int32, domain_info%n_x_local - 1_int32, 2_int32, domain_info%n_y_local - 1_int32, &
                domain_info%n_x_local, domain_info%n_y_local, &
                write_macro_fields, shear_wave_params%omega, &
                f_direct(:, :, :, active_buffer), f_direct(:, :, :, next_buffer), rho, u_x, u_y)
            call system_clock(clock_section_end)
            kernel_compute_seconds = kernel_compute_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

            call system_clock(clock_section_start)
            call finish_halos_direct_put(domain_info)
            call system_clock(clock_section_end)
            halo_exchange_seconds = halo_exchange_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

            call system_clock(clock_section_start)
            call fuzed_unrolled_pull_streaming_collision_range_SW( &
                1_int32, domain_info%n_x_local, 1_int32, 1_int32, &
                domain_info%n_x_local, domain_info%n_y_local, &
                write_macro_fields, shear_wave_params%omega, &
                f_direct(:, :, :, active_buffer), f_direct(:, :, :, next_buffer), rho, u_x, u_y)

            if (domain_info%n_y_local > 1) then
                call fuzed_unrolled_pull_streaming_collision_range_SW( &
                    1_int32, domain_info%n_x_local, domain_info%n_y_local, domain_info%n_y_local, &
                    domain_info%n_x_local, domain_info%n_y_local, &
                    write_macro_fields, shear_wave_params%omega, &
                    f_direct(:, :, :, active_buffer), f_direct(:, :, :, next_buffer), rho, u_x, u_y)
            end if

            if (domain_info%n_y_local > 2) then
                call fuzed_unrolled_pull_streaming_collision_range_SW( &
                    1_int32, 1_int32, 2_int32, domain_info%n_y_local - 1_int32, &
                    domain_info%n_x_local, domain_info%n_y_local, &
                    write_macro_fields, shear_wave_params%omega, &
                    f_direct(:, :, :, active_buffer), f_direct(:, :, :, next_buffer), rho, u_x, u_y)

                if (domain_info%n_x_local > 1) then
                    call fuzed_unrolled_pull_streaming_collision_range_SW( &
                        domain_info%n_x_local, domain_info%n_x_local, &
                        2_int32, domain_info%n_y_local - 1_int32, &
                        domain_info%n_x_local, domain_info%n_y_local, &
                        write_macro_fields, shear_wave_params%omega, &
                        f_direct(:, :, :, active_buffer), f_direct(:, :, :, next_buffer), rho, u_x, u_y)
                end if
            end if

            call system_clock(clock_section_end)
            kernel_compute_seconds = kernel_compute_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
        else if (use_distributed_shear_wave) then
            call system_clock(clock_section_start)
            call exchange_halos(domain_info, halo_buffers, domain_info%n_x_local, domain_info%n_y_local, f)
            call system_clock(clock_section_end)
            halo_exchange_seconds = halo_exchange_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

            call system_clock(clock_section_start)
            call fuzed_unrolled_pull_streaming_collision_local_SW( &
                domain_info%n_x_local, domain_info%n_y_local, &
                write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
            call system_clock(clock_section_end)
            kernel_compute_seconds = kernel_compute_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
        else
            call system_clock(clock_section_start)
            call execute_full_sim_step( &
                shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
                write_macro_fields, f, f_next, rho, u_x, u_y)
            call system_clock(clock_section_end)
            kernel_compute_seconds = kernel_compute_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
        end if

        call system_clock(clock_section_start)
        if (direct_coarray_halos_enabled) then
            active_buffer = 3_int32 - active_buffer
            next_buffer = 3_int32 - next_buffer
        else
            call swap_distribution_function_buffers(f, f_next)
        end if
        call system_clock(clock_section_end)
        buffer_swap_seconds = buffer_swap_seconds + &
            real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

        ! export selected field
        if (should_export_step(step, export_interval, &
            export_initial_state, export_final_state)) then
            if (use_distributed_shear_wave) then
                call system_clock(clock_section_start)
                call export_selected_data_distributed(domain_info, export_rho, export_u_x, export_u_y, export_u_mag, &
                    output_dir_name, export_num, step, rho, u_x, u_y)
                call system_clock(clock_section_end)
                export_seconds = export_seconds + &
                    real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
            else if (this_image() == 1) then
                call system_clock(clock_section_start)
                call export_selected_data(export_rho, export_u_x, export_u_y, export_u_mag, &
                    output_dir_name, export_num, step, rho, u_x, u_y)
                call system_clock(clock_section_end)
                export_seconds = export_seconds + &
                    real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
            end if
        end if
        
        ! print sim progress info
        if (this_image() == 1) then
            if (mod(step, progress_interval) == 0 .or. step == N_STEPS) then
                call system_clock(clock_section_start)
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

                call system_clock(clock_section_end)
                progress_seconds = progress_seconds + &
                    real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
            end if
        end if

    end do

    if (use_distributed_shear_wave) then
        sync all
    end if
    
    ! print sim finish timestamp
    if (this_image() == 1 .and. interactive_progress) then
        print *
    end if
    if (this_image() == 1) then
        call date_and_time(time=current_time)
        print '(A)', "[" // current_time(1:2) // ":" // current_time(3:4) // ":" // current_time(5:6) // "] &
            finished -------------------------------------------------------"
    end if

    ! finalize timing and print metrics
    call system_clock(clock_end)
    elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)
    seconds_per_step = elapsed_seconds / real(N_STEPS, real64)
    mlups = real(N_CELLS, real64) * real(N_STEPS, real64) / elapsed_seconds / 1.0e6_real64

    measured_seconds = kernel_compute_seconds + halo_exchange_seconds + &
        buffer_swap_seconds + export_seconds + progress_seconds
    other_seconds = max(0.0_real64, elapsed_seconds - measured_seconds)

    execution_time_values(1) = kernel_compute_seconds
    execution_time_values(2) = halo_exchange_seconds
    execution_time_values(3) = buffer_swap_seconds
    execution_time_values(4) = export_seconds
    execution_time_values(5) = progress_seconds
    execution_time_values(6) = other_seconds
    execution_time_values(7) = elapsed_seconds

    sync all

    if (this_image() == 1) then
        timing_image_id = 1
        do image_id = 2, domain_info%n_images
            if (execution_time_values(7)[image_id] > execution_time_values(7)[timing_image_id]) then
                timing_image_id = image_id
            end if
        end do

        print '(A)', ""
        print '(A,T42,A,T46,A,T59,A,T67,A)', "execution time", "|", "total [sec]", "|", "share [%]"
        print '(A)', "---------------------------------------------------------------------------"

        if (use_distributed_shear_wave) then
            call print_execution_time_row("kernel compute", &
                execution_time_values(1)[timing_image_id], execution_time_values(7)[timing_image_id])
            call print_execution_time_row("halo exchange", &
                execution_time_values(2)[timing_image_id], execution_time_values(7)[timing_image_id])
        else
            call print_execution_time_row("kernel compute", &
                execution_time_values(1)[timing_image_id], execution_time_values(7)[timing_image_id])
        end if

        call print_execution_time_row("buffer swap", &
            execution_time_values(3)[timing_image_id], execution_time_values(7)[timing_image_id])
        call print_execution_time_row("data export", &
            execution_time_values(4)[timing_image_id], execution_time_values(7)[timing_image_id])
        call print_execution_time_row("progress display", &
            execution_time_values(5)[timing_image_id], execution_time_values(7)[timing_image_id])
        call print_execution_time_row("other", &
            execution_time_values(6)[timing_image_id], execution_time_values(7)[timing_image_id])
        call print_execution_time_row("total", &
            execution_time_values(7)[timing_image_id], execution_time_values(7)[timing_image_id])

        print '(A)', ""
        print '(A)', "--- [ perf metrics ] ------------------------------------------------------"
        print '(A,I0,A,I0,A,I0,A)', "sim size [X/Y/N]:      [ ", N_X, " / ", N_Y, " / ", N_STEPS, " ]"
        print '(A,F12.3,A)',        "total time:     ", elapsed_seconds, " sec"
        print '(A,F12.3,A)',        "step time:      ", seconds_per_step * 1000.0_real64, " ms"
        print '(A,F12.3)',          "MLUPS:          ", mlups
    end if

contains

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

end program main
