program main
    ! imports
    use iso_fortran_env, only: int32, int64, real64
    use domain, only: domain_t, initialize_domain
    use exchange, only: halo_buffers_t, allocate_halo_buffers, exchange_halos
    use export, only: should_export_step, export_selected_data, export_selected_data_distributed, export_metadata
    use hardware_info, only: hardware_info_t, collect_hardware_info
    use initialization, only: initialize_sim_condition, apply_condition_shear_wave_local, &
        apply_condition_couette_flow_local, apply_condition_poiseuille_flow_local, apply_condition_sliding_lid_local
    use settings, only: N_X, N_Y, N_STEPS, N_CELLS, N_DIRS, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, &
        USE_PULL_SHIFT_KERNELS, EXPORT_RHO, EXPORT_U_X, EXPORT_U_Y, EXPORT_U_MAG, EXPORT_INTERVAL, &
        EXPORT_INITIAL_STATE, EXPORT_FINAL_STATE, EXPORT_NUM, INTERACTIVE_PROGRESS, &
        PROGRESS_INTERVAL, RHO_0, U_MAX, N_SIN
    use reporting, only: print_run_summary, print_launch_timestamp, print_progress_status, print_finish_timestamp, &
        print_execution_summary
    use simulation, only: execute_local_sim_step, swap_distribution_function_buffers
    implicit none

    ! misc
    integer(int32) :: step
    logical :: write_macro_fields
    logical :: use_distributed_domain
    type(domain_t) :: domain_info
    type(halo_buffers_t) :: halo_buffers
    type(hardware_info_t) :: machine_info

    ! metrics
    integer(int64) :: clock_start
    integer(int64) :: clock_end
    integer(int64) :: clock_rate
    integer(int64) :: clock_section_start
    integer(int64) :: clock_section_end
    integer(int64) :: bytes_fp
    integer(int64) :: dist_function_buffers_bytes
    integer(int64) :: macro_field_buffers_bytes
    integer(int64) :: total_buffer_bytes
    real(real64) :: elapsed_seconds
    real(real64) :: seconds_per_step
    real(real64) :: mlups
    real(real64) :: total_bytes_per_cell
    real(real64) :: kernel_compute_seconds
    real(real64) :: halo_exchange_seconds
    real(real64) :: buffer_swap_seconds
    real(real64) :: export_seconds
    real(real64) :: progress_seconds
    real(real64) :: measured_seconds
    real(real64) :: other_seconds
    real(real64) :: execution_time_values(7)[*]
    integer(int32) :: image_id
    integer(int32) :: timing_image_id

    ! allocate sim data structures (double-buffered distribution functions)
    real(FP), allocatable :: f(:, :, :) ! read-version of distribution functions f(x, y, dir)
    real(FP), allocatable :: f_next(:, :, :) ! write-version version of f(x, y, dir)
    real(FP), allocatable :: rho(:,:)
    real(FP), allocatable :: u_x(:,:)
    real(FP), allocatable :: u_y(:,:)

    ! setup domain decomposition
    call initialize_domain(domain_info)

    use_distributed_domain = SIM_MODE == SIM_SHEAR_WAVE .or. &
        SIM_MODE == SIM_COUETTE_FLOW .or. &
        SIM_MODE == SIM_POISEUILLE_FLOW .or. &
        SIM_MODE == SIM_SLIDING_LID

    if (domain_info%n_images > 1 .and. .not. use_distributed_domain) then
        error stop "error: distributed coarray execution is only implemented for known simulation modes yet"
    end if

    if ((SIM_MODE == SIM_COUETTE_FLOW .or. SIM_MODE == SIM_POISEUILLE_FLOW .or. &
        SIM_MODE == SIM_SLIDING_LID) .and. USE_PULL_SHIFT_KERNELS) then
        error stop "error: distributed pull-shift is not implemented for this simulation mode yet"
    end if

    if (this_image() == 1) then
        call collect_hardware_info(machine_info)
    end if

    if (use_distributed_domain) then
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
    dist_function_buffers_bytes = (size(f, kind=int64) + size(f_next, kind=int64)) * bytes_fp * &
        int(domain_info%n_images, int64)
    macro_field_buffers_bytes = (size(rho, kind=int64) + size(u_x, kind=int64) + size(u_y, kind=int64)) * bytes_fp * &
        int(domain_info%n_images, int64)
    total_buffer_bytes = dist_function_buffers_bytes + macro_field_buffers_bytes
    total_bytes_per_cell = real(total_buffer_bytes, real64) / real(N_CELLS, real64)

    ! inital condition
    if (SIM_MODE == SIM_SHEAR_WAVE) then
        call apply_condition_shear_wave_local( &
            RHO_0, U_MAX, N_SIN, &
            domain_info%n_x_local, domain_info%n_y_local, domain_info%y_global_start, f, rho, u_x, u_y)
    else if (SIM_MODE == SIM_COUETTE_FLOW) then
        call apply_condition_couette_flow_local( &
            RHO_0, domain_info%n_x_local, domain_info%n_y_local, f, rho, u_x, u_y)
    else if (SIM_MODE == SIM_POISEUILLE_FLOW) then
        call apply_condition_poiseuille_flow_local( &
            RHO_0, domain_info%n_x_local, domain_info%n_y_local, f, rho, u_x, u_y)
    else if (SIM_MODE == SIM_SLIDING_LID) then
        call apply_condition_sliding_lid_local( &
            RHO_0, domain_info%n_x_local, domain_info%n_y_local, f, rho, u_x, u_y)
    else
        call initialize_sim_condition(f, rho, u_x, u_y)
    end if

    if (use_distributed_domain .and. SIM_MODE == SIM_POISEUILLE_FLOW) then
        halo_buffers%send_macro_left(:, 1) = RHO_0
        halo_buffers%send_macro_left(:, 2) = 0.0_FP
        halo_buffers%send_macro_left(:, 3) = 0.0_FP

        halo_buffers%send_macro_right(:, 1) = RHO_0
        halo_buffers%send_macro_right(:, 2) = 0.0_FP
        halo_buffers%send_macro_right(:, 3) = 0.0_FP

        halo_buffers%recv_macro_left(:, 1) = RHO_0
        halo_buffers%recv_macro_left(:, 2) = 0.0_FP
        halo_buffers%recv_macro_left(:, 3) = 0.0_FP

        halo_buffers%recv_macro_right(:, 1) = RHO_0
        halo_buffers%recv_macro_right(:, 2) = 0.0_FP
        halo_buffers%recv_macro_right(:, 3) = 0.0_FP
    end if

    ! print sim info
    if (this_image() == 1) then
        call print_run_summary( &
            machine_info, domain_info, SIM_MODE, &
            EXPORT_RHO, EXPORT_U_X, EXPORT_U_Y, EXPORT_U_MAG, EXPORT_INTERVAL, EXPORT_INITIAL_STATE, &
            EXPORT_FINAL_STATE, EXPORT_NUM, dist_function_buffers_bytes, macro_field_buffers_bytes, &
            total_buffer_bytes, total_bytes_per_cell)
    end if

    ! export metadata
    if (this_image() == 1) then
        call export_metadata(machine_info, domain_info, SIM_MODE, &
            EXPORT_RHO, EXPORT_U_X, EXPORT_U_Y, EXPORT_U_MAG, EXPORT_INTERVAL, &
            EXPORT_NUM, EXPORT_INITIAL_STATE, EXPORT_FINAL_STATE, dist_function_buffers_bytes, &
            macro_field_buffers_bytes, total_buffer_bytes, total_bytes_per_cell)
    end if

    ! export initial condition
    if (should_export_step(0_int32, EXPORT_INTERVAL, &
        EXPORT_INITIAL_STATE, EXPORT_FINAL_STATE)) then
        if (use_distributed_domain) then
            call export_selected_data_distributed(domain_info, EXPORT_RHO, EXPORT_U_X, EXPORT_U_Y, EXPORT_U_MAG, &
                EXPORT_NUM, 0_int32, rho, u_x, u_y)
        else if (this_image() == 1) then
            call export_selected_data(EXPORT_RHO, EXPORT_U_X, EXPORT_U_Y, EXPORT_U_MAG, &
                EXPORT_NUM, 0_int32, rho, u_x, u_y)
        end if
    end if

    ! print sim launch timestamp
    if (this_image() == 1) then
        call print_launch_timestamp()
    end if

    if (use_distributed_domain) then
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
        write_macro_fields = should_export_step(step, EXPORT_INTERVAL, EXPORT_INITIAL_STATE, EXPORT_FINAL_STATE) .and. &
            (EXPORT_RHO .or. EXPORT_U_X .or. EXPORT_U_Y .or. EXPORT_U_MAG)

        if (use_distributed_domain) then
            call system_clock(clock_section_start)
            call exchange_halos( &
                domain_info, halo_buffers, domain_info%n_x_local, domain_info%n_y_local, f, &
                SIM_MODE == SIM_POISEUILLE_FLOW)
            call system_clock(clock_section_end)
            halo_exchange_seconds = halo_exchange_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

            call system_clock(clock_section_start)
            call execute_local_sim_step( &
                domain_info, halo_buffers, domain_info%n_x_local, domain_info%n_y_local, &
                write_macro_fields, f, f_next, rho, u_x, u_y)
            call system_clock(clock_section_end)
            kernel_compute_seconds = kernel_compute_seconds + &
                real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
        else
            error stop "error: selected sim mode is not implemented for distributed local execution"
        end if

        call system_clock(clock_section_start)
        call swap_distribution_function_buffers(f, f_next)
        call system_clock(clock_section_end)
        buffer_swap_seconds = buffer_swap_seconds + &
            real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

        ! export selected field
        if (should_export_step(step, EXPORT_INTERVAL, &
            EXPORT_INITIAL_STATE, EXPORT_FINAL_STATE)) then
            if (use_distributed_domain) then
                call system_clock(clock_section_start)
                call export_selected_data_distributed(domain_info, EXPORT_RHO, EXPORT_U_X, EXPORT_U_Y, EXPORT_U_MAG, &
                    EXPORT_NUM, step, rho, u_x, u_y)
                call system_clock(clock_section_end)
                export_seconds = export_seconds + &
                    real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
            else if (this_image() == 1) then
                call system_clock(clock_section_start)
                call export_selected_data(EXPORT_RHO, EXPORT_U_X, EXPORT_U_Y, EXPORT_U_MAG, &
                    EXPORT_NUM, step, rho, u_x, u_y)
                call system_clock(clock_section_end)
                export_seconds = export_seconds + &
                    real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
            end if
        end if

        ! print sim progress info
        if (this_image() == 1) then
            if (mod(step, PROGRESS_INTERVAL) == 0 .or. step == N_STEPS) then
                call system_clock(clock_section_start)
                call print_progress_status(step, clock_start, clock_rate, INTERACTIVE_PROGRESS)
                call system_clock(clock_section_end)
                progress_seconds = progress_seconds + &
                    real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
            end if
        end if

    end do

    if (use_distributed_domain) then
        sync all
    end if

    ! print sim finish timestamp
    if (this_image() == 1) then
        call print_finish_timestamp(INTERACTIVE_PROGRESS)
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

        call print_execution_summary( &
            use_distributed_domain, &
            execution_time_values(1)[timing_image_id], execution_time_values(2)[timing_image_id], &
            execution_time_values(3)[timing_image_id], execution_time_values(4)[timing_image_id], &
            execution_time_values(5)[timing_image_id], execution_time_values(6)[timing_image_id], &
            execution_time_values(7)[timing_image_id], seconds_per_step, mlups)
    end if


end program main
