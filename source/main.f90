program main
    ! imports
    use iso_fortran_env, only: int32, int64, real64
    use domain, only: domain_t, initialize_domain
    use exchange, only: halo_buffers_t, exchange_plan_t, exchange_timing_t, &
        BUF_MACRO_LEFT, BUF_MACRO_RIGHT, allocate_halo_buffers, build_exchange_plan, exchange_halos
    use export, only: should_export_step, export_selected_data_distributed, export_metadata
    use hardware_info, only: hardware_info_t, collect_hardware_info
    use initialization, only: apply_condition_shear_wave_local, &
        apply_condition_couette_flow_local, apply_condition_poiseuille_flow_local, apply_condition_sliding_lid_local
    use settings, only: N_STEPS, N_CELLS, N_DIRS, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, &
        EXPORT_MACROS, EXPORT_ENDPOINT_STATES, EXPORT_INTERVAL, &
        EXPORT_NUM, INTERACTIVE_PROGRESS, &
        PROGRESS_INTERVAL, RHO_0, U_MAX, N_SIN, RHO_IN, RHO_OUT
    use reporting, only: print_run_summary, print_launch_timestamp, print_progress_status, print_finish_timestamp, &
        print_execution_summary
    use simulation, only: execute_local_sim_step
    implicit none

    ! misc
    integer(int32) :: step
    logical :: write_macro_fields
    type(domain_t) :: domain_info
    type(halo_buffers_t) :: halo_buffers
    type(exchange_plan_t) :: exchange_plan
    type(exchange_timing_t) :: exchange_timing
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
    real(real64) :: halo_sync_seconds
    real(real64) :: halo_transfer_seconds
    real(real64) :: macro_exchange_seconds
    real(real64) :: measured_seconds
    real(real64) :: other_seconds
    real(real64) :: execution_time_values(7)[*]
    integer(int32) :: image_id
    integer(int32) :: timing_image_id

    ! allocate sim data structures
    real(FP), allocatable :: f_a(:, :, :)[:] ! distribution function buffer A
    real(FP), allocatable :: f_b(:, :, :)[:] ! distribution function buffer B
    logical :: read_from_a
    real(FP), allocatable :: rho(:,:)
    real(FP), allocatable :: u_x(:,:)
    real(FP), allocatable :: u_y(:,:)

    ! setup domain decomposition
    call initialize_domain(domain_info)
    call build_exchange_plan(domain_info, SIM_MODE, exchange_plan)

    if (this_image() == 1) then
        call collect_hardware_info(machine_info)
    end if

    allocate(f_a(0:domain_info%n_x+1, 0:domain_info%n_y+1, N_DIRS)[*])
    allocate(f_b(0:domain_info%n_x+1, 0:domain_info%n_y+1, N_DIRS)[*])
    read_from_a = .true.
    allocate(rho(domain_info%n_x, domain_info%n_y))
    allocate(u_x(domain_info%n_x, domain_info%n_y))
    allocate(u_y(domain_info%n_x, domain_info%n_y))
    call allocate_halo_buffers(domain_info, halo_buffers)

    ! compute memory metrics for persistent main sim buffers
    bytes_fp = int(storage_size(0.0_FP), int64) / 8_int64
    dist_function_buffers_bytes = (size(f_a, kind=int64) + size(f_b, kind=int64)) * bytes_fp * &
        int(domain_info%n_images, int64)
    macro_field_buffers_bytes = (size(rho, kind=int64) + size(u_x, kind=int64) + size(u_y, kind=int64)) * bytes_fp * &
        int(domain_info%n_images, int64)
    total_buffer_bytes = dist_function_buffers_bytes + macro_field_buffers_bytes
    total_bytes_per_cell = real(total_buffer_bytes, real64) / real(N_CELLS, real64)

    ! initial condition
    if (SIM_MODE == SIM_SHEAR_WAVE) then
        call apply_condition_shear_wave_local( &
            RHO_0, U_MAX, N_SIN, &
            domain_info%n_x, domain_info%n_y, domain_info%y_global_start, f_a, rho, u_x, u_y)
    else if (SIM_MODE == SIM_COUETTE_FLOW) then
        call apply_condition_couette_flow_local( &
            RHO_0, domain_info%n_x, domain_info%n_y, f_a, rho, u_x, u_y)
    else if (SIM_MODE == SIM_POISEUILLE_FLOW) then
        call apply_condition_poiseuille_flow_local( &
            RHO_0, domain_info%n_x, domain_info%n_y, f_a, rho, u_x, u_y)
    else if (SIM_MODE == SIM_SLIDING_LID) then
        call apply_condition_sliding_lid_local( &
            RHO_0, domain_info%n_x, domain_info%n_y, f_a, rho, u_x, u_y)
    else
        error stop "error: unknown sim mode in main initial condition"
    end if

    if (SIM_MODE == SIM_POISEUILLE_FLOW) then
        halo_buffers%window(:, 1, BUF_MACRO_LEFT) = RHO_IN
        halo_buffers%window(:, 2, BUF_MACRO_LEFT) = 0.0_FP
        halo_buffers%window(:, 3, BUF_MACRO_LEFT) = 0.0_FP

        halo_buffers%window(:, 1, BUF_MACRO_RIGHT) = RHO_OUT
        halo_buffers%window(:, 2, BUF_MACRO_RIGHT) = 0.0_FP
        halo_buffers%window(:, 3, BUF_MACRO_RIGHT) = 0.0_FP

        halo_buffers%recv_macro_left(:, 1) = RHO_IN
        halo_buffers%recv_macro_left(:, 2) = 0.0_FP
        halo_buffers%recv_macro_left(:, 3) = 0.0_FP

        halo_buffers%recv_macro_right(:, 1) = RHO_OUT
        halo_buffers%recv_macro_right(:, 2) = 0.0_FP
        halo_buffers%recv_macro_right(:, 3) = 0.0_FP
    end if

    ! print sim info
    if (this_image() == 1) then
        call print_run_summary( &
            machine_info, domain_info, SIM_MODE, &
            EXPORT_MACROS, EXPORT_ENDPOINT_STATES, EXPORT_INTERVAL, &
            EXPORT_NUM, dist_function_buffers_bytes, macro_field_buffers_bytes, &
            total_buffer_bytes, total_bytes_per_cell)
    end if

    ! export metadata
    if (this_image() == 1) then
        call export_metadata(machine_info, domain_info, SIM_MODE, &
            EXPORT_MACROS, EXPORT_ENDPOINT_STATES, EXPORT_INTERVAL, &
            EXPORT_NUM, dist_function_buffers_bytes, &
            macro_field_buffers_bytes, total_buffer_bytes, total_bytes_per_cell)
    end if

    ! export initial condition
    if (EXPORT_MACROS .and. should_export_step(0_int32, EXPORT_ENDPOINT_STATES, &
        EXPORT_INTERVAL)) then
        call export_selected_data_distributed(domain_info, EXPORT_NUM, 0_int32, rho, u_x, u_y)
    end if

    ! print sim launch timestamp
    if (this_image() == 1) then
        call print_launch_timestamp()
    end if

    sync all

    kernel_compute_seconds = 0.0_real64
    halo_exchange_seconds = 0.0_real64
    halo_sync_seconds = 0.0_real64
    halo_transfer_seconds = 0.0_real64
    macro_exchange_seconds = 0.0_real64

    call system_clock(clock_start, clock_rate)

    ! simulation loop
    do step = 1, N_STEPS

        ! only store density and velocity if required in this step
        write_macro_fields = should_export_step(step, EXPORT_ENDPOINT_STATES, EXPORT_INTERVAL) .and. &
            EXPORT_MACROS

        call system_clock(clock_section_start)
        if (read_from_a) then
            call exchange_halos( &
                domain_info, halo_buffers, domain_info%n_x, domain_info%n_y, f_a, &
                exchange_plan, clock_rate, exchange_timing)
        else
            call exchange_halos( &
                domain_info, halo_buffers, domain_info%n_x, domain_info%n_y, f_b, &
                exchange_plan, clock_rate, exchange_timing)
        end if
        call system_clock(clock_section_end)
        halo_exchange_seconds = halo_exchange_seconds + &
            real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)
        halo_sync_seconds = halo_sync_seconds + exchange_timing%halo_sync_seconds
        halo_transfer_seconds = halo_transfer_seconds + exchange_timing%halo_transfer_seconds
        macro_exchange_seconds = macro_exchange_seconds + exchange_timing%macro_exchange_seconds

        call system_clock(clock_section_start)
        if (read_from_a) then
            call execute_local_sim_step( &
                domain_info, halo_buffers, domain_info%n_x, domain_info%n_y, &
                write_macro_fields, f_a, f_b, rho, u_x, u_y)
        else
            call execute_local_sim_step( &
                domain_info, halo_buffers, domain_info%n_x, domain_info%n_y, &
                write_macro_fields, f_b, f_a, rho, u_x, u_y)
        end if
        call system_clock(clock_section_end)
        kernel_compute_seconds = kernel_compute_seconds + &
            real(clock_section_end - clock_section_start, real64) / real(clock_rate, real64)

        read_from_a = .not. read_from_a

        ! export selected field
        if (EXPORT_MACROS .and. should_export_step(step, EXPORT_ENDPOINT_STATES, &
            EXPORT_INTERVAL)) then
            call export_selected_data_distributed(domain_info, EXPORT_NUM, step, rho, u_x, u_y)
        end if

        ! print sim progress info
        if (this_image() == 1) then
            if (mod(step, PROGRESS_INTERVAL) == 0 .or. step == N_STEPS) then
                call print_progress_status(step, clock_start, clock_rate, INTERACTIVE_PROGRESS)
            end if
        end if

    end do

    sync all

    ! print sim finish timestamp
    if (this_image() == 1) then
        call print_finish_timestamp(INTERACTIVE_PROGRESS)
    end if

    ! finalize timing and print metrics
    call system_clock(clock_end)
    elapsed_seconds = real(clock_end - clock_start, real64) / real(clock_rate, real64)

    measured_seconds = kernel_compute_seconds + halo_exchange_seconds
    other_seconds = max(0.0_real64, elapsed_seconds - measured_seconds)

    execution_time_values(1) = kernel_compute_seconds
    execution_time_values(2) = halo_exchange_seconds
    execution_time_values(3) = halo_sync_seconds
    execution_time_values(4) = halo_transfer_seconds
    execution_time_values(5) = macro_exchange_seconds
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

        elapsed_seconds = execution_time_values(7)[timing_image_id]
        seconds_per_step = elapsed_seconds / real(N_STEPS, real64)
        mlups = real(N_CELLS, real64) * real(N_STEPS, real64) / elapsed_seconds / 1.0e6_real64

        call print_execution_summary( &
            execution_time_values(1)[timing_image_id], execution_time_values(2)[timing_image_id], &
            execution_time_values(3)[timing_image_id], execution_time_values(4)[timing_image_id], &
            execution_time_values(5)[timing_image_id], execution_time_values(6)[timing_image_id], &
            execution_time_values(7)[timing_image_id], &
            seconds_per_step, mlups)
    end if


end program main
