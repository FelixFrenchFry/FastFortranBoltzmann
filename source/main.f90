program main
    ! imports
    use iso_fortran_env, only: int32, int64, real32, real64, output_unit
    use export, only: should_export_step, export_selected_data, export_metadata
    use initialization, only: initialize_sim_condition
    use settings, only: N_X, N_Y, N_STEPS, N_CELLS, N_DIRS, SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, &
        sim_mode, shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t, sim_mode_to_string
    use simulation, only: execute_full_sim_step, swap_distribution_function_buffers
    implicit none

    ! misc
    integer(int32) :: step

    ! D2Q9 lattice velocities and weights
    integer(int32), parameter :: c_x(N_DIRS) = [ 0,  1,  0, -1,  0,  1, -1, -1,  1 ]
    integer(int32), parameter :: c_y(N_DIRS) = [ 0,  0,  1,  0, -1,  1,  1, -1, -1 ]
    real(real32), parameter :: c_x_fp(N_DIRS) = real(c_x, real32) ! fp-version for compute
    real(real32), parameter :: c_y_fp(N_DIRS) = real(c_y, real32) ! fp-version for compute
    ! ---------
    ! | 7 3 6 |
    ! | 4 1 2 |
    ! | 8 5 9 |
    ! ---------
    ! 1: ( 0,  0) = rest
    ! 2: ( 1,  0) = east
    ! 3: ( 0,  1) = north
    ! 4: (-1,  0) = west
    ! 5: ( 0, -1) = south
    ! 6: ( 1,  1) = north-east
    ! 7: (-1,  1) = north-west
    ! 8: (-1, -1) = south-west
    ! 9: ( 1, -1) = south-east
    real(real32), parameter :: w(N_DIRS) = [ &
        4.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/36.0_real32, &
        1.0_real32/36.0_real32, &
        1.0_real32/36.0_real32, &
        1.0_real32/36.0_real32]

    ! parameter set for shear wave
    type(shear_wave_params_t), parameter :: shear_wave_params = shear_wave_params_t( &
        rho_0 = 1.0_real32, &
        omega = 1.5_real32, &
        u_max = 0.1_real32, &
        n_sin = 2.0_real32 &
    )

    ! parameter set for couette flow
    type(couette_flow_params_t), parameter :: couette_flow_params = couette_flow_params_t( &
        rho_0 = 1.0_real32, &
        omega = 1.5_real32, &
        u_wall = 0.1_real32 &
    )

    ! parameter set for poiseuille flow
    type(poiseuille_flow_params_t), parameter :: poiseuille_flow_params = poiseuille_flow_params_t( &
        rho_0 = 1.0_real32, &
        omega = 1.5_real32, &
        rho_in = 1.001_real32, &
        rho_out = 0.999_real32 &
    )

    ! parameter set for sliding lid
    type(sliding_lid_params_t), parameter :: sliding_lid_params = sliding_lid_params_t( &
        rho_0 = 1.0_real32, &
        omega = 1.5_real32, &
        u_wall = 0.1_real32 &
    )

    ! general settings
    logical, parameter :: write_rho = .true.
    logical, parameter :: write_u_x = .true.
    logical, parameter :: write_u_y = .true.

    ! export settings
    logical, parameter :: export_rho = .true.
    logical, parameter :: export_u_x = .true.
    logical, parameter :: export_u_y = .true.
    logical, parameter :: export_u_mag = .true.
    integer(int32), parameter :: export_interval = 10000
    logical, parameter :: export_initial_state = .true.
    logical, parameter :: export_final_state = .true.
    character(len=*), parameter :: output_dir_name = "output"
    character(len=*), parameter :: export_num = "run_001"

    ! progress display settings
    logical, parameter :: interactive_progress = .true.
    integer(int32), parameter :: progress_interval = 1

    ! metrics
    integer(int64) :: clock_start
    integer(int64) :: clock_end
    integer(int64) :: clock_rate
    integer(int64) :: clock_now
    integer(int64) :: eta_total_seconds
    integer(int64) :: eta_hours
    integer(int64) :: eta_minutes
    integer(int64) :: eta_seconds_int
    real(real64) :: elapsed_seconds
    real(real64) :: seconds_per_step
    real(real64) :: avg_millisec_per_step
    real(real64) :: elapsed_now
    real(real64) :: eta_seconds
    real(real64) :: sim_percent
    real(real64) :: mlups
    character(len=10) :: current_time

    ! allocate sim data structures (double-buffered distribution functions)
    real(real32), allocatable :: f(:, :, :) ! read-version of distribution functions f(dir, x, y)
    real(real32), allocatable :: f_next(:, :, :) ! write-version version of f(dir, x, y)
    real(real32), allocatable :: rho(:,:)
    real(real32), allocatable :: u_x(:,:)
    real(real32), allocatable :: u_y(:,:)
    allocate(f(N_DIRS, N_X, N_Y))
    allocate(f_next(N_DIRS, N_X, N_Y))
    allocate(rho(N_X, N_Y))
    allocate(u_x(N_X, N_Y))
    allocate(u_y(N_X, N_Y))

    ! inital condition
    call initialize_sim_condition(sim_mode, shear_wave_params, couette_flow_params, poiseuille_flow_params, &
        sliding_lid_params, c_x_fp, c_y_fp, w, f, rho, u_x, u_y)

    ! print sim info
    if (this_image() == 1) then
        print '(A)', ""
        print '(A)', "--- [ simulation parameters ] ---------------------------------------------"
        print '(A,A)',     "sim_mode             = ", trim(sim_mode_to_string(sim_mode))

        select case (sim_mode)
        case (SIM_SHEAR_WAVE)
            print '(A,F8.6)', "rho_0                = ", shear_wave_params%rho_0
            print '(A,F8.6)', "omega                = ", shear_wave_params%omega
            print '(A,F8.6)', "u_max                = ", shear_wave_params%u_max
            print '(A,F8.6)', "n_sin                = ", shear_wave_params%n_sin
        case (SIM_COUETTE_FLOW)
            print '(A,F8.6)', "rho_0                = ", couette_flow_params%rho_0
            print '(A,F8.6)', "omega                = ", couette_flow_params%omega
            print '(A,F8.6)', "u_wall               = ", couette_flow_params%u_wall
        case (SIM_POISEUILLE_FLOW)
            print '(A,F8.6)', "rho_0                = ", poiseuille_flow_params%rho_0
            print '(A,F8.6)', "omega                = ", poiseuille_flow_params%omega
            print '(A,F8.6)', "rho_in               = ", poiseuille_flow_params%rho_in
            print '(A,F8.6)', "rho_out              = ", poiseuille_flow_params%rho_out
        case (SIM_SLIDING_LID)
            print '(A,F8.6)', "rho_0                = ", sliding_lid_params%rho_0
            print '(A,F8.6)', "omega                = ", sliding_lid_params%omega
            print '(A,F8.6)', "u_wall               = ", sliding_lid_params%u_wall
        case default
            error stop "error: unknown sim mode in main print block"
        end select

        print '(A)', ""
        print '(A)', "--- [ other parameters ] --------------------------------------------------"
        print '(A,I0)',    "N_X_TOTAL            = ", N_X
        print '(A,I0)',    "N_Y_TOTAL            = ", N_Y
        print '(A,I0)',    "N_STEPS              = ", N_STEPS
        print '(A,L1)',    "write_rho            = ", write_rho
        print '(A,L1)',    "write_u_x            = ", write_u_x
        print '(A,L1)',    "write_u_y            = ", write_u_y
        print '(A,L1)',    "export_rho           = ", export_rho
        print '(A,L1)',    "export_u_x           = ", export_u_x
        print '(A,L1)',    "export_u_y           = ", export_u_y
        print '(A,L1)',    "export_u_mag         = ", export_u_mag
        print '(A,I0)',    "export_interval      = ", export_interval
        print '(A,L1)',    "export_initial_state = ", export_initial_state
        print '(A,L1)',    "export_final_state   = ", export_final_state
        print '(A,A)',     "output_dir_name      = ", output_dir_name
        print '(A,A)',     "export_num           = ", export_num
        print *
    end if

    ! export metadata
    if (this_image() == 1) then
        call export_metadata(sim_mode, shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
            export_rho, export_u_x, export_u_y, export_u_mag, export_interval, &
            output_dir_name, export_num, export_initial_state, export_final_state)
    end if

    ! export initial condition
    if (this_image() == 1) then
        if (should_export_step(0_int32, export_interval, &
            export_initial_state, export_final_state)) then
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

    call system_clock(clock_start, clock_rate)

    ! simulation loop
    do step = 1, N_STEPS

        call execute_full_sim_step( &
            sim_mode, shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
            c_x, c_y, c_x_fp, c_y_fp, w, f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y)

        call swap_distribution_function_buffers(f, f_next)

        ! export selected field
        if (this_image() == 1) then
            if (should_export_step(step, export_interval, &
                export_initial_state, export_final_state)) then
                call export_selected_data(export_rho, export_u_x, export_u_y, export_u_mag, &
                    output_dir_name, export_num, step, rho, u_x, u_y)
            end if
        end if
        
        ! print sim progress info
        if (this_image() == 1) then
            if (mod(step, progress_interval) == 0 .or. step == N_STEPS) then
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
            end if
        end if

    end do
    
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

    if (this_image() == 1) then
        print '(A)', ""
        print '(A)', "--- [ test stats ] --------------------------------------------------------"
        print '(A,2(F15.8,1X))', "rho min/max:      ", minval(rho), maxval(rho)
        print '(A,2(F15.8,1X))', "u_x min/max:      ", minval(u_x), maxval(u_x)
        print '(A,2(F15.8,1X))', "u_y min/max:      ", minval(u_y), maxval(u_y)

        print '(A)', ""
        print '(A)', "--- [ perf metrics ] ------------------------------------------------------"
        print '(A,I0,A,I0,A,I0,A)', "sim size [X/Y/N]:      [ ", N_X, " / ", N_Y, " / ", N_STEPS, " ]"
        print '(A,F12.3,A)',        "total time:     ", elapsed_seconds, " sec"
        print '(A,F12.3,A)',        "step time:      ", seconds_per_step * 1000.0_real64, " ms"
        print '(A,F12.3)',          "MLUPS:          ", mlups
    end if

end program main
