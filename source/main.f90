program main
    ! imports
    use iso_fortran_env, only: int32, int64, real32, real64, output_unit
    use export, only: EXPORT_NONE, EXPORT_DENSITY, EXPORT_VELOCITY_X, EXPORT_VELOCITY_Y, &
        EXPORT_VELOCITY_MAG, export_mode_to_string, should_export_step, export_selected_data
    use initialization, only: apply_condition_shear_wave_decay
    use simulation, only: fuzed_pull_streaming_collision_shear_wave_decay, swap_distribution_function_buffers
    implicit none

    ! misc
    real(real32), parameter :: pi = 3.1415927410125732421875_real32
    integer(int32) :: step

    ! simulation size and duration
    integer(int32), parameter :: N_X = 600
    integer(int32), parameter :: N_Y = 400
    integer(int32), parameter :: N_STEPS = 500
    integer(int64), parameter :: N_CELLS = int(N_X, int64) * int(N_Y, int64)

    ! D2Q9 lattice velocities and weights
    integer(int32), parameter :: N_DIRS = 9
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

    ! general params
    real(real32), parameter :: rho_0 = 1.0_real32 ! rest density
    real(real32), parameter :: omega = 1.5_real32 ! relaxation factor

    ! specific params for shear wave decay
    real(real32), parameter :: u_max = 0.1_real32 ! initial velocity
    real(real32), parameter :: n_sin = 2.0_real32 ! num sin periods
    real(real32), parameter :: k = (2.0_real32 * pi * n_sin) / real(N_Y, real32) ! wave number

    ! general settings
    logical, parameter :: write_rho = .true.
    logical, parameter :: write_u_x = .true.
    logical, parameter :: write_u_y = .true.

    ! export settings
    integer(int32), parameter :: export_mode = EXPORT_VELOCITY_MAG
    integer(int32), parameter :: export_interval = 100
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
    call apply_condition_shear_wave_decay(N_X, N_Y, N_DIRS, c_x_fp, c_y_fp, w, rho_0, u_max, k, f, rho, u_x, u_y)

    ! print sim info
    if (this_image() == 1) then
        print '(A)', "--- [ simulation parameters ] ---------------------------------------------"
        print '(A,I0)',    "N_X_TOTAL            = ", N_X
        print '(A,I0)',    "N_Y_TOTAL            = ", N_Y
        print '(A,I0)',    "N_STEPS              = ", N_STEPS
        print '(A,F5.3)',  "omega                = ", omega
        print '(A,F5.3)',  "rho_0                = ", rho_0
        print '(A,F5.3)',  "u_max                = ", u_max
        print '(A,F5.3)',  "n_sin                = ", n_sin
        print '(A,L1)',    "write_rho            = ", write_rho
        print '(A,L1)',    "write_u_x            = ", write_u_x
        print '(A,L1)',    "write_u_y            = ", write_u_y
        print '(A,A)',     "export_mode          = ", trim(export_mode_to_string(export_mode))
        print '(A,I0)',    "export_interval      = ", export_interval
        print '(A,L1)',    "export_initial_state = ", export_initial_state
        print '(A,L1)',    "export_final_state   = ", export_final_state
        print '(A,A)',     "output_dir_name      = ", output_dir_name
        print '(A,A)',     "export_num           = ", export_num
        print '(A,L1)',    "shear_wave_decay     = ", .true.
        print *
    end if

    ! export initial condition
    if (this_image() == 1) then
        if (should_export_step(N_STEPS, 0_int32, export_mode, export_interval, &
            export_initial_state, export_final_state)) then
            call export_selected_data(export_mode, output_dir_name, export_num, 0_int32, rho, u_x, u_y)
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

        call fuzed_pull_streaming_collision_shear_wave_decay( &
            N_X, N_Y, N_DIRS, c_x, c_y, c_x_fp, c_y_fp, w, omega, &
            f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y)

        call swap_distribution_function_buffers(f, f_next)

        ! export selected field
        if (this_image() == 1) then
            if (should_export_step(N_STEPS, step, export_mode, export_interval, &
                export_initial_state, export_final_state)) then
                call export_selected_data(export_mode, output_dir_name, export_num, step, rho, u_x, u_y)
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
