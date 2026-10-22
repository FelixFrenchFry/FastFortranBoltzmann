program test_collision
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: FP, N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W
    use shear_wave, only: prepare_shear_wave_halos_SW
    use simulation, only: fuzed_pull_streaming_collision_local_universal, &
        fuzed_pull_streaming_collision_local_unrolled_universal
    use initialization, only: apply_condition_couette_flow_local
    use assertions, only: check_close
    implicit none

    ! sim size
    integer(int32), parameter :: N_X_LOCAL = 5
    integer(int32), parameter :: N_Y_LOCAL = 3

    ! distribution functions
    real(FP) :: f(0:N_X_LOCAL+1, 0:N_Y_LOCAL+1, N_DIRS)
    real(FP) :: f_next(0:N_X_LOCAL+1, 0:N_Y_LOCAL+1, N_DIRS)
    real(FP) :: f_next_unrolled(0:N_X_LOCAL+1, 0:N_Y_LOCAL+1, N_DIRS)

    ! macro fields
    real(FP) :: rho(N_X_LOCAL, N_Y_LOCAL)
    real(FP) :: u_x(N_X_LOCAL, N_Y_LOCAL)
    real(FP) :: u_y(N_X_LOCAL, N_Y_LOCAL)
    real(FP) :: rho_unrolled(N_X_LOCAL, N_Y_LOCAL)
    real(FP) :: u_x_unrolled(N_X_LOCAL, N_Y_LOCAL)
    real(FP) :: u_y_unrolled(N_X_LOCAL, N_Y_LOCAL)

    ! temp
    integer(int32) :: x, y, i
    real(FP) :: mass
    real(FP) :: momentum_x
    real(FP) :: momentum_y
    real(FP) :: f_pulled(N_DIRS)

    ! verify that a resting equilibrium stays stationary
    call apply_condition_couette_flow_local( &
        1.2_FP, N_X_LOCAL, N_Y_LOCAL, f, rho, u_x, u_y)
    call prepare_shear_wave_halos_SW( &
        1, 1, N_X_LOCAL, N_Y_LOCAL, f)
    call fuzed_pull_streaming_collision_local_universal( &
        N_X_LOCAL, N_Y_LOCAL, .true., 1.1_FP, f, f_next, rho, u_x, u_y)

    do i = 1, N_DIRS
        do y = 1, N_Y_LOCAL
            do x = 1, N_X_LOCAL
                call check_close( &
                    f_next(x, y, i), 1.2_FP * W(i), "stationary equilibrium distribution")
                call check_close(rho(x, y), 1.2_FP, "stationary equilibrium density")
                call check_close(u_x(x, y), 0.0_FP, "stationary equilibrium x velocity")
                call check_close(u_y(x, y), 0.0_FP, "stationary equilibrium y velocity")
            end do
        end do
    end do

    ! nonuniform positive distributions exercise streaming and relaxation together
    do i = 1, N_DIRS
        do y = 1, N_Y_LOCAL
            do x = 1, N_X_LOCAL
                f(x, y, i) = W(i) * (1.0_FP + 0.01_FP * sin(real(x + 3 * y + 7 * i, FP)))
            end do
        end do
    end do
    call prepare_shear_wave_halos_SW( &
        1, 1, N_X_LOCAL, N_Y_LOCAL, f)

    ! record mass and momentum before the collision step
    mass = sum(f(1:N_X_LOCAL, 1:N_Y_LOCAL, :))
    momentum_x = 0.0_FP
    momentum_y = 0.0_FP
    do i = 1, N_DIRS
        momentum_x = momentum_x + sum(f(1:N_X_LOCAL, 1:N_Y_LOCAL, i)) * C_X_FP(i)
        momentum_y = momentum_y + sum(f(1:N_X_LOCAL, 1:N_Y_LOCAL, i)) * C_Y_FP(i)
    end do

    ! execute both kernel variants with the same input
    call fuzed_pull_streaming_collision_local_universal( &
        N_X_LOCAL, N_Y_LOCAL, .true., 1.1_FP, f, f_next, rho, u_x, u_y)
    call fuzed_pull_streaming_collision_local_unrolled_universal( &
        N_X_LOCAL, N_Y_LOCAL, .true., 1.1_FP, f, f_next_unrolled, rho_unrolled, u_x_unrolled, u_y_unrolled)

    ! verify mass and momentum conservation
    call check_close( &
        sum(f_next(1:N_X_LOCAL, 1:N_Y_LOCAL, :)), mass, "collision conserves mass")
    do i = 1, N_DIRS
        momentum_x = momentum_x - sum(f_next(1:N_X_LOCAL, 1:N_Y_LOCAL, i)) * C_X_FP(i)
        momentum_y = momentum_y - sum(f_next(1:N_X_LOCAL, 1:N_Y_LOCAL, i)) * C_Y_FP(i)
    end do
    call check_close(momentum_x, 0.0_FP, "collision conserves x momentum")
    call check_close(momentum_y, 0.0_FP, "collision conserves y momentum")

    ! compare distributions and macro fields from both kernel variants
    do y = 1, N_Y_LOCAL
        do x = 1, N_X_LOCAL
            do i = 1, N_DIRS
                f_pulled(i) = f(x - C_X(i), y - C_Y(i), i)
                call check_close( &
                    f_next_unrolled(x, y, i), f_next(x, y, i), "unrolled distribution")
            end do
            call check_close(rho(x, y), sum(f_pulled), "streamed density")
            call check_close( &
                u_x(x, y), sum(f_pulled * C_X_FP) / sum(f_pulled), "streamed x velocity")
            call check_close( &
                u_y(x, y), sum(f_pulled * C_Y_FP) / sum(f_pulled), "streamed y velocity")
            call check_close(rho_unrolled(x, y), rho(x, y), "unrolled density")
            call check_close(u_x_unrolled(x, y), u_x(x, y), "unrolled x velocity")
            call check_close(u_y_unrolled(x, y), u_y(x, y), "unrolled y velocity")
        end do
    end do
end program test_collision
