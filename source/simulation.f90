module simulation
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_X, N_Y, N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, &
        shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t
    implicit none

contains

    subroutine execute_full_sim_step( &
        shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
        f, f_next, rho, u_x, u_y &
        )
        ! inputs
        type(shear_wave_params_t), intent(in) :: shear_wave_params
        type(couette_flow_params_t), intent(in) :: couette_flow_params
        type(poiseuille_flow_params_t), intent(in) :: poiseuille_flow_params
        type(sliding_lid_params_t), intent(in) :: sliding_lid_params
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(out) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! execute single sim step based on selected sim mode
        select case (SIM_MODE)
        case (SIM_SHEAR_WAVE)
            call push_shift_streaming_collision_shear_wave( &
                shear_wave_params%omega, f, f_next, rho, u_x, u_y)
        case (SIM_COUETTE_FLOW)
            call fuzed_pull_streaming_collision_inner_universal( &
                couette_flow_params%omega, f, f_next, rho, u_x, u_y)
            call fuzed_pull_streaming_collision_outer_couette_flow( &
                couette_flow_params%rho_0, couette_flow_params%omega, couette_flow_params%u_wall, &
                f, f_next, rho, u_x, u_y)
        case (SIM_POISEUILLE_FLOW)
            call fuzed_pull_streaming_collision_inner_universal( &
                poiseuille_flow_params%omega, f, f_next, rho, u_x, u_y)
            call fuzed_pull_streaming_collision_outer_poiseuille_flow( &
                poiseuille_flow_params%omega, poiseuille_flow_params%rho_in, poiseuille_flow_params%rho_out, &
                f, f_next, rho, u_x, u_y)
        case (SIM_SLIDING_LID)
            call fuzed_pull_streaming_collision_inner_universal( &
                sliding_lid_params%omega, f, f_next, rho, u_x, u_y)
            call fuzed_pull_streaming_collision_outer_sliding_lid( &
                sliding_lid_params%rho_0, sliding_lid_params%omega, sliding_lid_params%u_wall, &
                f, f_next, rho, u_x, u_y)
        case default
            error stop "error: unknown sim mode in execute_full_sim_step()"
        end select
    end subroutine execute_full_sim_step


    subroutine fuzed_pull_streaming_collision_inner_universal( &
        omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(out) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_1
        real(FP) :: f_2
        real(FP) :: f_3
        real(FP) :: f_4
        real(FP) :: f_5
        real(FP) :: f_6
        real(FP) :: f_7
        real(FP) :: f_8
        real(FP) :: f_9
        real(FP) :: f_pulled(N_DIRS)
        real(FP) :: rho_val
        real(FP) :: u_x_val
        real(FP) :: u_y_val
        real(FP) :: u_squ
        real(FP) :: c_dot_u
        real(FP) :: f_eq_val
        real(FP) :: f_next_val

        ! loop over rows and cols of inner cells only
        do y = 2, N_Y - 1
            do x = 2, N_X - 1

                ! 1: ( 0,  0) = rest
                ! 2: ( 1,  0) = east
                ! 3: ( 0,  1) = north
                ! 4: (-1,  0) = west
                ! 5: ( 0, -1) = south
                ! 6: ( 1,  1) = north-east
                ! 7: (-1,  1) = north-west
                ! 8: (-1, -1) = south-west
                ! 9: ( 1, -1) = south-east
                ! ---------
                ! | 7 3 6 |
                ! | 4 1 2 |
                ! | 8 5 9 |
                ! ---------
                ! pull streamed distribution functions from source cells in all channels
                ! (no boundary handling for inner cells, manually unrolled)
                f_1 = f(x, y, 1)
                f_2 = f(x - 1, y, 2)
                f_3 = f(x, y - 1, 3)
                f_4 = f(x + 1, y, 4)
                f_5 = f(x, y + 1, 5)
                f_6 = f(x - 1, y - 1, 6)
                f_7 = f(x + 1, y - 1, 7)
                f_8 = f(x + 1, y + 1, 8)
                f_9 = f(x - 1, y + 1, 9)

                f_pulled = [f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9]

                rho_val = f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8 + f_9
                u_x_val = f_2 - f_4 + f_6 - f_7 - f_8 + f_9
                u_y_val = f_3 - f_5 + f_6 + f_7 - f_8 - f_9

                ! safety check to avoid division by zero in case of wrong density
            #ifdef FFB_DENSITY_CHECKS
                if (rho_val <= 0.0_FP) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
                end if
            #endif

                ! finalize and store density and velocity
                u_x_val = u_x_val / rho_val
                u_y_val = u_y_val / rho_val
                u_squ = u_x_val * u_x_val + u_y_val * u_y_val
                rho(x, y) = rho_val
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val

                ! collide and stream to destination cells in all channels
                !DIR$ UNROLL(9)
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for channel i
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_val * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)

                    ! relax towards equilibrium and write to destination channel in this cell
                    f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)
                    f_next(x, y, i) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_inner_universal


    subroutine fuzed_pull_streaming_collision_outer_shear_wave( &
        omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y

        ! bottom row
        y = 1
        do x = 1, N_X
            call collide_stream_outer_cell_shear_wave(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call collide_stream_outer_cell_shear_wave(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_shear_wave(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_shear_wave(x, y)
        end do

    contains ! helper subroutine

        subroutine collide_stream_outer_cell_shear_wave( &
            x, y &
            )
            ! inputs
            integer(int32), intent(in) :: x
            integer(int32), intent(in) :: y

            ! temp
            integer(int32) :: i
            integer(int32) :: src_x, src_y
            real(FP) :: f_pulled(N_DIRS)
            real(FP) :: rho_val
            real(FP) :: u_x_val
            real(FP) :: u_y_val
            real(FP) :: u_squ
            real(FP) :: c_dot_u
            real(FP) :: f_eq_val
            real(FP) :: f_next_val

            rho_val = 0.0_FP
            u_x_val = 0.0_FP
            u_y_val = 0.0_FP

            ! 1: ( 0,  0) = rest
            ! 2: ( 1,  0) = east
            ! 3: ( 0,  1) = north
            ! 4: (-1,  0) = west
            ! 5: ( 0, -1) = south
            ! 6: ( 1,  1) = north-east
            ! 7: (-1,  1) = north-west
            ! 8: (-1, -1) = south-west
            ! 9: ( 1, -1) = south-east
            ! ---------
            ! | 7 3 6 |
            ! | 4 1 2 |
            ! | 8 5 9 |
            ! ---------
            ! pull streamed distribution functions from source cells in all channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                src_x = x - C_X(i)
                src_y = y - C_Y(i)

                ! periodic for left/right boundary
                if (src_x < 1) then
                    src_x = N_X
                else if (src_x > N_X) then
                    src_x = 1
                end if

                ! periodic for bottom/top boundary
                if (src_y < 1) then
                    src_y = N_Y
                else if (src_y > N_Y) then
                    src_y = 1
                end if

                f_pulled(i) = f(src_x, src_y, i)

                rho_val = rho_val + f_pulled(i)
                u_x_val = u_x_val + f_pulled(i) * C_X_FP(i)
                u_y_val = u_y_val + f_pulled(i) * C_Y_FP(i)
            end do

            ! safety check to avoid division by zero in case of wrong density
        #ifdef FFB_DENSITY_CHECKS
            if (rho_val <= 0.0_FP) then
                error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
            end if
        #endif

            ! finalize and store density and velocity
            u_x_val = u_x_val / rho_val
            u_y_val = u_y_val / rho_val
            u_squ = u_x_val * u_x_val + u_y_val * u_y_val
            rho(x, y) = rho_val
            u_x(x, y) = u_x_val
            u_y(x, y) = u_y_val

            ! collide and stream locally to destination channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                ! compute equilibrium distribution function for channel i
                c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                f_eq_val = W(i) * rho_val * ( &
                    1.0_FP + &
                    3.0_FP * c_dot_u + &
                    4.5_FP * c_dot_u * c_dot_u - &
                    1.5_FP * u_squ)

                ! relax towards equilibrium and write to destination channel in this cell
                f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)
                f_next(x, y, i) = f_next_val
            end do
        end subroutine collide_stream_outer_cell_shear_wave
    end subroutine fuzed_pull_streaming_collision_outer_shear_wave


    subroutine push_shift_streaming_collision_shear_wave( &
        omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_1
        real(FP) :: f_2
        real(FP) :: f_3
        real(FP) :: f_4
        real(FP) :: f_5
        real(FP) :: f_6
        real(FP) :: f_7
        real(FP) :: f_8
        real(FP) :: f_9
        real(FP) :: f_pulled(N_DIRS)
        real(FP) :: rho_val
        real(FP) :: u_x_val
        real(FP) :: u_y_val
        real(FP) :: u_squ
        real(FP) :: c_dot_u
        real(FP) :: f_eq_val
        real(FP) :: f_next_val

        ! 1: ( 0,  0) = rest
        ! 2: ( 1,  0) = east
        ! 3: ( 0,  1) = north
        ! 4: (-1,  0) = west
        ! 5: ( 0, -1) = south
        ! 6: ( 1,  1) = north-east
        ! 7: (-1,  1) = north-west
        ! 8: (-1, -1) = south-west
        ! 9: ( 1, -1) = south-east
        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! periodic push-streaming of f into f_next (temporary storage)
        f_next(:, :, 2) = cshift(f(:, :, 2), shift=-1, dim=1)
        f_next(:, :, 3) = cshift(f(:, :, 3), shift=-1, dim=2)
        f_next(:, :, 4) = cshift(f(:, :, 4), shift=1, dim=1)
        f_next(:, :, 5) = cshift(f(:, :, 5), shift=1, dim=2)
        f_next(:, :, 6) = cshift(cshift(f(:, :, 6), shift=-1, dim=1), shift=-1, dim=2)
        f_next(:, :, 7) = cshift(cshift(f(:, :, 7), shift=1, dim=1), shift=-1, dim=2)
        f_next(:, :, 8) = cshift(cshift(f(:, :, 8), shift=1, dim=1), shift=1, dim=2)
        f_next(:, :, 9) = cshift(cshift(f(:, :, 9), shift=-1, dim=1), shift=1, dim=2)

        ! loop over rows and cols of all cells
        do y = 1, N_Y
            do x = 1, N_X

                ! read streamed distribution functions from f_next in all moving channels
                ! (rest channel does not move)
                f_1 = f(x, y, 1)
                f_2 = f_next(x, y, 2)
                f_3 = f_next(x, y, 3)
                f_4 = f_next(x, y, 4)
                f_5 = f_next(x, y, 5)
                f_6 = f_next(x, y, 6)
                f_7 = f_next(x, y, 7)
                f_8 = f_next(x, y, 8)
                f_9 = f_next(x, y, 9)

                f_pulled = [f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9]

                rho_val = f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8 + f_9
                u_x_val = f_2 - f_4 + f_6 - f_7 - f_8 + f_9
                u_y_val = f_3 - f_5 + f_6 + f_7 - f_8 - f_9

                ! safety check to avoid division by zero in case of wrong density
            #ifdef FFB_DENSITY_CHECKS
                if (rho_val <= 0.0_FP) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
                end if
            #endif

                ! finalize and store density and velocity
                u_x_val = u_x_val / rho_val
                u_y_val = u_y_val / rho_val
                u_squ = u_x_val * u_x_val + u_y_val * u_y_val
                rho(x, y) = rho_val
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val

                ! collide and stream locally to destination channels
                !DIR$ UNROLL(9)
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for channel i
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_val * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)

                    ! relax towards equilibrium and write to destination channel in this cell
                    f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)
                    f_next(x, y, i) = f_next_val
                end do
            end do
        end do
    end subroutine push_shift_streaming_collision_shear_wave


    subroutine fuzed_pull_streaming_collision_outer_couette_flow( &
        rho_0, omega, u_wall, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: u_wall
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y

        ! bottom row
        y = 1
        do x = 1, N_X
            call collide_stream_outer_cell_couette_flow(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call collide_stream_outer_cell_couette_flow(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_couette_flow(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_couette_flow(x, y)
        end do

    contains ! helper subroutine

        subroutine collide_stream_outer_cell_couette_flow( &
            x, y &
            )
            ! inputs
            integer(int32), intent(in) :: x
            integer(int32), intent(in) :: y

            ! temp
            integer(int32) :: i
            integer(int32) :: src_x, src_y
            real(FP) :: f_pulled(N_DIRS)
            real(FP) :: rho_val
            real(FP) :: u_x_val
            real(FP) :: u_y_val
            real(FP) :: u_squ
            real(FP) :: c_dot_u
            real(FP) :: f_eq_val
            real(FP) :: f_next_val

            rho_val = 0.0_FP
            u_x_val = 0.0_FP
            u_y_val = 0.0_FP

            ! 1: ( 0,  0) = rest
            ! 2: ( 1,  0) = east
            ! 3: ( 0,  1) = north
            ! 4: (-1,  0) = west
            ! 5: ( 0, -1) = south
            ! 6: ( 1,  1) = north-east
            ! 7: (-1,  1) = north-west
            ! 8: (-1, -1) = south-west
            ! 9: ( 1, -1) = south-east
            ! ---------
            ! | 7 3 6 |
            ! | 4 1 2 |
            ! | 8 5 9 |
            ! ---------
            ! pull streamed distribution functions from source cells in all channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                src_x = x - C_X(i)
                src_y = y - C_Y(i)

                ! periodic for left/right boundary
                if (src_x < 1) then
                    src_x = N_X
                else if (src_x > N_X) then
                    src_x = 1
                end if

                ! no period or bounce-back
                if (src_y >= 1 .and. src_y <= N_Y) then
                    f_pulled(i) = f(src_x, src_y, i)

                ! bounce-back for bottom boundary (static)
                else if (src_y < 1) then
                    select case (i)
                    case (3)
                        f_pulled(i) = f(x, y, 5)
                    case (6)
                        f_pulled(i) = f(x, y, 8)
                    case (7)
                        f_pulled(i) = f(x, y, 9)
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid bottom boundary channel in couette flow"
                #endif
                    end select
                
                ! bounce-back for top boundary (moving)
                else if (src_y > N_Y) then
                    select case(i)
                    case (5)
                        f_pulled(i) = f(x, y, 3)
                    case (8)
                        f_pulled(i) = f(x, y, 6) - 6.0_FP * W(6) * rho_0 * u_wall
                    case (9)
                        f_pulled(i) = f(x, y, 7) + 6.0_FP * W(7) * rho_0 * u_wall
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid top boundary channel in couette flow"
                #endif
                    end select
                end if

                rho_val = rho_val + f_pulled(i)
                u_x_val = u_x_val + f_pulled(i) * C_X_FP(i)
                u_y_val = u_y_val + f_pulled(i) * C_Y_FP(i)
            end do

            ! safety check to avoid division by zero in case of wrong density
        #ifdef FFB_DENSITY_CHECKS
            if (rho_val <= 0.0_FP) then
                error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
            end if
        #endif

            ! finalize and store density and velocity
            u_x_val = u_x_val / rho_val
            u_y_val = u_y_val / rho_val
            u_squ = u_x_val * u_x_val + u_y_val * u_y_val
            rho(x, y) = rho_val
            u_x(x, y) = u_x_val
            u_y(x, y) = u_y_val

            ! collide and stream locally to destination channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                ! compute equilibrium distribution function for channel i
                c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                f_eq_val = W(i) * rho_val * ( &
                    1.0_FP + &
                    3.0_FP * c_dot_u + &
                    4.5_FP * c_dot_u * c_dot_u - &
                    1.5_FP * u_squ)

                ! relax towards equilibrium and write to destination channel in this cell
                f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)
                f_next(x, y, i) = f_next_val
            end do
        end subroutine collide_stream_outer_cell_couette_flow
    end subroutine fuzed_pull_streaming_collision_outer_couette_flow


    subroutine fuzed_pull_streaming_collision_outer_poiseuille_flow( &
        omega, rho_in, rho_out, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: rho_in
        real(FP), intent(in) :: rho_out
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y

        ! copied boundary densities and velocities for streaming step
        real(FP) :: rho_left(N_Y), rho_right(N_Y)
        real(FP) :: u_x_left(N_Y), u_x_right(N_Y)
        real(FP) :: u_y_left(N_Y), u_y_right(N_Y)
        rho_left(:) = rho(1, :)
        rho_right(:) = rho(N_X, :)
        u_x_left(:) = u_x(1, :)
        u_x_right(:) = u_x(N_X, :)
        u_y_left(:) = u_y(1, :)
        u_y_right(:) = u_y(N_X, :)

        ! bottom row
        y = 1
        do x = 1, N_X
            call collide_stream_outer_cell_poiseuille_flow(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call collide_stream_outer_cell_poiseuille_flow(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_poiseuille_flow(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_poiseuille_flow(x, y)
        end do

    contains ! helper subroutine

        subroutine collide_stream_outer_cell_poiseuille_flow( &
            x, y &
            )
            ! inputs
            integer(int32), intent(in) :: x
            integer(int32), intent(in) :: y

            ! temp
            integer(int32) :: i
            integer(int32) :: src_x, src_y
            real(FP) :: f_pulled(N_DIRS)
            real(FP) :: rho_val
            real(FP) :: u_x_val
            real(FP) :: u_y_val
            real(FP) :: u_squ
            real(FP) :: c_dot_u
            real(FP) :: f_eq_val
            real(FP) :: f_next_val
            real(FP) :: u_squ_src
            real(FP) :: c_dot_u_src
            real(FP) :: f_eq_src
            real(FP) :: f_eq_boundary

            rho_val = 0.0_FP
            u_x_val = 0.0_FP
            u_y_val = 0.0_FP

            ! 1: ( 0,  0) = rest
            ! 2: ( 1,  0) = east
            ! 3: ( 0,  1) = north
            ! 4: (-1,  0) = west
            ! 5: ( 0, -1) = south
            ! 6: ( 1,  1) = north-east
            ! 7: (-1,  1) = north-west
            ! 8: (-1, -1) = south-west
            ! 9: ( 1, -1) = south-east
            ! ---------
            ! | 7 3 6 |
            ! | 4 1 2 |
            ! | 8 5 9 |
            ! ---------
            ! pull streamed distribution functions from source cells in all channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                src_x = x - C_X(i)
                src_y = y - C_Y(i)

                ! no period or bounce-back
                if (src_x >= 1 .and. src_x <= N_X .and. &
                    src_y >= 1 .and. src_y <= N_Y) then
                    f_pulled(i) = f(src_x, src_y, i)
                
                ! bounce-back for bottom boundary (static)
                else if (src_y < 1) then
                    select case (i)
                    case (3)
                        f_pulled(i) = f(x, y, 5)
                    case (6)
                        f_pulled(i) = f(x, y, 8)
                    case (7)
                        f_pulled(i) = f(x, y, 9)
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid bottom boundary channel in poiseuille flow"
                #endif
                    end select

                ! bounce-back for top boundary (static)
                else if (src_y > N_Y) then
                    select case (i)
                    case (5)
                        f_pulled(i) = f(x, y, 3)
                    case (8)
                        f_pulled(i) = f(x, y, 6)
                    case (9)
                        f_pulled(i) = f(x, y, 7)
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid top boundary channel in poiseuille flow"
                #endif
                    end select
                
                ! pressure-periodic inlet for left boundary
                else if (src_x < 1) then

                    u_squ_src = u_x_right(y) * u_x_right(y) + u_y_right(y) * u_y_right(y)
                    c_dot_u_src = C_X_FP(i) * u_x_right(y) + C_Y_FP(i) * u_y_right(y)

                    ! equilibrium distribution function for the source cell at the opposite (right) boundary
                    f_eq_src = W(i) * rho_right(y) * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - &
                        1.5_FP * u_squ_src)

                    ! equilibrium distribution function at this cell
                    f_eq_boundary = W(i) * rho_in * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - &
                        1.5_FP * u_squ_src)

                    f_pulled(i) = f(N_X, y, i) - f_eq_src + f_eq_boundary

                ! pressure-periodic outlet for right boundary
                else if (src_x > N_X) then

                    u_squ_src = u_x_left(y) * u_x_left(y) + u_y_left(y) * u_y_left(y)
                    c_dot_u_src = C_X_FP(i) * u_x_left(y) + C_Y_FP(i) * u_y_left(y)

                    ! equilibrium distribution function for the source cell at the opposite (left) boundary
                    f_eq_src = W(i) * rho_left(y) * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - &
                        1.5_FP * u_squ_src)

                    ! equilibrium distribution function at this cell
                    f_eq_boundary = W(i) * rho_out * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - &
                        1.5_FP * u_squ_src)

                    f_pulled(i) = f(1, y, i) - f_eq_src + f_eq_boundary
                end if

                rho_val = rho_val + f_pulled(i)
                u_x_val = u_x_val + f_pulled(i) * C_X_FP(i)
                u_y_val = u_y_val + f_pulled(i) * C_Y_FP(i)
            end do

            ! safety check to avoid division by zero in case of wrong density
        #ifdef FFB_DENSITY_CHECKS
            if (rho_val <= 0.0_FP) then
                error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
            end if
        #endif

            ! finalize and store density and velocity
            u_x_val = u_x_val / rho_val
            u_y_val = u_y_val / rho_val
            u_squ = u_x_val * u_x_val + u_y_val * u_y_val
            rho(x, y) = rho_val
            u_x(x, y) = u_x_val
            u_y(x, y) = u_y_val

            ! collide and stream locally to destination channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                ! compute equilibrium distribution function for channels i
                c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                f_eq_val = W(i) * rho_val * ( &
                    1.0_FP + &
                    3.0_FP * c_dot_u + &
                    4.5_FP * c_dot_u * c_dot_u - &
                    1.5_FP * u_squ)

                ! relax towards equilibrium and write to destination channel in this cell
                f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)
                f_next(x, y, i) = f_next_val
            end do
        end subroutine collide_stream_outer_cell_poiseuille_flow
    end subroutine fuzed_pull_streaming_collision_outer_poiseuille_flow


    subroutine fuzed_pull_streaming_collision_outer_sliding_lid( &
        rho_0, omega, u_wall, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: u_wall
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y

        ! bottom row
        y = 1
        do x = 1, N_X
            call collide_stream_outer_cell_sliding_lid(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call collide_stream_outer_cell_sliding_lid(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_sliding_lid(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call collide_stream_outer_cell_sliding_lid(x, y)
        end do

    contains ! helper subroutine

        subroutine collide_stream_outer_cell_sliding_lid( &
            x, y &
            )
            ! inputs
            integer(int32), intent(in) :: x
            integer(int32), intent(in) :: y

            ! temp
            integer(int32) :: i
            integer(int32) :: src_x, src_y
            real(FP) :: f_pulled(N_DIRS)
            real(FP) :: rho_val
            real(FP) :: u_x_val
            real(FP) :: u_y_val
            real(FP) :: u_squ
            real(FP) :: c_dot_u
            real(FP) :: f_eq_val
            real(FP) :: f_next_val

            rho_val = 0.0_FP
            u_x_val = 0.0_FP
            u_y_val = 0.0_FP

            ! 1: ( 0,  0) = rest
            ! 2: ( 1,  0) = east
            ! 3: ( 0,  1) = north
            ! 4: (-1,  0) = west
            ! 5: ( 0, -1) = south
            ! 6: ( 1,  1) = north-east
            ! 7: (-1,  1) = north-west
            ! 8: (-1, -1) = south-west
            ! 9: ( 1, -1) = south-east
            ! ---------
            ! | 7 3 6 |
            ! | 4 1 2 |
            ! | 8 5 9 |
            ! ---------
            ! pull streamed distribution functions from source cells in all channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS
                
                src_x = x - C_X(i)
                src_y = y - C_Y(i)

                ! no bounce-back
                if (src_x >= 1 .and. src_x <= N_X .and. &
                    src_y >= 1 .and. src_y <= N_Y) then
                    f_pulled(i) = f(src_x, src_y, i)
                
                ! bounce-back for bottom boundary (static)
                else if (src_y < 1) then
                    select case (i)
                    case (3)
                        f_pulled(i) = f(x, y, 5)
                    case (6)
                        f_pulled(i) = f(x, y, 8)
                    case (7)
                        f_pulled(i) = f(x, y, 9)
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid bottom boundary channel in sliding lid"
                #endif
                    end select

                ! bounce-back for top boundary (moving)
                else if (src_y > N_Y) then
                    select case (i)
                    case (5)
                        f_pulled(i) = f(x, y, 3)
                    case (8)
                        f_pulled(i) = f(x, y, 6) - 6.0_FP * W(6) * rho_0 * u_wall
                    case (9)
                        f_pulled(i) = f(x, y, 7) + 6.0_FP * W(7) * rho_0 * u_wall
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid top boundary channel in sliding lid"
                #endif
                    end select

                ! bounce-back for left boundary (static)
                else if (src_x < 1) then
                    select case (i)
                    case (2)
                        f_pulled(i) = f(x, y, 4)
                    case (6)
                        f_pulled(i) = f(x, y, 8)
                    case (9)
                        f_pulled(i) = f(x, y, 7)
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid left boundary channel in sliding lid"
                #endif
                    end select
                
                ! bounce-back for right boundary (static)
                else if (src_x > N_X) then
                    select case (i)
                    case (4)
                        f_pulled(i) = f(x, y, 2)
                    case (7)
                        f_pulled(i) = f(x, y, 9)
                    case (8)
                        f_pulled(i) = f(x, y, 6)
                #ifdef FFB_BOUNDARY_CHECKS
                    case default
                        error stop "error: invalid right boundary channel in sliding lid"
                #endif
                    end select
                end if

                rho_val = rho_val + f_pulled(i)
                u_x_val = u_x_val + f_pulled(i) * C_X_FP(i)
                u_y_val = u_y_val + f_pulled(i) * C_Y_FP(i)
            end do

            ! safety check to avoid division by zero in case of wrong density
        #ifdef FFB_DENSITY_CHECKS
            if (rho_val <= 0.0_FP) then
                error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
            end if
        #endif

            ! finalize and store density and velocity
            u_x_val = u_x_val / rho_val
            u_y_val = u_y_val / rho_val
            u_squ = u_x_val * u_x_val + u_y_val * u_y_val
            rho(x, y) = rho_val
            u_x(x, y) = u_x_val
            u_y(x, y) = u_y_val

            ! collide and stream locally to destination channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                ! compute equilibrium distribution function for channels i
                c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                f_eq_val = W(i) * rho_val * ( &
                    1.0_FP + &
                    3.0_FP * c_dot_u + &
                    4.5_FP * c_dot_u * c_dot_u - &
                    1.5_FP * u_squ)

                ! relax towards equilibrium and write to destination channel in this cell
                f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)
                f_next(x, y, i) = f_next_val
            end do
         end subroutine collide_stream_outer_cell_sliding_lid
    end subroutine fuzed_pull_streaming_collision_outer_sliding_lid


    subroutine swap_distribution_function_buffers( &
        f, f_next &
        )
        ! read/write inputs
        real(FP), allocatable, intent(inout) :: f(:, :, :)
        real(FP), allocatable, intent(inout) :: f_next(:, :, :)
        real(FP), allocatable :: temp(:, :, :)

        ! swap ownership
        call move_alloc(f, temp)
        call move_alloc(f_next, f)
        call move_alloc(temp, f_next)
    end subroutine swap_distribution_function_buffers

end module simulation
