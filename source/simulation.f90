module simulation
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_X, N_Y, N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, &
        USE_UNROLLED_KERNELS, USE_PUSH_SHIFT_KERNELS, &
        shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t
    use shear_wave, only: fuzed_pull_streaming_collision_outer_SW, fuzed_unrolled_pull_streaming_collision_outer_SW, &
        fuzed_push_shift_streaming_collision_full_SW, fuzed_unrolled_push_shift_streaming_collision_full_SW
    use couette_flow, only: fuzed_pull_streaming_collision_outer_CF
    use poiseuille_flow, only: fuzed_pull_streaming_collision_outer_PF
    use sliding_lid, only: fuzed_pull_streaming_collision_outer_SL
    implicit none

contains

    subroutine execute_full_sim_step( &
        shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
        write_macro_fields, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        type(shear_wave_params_t), intent(in) :: shear_wave_params
        type(couette_flow_params_t), intent(in) :: couette_flow_params
        type(poiseuille_flow_params_t), intent(in) :: poiseuille_flow_params
        type(sliding_lid_params_t), intent(in) :: sliding_lid_params
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(out) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! execute single sim step based on selected sim mode
        select case (SIM_MODE)
        case (SIM_SHEAR_WAVE) ! shear wave
            if (USE_PUSH_SHIFT_KERNELS) then
                if (USE_UNROLLED_KERNELS) then
                    call fuzed_unrolled_push_shift_streaming_collision_full_SW( &
                        write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
                else
                    call fuzed_push_shift_streaming_collision_full_SW( &
                        write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
                end if
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_unrolled_pull_streaming_collision_inner_universal( &
                    write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
                call fuzed_unrolled_pull_streaming_collision_outer_SW( &
                    write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_inner_universal( &
                    write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
                call fuzed_pull_streaming_collision_outer_SW( &
                    write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
            end if
        case (SIM_COUETTE_FLOW) ! couette flow
            if (USE_PUSH_SHIFT_KERNELS) then
                error stop "error: push-shift couette flow is not implemented yet"
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_unrolled_pull_streaming_collision_inner_universal( &
                    write_macro_fields, couette_flow_params%omega, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_inner_universal( &
                    write_macro_fields, couette_flow_params%omega, f, f_next, rho, u_x, u_y)
            end if
            call fuzed_pull_streaming_collision_outer_CF( &
                write_macro_fields, couette_flow_params%rho_0, couette_flow_params%omega, couette_flow_params%u_wall, &
                f, f_next, rho, u_x, u_y)
        case (SIM_POISEUILLE_FLOW) ! poiseuille flow
            if (USE_PUSH_SHIFT_KERNELS) then
                error stop "error: push-shift poiseuille flow is not implemented yet"
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_unrolled_pull_streaming_collision_inner_universal( &
                    write_macro_fields, poiseuille_flow_params%omega, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_inner_universal( &
                    write_macro_fields, poiseuille_flow_params%omega, f, f_next, rho, u_x, u_y)
            end if
            call fuzed_pull_streaming_collision_outer_PF( &
                write_macro_fields, poiseuille_flow_params%omega, &
                poiseuille_flow_params%rho_in, poiseuille_flow_params%rho_out, &
                f, f_next, rho, u_x, u_y)
        case (SIM_SLIDING_LID) ! sliding lid
            if (USE_PUSH_SHIFT_KERNELS) then
                error stop "error: push-shift sliding lid is not implemented yet"
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_unrolled_pull_streaming_collision_inner_universal( &
                    write_macro_fields, sliding_lid_params%omega, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_inner_universal( &
                    write_macro_fields, sliding_lid_params%omega, f, f_next, rho, u_x, u_y)
            end if
            call fuzed_pull_streaming_collision_outer_SL( &
                write_macro_fields, sliding_lid_params%rho_0, sliding_lid_params%omega, sliding_lid_params%u_wall, &
                f, f_next, rho, u_x, u_y)
        case default
            error stop "error: unknown sim mode in execute_full_sim_step()"
        end select
    end subroutine execute_full_sim_step


    subroutine fuzed_pull_streaming_collision_inner_universal( &
        write_macro_fields, omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        logical, intent(in) :: write_macro_fields
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

                ! finalize density and velocity
                u_x_val = u_x_val / rho_val
                u_y_val = u_y_val / rho_val
                u_squ = u_x_val * u_x_val + u_y_val * u_y_val

                if (write_macro_fields) then
                    rho(x, y) = rho_val
                    u_x(x, y) = u_x_val
                    u_y(x, y) = u_y_val
                end if

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


    subroutine fuzed_unrolled_pull_streaming_collision_inner_universal( &
        write_macro_fields, omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(out) :: f_next(N_X, N_Y, N_DIRS)
        real(FP), intent(inout) :: rho(N_X, N_Y)
        real(FP), intent(inout) :: u_x(N_X, N_Y)
        real(FP), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y
        real(FP) :: f_1
        real(FP) :: f_2
        real(FP) :: f_3
        real(FP) :: f_4
        real(FP) :: f_5
        real(FP) :: f_6
        real(FP) :: f_7
        real(FP) :: f_8
        real(FP) :: f_9
        real(FP) :: rho_val
        real(FP) :: u_x_val
        real(FP) :: u_y_val
        real(FP) :: u_squ

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

                rho_val = f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8 + f_9
                u_x_val = f_2 - f_4 + f_6 - f_7 - f_8 + f_9
                u_y_val = f_3 - f_5 + f_6 + f_7 - f_8 - f_9

                ! safety check to avoid division by zero in case of wrong density
            #ifdef FFB_DENSITY_CHECKS
                if (rho_val <= 0.0_FP) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
                end if
            #endif

                ! finalize density and velocity
                u_x_val = u_x_val / rho_val
                u_y_val = u_y_val / rho_val
                u_squ = u_x_val * u_x_val + u_y_val * u_y_val

                if (write_macro_fields) then
                    rho(x, y) = rho_val
                    u_x(x, y) = u_x_val
                    u_y(x, y) = u_y_val
                end if

                ! collide and stream to destination cells in all channels
                ! (manually unrolled)
                ! 1: (0, 0)
                f_next(x, y, 1) = f_1 - omega * (f_1 - (4.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP - 1.5_FP * u_squ))
                
                ! 2: (1, 0)
                f_next(x, y, 2) = f_2 - omega * (f_2 - (1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * u_x_val + 4.5_FP * u_x_val * u_x_val - &
                    1.5_FP * u_squ))

                ! 3: (0, 1)
                f_next(x, y, 3) = f_3 - omega * (f_3 - (1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * u_y_val + 4.5_FP * u_y_val * u_y_val - &
                    1.5_FP * u_squ))

                ! 4: (-1, 0)
                f_next(x, y, 4) = f_4 - omega * (f_4 - (1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP - 3.0_FP * u_x_val + 4.5_FP * u_x_val * u_x_val - &
                    1.5_FP * u_squ))

                ! 5: (0, -1)
                f_next(x, y, 5) = f_5 - omega * (f_5 - (1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP - 3.0_FP * u_y_val + 4.5_FP * u_y_val * u_y_val - &
                    1.5_FP * u_squ))

                ! 6: (1, 1)
                f_next(x, y, 6) = f_6 - omega * (f_6 - (1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * (u_x_val + u_y_val) + &
                    4.5_FP * (u_x_val + u_y_val) * (u_x_val + u_y_val) - &
                    1.5_FP * u_squ))

                ! 7: (-1, 1)
                f_next(x, y, 7) = f_7 - omega * (f_7 - (1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * (-u_x_val + u_y_val) + &
                    4.5_FP * (-u_x_val + u_y_val) * (-u_x_val + u_y_val) - &
                    1.5_FP * u_squ))

                ! 8: (-1, -1)
                f_next(x, y, 8) = f_8 - omega * (f_8 - (1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP - 3.0_FP * (u_x_val + u_y_val) + &
                    4.5_FP * (u_x_val + u_y_val) * (u_x_val + u_y_val) - &
                    1.5_FP * u_squ))
                
                ! 9: (1, -1)
                f_next(x, y, 9) = f_9 - omega * (f_9 - (1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * (u_x_val - u_y_val) + &
                    4.5_FP * (u_x_val - u_y_val) * (u_x_val - u_y_val) - &
                    1.5_FP * u_squ))
            end do
        end do
    end subroutine fuzed_unrolled_pull_streaming_collision_inner_universal


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
