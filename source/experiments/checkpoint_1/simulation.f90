module simulation
    ! imports
    use iso_fortran_env, only: int32
    use domain, only: domain_t
    use exchange, only: halo_buffers_t
    use settings, only: N_DIRS, C_X_FP, C_Y_FP, W, &
        SIM_SLIDING_LID, SIM_MODE, FP, &
        RHO_0, OMEGA, U_LID
    implicit none

contains

    subroutine execute_local_sim_step( &
        domain_info, halo_buffers, n_x_local, n_y_local, &
        write_macro_fields, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: write_macro_fields
        real(FP), intent(inout) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)

        ! write destinations
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f_next(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)

        if (SIM_MODE /= SIM_SLIDING_LID) then
            error stop "error: checkpoint 1 only supports sliding lid"
        end if

        call fuzed_pull_streaming_collision_local_branchy( &
            domain_info, n_x_local, n_y_local, &
            write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
    end subroutine execute_local_sim_step


    subroutine fuzed_pull_streaming_collision_local_branchy( &
        domain_info, n_x_local, n_y_local, write_macro_fields, omega, &
        f, f_next, rho, u_x, u_y &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_pulled(N_DIRS)
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
        real(FP) :: c_dot_u
        real(FP) :: f_eq_val

        do y = 1, n_y_local
            do x = 1, n_x_local
                call pull_sliding_lid_populations_branchy( &
                    n_x_local, n_y_local, x, y, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    RHO_0, U_LID, f, &
                    f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9)
                f_pulled = [f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9]

                rho_val = 0.0_FP
                u_x_val = 0.0_FP
                u_y_val = 0.0_FP

                do i = 1, N_DIRS
                    rho_val = rho_val + f_pulled(i)
                    u_x_val = u_x_val + f_pulled(i) * C_X_FP(i)
                    u_y_val = u_y_val + f_pulled(i) * C_Y_FP(i)
                end do

            #ifdef FFB_DENSITY_CHECKS
                if (rho_val <= 0.0_FP) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
                end if
            #endif

                u_x_val = u_x_val / rho_val
                u_y_val = u_y_val / rho_val
                u_squ = u_x_val * u_x_val + u_y_val * u_y_val

                if (write_macro_fields) then
                    rho(x, y) = rho_val
                    u_x(x, y) = u_x_val
                    u_y(x, y) = u_y_val
                end if

                do i = 1, N_DIRS
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_val * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)
                    f_next(i, x, y) = f_pulled(i) + omega * (f_eq_val - f_pulled(i))
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_local_branchy


    pure subroutine pull_sliding_lid_populations_branchy( &
        n_x_local, n_y_local, x, y, &
        at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        rho_0, u_wall, f, &
        f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9 &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        integer(int32), intent(in) :: x
        integer(int32), intent(in) :: y
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: u_wall
        real(FP), intent(in) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)

        ! output
        real(FP), intent(out) :: f_1
        real(FP), intent(out) :: f_2
        real(FP), intent(out) :: f_3
        real(FP), intent(out) :: f_4
        real(FP), intent(out) :: f_5
        real(FP), intent(out) :: f_6
        real(FP), intent(out) :: f_7
        real(FP), intent(out) :: f_8
        real(FP), intent(out) :: f_9

        ! temp
        real(FP) :: moving_wall_correction_8
        real(FP) :: moving_wall_correction_9

        ! pull streamed distribution functions from source cells
        f_1 = f(1, x, y)
        f_2 = f(2, x-1, y)
        f_3 = f(3, x, y-1)
        f_4 = f(4, x+1, y)
        f_5 = f(5, x, y+1)
        f_6 = f(6, x-1, y-1)
        f_7 = f(7, x+1, y-1)
        f_8 = f(8, x+1, y+1)
        f_9 = f(9, x-1, y+1)

        ! left bounce-back boundary
        if (at_left_boundary .and. x == 1) then
            f_2 = f(4, 1, y)
            f_6 = f(8, 1, y)
            f_9 = f(7, 1, y)
        end if

        ! right bounce-back boundary
        if (at_right_boundary .and. x == n_x_local) then
            f_4 = f(2, n_x_local, y)
            f_7 = f(9, n_x_local, y)
            f_8 = f(6, n_x_local, y)
        end if

        ! bottom bounce-back boundary
        if (at_bottom_boundary .and. y == 1) then
            f_3 = f(5, x, 1)
            f_6 = f(8, x, 1)
            f_7 = f(9, x, 1)
        end if

        ! top moving-wall bounce-back boundary
        if (at_top_boundary .and. y == n_y_local) then
            moving_wall_correction_8 = 6.0_FP * W(6) * rho_0 * u_wall
            moving_wall_correction_9 = 6.0_FP * W(7) * rho_0 * u_wall

            f_5 = f(3, x, n_y_local)
            f_8 = f(6, x, n_y_local) - moving_wall_correction_8
            f_9 = f(7, x, n_y_local) + moving_wall_correction_9
        end if
    end subroutine pull_sliding_lid_populations_branchy


end module simulation
