module simulation
    ! imports
    use iso_fortran_env, only: int32
    use domain, only: domain_t
    use exchange, only: halo_buffers_t
    use settings, only: N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, &
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
            error stop "error: checkpoint 2 only supports sliding lid"
        end if

        ! prepare sliding-lid boundary halos before universal pull streaming
        call prepare_sliding_lid_halos_SL( &
            n_x_local, n_y_local, &
            domain_info%at_left_boundary, domain_info%at_right_boundary, &
            domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
            RHO_0, U_LID, f)

        ! universal pull streaming + collision kernel
        call fuzed_pull_streaming_collision_local_universal( &
            n_x_local, n_y_local, &
            write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
    end subroutine execute_local_sim_step


    subroutine fuzed_pull_streaming_collision_local_universal( &
        n_x_local, n_y_local, write_macro_fields, omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
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
        integer(int32) :: src_x, src_y
        real(FP) :: f_pulled(N_DIRS)
        real(FP) :: rho_val
        real(FP) :: u_x_val
        real(FP) :: u_y_val
        real(FP) :: u_squ
        real(FP) :: c_dot_u
        real(FP) :: f_eq_val
        real(FP) :: f_next_val

        ! loop over all image-owned cells
        do y = 1, n_y_local
            do x = 1, n_x_local

                rho_val = 0.0_FP
                u_x_val = 0.0_FP
                u_y_val = 0.0_FP

                ! ---------
                ! | 7 3 6 |
                ! | 4 1 2 |
                ! | 8 5 9 |
                ! ---------
                ! pull streamed distribution functions from source cells
                ! (boundaries handled separately in sim-mode-specific halo preparation step)
                do i = 1, N_DIRS

                    src_x = x - C_X(i)
                    src_y = y - C_Y(i)

                    f_pulled(i) = f(i, src_x, src_y)

                    rho_val = rho_val + f_pulled(i)
                    u_x_val = u_x_val + f_pulled(i) * C_X_FP(i)
                    u_y_val = u_y_val + f_pulled(i) * C_Y_FP(i)
                end do

                ! debug check
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

                ! collide and stream locally
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for channel i
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_val * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)

                    ! relax towards equilibrium and write to destination channel
                    f_next_val = f_pulled(i) + omega * (f_eq_val - f_pulled(i))
                    f_next(i, x, y) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_local_universal


    subroutine prepare_sliding_lid_halos_SL( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        rho_0, u_wall, f &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: u_wall

        ! read/write inputs
        real(FP), intent(inout) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! temp
        integer(int32) :: x, y
        real(FP) :: moving_wall_correction_8
        real(FP) :: moving_wall_correction_9

        ! left bounce-back boundary, written into the halo column used by pull streaming
        if (at_left_boundary) then
            do y = 1, n_y_local
                f(2, 0, y) = f(4, 1, y)
                f(6, 0, y-1) = f(8, 1, y)
                f(9, 0, y+1) = f(7, 1, y)
            end do
        end if

        ! right bounce-back boundary, written into the halo column used by pull streaming
        if (at_right_boundary) then
            do y = 1, n_y_local
                f(4, n_x_local+1, y) = f(2, n_x_local, y)
                f(7, n_x_local+1, y-1) = f(9, n_x_local, y)
                f(8, n_x_local+1, y+1) = f(6, n_x_local, y)
            end do
        end if

        ! bottom bounce-back boundary, written into the halo row used by pull streaming
        if (at_bottom_boundary) then
            do x = 1, n_x_local
                f(3, x, 0) = f(5, x, 1)
                f(6, x-1, 0) = f(8, x, 1)
                f(7, x+1, 0) = f(9, x, 1)
            end do
        end if

        ! top bounce-back boundary, written into the halo row used by pull streaming
        if (at_top_boundary) then
            moving_wall_correction_8 = 6.0_FP * W(6) * rho_0 * u_wall
            moving_wall_correction_9 = 6.0_FP * W(7) * rho_0 * u_wall

            do x = 1, n_x_local
                f(5, x, n_y_local+1) = f(3, x, n_y_local)
                f(8, x+1, n_y_local+1) = f(6, x, n_y_local) - moving_wall_correction_8
                f(9, x-1, n_y_local+1) = f(7, x, n_y_local) + moving_wall_correction_9
            end do
        end if
    end subroutine prepare_sliding_lid_halos_SL


end module simulation
