module sliding_lid
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, FP
    implicit none

contains

    subroutine fuzed_pull_streaming_collision_local_SL( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        write_macro_fields, rho_0, omega, u_wall, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: u_wall
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
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
                !DIR$ UNROLL(9)
                do i = 1, N_DIRS

                    src_x = x - C_X(i)
                    src_y = y - C_Y(i)

                    ! bounce-back for global left boundary (static)
                    if (src_x < 1 .and. at_left_boundary) then
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

                    ! bounce-back for global right boundary (static)
                    else if (src_x > n_x_local .and. at_right_boundary) then
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

                    ! bounce-back for global bottom boundary (static)
                    else if (src_y < 1 .and. at_bottom_boundary) then
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

                    ! bounce-back for global top boundary (moving)
                    else if (src_y > n_y_local .and. at_top_boundary) then
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

                    ! no boundary
                    else
                        f_pulled(i) = f(src_x, src_y, i)
                    end if

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
                !DIR$ UNROLL(9)
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for channels i
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_val * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)

                    ! relax towards equilibrium and write to destination channel
                    f_next_val = f_pulled(i) + omega * (f_eq_val - f_pulled(i))
                    f_next(x, y, i) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_local_SL


    subroutine fuzed_pull_streaming_collision_local_unrolled_SL( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        write_macro_fields, rho_0, omega, u_wall, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: u_wall
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y
        logical :: bottom_wall_row
        logical :: top_wall_row
        logical :: left_wall_col
        logical :: right_wall_col
        real(FP) :: moving_wall_correction_8
        real(FP) :: moving_wall_correction_9
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

        moving_wall_correction_8 = 6.0_FP * W(6) * rho_0 * u_wall
        moving_wall_correction_9 = 6.0_FP * W(7) * rho_0 * u_wall

        ! loop over all image-owned cells
        do y = 1, n_y_local

            bottom_wall_row = at_bottom_boundary .and. y == 1
            top_wall_row = at_top_boundary .and. y == n_y_local

            do x = 1, n_x_local

                left_wall_col = at_left_boundary .and. x == 1
                right_wall_col = at_right_boundary .and. x == n_x_local

                ! ---------
                ! | 7 3 6 |
                ! | 4 1 2 |
                ! | 8 5 9 |
                ! ---------
                ! pull streamed distribution functions from source cells
                f_1 = f(x, y, 1)

                if (left_wall_col) then
                    ! bounce-back for global left boundary (static)
                    f_2 = f(x, y, 4)
                    f_6 = f(x, y, 8)
                    f_9 = f(x, y, 7)
                else
                    f_2 = f(x - 1, y, 2)

                    if (bottom_wall_row) then
                        ! bounce-back for global bottom boundary (static)
                        f_6 = f(x, y, 8)
                    else
                        f_6 = f(x - 1, y - 1, 6)
                    end if

                    if (top_wall_row) then
                        ! bounce-back for global top boundary (moving)
                        f_9 = f(x, y, 7) + moving_wall_correction_9
                    else
                        f_9 = f(x - 1, y + 1, 9)
                    end if
                end if

                if (right_wall_col) then
                    ! bounce-back for global right boundary (static)
                    f_4 = f(x, y, 2)
                    f_7 = f(x, y, 9)
                    f_8 = f(x, y, 6)
                else
                    f_4 = f(x + 1, y, 4)

                    if (bottom_wall_row) then
                        ! bounce-back for global bottom boundary (static)
                        f_7 = f(x, y, 9)
                    else
                        f_7 = f(x + 1, y - 1, 7)
                    end if

                    if (top_wall_row) then
                        ! bounce-back for global top boundary (moving)
                        f_8 = f(x, y, 6) - moving_wall_correction_8
                    else
                        f_8 = f(x + 1, y + 1, 8)
                    end if
                end if

                if (bottom_wall_row) then
                    ! bounce-back for global bottom boundary (static)
                    f_3 = f(x, y, 5)
                else
                    f_3 = f(x, y - 1, 3)
                end if

                if (top_wall_row) then
                    ! bounce-back for global top boundary (moving)
                    f_5 = f(x, y, 3)
                else
                    f_5 = f(x, y + 1, 5)
                end if

                rho_val = f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8 + f_9
                u_x_val = f_2 - f_4 + f_6 - f_7 - f_8 + f_9
                u_y_val = f_3 - f_5 + f_6 + f_7 - f_8 - f_9

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
                ! 1: (0, 0)
                f_next(x, y, 1) = f_1 + omega * ((4.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP - 1.5_FP * u_squ) - f_1)

                ! 2: (1, 0)
                f_next(x, y, 2) = f_2 + omega * ((1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * u_x_val + 4.5_FP * u_x_val * u_x_val - &
                    1.5_FP * u_squ) - f_2)

                ! 3: (0, 1)
                f_next(x, y, 3) = f_3 + omega * ((1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * u_y_val + 4.5_FP * u_y_val * u_y_val - &
                    1.5_FP * u_squ) - f_3)

                ! 4: (-1, 0)
                f_next(x, y, 4) = f_4 + omega * ((1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP - 3.0_FP * u_x_val + 4.5_FP * u_x_val * u_x_val - &
                    1.5_FP * u_squ) - f_4)

                ! 5: (0, -1)
                f_next(x, y, 5) = f_5 + omega * ((1.0_FP/9.0_FP) * rho_val * ( &
                    1.0_FP - 3.0_FP * u_y_val + 4.5_FP * u_y_val * u_y_val - &
                    1.5_FP * u_squ) - f_5)

                ! 6: (1, 1)
                f_next(x, y, 6) = f_6 + omega * ((1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * (u_x_val + u_y_val) + &
                    4.5_FP * (u_x_val + u_y_val) * (u_x_val + u_y_val) - &
                    1.5_FP * u_squ) - f_6)

                ! 7: (-1, 1)
                f_next(x, y, 7) = f_7 + omega * ((1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * (-u_x_val + u_y_val) + &
                    4.5_FP * (-u_x_val + u_y_val) * (-u_x_val + u_y_val) - &
                    1.5_FP * u_squ) - f_7)

                ! 8: (-1, -1)
                f_next(x, y, 8) = f_8 + omega * ((1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP - 3.0_FP * (u_x_val + u_y_val) + &
                    4.5_FP * (u_x_val + u_y_val) * (u_x_val + u_y_val) - &
                    1.5_FP * u_squ) - f_8)

                ! 9: (1, -1)
                f_next(x, y, 9) = f_9 + omega * ((1.0_FP/36.0_FP) * rho_val * ( &
                    1.0_FP + 3.0_FP * (u_x_val - u_y_val) + &
                    4.5_FP * (u_x_val - u_y_val) * (u_x_val - u_y_val) - &
                    1.5_FP * u_squ) - f_9)
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_local_unrolled_SL


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
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

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
                f(0, y, 2) = f(1, y, 4)
                f(0, y-1, 6) = f(1, y, 8)
                f(0, y+1, 9) = f(1, y, 7)
            end do
        end if

        ! right bounce-back boundary, written into the halo column used by pull streaming
        if (at_right_boundary) then
            do y = 1, n_y_local
                f(n_x_local+1, y, 4) = f(n_x_local, y, 2)
                f(n_x_local+1, y-1, 7) = f(n_x_local, y, 9)
                f(n_x_local+1, y+1, 8) = f(n_x_local, y, 6)
            end do
        end if

        ! bottom bounce-back boundary, written into the halo row used by pull streaming
        if (at_bottom_boundary) then
            do x = 1, n_x_local
                f(x, 0, 3) = f(x, 1, 5)
                f(x-1, 0, 6) = f(x, 1, 8)
                f(x+1, 0, 7) = f(x, 1, 9)
            end do
        end if

        ! top bounce-back boundary, written into the halo row used by pull streaming
        if (at_top_boundary) then
            moving_wall_correction_8 = 6.0_FP * W(6) * rho_0 * u_wall
            moving_wall_correction_9 = 6.0_FP * W(7) * rho_0 * u_wall

            do x = 1, n_x_local
                f(x, n_y_local+1, 5) = f(x, n_y_local, 3)
                f(x+1, n_y_local+1, 8) = f(x, n_y_local, 6) - moving_wall_correction_8
                f(x-1, n_y_local+1, 9) = f(x, n_y_local, 7) + moving_wall_correction_9
            end do
        end if
    end subroutine prepare_sliding_lid_halos_SL


end module sliding_lid
