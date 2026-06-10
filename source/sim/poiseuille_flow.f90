module poiseuille_flow
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, FP
    implicit none

contains

    subroutine fuzed_pull_streaming_collision_local_PF( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        write_macro_fields, omega, rho_in, rho_out, f, f_next, rho, u_x, u_y, &
        macro_left, macro_right, macro_send_left, macro_send_right &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: rho_in
        real(FP), intent(in) :: rho_out
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(in) :: macro_left(n_y_local, 3)
        real(FP), intent(in) :: macro_right(n_y_local, 3)

        ! write destinations
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)
        real(FP), intent(inout) :: macro_send_left(n_y_local, 3)
        real(FP), intent(inout) :: macro_send_right(n_y_local, 3)

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
                ! (periodic left/right boundaries handled by wrapped halo exchange)
                !DIR$ UNROLL(9)
                do i = 1, N_DIRS

                    src_x = x - C_X(i)
                    src_y = y - C_Y(i)

                    ! no boundary
                    if (src_x >= 1 .and. src_x <= n_x_local .and. &
                        src_y >= 1 .and. src_y <= n_y_local) then
                        f_pulled(i) = f(src_x, src_y, i)

                    ! pressure-periodic inlet for left boundary
                    else if (src_x < 1 .and. at_left_boundary) then
                        f_pulled(i) = pressure_periodic_distribution( &
                            i, f(0, y, i), macro_left(y, 1), macro_left(y, 2), macro_left(y, 3), rho_in)

                    ! pressure-periodic outlet for right boundary
                    else if (src_x > n_x_local .and. at_right_boundary) then
                        f_pulled(i) = pressure_periodic_distribution( &
                            i, f(n_x_local+1, y, i), macro_right(y, 1), macro_right(y, 2), macro_right(y, 3), rho_out)

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
                            error stop "error: invalid bottom boundary channel in poiseuille flow"
                    #endif
                        end select

                    ! bounce-back for global top boundary (static)
                    else if (src_y > n_y_local .and. at_top_boundary) then
                        select case(i)
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

                if (at_left_boundary .and. x == 1) then
                    macro_send_left(y, 1) = rho_val
                    macro_send_left(y, 2) = u_x_val
                    macro_send_left(y, 3) = u_y_val
                end if

                if (at_right_boundary .and. x == n_x_local) then
                    macro_send_right(y, 1) = rho_val
                    macro_send_right(y, 2) = u_x_val
                    macro_send_right(y, 3) = u_y_val
                end if

                if (write_macro_fields) then
                    rho(x, y) = rho_val
                    u_x(x, y) = u_x_val
                    u_y(x, y) = u_y_val
                end if

                ! collide and stream locally
                !DIR$ UNROLL(9)
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
                    f_next(x, y, i) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_local_PF


    subroutine fuzed_pull_streaming_collision_local_unrolled_PF( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        write_macro_fields, omega, rho_in, rho_out, f, f_next, rho, u_x, u_y, &
        macro_left, macro_right, macro_send_left, macro_send_right &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: rho_in
        real(FP), intent(in) :: rho_out
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(in) :: macro_left(n_y_local, 3)
        real(FP), intent(in) :: macro_right(n_y_local, 3)

        ! write destinations
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)
        real(FP), intent(inout) :: macro_send_left(n_y_local, 3)
        real(FP), intent(inout) :: macro_send_right(n_y_local, 3)

        ! temp
        integer(int32) :: x, y
        logical :: left_pressure_cell
        logical :: right_pressure_cell
        logical :: bottom_wall_row
        logical :: top_wall_row
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
        real(FP) :: rho_src
        real(FP) :: u_x_src
        real(FP) :: u_y_src
        real(FP) :: u_squ_src
        real(FP) :: c_dot_u_src
        real(FP) :: pressure_factor

        ! loop over all image-owned cells
        do y = 1, n_y_local

            bottom_wall_row = at_bottom_boundary .and. y == 1
            top_wall_row = at_top_boundary .and. y == n_y_local

            do x = 1, n_x_local

                left_pressure_cell = at_left_boundary .and. x == 1
                right_pressure_cell = at_right_boundary .and. x == n_x_local

                ! ---------
                ! | 7 3 6 |
                ! | 4 1 2 |
                ! | 8 5 9 |
                ! ---------
                ! pull streamed distribution functions from source cells
                ! (periodic left/right boundaries handled by wrapped halo exchange)
                f_1 = f(x, y, 1)

                if (left_pressure_cell) then
                    rho_src = macro_left(y, 1)
                    u_x_src = macro_left(y, 2)
                    u_y_src = macro_left(y, 3)
                    u_squ_src = u_x_src * u_x_src + u_y_src * u_y_src

                    ! 2: (1, 0)
                    c_dot_u_src = u_x_src
                    pressure_factor = 1.0_FP + 3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - 1.5_FP * u_squ_src
                    f_2 = f(0, y, 2) + W(2) * (rho_in - rho_src) * pressure_factor

                    ! 6: (1, 1)
                    c_dot_u_src = u_x_src + u_y_src
                    pressure_factor = 1.0_FP + 3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - 1.5_FP * u_squ_src
                    f_6 = f(0, y, 6) + W(6) * (rho_in - rho_src) * pressure_factor

                    ! 9: (1, -1)
                    c_dot_u_src = u_x_src - u_y_src
                    pressure_factor = 1.0_FP + 3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - 1.5_FP * u_squ_src
                    f_9 = f(0, y, 9) + W(9) * (rho_in - rho_src) * pressure_factor
                else
                    f_2 = f(x - 1, y, 2)

                    if (bottom_wall_row) then
                        ! bounce-back for global bottom boundary (static)
                        f_6 = f(x, y, 8)
                    else
                        f_6 = f(x - 1, y - 1, 6)
                    end if

                    if (top_wall_row) then
                        ! bounce-back for global top boundary (static)
                        f_9 = f(x, y, 7)
                    else
                        f_9 = f(x - 1, y + 1, 9)
                    end if
                end if

                if (right_pressure_cell) then
                    rho_src = macro_right(y, 1)
                    u_x_src = macro_right(y, 2)
                    u_y_src = macro_right(y, 3)
                    u_squ_src = u_x_src * u_x_src + u_y_src * u_y_src

                    ! 4: (-1, 0)
                    c_dot_u_src = -u_x_src
                    pressure_factor = 1.0_FP + 3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - 1.5_FP * u_squ_src
                    f_4 = f(n_x_local+1, y, 4) + W(4) * (rho_out - rho_src) * pressure_factor

                    ! 7: (-1, 1)
                    c_dot_u_src = -u_x_src + u_y_src
                    pressure_factor = 1.0_FP + 3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - 1.5_FP * u_squ_src
                    f_7 = f(n_x_local+1, y, 7) + W(7) * (rho_out - rho_src) * pressure_factor

                    ! 8: (-1, -1)
                    c_dot_u_src = -u_x_src - u_y_src
                    pressure_factor = 1.0_FP + 3.0_FP * c_dot_u_src + &
                        4.5_FP * c_dot_u_src * c_dot_u_src - 1.5_FP * u_squ_src
                    f_8 = f(n_x_local+1, y, 8) + W(8) * (rho_out - rho_src) * pressure_factor
                else
                    f_4 = f(x + 1, y, 4)

                    if (bottom_wall_row) then
                        ! bounce-back for global bottom boundary (static)
                        f_7 = f(x, y, 9)
                    else
                        f_7 = f(x + 1, y - 1, 7)
                    end if

                    if (top_wall_row) then
                        ! bounce-back for global top boundary (static)
                        f_8 = f(x, y, 6)
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
                    ! bounce-back for global top boundary (static)
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

                if (left_pressure_cell) then
                    macro_send_left(y, 1) = rho_val
                    macro_send_left(y, 2) = u_x_val
                    macro_send_left(y, 3) = u_y_val
                end if

                if (right_pressure_cell) then
                    macro_send_right(y, 1) = rho_val
                    macro_send_right(y, 2) = u_x_val
                    macro_send_right(y, 3) = u_y_val
                end if

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
    end subroutine fuzed_pull_streaming_collision_local_unrolled_PF


    subroutine prepare_poiseuille_flow_halos_PF( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        rho_in, rho_out, f, macro_left, macro_right &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        real(FP), intent(in) :: rho_in
        real(FP), intent(in) :: rho_out
        real(FP), intent(in) :: macro_left(n_y_local, 3)
        real(FP), intent(in) :: macro_right(n_y_local, 3)

        ! read/write inputs
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! temp
        integer(int32) :: x, y

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
            do x = 1, n_x_local
                f(x, n_y_local+1, 5) = f(x, n_y_local, 3)
                f(x+1, n_y_local+1, 8) = f(x, n_y_local, 6)
                f(x-1, n_y_local+1, 9) = f(x, n_y_local, 7)
            end do
        end if

        ! pressure-periodic inlet, written into the halo column used by pull streaming
        if (at_left_boundary) then
            do y = 1, n_y_local
                f(0, y, 2) = pressure_periodic_distribution( &
                    2_int32, f(0, y, 2), macro_left(y, 1), macro_left(y, 2), macro_left(y, 3), rho_in)
                f(0, y, 6) = pressure_periodic_distribution( &
                    6_int32, f(0, y, 6), macro_left(y, 1), macro_left(y, 2), macro_left(y, 3), rho_in)
                f(0, y, 9) = pressure_periodic_distribution( &
                    9_int32, f(0, y, 9), macro_left(y, 1), macro_left(y, 2), macro_left(y, 3), rho_in)
            end do

            f(0, 0, 6) = pressure_periodic_distribution( &
                6_int32, f(0, 0, 6), macro_left(1, 1), macro_left(1, 2), macro_left(1, 3), rho_in)
            f(0, n_y_local+1, 9) = pressure_periodic_distribution( &
                9_int32, f(0, n_y_local+1, 9), macro_left(n_y_local, 1), &
                macro_left(n_y_local, 2), macro_left(n_y_local, 3), rho_in)
        end if

        ! pressure-periodic outlet, written into the halo column used by pull streaming
        if (at_right_boundary) then
            do y = 1, n_y_local
                f(n_x_local+1, y, 4) = pressure_periodic_distribution( &
                    4_int32, f(n_x_local+1, y, 4), macro_right(y, 1), macro_right(y, 2), macro_right(y, 3), rho_out)
                f(n_x_local+1, y, 7) = pressure_periodic_distribution( &
                    7_int32, f(n_x_local+1, y, 7), macro_right(y, 1), macro_right(y, 2), macro_right(y, 3), rho_out)
                f(n_x_local+1, y, 8) = pressure_periodic_distribution( &
                    8_int32, f(n_x_local+1, y, 8), macro_right(y, 1), macro_right(y, 2), macro_right(y, 3), rho_out)
            end do

            f(n_x_local+1, 0, 7) = pressure_periodic_distribution( &
                7_int32, f(n_x_local+1, 0, 7), macro_right(1, 1), macro_right(1, 2), macro_right(1, 3), rho_out)
            f(n_x_local+1, n_y_local+1, 8) = pressure_periodic_distribution( &
                8_int32, f(n_x_local+1, n_y_local+1, 8), macro_right(n_y_local, 1), &
                macro_right(n_y_local, 2), macro_right(n_y_local, 3), rho_out)
        end if
    end subroutine prepare_poiseuille_flow_halos_PF


    subroutine update_poiseuille_flow_macro_strips_PF( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, f, macro_send_left, macro_send_right &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: macro_send_left(n_y_local, 3)
        real(FP), intent(inout) :: macro_send_right(n_y_local, 3)

        if (at_left_boundary) then
            call update_macro_strip(1_int32, macro_send_left)
        end if

        if (at_right_boundary) then
            call update_macro_strip(n_x_local, macro_send_right)
        end if

    contains

        subroutine update_macro_strip( &
            x, macro_strip &
            )
            ! inputs
            integer(int32), intent(in) :: x

            ! write destinations
            real(FP), intent(inout) :: macro_strip(n_y_local, 3)

            ! temp
            integer(int32) :: y
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

            do y = 1, n_y_local
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

                ! debug check
            #ifdef FFB_DENSITY_CHECKS
                if (rho_val <= 0.0_FP) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
                end if
            #endif

                macro_strip(y, 1) = rho_val
                macro_strip(y, 2) = u_x_val / rho_val
                macro_strip(y, 3) = u_y_val / rho_val
            end do
        end subroutine update_macro_strip
    end subroutine update_poiseuille_flow_macro_strips_PF


    pure function pressure_periodic_distribution( &
        i, f_src, rho_src, u_x_src, u_y_src, rho_boundary &
        ) result(f_boundary)
        ! inputs
        integer(int32), intent(in) :: i
        real(FP), intent(in) :: f_src
        real(FP), intent(in) :: rho_src
        real(FP), intent(in) :: u_x_src
        real(FP), intent(in) :: u_y_src
        real(FP), intent(in) :: rho_boundary

        ! output
        real(FP) :: f_boundary

        ! temp
        real(FP) :: u_squ_src
        real(FP) :: c_dot_u_src
        real(FP) :: pressure_factor

        u_squ_src = u_x_src * u_x_src + u_y_src * u_y_src
        c_dot_u_src = C_X_FP(i) * u_x_src + C_Y_FP(i) * u_y_src
        pressure_factor = 1.0_FP + 3.0_FP * c_dot_u_src + &
            4.5_FP * c_dot_u_src * c_dot_u_src - 1.5_FP * u_squ_src

        f_boundary = f_src + W(i) * (rho_boundary - rho_src) * pressure_factor
    end function pressure_periodic_distribution


end module poiseuille_flow
