module poiseuille_flow
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_DIRS, C_X_FP, C_Y_FP, W, FP
    implicit none

contains

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
