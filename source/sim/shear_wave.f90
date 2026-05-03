module shear_wave
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_X, N_Y, N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, FP
    implicit none

contains

    subroutine fuzed_pull_streaming_collision_outer_SW( &
        write_macro_fields, omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        logical, intent(in) :: write_macro_fields
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
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

    contains ! helper subroutine

        subroutine stream_collide_outer_cell_SW( &
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

            ! finalize density and velocity
            u_x_val = u_x_val / rho_val
            u_y_val = u_y_val / rho_val
            u_squ = u_x_val * u_x_val + u_y_val * u_y_val

            if (write_macro_fields) then
                rho(x, y) = rho_val
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val
            end if

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
        end subroutine stream_collide_outer_cell_SW
    end subroutine fuzed_pull_streaming_collision_outer_SW


    subroutine fuzed_unrolled_pull_streaming_collision_outer_SW( &
        write_macro_fields, omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        logical, intent(in) :: write_macro_fields
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
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

    contains ! helper subroutine

        subroutine stream_collide_outer_cell_SW( &
            x, y &
            )
            ! inputs
            integer(int32), intent(in) :: x
            integer(int32), intent(in) :: y

            ! temp
            integer(int32) :: src_x_minus
            integer(int32) :: src_x_plus
            integer(int32) :: src_y_minus
            integer(int32) :: src_y_plus
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
            ! periodic source cell indices for left/right boundary
            src_x_minus = x - 1
            src_x_plus = x + 1

            if (src_x_minus < 1) then
                src_x_minus = N_X
            end if

            if (src_x_plus > N_X) then
                src_x_plus = 1
            end if

            ! periodic source cell indices for bottom/top boundary
            src_y_minus = y - 1
            src_y_plus = y + 1

            if (src_y_minus < 1) then
                src_y_minus = N_Y
            end if

            if (src_y_plus > N_Y) then
                src_y_plus = 1
            end if

            ! pull streamed distribution functions from source cells in all channels
            ! (periodic boundary handling, manually unrolled)
            f_1 = f(x, y, 1)
            f_2 = f(src_x_minus, y, 2)
            f_3 = f(x, src_y_minus, 3)
            f_4 = f(src_x_plus, y, 4)
            f_5 = f(x, src_y_plus, 5)
            f_6 = f(src_x_minus, src_y_minus, 6)
            f_7 = f(src_x_plus, src_y_minus, 7)
            f_8 = f(src_x_plus, src_y_plus, 8)
            f_9 = f(src_x_minus, src_y_plus, 9)

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

            ! collide and stream locally to destination channels
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
        end subroutine stream_collide_outer_cell_SW
    end subroutine fuzed_unrolled_pull_streaming_collision_outer_SW


    subroutine fuzed_push_streaming_collision_outer_SW( &
        omega, f, f_next &
        )
        ! inputs
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_X, N_Y, N_DIRS)

        ! temp
        integer(int32) :: x, y

        ! bottom row
        y = 1
        do x = 1, N_X
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

    contains ! helper subroutine

        subroutine stream_collide_outer_cell_SW( &
            x, y &
            )
            ! inputs
            integer(int32), intent(in) :: x
            integer(int32), intent(in) :: y

            ! temp
            integer(int32) :: i
            integer(int32) :: dst_x, dst_y
            real(FP) :: f_local(N_DIRS)
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
            ! read local distribution functions from this source cell in all channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                f_local(i) = f(x, y, i)

                rho_val = rho_val + f_local(i)
                u_x_val = u_x_val + f_local(i) * C_X_FP(i)
                u_y_val = u_y_val + f_local(i) * C_Y_FP(i)
            end do

            ! safety check to avoid division by zero in case of wrong density
        #ifdef FFB_DENSITY_CHECKS
            if (rho_val <= 0.0_FP) then
                error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
            end if
        #endif

            ! finalize velocity for equilibrium computation
            u_x_val = u_x_val / rho_val
            u_y_val = u_y_val / rho_val
            u_squ = u_x_val * u_x_val + u_y_val * u_y_val

            ! collide locally and push-stream to destination cells in all channels
            !DIR$ UNROLL(9)
            do i = 1, N_DIRS

                ! compute equilibrium distribution function for channel i
                c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                f_eq_val = W(i) * rho_val * ( &
                    1.0_FP + &
                    3.0_FP * c_dot_u + &
                    4.5_FP * c_dot_u * c_dot_u - &
                    1.5_FP * u_squ)

                ! periodic push-streaming for left/right boundary
                dst_x = x + C_X(i)
                if (dst_x < 1) then
                    dst_x = N_X
                else if (dst_x > N_X) then
                    dst_x = 1
                end if

                ! periodic push-streaming for bottom/top boundary
                dst_y = y + C_Y(i)
                if (dst_y < 1) then
                    dst_y = N_Y
                else if (dst_y > N_Y) then
                    dst_y = 1
                end if

                ! relax towards equilibrium and write to destination channel in destination cell
                f_next_val = f_local(i) - omega * (f_local(i) - f_eq_val)
                f_next(dst_x, dst_y, i) = f_next_val
            end do
        end subroutine stream_collide_outer_cell_SW
    end subroutine fuzed_push_streaming_collision_outer_SW


    subroutine fuzed_unrolled_push_streaming_collision_outer_SW( &
        omega, f, f_next &
        )
        ! inputs
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(N_X, N_Y, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(N_X, N_Y, N_DIRS)

        ! temp
        integer(int32) :: x, y

        ! bottom row
        y = 1
        do x = 1, N_X
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SW(x, y)
        end do

    contains ! helper subroutine

        subroutine stream_collide_outer_cell_SW( &
            x, y &
            )
            ! inputs
            integer(int32), intent(in) :: x
            integer(int32), intent(in) :: y

            ! temp
            integer(int32) :: dst_x_minus
            integer(int32) :: dst_x_plus
            integer(int32) :: dst_y_minus
            integer(int32) :: dst_y_plus
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
            ! read local distribution functions from this source cell in all channels
            ! (manually unrolled)
            f_1 = f(x, y, 1)
            f_2 = f(x, y, 2)
            f_3 = f(x, y, 3)
            f_4 = f(x, y, 4)
            f_5 = f(x, y, 5)
            f_6 = f(x, y, 6)
            f_7 = f(x, y, 7)
            f_8 = f(x, y, 8)
            f_9 = f(x, y, 9)

            rho_val = f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8 + f_9
            u_x_val = f_2 - f_4 + f_6 - f_7 - f_8 + f_9
            u_y_val = f_3 - f_5 + f_6 + f_7 - f_8 - f_9

            ! safety check to avoid division by zero in case of wrong density
        #ifdef FFB_DENSITY_CHECKS
            if (rho_val <= 0.0_FP) then
                error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
            end if
        #endif

            ! finalize velocity for equilibrium computation
            u_x_val = u_x_val / rho_val
            u_y_val = u_y_val / rho_val
            u_squ = u_x_val * u_x_val + u_y_val * u_y_val

            ! periodic destination cell indices for left/right boundary
            dst_x_minus = x - 1
            dst_x_plus = x + 1

            if (dst_x_minus < 1) then
                dst_x_minus = N_X
            end if

            if (dst_x_plus > N_X) then
                dst_x_plus = 1
            end if

            ! periodic destination cell indices for bottom/top boundary
            dst_y_minus = y - 1
            dst_y_plus = y + 1

            if (dst_y_minus < 1) then
                dst_y_minus = N_Y
            end if

            if (dst_y_plus > N_Y) then
                dst_y_plus = 1
            end if

            ! collide locally and push-stream to destination cells in all channels
            ! (manually unrolled)
            ! 1: (0, 0)
            f_next(x, y, 1) = f_1 - omega * (f_1 - (4.0_FP/9.0_FP) * rho_val * ( &
                1.0_FP - 1.5_FP * u_squ))

            ! 2: (1, 0)
            f_next(dst_x_plus, y, 2) = f_2 - omega * (f_2 - (1.0_FP/9.0_FP) * rho_val * ( &
                1.0_FP + 3.0_FP * u_x_val + 4.5_FP * u_x_val * u_x_val - &
                1.5_FP * u_squ))

            ! 3: (0, 1)
            f_next(x, dst_y_plus, 3) = f_3 - omega * (f_3 - (1.0_FP/9.0_FP) * rho_val * ( &
                1.0_FP + 3.0_FP * u_y_val + 4.5_FP * u_y_val * u_y_val - &
                1.5_FP * u_squ))

            ! 4: (-1, 0)
            f_next(dst_x_minus, y, 4) = f_4 - omega * (f_4 - (1.0_FP/9.0_FP) * rho_val * ( &
                1.0_FP - 3.0_FP * u_x_val + 4.5_FP * u_x_val * u_x_val - &
                1.5_FP * u_squ))

            ! 5: (0, -1)
            f_next(x, dst_y_minus, 5) = f_5 - omega * (f_5 - (1.0_FP/9.0_FP) * rho_val * ( &
                1.0_FP - 3.0_FP * u_y_val + 4.5_FP * u_y_val * u_y_val - &
                1.5_FP * u_squ))

            ! 6: (1, 1)
            f_next(dst_x_plus, dst_y_plus, 6) = f_6 - omega * (f_6 - (1.0_FP/36.0_FP) * rho_val * ( &
                1.0_FP + 3.0_FP * (u_x_val + u_y_val) + &
                4.5_FP * (u_x_val + u_y_val) * (u_x_val + u_y_val) - &
                1.5_FP * u_squ))

            ! 7: (-1, 1)
            f_next(dst_x_minus, dst_y_plus, 7) = f_7 - omega * (f_7 - (1.0_FP/36.0_FP) * rho_val * ( &
                1.0_FP + 3.0_FP * (-u_x_val + u_y_val) + &
                4.5_FP * (-u_x_val + u_y_val) * (-u_x_val + u_y_val) - &
                1.5_FP * u_squ))

            ! 8: (-1, -1)
            f_next(dst_x_minus, dst_y_minus, 8) = f_8 - omega * (f_8 - (1.0_FP/36.0_FP) * rho_val * ( &
                1.0_FP - 3.0_FP * (u_x_val + u_y_val) + &
                4.5_FP * (u_x_val + u_y_val) * (u_x_val + u_y_val) - &
                1.5_FP * u_squ))

            ! 9: (1, -1)
            f_next(dst_x_plus, dst_y_minus, 9) = f_9 - omega * (f_9 - (1.0_FP/36.0_FP) * rho_val * ( &
                1.0_FP + 3.0_FP * (u_x_val - u_y_val) + &
                4.5_FP * (u_x_val - u_y_val) * (u_x_val - u_y_val) - &
                1.5_FP * u_squ))
        end subroutine stream_collide_outer_cell_SW
    end subroutine fuzed_unrolled_push_streaming_collision_outer_SW


    ! TODO: add launcher path that uses this version, controlled by flag?
    subroutine fuzed_push_shift_streaming_collision_full_SW( &
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
    end subroutine fuzed_push_shift_streaming_collision_full_SW

end module shear_wave
