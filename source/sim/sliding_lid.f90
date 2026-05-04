module sliding_lid
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_X, N_Y, N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, FP
    implicit none

contains

    subroutine fuzed_pull_streaming_collision_outer_SL( &
        write_macro_fields, rho_0, omega, u_wall, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        logical, intent(in) :: write_macro_fields
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
            call stream_collide_outer_cell_SL(x, y)
        end do

        ! top row
        y = N_Y
        do x = 1, N_X
            call stream_collide_outer_cell_SL(x, y)
        end do

        ! left col (no corners)
        x = 1
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SL(x, y)
        end do

        ! right col (no corners)
        x = N_X
        do y = 2, N_Y - 1
            call stream_collide_outer_cell_SL(x, y)
        end do

    contains ! helper subroutine

        subroutine stream_collide_outer_cell_SL( &
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

                ! compute equilibrium distribution function for channels i
                c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                f_eq_val = W(i) * rho_val * ( &
                    1.0_FP + &
                    3.0_FP * c_dot_u + &
                    4.5_FP * c_dot_u * c_dot_u - &
                    1.5_FP * u_squ)

                ! relax towards equilibrium and write to destination channel in this cell
                f_next_val = f_pulled(i) + omega * (f_eq_val - f_pulled(i))
                f_next(x, y, i) = f_next_val
            end do
         end subroutine stream_collide_outer_cell_SL
    end subroutine fuzed_pull_streaming_collision_outer_SL

end module sliding_lid
