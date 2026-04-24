module simulation
    ! imports
    use iso_fortran_env, only: int32, real32
    implicit none

contains

    subroutine fuzed_pull_streaming_collision_shear_wave_decay( &
        N_X, N_Y, N_DIRS, c_x, c_y, c_x_fp, c_y_fp, w, omega, &
        f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: N_X
        integer(int32), intent(in) :: N_Y
        integer(int32), intent(in) :: N_DIRS
        integer(int32), intent(in) :: c_x(:)
        integer(int32), intent(in) :: c_y(:)
        real(real32), intent(in) :: c_x_fp(:)
        real(real32), intent(in) :: c_y_fp(:)
        real(real32), intent(in) :: w(:)
        real(real32), intent(in) :: omega
        real(real32), intent(in) :: f(:, :, :)
        logical, intent(in) :: write_rho
        logical, intent(in) :: write_u_x
        logical, intent(in) :: write_u_y

        ! write destinations
        real(real32), intent(out) :: f_next(:, :, :)

        ! optional write densitions
        real(real32), intent(inout) :: rho(:,:)
        real(real32), intent(inout) :: u_x(:,:)
        real(real32), intent(inout) :: u_y(:,:)

        ! temp
        integer(int32) :: x, y, i
        integer(int32) :: src_x, src_y
        real(real32) :: f_pulled(N_DIRS)
        real(real32) :: rho_val
        real(real32) :: u_x_val
        real(real32) :: u_y_val
        real(real32) :: u_squ
        real(real32) :: c_dot_u
        real(real32) :: f_eq
        real(real32) :: f_next_val

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho_val = 0.0_real32
                u_x_val = 0.0_real32
                u_y_val = 0.0_real32

                ! ---------
                ! | 7 3 6 |
                ! | 4 1 2 |
                ! | 8 5 9 |
                ! ---------
                ! pull streamed distribution functions from source cells in all dirs
                do i = 1, N_DIRS
                    src_x = modulo((x - 1) - c_x(i), N_X) + 1
                    src_y = modulo((y - 1) - c_y(i), N_Y) + 1

                    f_pulled(i) = f(i, src_x, src_y)

                    rho_val = rho_val + f_pulled(i)
                    u_x_val = u_x_val + f_pulled(i) * c_x_fp(i)
                    u_y_val = u_y_val + f_pulled(i) * c_y_fp(i)
                end do

                ! safety check
                if (rho_val <= 0.0_real32) then
                    error stop "density is zero in collision/streaming step (rho_val <= 0)"
                end if

                ! finalize velocity
                u_x_val = u_x_val / rho_val
                u_y_val = u_y_val / rho_val
                u_squ = u_x_val * u_x_val + u_y_val * u_y_val

                ! store density and velocity values only if required
                if (write_rho) then
                    rho(x, y) = rho_val
                end if
                if (write_u_x) then
                    u_x(x, y) = u_x_val
                end if
                if (write_u_y) then
                    u_y(x, y) = u_y_val
                end if

                ! collide and stream to destination cells in all dirs
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for dir i
                    c_dot_u = c_x_fp(i) * u_x_val + c_y_fp(i) * u_y_val
                    f_eq = w(i) * rho_val * ( &
                        1.0_real32 + &
                        3.0_real32 * c_dot_u + &
                        4.5_real32 * c_dot_u * c_dot_u - &
                        1.5_real32 * u_squ)

                    ! relax towards equilibrium
                    f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq)

                    ! write to destination dir i of this cell in next distribution function buffer
                    f_next(i, x, y) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_shear_wave_decay


    subroutine swap_distribution_function_buffers( &
        f, f_next &
        )
        ! read/write inputs
        real(real32), allocatable, intent(inout) :: f(:, :, :)
        real(real32), allocatable, intent(inout) :: f_next(:, :, :)
        real(real32), allocatable :: temp(:, :, :)

        ! swap ownership
        call move_alloc(f, temp)
        call move_alloc(f_next, f)
        call move_alloc(temp, f_next)
    end subroutine swap_distribution_function_buffers

end module simulation
