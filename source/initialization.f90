module initialization
    ! imports
    use iso_fortran_env, only: int32, int64, real32, real64
    implicit none

contains

    subroutine apply_condition_shear_wave_decay(N_X, N_Y, N_DIR, c_x_fp, c_y_fp, w, rho_0, u_max, k, f, rho, u_x, u_y)
        ! read-only inputs
        integer(int32), intent(in) :: N_X
        integer(int32), intent(in) :: N_Y
        integer(int32), intent(in) :: N_DIR
        real(real32), intent(in) :: c_x_fp(:)
        real(real32), intent(in) :: c_y_fp(:)
        real(real32), intent(in) :: w(:)
        real(real32), intent(in) :: rho_0
        real(real32), intent(in) :: u_max
        real(real32), intent(in) :: k

        ! write destinations
        real(real32), intent(out) :: f(:, :, :)
        real(real32), intent(out) :: rho(:,:)
        real(real32), intent(out) :: u_x(:,:)
        real(real32), intent(out) :: u_y(:,:)

        ! temp
        integer(int32) :: x, y, i
        real(real32) :: u_x_val
        real(real32) :: u_y_val
        real(real32) :: u_squ
        real(real32) :: c_dot_u

        ! outer loop over rows
        do y = 1, N_Y
            ! shear wave velocity for this row
            u_x_val = u_max * sin(k * real(y - 1, real32))
            u_y_val = 0.0_real32
            u_squ = u_x_val**2

            ! inner loop over columns
            do x = 1, N_X
                rho(x, y) = rho_0
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val

                ! init distribution functions based on equilibrium distribution in each direction
                do i = 1, N_DIR
                    c_dot_u = c_x_fp(i) * u_x_val + c_y_fp(i) * u_y_val

                    f(i, x, y) = w(i) * rho_0 * ( &
                        1.0_real32 + &
                        3.0_real32 * c_dot_u + &
                        4.5_real32 * c_dot_u**2 - &
                        1.5_real32 * u_squ &
                    )
                end do
            end do
        end do
    end subroutine apply_condition_shear_wave_decay

end module initialization
