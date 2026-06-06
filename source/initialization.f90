module initialization
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_Y, N_DIRS, C_X_FP, C_Y_FP, W, FP, PI
    implicit none

contains

    subroutine apply_condition_shear_wave_local( &
        rho_0, u_max, n_sin, n_x_local, n_y_local, y_global_start, f, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: u_max
        real(FP), intent(in) :: n_sin
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        integer(int32), intent(in) :: y_global_start

        ! write destinations
        real(FP), intent(out) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(out) :: rho(n_x_local, n_y_local)
        real(FP), intent(out) :: u_x(n_x_local, n_y_local)
        real(FP), intent(out) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y, i
        integer(int32) :: y_global
        real(FP) :: u_x_val
        real(FP) :: u_y_val
        real(FP) :: u_squ
        real(FP) :: c_dot_u
        real(FP) :: f_eq_val

        f = 0.0_FP

        ! loop over all image-owned rows
        do y = 1, n_y_local

            ! shear wave velocity for this global row
            y_global = y_global_start + y - 1
            u_x_val = u_max * sin((2.0_FP * PI * real(n_sin, FP) * real(y_global - 1, FP)) / real(N_Y, FP))
            u_y_val = 0.0_FP
            u_squ = u_x_val * u_x_val

            ! loop over all image-owned cols
            !DIR$ SIMD
            do x = 1, n_x_local

                rho(x, y) = rho_0
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val

                ! init distribution functions
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for channel i
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_0 * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)

                    ! write to destination channel
                    f(x, y, i) = f_eq_val
                end do
            end do
        end do
    end subroutine apply_condition_shear_wave_local


    subroutine apply_condition_couette_flow_local( &
        rho_0, n_x_local, n_y_local, f, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: rho_0
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local

        ! write destinations
        real(FP), intent(out) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(out) :: rho(n_x_local, n_y_local)
        real(FP), intent(out) :: u_x(n_x_local, n_y_local)
        real(FP), intent(out) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_eq(N_DIRS)

        f = 0.0_FP

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = W(i) * rho_0
        end do

        ! loop over all image-owned cells
        do y = 1, n_y_local
            !DIR$ SIMD
            do x = 1, n_x_local

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_FP
                u_y(x, y) = 0.0_FP

                ! init distribution functions
                do i = 1, N_DIRS
                    f(x, y, i) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_couette_flow_local


    subroutine apply_condition_poiseuille_flow_local( &
        rho_0, n_x_local, n_y_local, f, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: rho_0
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local

        ! write destinations
        real(FP), intent(out) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(out) :: rho(n_x_local, n_y_local)
        real(FP), intent(out) :: u_x(n_x_local, n_y_local)
        real(FP), intent(out) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_eq(N_DIRS)

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = W(i) * rho_0
        end do

        ! init pressure halo values
        do i = 1, N_DIRS
            f(:, :, i) = f_eq(i)
        end do

        ! loop over all image-owned cells
        do y = 1, n_y_local
            !DIR$ SIMD
            do x = 1, n_x_local

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_FP
                u_y(x, y) = 0.0_FP

                ! init distribution functions
                do i = 1, N_DIRS
                    f(x, y, i) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_poiseuille_flow_local


    subroutine apply_condition_sliding_lid_local( &
        rho_0, n_x_local, n_y_local, f, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: rho_0
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local

        ! write destinations
        real(FP), intent(out) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(out) :: rho(n_x_local, n_y_local)
        real(FP), intent(out) :: u_x(n_x_local, n_y_local)
        real(FP), intent(out) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_eq(N_DIRS)

        f = 0.0_FP

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = W(i) * rho_0
        end do

        ! loop over all image-owned cells
        do y = 1, n_y_local
            !DIR$ SIMD
            do x = 1, n_x_local

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_FP
                u_y(x, y) = 0.0_FP

                ! init distribution functions
                do i = 1, N_DIRS
                    f(x, y, i) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_sliding_lid_local


end module initialization
