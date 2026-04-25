module initialization
    ! imports
    use iso_fortran_env, only: int32, real32
    use settings, only: SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, &
        shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t
    implicit none

contains

    subroutine initialize_sim_condition( &
        sim_mode, shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
        N_X, N_Y, N_DIRS, pi, c_x_fp, c_y_fp, w, f, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: sim_mode
        type(shear_wave_params_t), intent(in) :: shear_wave_params
        type(couette_flow_params_t), intent(in) :: couette_flow_params
        type(poiseuille_flow_params_t), intent(in) :: poiseuille_flow_params
        type(sliding_lid_params_t), intent(in) :: sliding_lid_params
        integer(int32), intent(in) :: N_X
        integer(int32), intent(in) :: N_Y
        integer(int32), intent(in) :: N_DIRS
        real(real32), intent(in) :: pi
        real(real32), intent(in) :: c_x_fp(:)
        real(real32), intent(in) :: c_y_fp(:)
        real(real32), intent(in) :: w(:)

        ! write destinations
        real(real32), intent(out) :: f(:, :, :)
        real(real32), intent(out) :: rho(:,:)
        real(real32), intent(out) :: u_x(:,:)
        real(real32), intent(out) :: u_y(:,:)

        ! apply initial condition based on selected sim mode
        select case (sim_mode)
        case (SIM_SHEAR_WAVE)
            call apply_condition_shear_wave(N_X, N_Y, N_DIRS, pi, c_x_fp, c_y_fp, w, &
                shear_wave_params%rho_0, shear_wave_params%u_max, shear_wave_params%n_sin, f, rho, u_x, u_y)
        case (SIM_COUETTE_FLOW)
            call apply_condition_couette_flow(N_X, N_Y, N_DIRS, w, &
                couette_flow_params%rho_0, f, rho, u_x, u_y)
        case (SIM_POISEUILLE_FLOW)
            error stop "error: initial condition for poiseuille flow not implemented yet"
            ! TODO: implement
        case (SIM_SLIDING_LID)
            call apply_condition_sliding_lid(N_X, N_Y, N_DIRS, w, &
                sliding_lid_params%rho_0, f, rho, u_x, u_y)
        case default
            error stop "error: unknown sim mode in initialize_sim_condition()"
        end select
    end subroutine initialize_sim_condition


    subroutine apply_condition_shear_wave( &
        N_X, N_Y, N_DIRS, pi, c_x_fp, c_y_fp, w, rho_0, u_max, n_sin, f, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: N_X
        integer(int32), intent(in) :: N_Y
        integer(int32), intent(in) :: N_DIRS
        real(real32), intent(in) :: pi
        real(real32), intent(in) :: c_x_fp(:)
        real(real32), intent(in) :: c_y_fp(:)
        real(real32), intent(in) :: w(:)
        real(real32), intent(in) :: rho_0
        real(real32), intent(in) :: u_max
        real(real32), intent(in) :: n_sin

        ! write destinations
        real(real32), intent(out) :: f(:, :, :)
        real(real32), intent(out) :: rho(:,:)
        real(real32), intent(out) :: u_x(:,:)
        real(real32), intent(out) :: u_y(:,:)

        ! temp
        integer(int32) :: x, y, i
        real(real32) :: k ! wave number
        real(real32) :: u_x_val
        real(real32) :: u_y_val
        real(real32) :: u_squ
        real(real32) :: c_dot_u
        real(real32) :: f_eq_val

        ! wave number
        k = (2.0_real32 * pi * n_sin) / real(N_Y, real32)

        ! loop over rows
        do y = 1, N_Y

            ! shear wave velocity for this row
            u_x_val = u_max * sin(k * real(y - 1, real32))
            u_y_val = 0.0_real32
            u_squ = u_x_val * u_x_val

            ! loop over cols
            do x = 1, N_X

                rho(x, y) = rho_0
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val

                ! init distribution functions in all dirs
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for dir i
                    c_dot_u = c_x_fp(i) * u_x_val + c_y_fp(i) * u_y_val
                    f_eq_val = w(i) * rho_0 * ( &
                        1.0_real32 + &
                        3.0_real32 * c_dot_u + &
                        4.5_real32 * c_dot_u * c_dot_u - &
                        1.5_real32 * u_squ)

                    ! write to destination dir i of this cell
                    f(i, x, y) = f_eq_val
                end do
            end do
        end do
    end subroutine apply_condition_shear_wave


    subroutine apply_condition_couette_flow( &
        N_X, N_Y, N_DIRS, w, rho_0, f, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: N_X
        integer(int32), intent(in) :: N_Y
        integer(int32), intent(in) :: N_DIRS
        real(real32), intent(in) :: w(:)
        real(real32), intent(in) :: rho_0

        ! write destinations
        real(real32), intent(out) :: f(:, :, :)
        real(real32), intent(out) :: rho(:,:)
        real(real32), intent(out) :: u_x(:,:)
        real(real32), intent(out) :: u_y(:,:)

        ! temp
        integer(int32) :: x, y, i
        real(real32) :: f_eq(N_DIRS)

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = w(i) * rho_0
        end do

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_real32
                u_y(x, y) = 0.0_real32

                ! init distribution functions in all dirs
                do i = 1, N_DIRS
                    f(i, x, y) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_couette_flow


    subroutine apply_condition_sliding_lid( &
        N_X, N_Y, N_DIRS, w, rho_0, f, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: N_X
        integer(int32), intent(in) :: N_Y
        integer(int32), intent(in) :: N_DIRS
        real(real32), intent(in) :: w(:)
        real(real32), intent(in) :: rho_0

        ! write destinations
        real(real32), intent(out) :: f(:, :, :)
        real(real32), intent(out) :: rho(:,:)
        real(real32), intent(out) :: u_x(:,:)
        real(real32), intent(out) :: u_y(:,:)

        ! temp
        integer(int32) :: x, y, i
        real(real32) :: f_eq(N_DIRS)

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = w(i) * rho_0
        end do

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_real32
                u_y(x, y) = 0.0_real32

                ! init distribution functions in all dirs
                do i = 1, N_DIRS
                    f(i, x, y) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_sliding_lid

end module initialization
