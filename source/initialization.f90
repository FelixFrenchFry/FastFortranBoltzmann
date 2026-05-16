module initialization
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_X, N_Y, N_DIRS, C_X_FP, C_Y_FP, W, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, PI, &
        shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t
    implicit none

contains

    subroutine initialize_sim_condition( &
        shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
        f, rho, u_x, u_y &
        )
        ! read-only inputs
        type(shear_wave_params_t), intent(in) :: shear_wave_params
        type(couette_flow_params_t), intent(in) :: couette_flow_params
        type(poiseuille_flow_params_t), intent(in) :: poiseuille_flow_params
        type(sliding_lid_params_t), intent(in) :: sliding_lid_params

        ! write destinations
        real(FP), intent(out) :: f(N_X, N_Y, N_DIRS)
        real(FP), intent(out) :: rho(N_X, N_Y)
        real(FP), intent(out) :: u_x(N_X, N_Y)
        real(FP), intent(out) :: u_y(N_X, N_Y)

        ! apply initial condition based on selected sim mode
        select case (SIM_MODE)
        case (SIM_SHEAR_WAVE)
            call apply_condition_shear_wave(shear_wave_params%rho_0, &
                shear_wave_params%u_max, shear_wave_params%n_sin, f, rho, u_x, u_y)
        case (SIM_COUETTE_FLOW)
            call apply_condition_couette_flow(couette_flow_params%rho_0, f, rho, u_x, u_y)
        case (SIM_POISEUILLE_FLOW)
            call apply_condition_poiseuille_flow(poiseuille_flow_params%rho_0, f, rho, u_x, u_y)
        case (SIM_SLIDING_LID)
            call apply_condition_sliding_lid(sliding_lid_params%rho_0, f, rho, u_x, u_y)
        case default
            error stop "error: unknown sim mode in initialize_sim_condition()"
        end select
    end subroutine initialize_sim_condition


    subroutine apply_condition_shear_wave( &
        rho_0, u_max, n_sin, f, rho, u_x, u_y &
        )
        ! read-only inputs
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: u_max
        real(FP), intent(in) :: n_sin

        ! write destinations
        real(FP), intent(out) :: f(N_X, N_Y, N_DIRS)
        real(FP), intent(out) :: rho(N_X, N_Y)
        real(FP), intent(out) :: u_x(N_X, N_Y)
        real(FP), intent(out) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: u_x_val
        real(FP) :: u_y_val
        real(FP) :: u_squ
        real(FP) :: c_dot_u
        real(FP) :: f_eq_val

        ! loop over rows
        do y = 1, N_Y

            ! shear wave velocity for this row
            u_x_val = u_max * sin((2.0_FP * PI * real(n_sin, FP) * real(y - 1, FP)) / real(N_Y, FP))
            u_y_val = 0.0_FP
            u_squ = u_x_val * u_x_val

            ! loop over cols
            do x = 1, N_X

                rho(x, y) = rho_0
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val

                ! init distribution functions in all dirs
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for dir i
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_0 * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)

                    ! write to destination dir i of this cell
                    f(x, y, i) = f_eq_val
                end do
            end do
        end do
    end subroutine apply_condition_shear_wave


    subroutine apply_condition_shear_wave_local( &
        rho_0, u_max, n_sin, n_x_local, n_y_local, y_global_start, f, rho, u_x, u_y &
        )
        ! read-only inputs
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

        ! loop over local rows
        do y = 1, n_y_local

            ! shear wave velocity for this global row
            y_global = y_global_start + y - 1
            u_x_val = u_max * sin((2.0_FP * PI * real(n_sin, FP) * real(y_global - 1, FP)) / real(N_Y, FP))
            u_y_val = 0.0_FP
            u_squ = u_x_val * u_x_val

            ! loop over local cols
            do x = 1, n_x_local

                rho(x, y) = rho_0
                u_x(x, y) = u_x_val
                u_y(x, y) = u_y_val

                ! init distribution functions in all dirs
                do i = 1, N_DIRS

                    ! compute equilibrium distribution function for dir i
                    c_dot_u = C_X_FP(i) * u_x_val + C_Y_FP(i) * u_y_val
                    f_eq_val = W(i) * rho_0 * ( &
                        1.0_FP + &
                        3.0_FP * c_dot_u + &
                        4.5_FP * c_dot_u * c_dot_u - &
                        1.5_FP * u_squ)

                    ! write to destination dir i of this local cell
                    f(x, y, i) = f_eq_val
                end do
            end do
        end do
    end subroutine apply_condition_shear_wave_local


    subroutine apply_condition_couette_flow( &
        rho_0, f, rho, u_x, u_y &
        )
        ! read-only inputs
        real(FP), intent(in) :: rho_0

        ! write destinations
        real(FP), intent(out) :: f(N_X, N_Y, N_DIRS)
        real(FP), intent(out) :: rho(N_X, N_Y)
        real(FP), intent(out) :: u_x(N_X, N_Y)
        real(FP), intent(out) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_eq(N_DIRS)

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = W(i) * rho_0
        end do

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_FP
                u_y(x, y) = 0.0_FP

                ! init distribution functions in all dirs
                do i = 1, N_DIRS
                    f(x, y, i) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_couette_flow


    subroutine apply_condition_poiseuille_flow( &
        rho_0, f, rho, u_x, u_y &
        )
        ! read-only inputs
        real(FP), intent(in) :: rho_0

        ! write destinations
        real(FP), intent(out) :: f(N_X, N_Y, N_DIRS)
        real(FP), intent(out) :: rho(N_X, N_Y)
        real(FP), intent(out) :: u_x(N_X, N_Y)
        real(FP), intent(out) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_eq(N_DIRS)

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = W(i) * rho_0
        end do

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_FP
                u_y(x, y) = 0.0_FP

                ! init distribution functions in all dirs
                do i = 1, N_DIRS
                    f(x, y, i) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_poiseuille_flow


    subroutine apply_condition_sliding_lid( &
        rho_0, f, rho, u_x, u_y &
        )
        ! read-only inputs
        real(FP), intent(in) :: rho_0

        ! write destinations
        real(FP), intent(out) :: f(N_X, N_Y, N_DIRS)
        real(FP), intent(out) :: rho(N_X, N_Y)
        real(FP), intent(out) :: u_x(N_X, N_Y)
        real(FP), intent(out) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_eq(N_DIRS)

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = W(i) * rho_0
        end do

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_FP
                u_y(x, y) = 0.0_FP

                ! init distribution functions in all dirs
                do i = 1, N_DIRS
                    f(x, y, i) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_sliding_lid

end module initialization
