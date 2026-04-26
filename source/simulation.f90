module simulation
    ! imports
    use iso_fortran_env, only: int32, real32
    use settings, only: N_X, N_Y, N_DIRS, SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, &
        PI, shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t
    implicit none

contains

    subroutine execute_full_sim_step( &
        shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
        c_x, c_y, c_x_fp, c_y_fp, w, f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y &
        )
        ! read-only inputs
        type(shear_wave_params_t), intent(in) :: shear_wave_params
        type(couette_flow_params_t), intent(in) :: couette_flow_params
        type(poiseuille_flow_params_t), intent(in) :: poiseuille_flow_params
        type(sliding_lid_params_t), intent(in) :: sliding_lid_params
        integer(int32), intent(in) :: c_x(N_DIRS)
        integer(int32), intent(in) :: c_y(N_DIRS)
        real(real32), intent(in) :: c_x_fp(N_DIRS)
        real(real32), intent(in) :: c_y_fp(N_DIRS)
        real(real32), intent(in) :: w(N_DIRS)
        real(real32), intent(in) :: f(N_DIRS, N_X, N_Y)
        logical, intent(in) :: write_rho
        logical, intent(in) :: write_u_x
        logical, intent(in) :: write_u_y

        ! write destinations
        real(real32), intent(out) :: f_next(N_DIRS, N_X, N_Y)

        ! optional write destinations
        real(real32), intent(inout) :: rho(N_X, N_Y)
        real(real32), intent(inout) :: u_x(N_X, N_Y)
        real(real32), intent(inout) :: u_y(N_X, N_Y)

        ! execute single sim step based on selected sim mode
        select case (SIM_MODE)
        case (SIM_SHEAR_WAVE)
            call fuzed_pull_streaming_collision_shear_wave( &
                c_x, c_y, c_x_fp, c_y_fp, w, shear_wave_params%omega, &
                f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y)
        case (SIM_COUETTE_FLOW)
            call fuzed_pull_streaming_collision_couette_flow( &
                c_x, c_y, c_x_fp, c_y_fp, w, couette_flow_params%rho_0, &
                couette_flow_params%omega, couette_flow_params%u_wall, &
                f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y)
        case (SIM_POISEUILLE_FLOW)
            call fuzed_pull_streaming_collision_poiseuille_flow( &
                c_x, c_y, c_x_fp, c_y_fp, w, poiseuille_flow_params%omega, &
                poiseuille_flow_params%rho_in, poiseuille_flow_params%rho_out, &
                f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y)
        case (SIM_SLIDING_LID)
            call fuzed_pull_streaming_collision_sliding_lid( &
                c_x, c_y, c_x_fp, c_y_fp, w, sliding_lid_params%rho_0, &
                sliding_lid_params%omega, sliding_lid_params%u_wall, &
                f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y)
        case default
            error stop "error: unknown sim mode in execute_full_sim_step()"
        end select
    end subroutine execute_full_sim_step


    subroutine fuzed_pull_streaming_collision_shear_wave( &
        c_x, c_y, c_x_fp, c_y_fp, w, omega, &
        f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: c_x(N_DIRS)
        integer(int32), intent(in) :: c_y(N_DIRS)
        real(real32), intent(in) :: c_x_fp(N_DIRS)
        real(real32), intent(in) :: c_y_fp(N_DIRS)
        real(real32), intent(in) :: w(N_DIRS)
        real(real32), intent(in) :: omega
        real(real32), intent(in) :: f(N_DIRS, N_X, N_Y)
        logical, intent(in) :: write_rho
        logical, intent(in) :: write_u_x
        logical, intent(in) :: write_u_y

        ! write destinations
        real(real32), intent(out) :: f_next(N_DIRS, N_X, N_Y)

        ! optional write destinations
        real(real32), intent(inout) :: rho(N_X, N_Y)
        real(real32), intent(inout) :: u_x(N_X, N_Y)
        real(real32), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        integer(int32) :: src_x, src_y
        real(real32) :: f_pulled(N_DIRS)
        real(real32) :: rho_val
        real(real32) :: u_x_val
        real(real32) :: u_y_val
        real(real32) :: u_squ
        real(real32) :: c_dot_u
        real(real32) :: f_eq_val
        real(real32) :: f_next_val

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho_val = 0.0_real32
                u_x_val = 0.0_real32
                u_y_val = 0.0_real32

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
                ! pull streamed distribution functions from source cells in all dirs
                do i = 1, N_DIRS

                    src_x = x - c_x(i)
                    src_y = y - c_y(i)

                    ! periodic boundary in x-dir
                    ! TODO: optimize by using a branchless method?
                    if (src_x < 1) then
                        src_x = N_X
                    else if (src_x > N_X) then
                        src_x = 1
                    end if

                    ! periodic boundary in y-dir
                    ! TODO: optimize by using a branchless method?
                    if (src_y < 1) then
                        src_y = N_Y
                    else if (src_y > N_Y) then
                        src_y = 1
                    end if

                    f_pulled(i) = f(i, src_x, src_y)

                    rho_val = rho_val + f_pulled(i)
                    u_x_val = u_x_val + f_pulled(i) * c_x_fp(i)
                    u_y_val = u_y_val + f_pulled(i) * c_y_fp(i)
                end do

                ! safety check
                if (rho_val <= 0.0_real32) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
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
                    f_eq_val = w(i) * rho_val * ( &
                        1.0_real32 + &
                        3.0_real32 * c_dot_u + &
                        4.5_real32 * c_dot_u * c_dot_u - &
                        1.5_real32 * u_squ)

                    ! relax towards equilibrium
                    f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)

                    ! write to destination dir i of this cell in next distribution function buffer
                    f_next(i, x, y) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_shear_wave


    subroutine fuzed_pull_streaming_collision_couette_flow( &
        c_x, c_y, c_x_fp, c_y_fp, w, rho_0, omega, u_wall, &
        f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: c_x(N_DIRS)
        integer(int32), intent(in) :: c_y(N_DIRS)
        real(real32), intent(in) :: c_x_fp(N_DIRS)
        real(real32), intent(in) :: c_y_fp(N_DIRS)
        real(real32), intent(in) :: w(N_DIRS)
        real(real32), intent(in) :: rho_0
        real(real32), intent(in) :: omega
        real(real32), intent(in) :: u_wall
        real(real32), intent(in) :: f(N_DIRS, N_X, N_Y)
        logical, intent(in) :: write_rho
        logical, intent(in) :: write_u_x
        logical, intent(in) :: write_u_y

        ! write destinations
        real(real32), intent(out) :: f_next(N_DIRS, N_X, N_Y)

        ! optional write destinations
        real(real32), intent(inout) :: rho(N_X, N_Y)
        real(real32), intent(inout) :: u_x(N_X, N_Y)
        real(real32), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        integer(int32) :: src_x, src_y
        real(real32) :: f_pulled(N_DIRS)
        real(real32) :: rho_val
        real(real32) :: u_x_val
        real(real32) :: u_y_val
        real(real32) :: u_squ
        real(real32) :: c_dot_u
        real(real32) :: f_eq_val
        real(real32) :: f_next_val

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho_val = 0.0_real32
                u_x_val = 0.0_real32
                u_y_val = 0.0_real32

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
                ! pull streamed distribution functions from source cells in all dirs
                do i = 1, N_DIRS

                    src_x = x - c_x(i)
                    src_y = y - c_y(i)

                    ! periodic boundary in x-dir
                    if (src_x < 1) then
                        src_x = N_X
                    else if (src_x > N_X) then
                        src_x = 1
                    end if

                    ! non-periodic boundary in y-dir with bounce-back
                    if (src_y >= 1 .and. src_y <= N_Y) then ! inner cell -> normal streaming
                        f_pulled(i) = f(i, src_x, src_y)
                    
                    else if (src_y < 1) then ! bottom wall -> static bounce-back
                        select case (i)
                        case (3)
                            f_pulled(i) = f(5, x, y)
                        case (6)
                            f_pulled(i) = f(8, x, y)
                        case (7)
                            f_pulled(i) = f(9, x, y)
                        case default
                            error stop "error: invalid bottom wall dir in couette flow"
                        end select

                    else if (src_y > N_Y) then ! top wall -> moving bounce-back
                        select case (i)
                        case (5)
                            f_pulled(i) = f(3, x, y)
                        case (8)
                            f_pulled(i) = f(6, x, y) - 6.0_real32 * w(6) * rho_0 * u_wall
                        case (9)
                            f_pulled(i) = f(7, x, y) + 6.0_real32 * w(7) * rho_0 * u_wall
                        case default
                            error stop "error: invalid top wall dir in couette flow"
                        end select
                    end if

                    rho_val = rho_val + f_pulled(i)
                    u_x_val = u_x_val + f_pulled(i) * c_x_fp(i)
                    u_y_val = u_y_val + f_pulled(i) * c_y_fp(i)
                end do

                ! safety check
                if (rho_val <= 0.0_real32) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
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
                    f_eq_val = w(i) * rho_val * ( &
                        1.0_real32 + &
                        3.0_real32 * c_dot_u + &
                        4.5_real32 * c_dot_u * c_dot_u - &
                        1.5_real32 * u_squ)

                    ! relax towards equilibrium
                    f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)

                    ! write to destination dir i of this cell in next distribution function buffer
                    f_next(i, x, y) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_couette_flow


    subroutine fuzed_pull_streaming_collision_poiseuille_flow( &
        c_x, c_y, c_x_fp, c_y_fp, w, omega, rho_in, rho_out, &
        f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: c_x(N_DIRS)
        integer(int32), intent(in) :: c_y(N_DIRS)
        real(real32), intent(in) :: c_x_fp(N_DIRS)
        real(real32), intent(in) :: c_y_fp(N_DIRS)
        real(real32), intent(in) :: w(N_DIRS)
        real(real32), intent(in) :: omega
        real(real32), intent(in) :: rho_in
        real(real32), intent(in) :: rho_out
        real(real32), intent(in) :: f(N_DIRS, N_X, N_Y)
        logical, intent(in) :: write_rho
        logical, intent(in) :: write_u_x
        logical, intent(in) :: write_u_y

        ! write destinations
        real(real32), intent(out) :: f_next(N_DIRS, N_X, N_Y)

        ! optional write destinations
        real(real32), intent(inout) :: rho(N_X, N_Y)
        real(real32), intent(inout) :: u_x(N_X, N_Y)
        real(real32), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        integer(int32) :: src_x, src_y
        real(real32) :: f_pulled(N_DIRS)
        real(real32) :: rho_val
        real(real32) :: u_x_val
        real(real32) :: u_y_val
        real(real32) :: u_squ
        real(real32) :: c_dot_u
        real(real32) :: f_eq_val
        real(real32) :: f_next_val
        real(real32) :: u_squ_src
        real(real32) :: c_dot_u_src
        real(real32) :: f_eq_src
        real(real32) :: f_eq_boundary

        ! copied boundary densities and velocities for streaming step
        real(real32) :: rho_left(N_Y), rho_right(N_Y)
        real(real32) :: u_x_left(N_Y), u_x_right(N_Y)
        real(real32) :: u_y_left(N_Y), u_y_right(N_Y)
        rho_left(:) = rho(1, :)
        rho_right(:) = rho(N_X, :)
        u_x_left(:) = u_x(1, :)
        u_x_right(:) = u_x(N_X, :)
        u_y_left(:) = u_y(1, :)
        u_y_right(:) = u_y(N_X, :)

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho_val = 0.0_real32
                u_x_val = 0.0_real32
                u_y_val = 0.0_real32

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
                ! pull streamed distribution functions from source cells in all dirs
                do i = 1, N_DIRS

                    src_x = x - c_x(i)
                    src_y = y - c_y(i)

                    ! non-periodic boundary in x-dir and y-dir with bounce-back
                    if (src_x >= 1 .and. src_x <= N_X .and. &
                        src_y >= 1 .and. src_y <= N_Y) then ! inner cell -> normal streaming
                        f_pulled(i) = f(i, src_x, src_y)
                    
                    else if (src_y > N_Y) then ! top wall -> static bounce-back
                        select case (i)
                        case (5)
                            f_pulled(i) = f(3, x, y)
                        case (8)
                            f_pulled(i) = f(6, x, y)
                        case (9)
                            f_pulled(i) = f(7, x, y)
                        case default
                            error stop "error: invalid top wall dir in poiseuille flow"
                        end select

                    else if (src_y < 1) then ! bottom wall -> static bounce-back
                        select case (i)
                        case (3)
                            f_pulled(i) = f(5, x, y)
                        case (6)
                            f_pulled(i) = f(8, x, y)
                        case (7)
                            f_pulled(i) = f(9, x, y)
                        case default
                            error stop "error: invalid bottom wall dir in poiseuille flow"
                        end select
                    
                    else if (src_x > N_X) then ! right boundary -> pressure-periodic outlet
                        u_squ_src = u_x_left(y) * u_x_left(y) + u_y_left(y) * u_y_left(y)
                        c_dot_u_src = c_x_fp(i) * u_x_left(y) + c_y_fp(i) * u_y_left(y)

                        ! equilibrium distribution function for the source cell at the opposite (left) boundary
                        f_eq_src = w(i) * rho_left(y) * ( &
                            1.0_real32 + &
                            3.0_real32 * c_dot_u_src + &
                            4.5_real32 * c_dot_u_src * c_dot_u_src - &
                            1.5_real32 * u_squ_src)

                        ! equilibrium distribution function at this cell
                        f_eq_boundary = w(i) * rho_out * ( &
                            1.0_real32 + &
                            3.0_real32 * c_dot_u_src + &
                            4.5_real32 * c_dot_u_src * c_dot_u_src - &
                            1.5_real32 * u_squ_src)

                        f_pulled(i) = f(i, 1, y) - f_eq_src + f_eq_boundary
                    
                    else if (src_x < 1) then ! left boundary -> pressure-periodic inlet
                        u_squ_src = u_x_right(y) * u_x_right(y) + u_y_right(y) * u_y_right(y)
                        c_dot_u_src = c_x_fp(i) * u_x_right(y) + c_y_fp(i) * u_y_right(y)

                        ! equilibrium distribution function for the source cell at the opposite (right) boundary
                        f_eq_src = w(i) * rho_right(y) * ( &
                            1.0_real32 + &
                            3.0_real32 * c_dot_u_src + &
                            4.5_real32 * c_dot_u_src * c_dot_u_src - &
                            1.5_real32 * u_squ_src)

                        ! equilibrium distribution function at this cell
                        f_eq_boundary = w(i) * rho_in * ( &
                            1.0_real32 + &
                            3.0_real32 * c_dot_u_src + &
                            4.5_real32 * c_dot_u_src * c_dot_u_src - &
                            1.5_real32 * u_squ_src)

                        f_pulled(i) = f(i, N_X, y) - f_eq_src + f_eq_boundary
                    end if

                    rho_val = rho_val + f_pulled(i)
                    u_x_val = u_x_val + f_pulled(i) * c_x_fp(i)
                    u_y_val = u_y_val + f_pulled(i) * c_y_fp(i)
                end do

                ! safety check
                if (rho_val <= 0.0_real32) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
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
                    f_eq_val = w(i) * rho_val * ( &
                        1.0_real32 + &
                        3.0_real32 * c_dot_u + &
                        4.5_real32 * c_dot_u * c_dot_u - &
                        1.5_real32 * u_squ)

                    ! relax towards equilibrium
                    f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)

                    ! write to destination dir i of this cell in next distribution function buffer
                    f_next(i, x, y) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_poiseuille_flow


    subroutine fuzed_pull_streaming_collision_sliding_lid( &
        c_x, c_y, c_x_fp, c_y_fp, w, rho_0, omega, u_wall, &
        f, write_rho, write_u_x, write_u_y, f_next, rho, u_x, u_y &
        )
        ! read-only inputs
        integer(int32), intent(in) :: c_x(N_DIRS)
        integer(int32), intent(in) :: c_y(N_DIRS)
        real(real32), intent(in) :: c_x_fp(N_DIRS)
        real(real32), intent(in) :: c_y_fp(N_DIRS)
        real(real32), intent(in) :: w(N_DIRS)
        real(real32), intent(in) :: rho_0
        real(real32), intent(in) :: omega
        real(real32), intent(in) :: u_wall
        real(real32), intent(in) :: f(N_DIRS, N_X, N_Y)
        logical, intent(in) :: write_rho
        logical, intent(in) :: write_u_x
        logical, intent(in) :: write_u_y

        ! write destinations
        real(real32), intent(out) :: f_next(N_DIRS, N_X, N_Y)

        ! optional write destinations
        real(real32), intent(inout) :: rho(N_X, N_Y)
        real(real32), intent(inout) :: u_x(N_X, N_Y)
        real(real32), intent(inout) :: u_y(N_X, N_Y)

        ! temp
        integer(int32) :: x, y, i
        integer(int32) :: src_x, src_y
        real(real32) :: f_pulled(N_DIRS)
        real(real32) :: rho_val
        real(real32) :: u_x_val
        real(real32) :: u_y_val
        real(real32) :: u_squ
        real(real32) :: c_dot_u
        real(real32) :: f_eq_val
        real(real32) :: f_next_val

        ! loop over rows and cols
        do y = 1, N_Y
            do x = 1, N_X

                rho_val = 0.0_real32
                u_x_val = 0.0_real32
                u_y_val = 0.0_real32

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
                ! pull streamed distribution functions from source cells in all dirs
                do i = 1, N_DIRS

                    src_x = x - c_x(i)
                    src_y = y - c_y(i)

                    ! non-periodic boundary in x-dir and y-dir with bounce-back
                    if (src_x >= 1 .and. src_x <= N_X .and. &
                        src_y >= 1 .and. src_y <= N_Y) then ! inner cell -> normal streaming
                        f_pulled(i) = f(i, src_x, src_y)
                    
                    else if (src_y > N_Y) then ! top wall -> moving bounce-back
                        select case (i)
                        case (5)
                            f_pulled(i) = f(3, x, y)
                        case (8)
                            f_pulled(i) = f(6, x, y) - 6.0_real32 * w(6) * rho_0 * u_wall
                        case (9)
                            f_pulled(i) = f(7, x, y) + 6.0_real32 * w(7) * rho_0 * u_wall
                        case default
                            error stop "error: invalid top wall dir in sliding lid"
                        end select

                    else if (src_y < 1) then ! bottom wall -> static bounce-back
                        select case (i)
                        case (3)
                            f_pulled(i) = f(5, x, y)
                        case (6)
                            f_pulled(i) = f(8, x, y)
                        case (7)
                            f_pulled(i) = f(9, x, y)
                        case default
                            error stop "error: invalid bottom wall dir in sliding lid"
                        end select
                    
                    else if (src_x > N_X) then ! right wall -> static bounce-back
                        select case (i)
                        case (4)
                            f_pulled(i) = f(2, x, y)
                        case (7)
                            f_pulled(i) = f(9, x, y)
                        case (8)
                            f_pulled(i) = f(6, x, y)
                        case default
                            error stop "error: invalid right wall dir in sliding lid"
                        end select
                    
                    else if (src_x < 1) then ! left wall -> static bounce-back
                        select case (i)
                        case (2)
                            f_pulled(i) = f(4, x, y)
                        case (6)
                            f_pulled(i) = f(8, x, y)
                        case (9)
                            f_pulled(i) = f(7, x, y)
                        case default
                            error stop "error: invalid left wall dir in sliding lid"
                        end select
                    end if

                    rho_val = rho_val + f_pulled(i)
                    u_x_val = u_x_val + f_pulled(i) * c_x_fp(i)
                    u_y_val = u_y_val + f_pulled(i) * c_y_fp(i)
                end do

                ! safety check
                if (rho_val <= 0.0_real32) then
                    error stop "error: density is zero in collision/streaming step (rho_val <= 0)"
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
                    f_eq_val = w(i) * rho_val * ( &
                        1.0_real32 + &
                        3.0_real32 * c_dot_u + &
                        4.5_real32 * c_dot_u * c_dot_u - &
                        1.5_real32 * u_squ)

                    ! relax towards equilibrium
                    f_next_val = f_pulled(i) - omega * (f_pulled(i) - f_eq_val)

                    ! write to destination dir i of this cell in next distribution function buffer
                    f_next(i, x, y) = f_next_val
                end do
            end do
        end do
    end subroutine fuzed_pull_streaming_collision_sliding_lid


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
