module simulation
    ! imports
    use iso_fortran_env, only: int32
    use domain, only: domain_t
    use exchange, only: halo_buffers_t
    use settings, only: N_DIRS, C_X, C_Y, C_X_FP, C_Y_FP, W, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, &
        USE_UNROLLED_KERNELS, USE_UNIVERSAL_KERNELS, &
        RHO_0, OMEGA, U_WALL, U_LID, RHO_IN, RHO_OUT
    use shear_wave, only: prepare_shear_wave_halos_SW, fuzed_pull_streaming_collision_local_SW, &
        fuzed_pull_streaming_collision_local_unrolled_SW
    use couette_flow, only: prepare_couette_flow_halos_CF, fuzed_pull_streaming_collision_local_CF, &
        fuzed_pull_streaming_collision_local_unrolled_CF
    use poiseuille_flow, only: prepare_poiseuille_flow_halos_PF, update_poiseuille_flow_macro_strips_PF, &
        fuzed_pull_streaming_collision_local_PF, &
        fuzed_pull_streaming_collision_local_unrolled_PF
    use sliding_lid, only: prepare_sliding_lid_halos_SL, fuzed_pull_streaming_collision_local_SL, &
        fuzed_pull_streaming_collision_local_unrolled_SL
    implicit none

contains

    subroutine execute_local_sim_step( &
        domain_info, halo_buffers, n_x_local, n_y_local, &
        write_macro_fields, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: write_macro_fields
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! write destinations
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)

        ! single sim step based on selected sim mode
        select case (SIM_MODE)
        case (SIM_SHEAR_WAVE) ! shear wave
            if (USE_UNIVERSAL_KERNELS) then
                call prepare_shear_wave_halos_SW( &
                    domain_info%n_images_x, domain_info%n_images_y, n_x_local, n_y_local, f)
                if (USE_UNROLLED_KERNELS) then
                    call fuzed_pull_streaming_collision_local_unrolled_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                else
                    call fuzed_pull_streaming_collision_local_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                end if
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_SW( &
                    n_x_local, n_y_local, &
                    write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_local_SW( &
                    n_x_local, n_y_local, &
                    write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
            end if
        case (SIM_COUETTE_FLOW) ! couette flow
            if (USE_UNIVERSAL_KERNELS) then
                call prepare_couette_flow_halos_CF( &
                    domain_info%n_images_x, n_x_local, n_y_local, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    RHO_0, U_WALL, f)
                if (USE_UNROLLED_KERNELS) then
                    call fuzed_pull_streaming_collision_local_unrolled_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                else
                    call fuzed_pull_streaming_collision_local_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                end if
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_CF( &
                    n_x_local, n_y_local, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, RHO_0, OMEGA, U_WALL, &
                    f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_local_CF( &
                    n_x_local, n_y_local, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, RHO_0, OMEGA, U_WALL, &
                    f, f_next, rho, u_x, u_y)
            end if
        case (SIM_POISEUILLE_FLOW) ! poiseuille flow
            if (USE_UNIVERSAL_KERNELS) then
                call prepare_poiseuille_flow_halos_PF( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    RHO_IN, RHO_OUT, &
                    f, halo_buffers%recv_macro_left, halo_buffers%recv_macro_right)
                call update_poiseuille_flow_macro_strips_PF( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, f, &
                    halo_buffers%send_macro_left, halo_buffers%send_macro_right)
                if (USE_UNROLLED_KERNELS) then
                    call fuzed_pull_streaming_collision_local_unrolled_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                else
                    call fuzed_pull_streaming_collision_local_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                end if
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_PF( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, OMEGA, &
                    RHO_IN, RHO_OUT, &
                    f, f_next, rho, u_x, u_y, &
                    halo_buffers%recv_macro_left, halo_buffers%recv_macro_right, &
                    halo_buffers%send_macro_left, halo_buffers%send_macro_right)
            else
                call fuzed_pull_streaming_collision_local_PF( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, OMEGA, &
                    RHO_IN, RHO_OUT, &
                    f, f_next, rho, u_x, u_y, &
                    halo_buffers%recv_macro_left, halo_buffers%recv_macro_right, &
                    halo_buffers%send_macro_left, halo_buffers%send_macro_right)
            end if
        case (SIM_SLIDING_LID) ! sliding lid
            if (USE_UNIVERSAL_KERNELS) then
                call prepare_sliding_lid_halos_SL( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    RHO_0, U_LID, f)
                if (USE_UNROLLED_KERNELS) then
                    call fuzed_pull_streaming_collision_local_unrolled_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                else
                    call fuzed_pull_streaming_collision_local_universal( &
                        n_x_local, n_y_local, &
                        write_macro_fields, OMEGA, f, f_next, rho, u_x, u_y)
                end if
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_SL( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, RHO_0, OMEGA, &
                    U_LID, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_local_SL( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, RHO_0, OMEGA, &
                    U_LID, f, f_next, rho, u_x, u_y)
            end if
        case default
            error stop "error: unknown sim mode in execute_local_sim_step()"
        end select
    end subroutine execute_local_sim_step


    subroutine fuzed_pull_streaming_collision_local_universal( &
        n_x_local, n_y_local, write_macro_fields, omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)

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
        !DIR$ ASSUME_ALIGNED f: 64, f_next: 64, rho: 64, u_x: 64, u_y: 64
        do y = 1, n_y_local
            !$OMP SIMD
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
                ! (boundaries handled separately in sim-mode-specific halo preparation step)
                !DIR$ UNROLL(9)
                do i = 1, N_DIRS

                    src_x = x - C_X(i)
                    src_y = y - C_Y(i)

                    f_pulled(i) = f(src_x, src_y, i)

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
    end subroutine fuzed_pull_streaming_collision_local_universal


    subroutine fuzed_pull_streaming_collision_local_unrolled_universal( &
        n_x_local, n_y_local, write_macro_fields, omega, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: omega
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! write destinations
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y
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

        ! loop over all image-owned cells
        !DIR$ ASSUME_ALIGNED f: 64, f_next: 64, rho: 64, u_x: 64, u_y: 64
        do y = 1, n_y_local
            !$OMP SIMD
            do x = 1, n_x_local

                ! ---------
                ! | 7 3 6 |
                ! | 4 1 2 |
                ! | 8 5 9 |
                ! ---------
                ! pull streamed distribution functions from source cells
                ! (boundaries handled separately in sim-mode-specific halo preparation step)
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

                ! finalize density and velocity
                u_x_val = u_x_val / rho_val
                u_y_val = u_y_val / rho_val
                u_squ = u_x_val * u_x_val + u_y_val * u_y_val

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
    end subroutine fuzed_pull_streaming_collision_local_unrolled_universal


end module simulation
