module simulation
    ! imports
    use iso_fortran_env, only: int32
    use domain, only: domain_t
    use exchange, only: halo_buffers_t
    use settings, only: N_DIRS, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID, SIM_MODE, FP, &
        USE_UNROLLED_KERNELS, USE_PULL_SHIFT_KERNELS, &
        shear_wave_params_t, couette_flow_params_t, poiseuille_flow_params_t, sliding_lid_params_t
    use shear_wave, only: fuzed_pull_streaming_collision_local_SW, fuzed_pull_shift_streaming_collision_local_SW, &
        fuzed_pull_shift_streaming_collision_local_unrolled_SW, fuzed_pull_streaming_collision_local_unrolled_SW
    use couette_flow, only: fuzed_pull_streaming_collision_local_CF, fuzed_pull_streaming_collision_local_unrolled_CF
    use poiseuille_flow, only: fuzed_pull_streaming_collision_local_PF, &
        fuzed_pull_streaming_collision_local_unrolled_PF
    use sliding_lid, only: fuzed_pull_streaming_collision_local_SL, fuzed_pull_streaming_collision_local_unrolled_SL
    implicit none

contains

    subroutine execute_local_sim_step( &
        domain_info, halo_buffers, n_x_local, n_y_local, &
        shear_wave_params, couette_flow_params, poiseuille_flow_params, sliding_lid_params, &
        write_macro_fields, f, f_next, rho, u_x, u_y &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        type(shear_wave_params_t), intent(in) :: shear_wave_params
        type(couette_flow_params_t), intent(in) :: couette_flow_params
        type(poiseuille_flow_params_t), intent(in) :: poiseuille_flow_params
        type(sliding_lid_params_t), intent(in) :: sliding_lid_params
        logical, intent(in) :: write_macro_fields
        real(FP), intent(in) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! write destinations
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f_next(0:n_x_local+1, 0:n_y_local+1, N_DIRS)
        real(FP), intent(inout) :: rho(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_x(n_x_local, n_y_local)
        real(FP), intent(inout) :: u_y(n_x_local, n_y_local)

        ! execute single local sim step based on selected sim mode
        select case (SIM_MODE)
        case (SIM_SHEAR_WAVE)
            if (USE_PULL_SHIFT_KERNELS) then
                if (USE_UNROLLED_KERNELS) then
                    call fuzed_pull_shift_streaming_collision_local_unrolled_SW( &
                        n_x_local, n_y_local, &
                        write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
                else
                    call fuzed_pull_shift_streaming_collision_local_SW( &
                        n_x_local, n_y_local, &
                        write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
                end if
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_SW( &
                    n_x_local, n_y_local, &
                    write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_local_SW( &
                    n_x_local, n_y_local, &
                    write_macro_fields, shear_wave_params%omega, f, f_next, rho, u_x, u_y)
            end if
        case (SIM_COUETTE_FLOW)
            if (USE_PULL_SHIFT_KERNELS) then
                error stop "error: distributed pull-shift is not implemented for this simulation mode yet"
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_CF( &
                    n_x_local, n_y_local, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, couette_flow_params%rho_0, couette_flow_params%omega, couette_flow_params%u_wall, &
                    f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_local_CF( &
                    n_x_local, n_y_local, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, couette_flow_params%rho_0, couette_flow_params%omega, couette_flow_params%u_wall, &
                    f, f_next, rho, u_x, u_y)
            end if
        case (SIM_POISEUILLE_FLOW)
            if (USE_PULL_SHIFT_KERNELS) then
                error stop "error: distributed pull-shift is not implemented for this simulation mode yet"
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_PF( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, poiseuille_flow_params%omega, &
                    poiseuille_flow_params%rho_in, poiseuille_flow_params%rho_out, &
                    f, f_next, rho, u_x, u_y, &
                    halo_buffers%recv_macro_left, halo_buffers%recv_macro_right, &
                    halo_buffers%send_macro_left, halo_buffers%send_macro_right)
            else
                call fuzed_pull_streaming_collision_local_PF( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, poiseuille_flow_params%omega, &
                    poiseuille_flow_params%rho_in, poiseuille_flow_params%rho_out, &
                    f, f_next, rho, u_x, u_y, &
                    halo_buffers%recv_macro_left, halo_buffers%recv_macro_right, &
                    halo_buffers%send_macro_left, halo_buffers%send_macro_right)
            end if
        case (SIM_SLIDING_LID)
            if (USE_PULL_SHIFT_KERNELS) then
                error stop "error: distributed pull-shift is not implemented for this simulation mode yet"
            else if (USE_UNROLLED_KERNELS) then
                call fuzed_pull_streaming_collision_local_unrolled_SL( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, sliding_lid_params%rho_0, sliding_lid_params%omega, &
                    sliding_lid_params%u_wall, f, f_next, rho, u_x, u_y)
            else
                call fuzed_pull_streaming_collision_local_SL( &
                    n_x_local, n_y_local, &
                    domain_info%at_left_boundary, domain_info%at_right_boundary, &
                    domain_info%at_bottom_boundary, domain_info%at_top_boundary, &
                    write_macro_fields, sliding_lid_params%rho_0, sliding_lid_params%omega, &
                    sliding_lid_params%u_wall, f, f_next, rho, u_x, u_y)
            end if
        case default
            error stop "error: unknown sim mode in execute_local_sim_step()"
        end select
    end subroutine execute_local_sim_step


    subroutine swap_distribution_function_buffers( &
        f, f_next &
        )
        ! read/write inputs
        real(FP), allocatable, intent(inout) :: f(:, :, :)
        real(FP), allocatable, intent(inout) :: f_next(:, :, :)
        real(FP), allocatable :: temp(:, :, :)

        ! swap ownership
        call move_alloc(f, temp)
        call move_alloc(f_next, f)
        call move_alloc(temp, f_next)
    end subroutine swap_distribution_function_buffers

end module simulation
