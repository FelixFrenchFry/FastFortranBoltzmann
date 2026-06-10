module sliding_lid
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_DIRS, W, FP
    implicit none

contains

    subroutine prepare_sliding_lid_halos_SL( &
        n_x_local, n_y_local, at_left_boundary, at_right_boundary, at_bottom_boundary, at_top_boundary, &
        rho_0, u_wall, f &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_left_boundary
        logical, intent(in) :: at_right_boundary
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: u_wall

        ! read/write inputs
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! temp
        integer(int32) :: x, y
        real(FP) :: moving_wall_correction_8
        real(FP) :: moving_wall_correction_9

        ! left bounce-back boundary, written into the halo column used by pull streaming
        if (at_left_boundary) then
            do y = 1, n_y_local
                f(0, y, 2) = f(1, y, 4)
                f(0, y-1, 6) = f(1, y, 8)
                f(0, y+1, 9) = f(1, y, 7)
            end do
        end if

        ! right bounce-back boundary, written into the halo column used by pull streaming
        if (at_right_boundary) then
            do y = 1, n_y_local
                f(n_x_local+1, y, 4) = f(n_x_local, y, 2)
                f(n_x_local+1, y-1, 7) = f(n_x_local, y, 9)
                f(n_x_local+1, y+1, 8) = f(n_x_local, y, 6)
            end do
        end if

        ! bottom bounce-back boundary, written into the halo row used by pull streaming
        if (at_bottom_boundary) then
            do x = 1, n_x_local
                f(x, 0, 3) = f(x, 1, 5)
                f(x-1, 0, 6) = f(x, 1, 8)
                f(x+1, 0, 7) = f(x, 1, 9)
            end do
        end if

        ! top bounce-back boundary, written into the halo row used by pull streaming
        if (at_top_boundary) then
            moving_wall_correction_8 = 6.0_FP * W(6) * rho_0 * u_wall
            moving_wall_correction_9 = 6.0_FP * W(7) * rho_0 * u_wall

            do x = 1, n_x_local
                f(x, n_y_local+1, 5) = f(x, n_y_local, 3)
                f(x+1, n_y_local+1, 8) = f(x, n_y_local, 6) - moving_wall_correction_8
                f(x-1, n_y_local+1, 9) = f(x, n_y_local, 7) + moving_wall_correction_9
            end do
        end if
    end subroutine prepare_sliding_lid_halos_SL


end module sliding_lid
