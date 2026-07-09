module couette_flow
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_DIRS, W, FP
    implicit none

contains

    subroutine prepare_couette_flow_halos_CF( &
        n_images_x, n_x_local, n_y_local, at_bottom_boundary, at_top_boundary, rho_0, u_wall, f &
        )
        ! inputs
        integer(int32), intent(in) :: n_images_x
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in) :: at_bottom_boundary
        logical, intent(in) :: at_top_boundary
        real(FP), intent(in) :: rho_0
        real(FP), intent(in) :: u_wall

        ! read/write inputs
        real(FP), intent(inout) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! temp
        integer(int32) :: x, y
        real(FP) :: moving_wall_correction_8
        real(FP) :: moving_wall_correction_9

        ! periodic boundary handling for left/right sides
        ! (only used in single-image decompositions, otherwise handled by wrapped halo exchange)
        if (n_images_x == 1) then
            do y = 1, n_y_local
                f(2, 0, y) = f(2, n_x_local, y)
                f(6, 0, y) = f(6, n_x_local, y)
                f(9, 0, y) = f(9, n_x_local, y)

                f(4, n_x_local+1, y) = f(4, 1, y)
                f(7, n_x_local+1, y) = f(7, 1, y)
                f(8, n_x_local+1, y) = f(8, 1, y)
            end do
        end if

        ! bottom bounce-back boundary, written into the halo row used by pull streaming
        if (at_bottom_boundary) then
            do x = 1, n_x_local
                f(3, x, 0) = f(5, x, 1)
                f(6, x-1, 0) = f(8, x, 1)
                f(7, x+1, 0) = f(9, x, 1)
            end do
        end if

        ! top bounce-back boundary, written into the halo row used by pull streaming
        if (at_top_boundary) then
            moving_wall_correction_8 = 6.0_FP * W(6) * rho_0 * u_wall
            moving_wall_correction_9 = 6.0_FP * W(7) * rho_0 * u_wall

            do x = 1, n_x_local
                f(5, x, n_y_local+1) = f(3, x, n_y_local)
                f(8, x+1, n_y_local+1) = f(6, x, n_y_local) - moving_wall_correction_8
                f(9, x-1, n_y_local+1) = f(7, x, n_y_local) + moving_wall_correction_9
            end do
        end if
    end subroutine prepare_couette_flow_halos_CF


end module couette_flow
