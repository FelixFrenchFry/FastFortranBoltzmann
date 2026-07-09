module shear_wave
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_DIRS, FP
    implicit none

contains

    subroutine prepare_shear_wave_halos_SW( &
        n_images_x, n_images_y, n_x_local, n_y_local, f &
        )
        ! inputs
        integer(int32), intent(in) :: n_images_x
        integer(int32), intent(in) :: n_images_y
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local

        ! read/write inputs
        real(FP), intent(inout) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! temp
        integer(int32) :: x, y

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

        ! periodic boundary handling for bottom/top sides
        ! (only used in single-image decompositions, otherwise handled by wrapped halo exchange)
        if (n_images_y == 1) then
            do x = 0, n_x_local + 1
                f(3, x, 0) = f(3, x, n_y_local)
                f(6, x, 0) = f(6, x, n_y_local)
                f(7, x, 0) = f(7, x, n_y_local)

                f(5, x, n_y_local+1) = f(5, x, 1)
                f(8, x, n_y_local+1) = f(8, x, 1)
                f(9, x, n_y_local+1) = f(9, x, 1)
            end do
        end if
    end subroutine prepare_shear_wave_halos_SW


end module shear_wave
