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
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

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
                f(0, y, 2) = f(n_x_local, y, 2)
                f(0, y, 6) = f(n_x_local, y, 6)
                f(0, y, 9) = f(n_x_local, y, 9)

                f(n_x_local+1, y, 4) = f(1, y, 4)
                f(n_x_local+1, y, 7) = f(1, y, 7)
                f(n_x_local+1, y, 8) = f(1, y, 8)
            end do
        end if

        ! periodic boundary handling for bottom/top sides
        ! (only used in single-image decompositions, otherwise handled by wrapped halo exchange)
        if (n_images_y == 1) then
            do x = 0, n_x_local + 1
                f(x, 0, 3) = f(x, n_y_local, 3)
                f(x, 0, 6) = f(x, n_y_local, 6)
                f(x, 0, 7) = f(x, n_y_local, 7)

                f(x, n_y_local+1, 5) = f(x, 1, 5)
                f(x, n_y_local+1, 8) = f(x, 1, 8)
                f(x, n_y_local+1, 9) = f(x, 1, 9)
            end do
        end if
    end subroutine prepare_shear_wave_halos_SW


end module shear_wave
