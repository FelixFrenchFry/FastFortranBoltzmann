program test_boundaries
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: FP, N_DIRS, C_X, C_Y
    use shear_wave, only: prepare_shear_wave_halos_SW
    use couette_flow, only: prepare_couette_flow_halos_CF
    use assertions, only: check, check_close
    implicit none

    ! sim size and parameters
    integer(int32), parameter :: N_X_LOCAL = 3
    integer(int32), parameter :: N_Y_LOCAL = 2
    real(FP), parameter :: RHO_0 = 1.2_FP
    real(FP), parameter :: U_WALL = 0.03_FP

    ! distribution functions
    real(FP) :: f(0:N_X_LOCAL+1, 0:N_Y_LOCAL+1, N_DIRS)
    real(FP) :: f_original(0:N_X_LOCAL+1, 0:N_Y_LOCAL+1, N_DIRS)

    ! temp
    integer(int32) :: x, y, i
    integer(int32) :: src_x, src_y

    ! distinct cell/channel values expose incorrect directions and corner indices
    do i = 1, N_DIRS
        do y = 0, N_Y_LOCAL + 1
            do x = 0, N_X_LOCAL + 1
                f_original(x, y, i) = real(100 * i + 10 * y + x, FP) / 1000.0_FP
            end do
        end do
    end do

    ! verify periodic incoming populations, including corners
    f = f_original
    call prepare_shear_wave_halos_SW( &
        1, 1, N_X_LOCAL, N_Y_LOCAL, f)

    do i = 1, N_DIRS
        do y = 1, N_Y_LOCAL
            do x = 1, N_X_LOCAL
                src_x = x - C_X(i)
                src_y = y - C_Y(i)
                call check_close( &
                    f(src_x, src_y, i), &
                    f_original(modulo(src_x - 1, N_X_LOCAL) + 1, modulo(src_y - 1, N_Y_LOCAL) + 1, i), &
                    "periodic pull source including corners")
            end do
        end do
    end do
    call check( &
        all(f(1:N_X_LOCAL, 1:N_Y_LOCAL, :) == f_original(1:N_X_LOCAL, 1:N_Y_LOCAL, :)), &
        "periodic halos preserve interior")

    ! interior subdomains must preserve the halos received from other images
    f = f_original
    call prepare_couette_flow_halos_CF( &
        2, N_X_LOCAL, N_Y_LOCAL, .false., .false., RHO_0, U_WALL, f)
    call check(all(f == f_original), "interior subdomains preserve received halos")

    ! verify stationary bottom and moving top wall bounce-back
    call prepare_couette_flow_halos_CF( &
        2, N_X_LOCAL, N_Y_LOCAL, .true., .true., RHO_0, U_WALL, f)

    do x = 1, N_X_LOCAL
        call check_close(f(x, 0, 3), f_original(x, 1, 5), "bottom normal bounce-back")
        call check_close(f(x-1, 0, 6), f_original(x, 1, 8), "bottom diagonal bounce-back")
        call check_close(f(x+1, 0, 7), f_original(x, 1, 9), "bottom opposite diagonal")
        call check_close(f(x, N_Y_LOCAL+1, 5), f_original(x, N_Y_LOCAL, 3), "top normal bounce-back")
        call check_close( &
            f(x+1, N_Y_LOCAL+1, 8), f_original(x, N_Y_LOCAL, 6) - RHO_0 * U_WALL / 6.0_FP, &
            "moving wall negative correction")
        call check_close( &
            f(x-1, N_Y_LOCAL+1, 9), f_original(x, N_Y_LOCAL, 7) + RHO_0 * U_WALL / 6.0_FP, &
            "moving wall positive correction")
    end do
end program test_boundaries
