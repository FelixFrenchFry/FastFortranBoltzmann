program main
    ! imports
    use iso_fortran_env, only: int32, int64, real32
    use initialization, only: apply_condition_shear_wave_decay
    implicit none

    ! misc constants
    real(real32), parameter :: pi = 3.1415927410125732421875

    ! simulation size and duration
    integer(int32), parameter :: N_X = 100
    integer(int32), parameter :: N_Y = 100
    integer(int32), parameter :: N_STEPS = 1000
    integer(int64), parameter :: N_CELLS = int(N_X, int64) * int(N_Y, int64)

    ! D2Q9 lattice velocities and weights
    integer(int32), parameter :: N_DIRS = 9
    integer(int32), parameter :: c_x(N_DIRS) = [ 0,  1,  0, -1,  0,  1, -1, -1,  1 ]
    integer(int32), parameter :: c_y(N_DIRS) = [ 0,  0,  1,  0, -1,  1,  1, -1, -1 ]
    real(real32), parameter :: c_x_fp(N_DIRS) = real(c_x, real32) ! fp-version for compute
    real(real32), parameter :: c_y_fp(N_DIRS) = real(c_y, real32) ! fp-version for compute
    ! ---------
    ! | 7 3 6 |
    ! | 4 1 2 |
    ! | 8 5 9 |
    ! ---------
    ! 1: ( 0,  0) = rest
    ! 2: ( 1,  0) = east
    ! 3: ( 0,  1) = north
    ! 4: (-1,  0) = west
    ! 5: ( 0, -1) = south
    ! 6: ( 1,  1) = north-east
    ! 7: (-1,  1) = north-west
    ! 8: (-1, -1) = south-west
    ! 9: ( 1, -1) = south-east
    real(real32), parameter :: w(N_DIRS) = [ &
        4.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/9.0_real32, &
        1.0_real32/36.0_real32, &
        1.0_real32/36.0_real32, &
        1.0_real32/36.0_real32, &
        1.0_real32/36.0_real32 &
    ]

    ! general params
    real(real32), parameter :: rho_0 = 1.0_real32 ! rest density

    ! specific params for shear wave decay
    real(real32), parameter :: u_max = 0.1_real32 ! initial velocity
    real(real32), parameter :: n = 2.0_real32 ! num waves
    real(real32), parameter :: k = (2.0_real32 * pi * n) / real(N_Y, real32)

    ! allocate sim data structures
    real(real32), allocatable :: f(:, :, :) ! distribution functions as f(dir, x, y)
    real(real32), allocatable :: rho(:,:)
    real(real32), allocatable :: u_x(:,:)
    real(real32), allocatable :: u_y(:,:)
    allocate(f(N_DIRS, N_X, N_Y))
    allocate(rho(N_X, N_Y))
    allocate(u_x(N_X, N_Y))
    allocate(u_y(N_X, N_Y))

    ! inital condition
    call apply_condition_shear_wave_decay(N_X, N_Y, N_DIRS, c_x_fp, c_y_fp, w, rho_0, u_max, k, f, rho, u_x, u_y)

end program main
