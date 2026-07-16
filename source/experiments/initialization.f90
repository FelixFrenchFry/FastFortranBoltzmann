module initialization
    ! imports
    use iso_fortran_env, only: int32
    use settings, only: N_DIRS, W, FP
    implicit none

contains

    subroutine apply_condition_sliding_lid_local( &
        rho_0, n_x_local, n_y_local, f, rho, u_x, u_y &
        )
        ! inputs
        real(FP), intent(in) :: rho_0
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local

        ! write destinations
        real(FP), intent(out) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)
        real(FP), intent(out) :: rho(n_x_local, n_y_local)
        real(FP), intent(out) :: u_x(n_x_local, n_y_local)
        real(FP), intent(out) :: u_y(n_x_local, n_y_local)

        ! temp
        integer(int32) :: x, y, i
        real(FP) :: f_eq(N_DIRS)

        f = 0.0_FP

        ! pre-computed equilibrium distribution at t=0
        do i = 1, N_DIRS
            f_eq(i) = W(i) * rho_0
        end do

        ! loop over all image-owned cells
        do y = 1, n_y_local
            do x = 1, n_x_local

                rho(x, y) = rho_0
                u_x(x, y) = 0.0_FP
                u_y(x, y) = 0.0_FP

                ! init distribution functions
                do i = 1, N_DIRS
                    f(i, x, y) = f_eq(i)
                end do
            end do
        end do
    end subroutine apply_condition_sliding_lid_local


end module initialization
