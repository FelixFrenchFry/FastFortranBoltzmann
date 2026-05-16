module exchange
    ! imports
    use iso_fortran_env, only: int32
    use domain, only: domain_t
    use settings, only: N_DIRS, FP
    implicit none
    private

    public :: halo_buffers_t
    public :: allocate_halo_buffers
    public :: exchange_halos

    type :: halo_buffers_t

        ! x-direction send buffers
        real(FP), allocatable :: send_left(:,:)[:]
        real(FP), allocatable :: send_right(:,:)[:]

        ! y-direction send buffers
        real(FP), allocatable :: send_bottom(:,:)[:]
        real(FP), allocatable :: send_top(:,:)[:]

    end type halo_buffers_t

contains

    subroutine allocate_halo_buffers( &
        domain_info, halo_buffers &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers

        if (allocated(halo_buffers%send_left)) then
            error stop "error: halo buffers are already allocated"
        end if

        ! bottom/top buffers include corners
        allocate(halo_buffers%send_left(domain_info%n_y_local, N_DIRS)[*])
        allocate(halo_buffers%send_right(domain_info%n_y_local, N_DIRS)[*])
        allocate(halo_buffers%send_bottom(0:domain_info%n_x_local+1, N_DIRS)[*])
        allocate(halo_buffers%send_top(0:domain_info%n_x_local+1, N_DIRS)[*])
    end subroutine allocate_halo_buffers


    subroutine exchange_halos( &
        domain_info, halo_buffers, n_x_local, n_y_local, f &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! pack owned left/right borders
        halo_buffers%send_left(:, :) = f(1, 1:n_y_local, :)
        halo_buffers%send_right(:, :) = f(n_x_local, 1:n_y_local, :)

        sync all

        ! unpack left/right halos from neighboring images
        f(0, 1:n_y_local, :) = halo_buffers%send_right(:, :)[domain_info%left_image_id]
        f(n_x_local+1, 1:n_y_local, :) = halo_buffers%send_left(:, :)[domain_info%right_image_id]

        sync all

        ! pack bottom/top borders, including updated x-halos carrying corner halo values
        halo_buffers%send_bottom(:, :) = f(0:n_x_local+1, 1, :)
        halo_buffers%send_top(:, :) = f(0:n_x_local+1, n_y_local, :)

        sync all

        ! unpack bottom/top halos from neighboring images
        f(0:n_x_local+1, 0, :) = halo_buffers%send_top(:, :)[domain_info%bottom_image_id]
        f(0:n_x_local+1, n_y_local+1, :) = halo_buffers%send_bottom(:, :)[domain_info%top_image_id]

        sync all
    end subroutine exchange_halos

end module exchange
