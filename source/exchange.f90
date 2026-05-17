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

    integer(int32), parameter :: N_HALO_DIRS = 3_int32

    type :: halo_buffers_t

        ! x-direction send buffers
        real(FP), allocatable :: send_left(:,:)[:]
        real(FP), allocatable :: send_right(:,:)[:]

        ! y-direction send buffers
        real(FP), allocatable :: send_bottom(:,:)[:]
        real(FP), allocatable :: send_top(:,:)[:]

        ! x-direction receive buffers
        real(FP), allocatable :: recv_left(:,:)
        real(FP), allocatable :: recv_right(:,:)

        ! y-direction receive buffers
        real(FP), allocatable :: recv_bottom(:,:)
        real(FP), allocatable :: recv_top(:,:)

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
        allocate(halo_buffers%send_left(domain_info%n_y_local, N_HALO_DIRS)[*])
        allocate(halo_buffers%send_right(domain_info%n_y_local, N_HALO_DIRS)[*])
        allocate(halo_buffers%send_bottom(0:domain_info%n_x_local+1, N_HALO_DIRS)[*])
        allocate(halo_buffers%send_top(0:domain_info%n_x_local+1, N_HALO_DIRS)[*])

        allocate(halo_buffers%recv_left(domain_info%n_y_local, N_HALO_DIRS))
        allocate(halo_buffers%recv_right(domain_info%n_y_local, N_HALO_DIRS))
        allocate(halo_buffers%recv_bottom(0:domain_info%n_x_local+1, N_HALO_DIRS))
        allocate(halo_buffers%recv_top(0:domain_info%n_x_local+1, N_HALO_DIRS))
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

        ! locals
        integer(int32) :: n_x_neighbor_images
        integer(int32) :: n_y_neighbor_images
        integer(int32) :: x_neighbor_images(2)
        integer(int32) :: y_neighbor_images(2)

        ! left/right or bottom/top can be the same image for 2-image axes
        x_neighbor_images(1) = domain_info%left_image_id
        if (domain_info%right_image_id == domain_info%left_image_id) then
            n_x_neighbor_images = 1
        else
            x_neighbor_images(2) = domain_info%right_image_id
            n_x_neighbor_images = 2
        end if

        y_neighbor_images(1) = domain_info%bottom_image_id
        if (domain_info%top_image_id == domain_info%bottom_image_id) then
            n_y_neighbor_images = 1
        else
            y_neighbor_images(2) = domain_info%top_image_id
            n_y_neighbor_images = 2
        end if

        ! pack owned left/right borders
        halo_buffers%send_left(:, 1) = f(1, 1:n_y_local, 4)
        halo_buffers%send_left(:, 2) = f(1, 1:n_y_local, 7)
        halo_buffers%send_left(:, 3) = f(1, 1:n_y_local, 8)

        halo_buffers%send_right(:, 1) = f(n_x_local, 1:n_y_local, 2)
        halo_buffers%send_right(:, 2) = f(n_x_local, 1:n_y_local, 6)
        halo_buffers%send_right(:, 3) = f(n_x_local, 1:n_y_local, 9)

        call sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! unpack left/right halos from neighboring images
        halo_buffers%recv_left(:, :) = halo_buffers%send_right(:, :)[domain_info%left_image_id]
        halo_buffers%recv_right(:, :) = halo_buffers%send_left(:, :)[domain_info%right_image_id]

        f(0, 1:n_y_local, 2) = halo_buffers%recv_left(:, 1)
        f(0, 1:n_y_local, 6) = halo_buffers%recv_left(:, 2)
        f(0, 1:n_y_local, 9) = halo_buffers%recv_left(:, 3)

        f(n_x_local+1, 1:n_y_local, 4) = halo_buffers%recv_right(:, 1)
        f(n_x_local+1, 1:n_y_local, 7) = halo_buffers%recv_right(:, 2)
        f(n_x_local+1, 1:n_y_local, 8) = halo_buffers%recv_right(:, 3)

        call sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! pack bottom/top borders, including updated x-halos carrying corner halo values
        halo_buffers%send_bottom(:, 1) = f(0:n_x_local+1, 1, 5)
        halo_buffers%send_bottom(:, 2) = f(0:n_x_local+1, 1, 8)
        halo_buffers%send_bottom(:, 3) = f(0:n_x_local+1, 1, 9)

        halo_buffers%send_top(:, 1) = f(0:n_x_local+1, n_y_local, 3)
        halo_buffers%send_top(:, 2) = f(0:n_x_local+1, n_y_local, 6)
        halo_buffers%send_top(:, 3) = f(0:n_x_local+1, n_y_local, 7)

        call sync_neighbor_images(y_neighbor_images, n_y_neighbor_images)

        ! unpack bottom/top halos from neighboring images
        halo_buffers%recv_bottom(:, :) = halo_buffers%send_top(:, :)[domain_info%bottom_image_id]
        halo_buffers%recv_top(:, :) = halo_buffers%send_bottom(:, :)[domain_info%top_image_id]

        f(0:n_x_local+1, 0, 3) = halo_buffers%recv_bottom(:, 1)
        f(0:n_x_local+1, 0, 6) = halo_buffers%recv_bottom(:, 2)
        f(0:n_x_local+1, 0, 7) = halo_buffers%recv_bottom(:, 3)

        f(0:n_x_local+1, n_y_local+1, 5) = halo_buffers%recv_top(:, 1)
        f(0:n_x_local+1, n_y_local+1, 8) = halo_buffers%recv_top(:, 2)
        f(0:n_x_local+1, n_y_local+1, 9) = halo_buffers%recv_top(:, 3)

        call sync_neighbor_images(y_neighbor_images, n_y_neighbor_images)
    contains

        subroutine sync_neighbor_images( &
            neighbor_images, n_neighbor_images &
            )
            ! inputs
            integer(int32), intent(in) :: neighbor_images(2)
            integer(int32), intent(in) :: n_neighbor_images

            if (n_neighbor_images == 1) then
                sync images(neighbor_images(1))
            else
                sync images(neighbor_images)
            end if
        end subroutine sync_neighbor_images
    end subroutine exchange_halos

end module exchange
