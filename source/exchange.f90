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

        ! corner send buffers
        real(FP), allocatable :: send_corners(:)[:]

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

        allocate(halo_buffers%send_left(domain_info%n_y_local, N_HALO_DIRS)[*])
        allocate(halo_buffers%send_right(domain_info%n_y_local, N_HALO_DIRS)[*])
        allocate(halo_buffers%send_bottom(domain_info%n_x_local, N_HALO_DIRS)[*])
        allocate(halo_buffers%send_top(domain_info%n_x_local, N_HALO_DIRS)[*])
        allocate(halo_buffers%send_corners(4)[*])

        allocate(halo_buffers%recv_left(domain_info%n_y_local, N_HALO_DIRS))
        allocate(halo_buffers%recv_right(domain_info%n_y_local, N_HALO_DIRS))
        allocate(halo_buffers%recv_bottom(domain_info%n_x_local, N_HALO_DIRS))
        allocate(halo_buffers%recv_top(domain_info%n_x_local, N_HALO_DIRS))
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
        integer(int32) :: n_y_corner_neighbor_images
        integer(int32) :: x_neighbor_images(2)
        integer(int32) :: y_corner_neighbor_images(6)
        integer(int32) :: x
        integer(int32) :: y

        ! left/right can be the same image for 2-image axes
        x_neighbor_images(1) = domain_info%left_image_id
        if (domain_info%right_image_id == domain_info%left_image_id) then
            n_x_neighbor_images = 1
        else
            x_neighbor_images(2) = domain_info%right_image_id
            n_x_neighbor_images = 2
        end if

        ! pack owned left/right borders
        do y = 1, n_y_local
            halo_buffers%send_left(y, 1) = f(1, y, 4)
            halo_buffers%send_left(y, 2) = f(1, y, 7)
            halo_buffers%send_left(y, 3) = f(1, y, 8)

            halo_buffers%send_right(y, 1) = f(n_x_local, y, 2)
            halo_buffers%send_right(y, 2) = f(n_x_local, y, 6)
            halo_buffers%send_right(y, 3) = f(n_x_local, y, 9)
        end do

        ! pack owned bottom/top borders
        do x = 1, n_x_local
            halo_buffers%send_bottom(x, 1) = f(x, 1, 5)
            halo_buffers%send_bottom(x, 2) = f(x, 1, 8)
            halo_buffers%send_bottom(x, 3) = f(x, 1, 9)

            halo_buffers%send_top(x, 1) = f(x, n_y_local, 3)
            halo_buffers%send_top(x, 2) = f(x, n_y_local, 6)
            halo_buffers%send_top(x, 3) = f(x, n_y_local, 7)
        end do

        ! pack owned diagonal corner values
        halo_buffers%send_corners(1) = f(n_x_local, n_y_local, 6)
        halo_buffers%send_corners(2) = f(1, n_y_local, 7)
        halo_buffers%send_corners(3) = f(1, 1, 8)
        halo_buffers%send_corners(4) = f(n_x_local, 1, 9)

        n_y_corner_neighbor_images = 0
        call append_unique_image(domain_info%bottom_image_id, y_corner_neighbor_images, n_y_corner_neighbor_images)
        call append_unique_image(domain_info%top_image_id, y_corner_neighbor_images, n_y_corner_neighbor_images)
        call append_unique_image(domain_info%bottom_left_image_id, y_corner_neighbor_images, n_y_corner_neighbor_images)
        call append_unique_image(domain_info%bottom_right_image_id, y_corner_neighbor_images, n_y_corner_neighbor_images)
        call append_unique_image(domain_info%top_left_image_id, y_corner_neighbor_images, n_y_corner_neighbor_images)
        call append_unique_image(domain_info%top_right_image_id, y_corner_neighbor_images, n_y_corner_neighbor_images)

        call sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! unpack left/right halos from neighboring images
        halo_buffers%recv_left(:, :) = halo_buffers%send_right(:, :)[domain_info%left_image_id]
        halo_buffers%recv_right(:, :) = halo_buffers%send_left(:, :)[domain_info%right_image_id]

        do y = 1, n_y_local
            f(0, y, 2) = halo_buffers%recv_left(y, 1)
            f(0, y, 6) = halo_buffers%recv_left(y, 2)
            f(0, y, 9) = halo_buffers%recv_left(y, 3)

            f(n_x_local+1, y, 4) = halo_buffers%recv_right(y, 1)
            f(n_x_local+1, y, 7) = halo_buffers%recv_right(y, 2)
            f(n_x_local+1, y, 8) = halo_buffers%recv_right(y, 3)
        end do

        call sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)
        call sync_neighbor_images(y_corner_neighbor_images, n_y_corner_neighbor_images)

        ! unpack bottom/top halos from neighboring images
        halo_buffers%recv_bottom(:, :) = halo_buffers%send_top(:, :)[domain_info%bottom_image_id]
        halo_buffers%recv_top(:, :) = halo_buffers%send_bottom(:, :)[domain_info%top_image_id]
 
        do x = 1, n_x_local
            f(x, 0, 3) = halo_buffers%recv_bottom(x, 1)
            f(x, 0, 6) = halo_buffers%recv_bottom(x, 2)
            f(x, 0, 7) = halo_buffers%recv_bottom(x, 3)

            f(x, n_y_local+1, 5) = halo_buffers%recv_top(x, 1)
            f(x, n_y_local+1, 8) = halo_buffers%recv_top(x, 2)
            f(x, n_y_local+1, 9) = halo_buffers%recv_top(x, 3)
        end do

        ! unpack diagonal corner halos from neighboring images
        f(0, 0, 6) = halo_buffers%send_corners(1)[domain_info%bottom_left_image_id]
        f(n_x_local+1, 0, 7) = halo_buffers%send_corners(2)[domain_info%bottom_right_image_id]
        f(n_x_local+1, n_y_local+1, 8) = halo_buffers%send_corners(3)[domain_info%top_right_image_id]
        f(0, n_y_local+1, 9) = halo_buffers%send_corners(4)[domain_info%top_left_image_id]

        call sync_neighbor_images(y_corner_neighbor_images, n_y_corner_neighbor_images)
    end subroutine exchange_halos


    subroutine append_unique_image( &
        image_id, neighbor_images, n_neighbor_images &
        )
        ! inputs
        integer(int32), intent(in) :: image_id

        ! read/write inputs
        integer(int32), intent(inout) :: neighbor_images(:)
        integer(int32), intent(inout) :: n_neighbor_images

        ! locals
        integer(int32) :: i

        do i = 1, n_neighbor_images
            if (neighbor_images(i) == image_id) then
                return
            end if
        end do

        n_neighbor_images = n_neighbor_images + 1
        neighbor_images(n_neighbor_images) = image_id
    end subroutine append_unique_image


    subroutine sync_neighbor_images( &
        neighbor_images, n_neighbor_images &
        )
        ! inputs
        integer(int32), intent(in) :: neighbor_images(:)
        integer(int32), intent(in) :: n_neighbor_images

        if (n_neighbor_images == 1) then
            sync images(neighbor_images(1))
        else
            sync images(neighbor_images(1:n_neighbor_images))
        end if
    end subroutine sync_neighbor_images

end module exchange
