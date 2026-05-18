module exchange
    ! imports
    use iso_fortran_env, only: int32
    use domain, only: domain_t
    use settings, only: N_DIRS, FP, USE_SCALAR_COLUMN_REMOTE_PUTS
    implicit none
    private

    public :: halo_buffers_t
    public :: allocate_halo_buffers
    public :: exchange_halos
    public :: exchange_halos_direct_put
    public :: finish_halos_direct_put
    public :: exchange_halos_hybrid_start
    public :: exchange_halos_hybrid_finish

    integer(int32), parameter :: N_HALO_DIRS = 3_int32

    type :: halo_buffers_t

        ! x-direction send buffers
        real(FP), allocatable :: send_left(:,:)[:]
        real(FP), allocatable :: send_right(:,:)[:]
        real(FP), allocatable :: send_left_hybrid(:,:,:)[:]
        real(FP), allocatable :: send_right_hybrid(:,:,:)[:]

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
        allocate(halo_buffers%send_left_hybrid(domain_info%n_y_local, N_HALO_DIRS, 2)[*])
        allocate(halo_buffers%send_right_hybrid(domain_info%n_y_local, N_HALO_DIRS, 2)[*])
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
        integer(int32) :: x
        integer(int32) :: y

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
        do y = 1, n_y_local
            halo_buffers%send_left(y, 1) = f(1, y, 4)
            halo_buffers%send_left(y, 2) = f(1, y, 7)
            halo_buffers%send_left(y, 3) = f(1, y, 8)

            halo_buffers%send_right(y, 1) = f(n_x_local, y, 2)
            halo_buffers%send_right(y, 2) = f(n_x_local, y, 6)
            halo_buffers%send_right(y, 3) = f(n_x_local, y, 9)
        end do

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

        ! pack bottom/top borders, including updated x-halos carrying corner halo values
        do x = 0, n_x_local + 1
            halo_buffers%send_bottom(x, 1) = f(x, 1, 5)
            halo_buffers%send_bottom(x, 2) = f(x, 1, 8)
            halo_buffers%send_bottom(x, 3) = f(x, 1, 9)

            halo_buffers%send_top(x, 1) = f(x, n_y_local, 3)
            halo_buffers%send_top(x, 2) = f(x, n_y_local, 6)
            halo_buffers%send_top(x, 3) = f(x, n_y_local, 7)
        end do

        call sync_neighbor_images(y_neighbor_images, n_y_neighbor_images)

        ! unpack bottom/top halos from neighboring images
        halo_buffers%recv_bottom(:, :) = halo_buffers%send_top(:, :)[domain_info%bottom_image_id]
        halo_buffers%recv_top(:, :) = halo_buffers%send_bottom(:, :)[domain_info%top_image_id]
 
        do x = 0, n_x_local + 1
            f(x, 0, 3) = halo_buffers%recv_bottom(x, 1)
            f(x, 0, 6) = halo_buffers%recv_bottom(x, 2)
            f(x, 0, 7) = halo_buffers%recv_bottom(x, 3)

            f(x, n_y_local+1, 5) = halo_buffers%recv_top(x, 1)
            f(x, n_y_local+1, 8) = halo_buffers%recv_top(x, 2)
            f(x, n_y_local+1, 9) = halo_buffers%recv_top(x, 3)
        end do

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


    subroutine exchange_halos_direct_put( &
        domain_info, n_x_local, n_y_local, active_buffer, f &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        integer(int32), intent(in) :: active_buffer

        ! read/write inputs
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS, 2)[*]

        if (domain_info%n_images == 1) then
            call exchange_halos_direct_local(n_x_local, n_y_local, f(:, :, :, active_buffer))
            return
        end if

        if (USE_SCALAR_COLUMN_REMOTE_PUTS) then
            call exchange_x_halos_direct_scalar_put( &
                domain_info, n_x_local, n_y_local, active_buffer, f)
        else
            ! write left/right halos on neighboring images directly from owned borders
            f(n_x_local+1, 1:n_y_local, 4, active_buffer)[domain_info%left_image_id] = &
                f(1, 1:n_y_local, 4, active_buffer)
            f(0, 1:n_y_local, 2, active_buffer)[domain_info%right_image_id] = &
                f(n_x_local, 1:n_y_local, 2, active_buffer)

            if (n_y_local > 1) then
                f(n_x_local+1, 1:n_y_local-1, 7, active_buffer)[domain_info%left_image_id] = &
                    f(1, 1:n_y_local-1, 7, active_buffer)
                f(n_x_local+1, 2:n_y_local, 8, active_buffer)[domain_info%left_image_id] = &
                    f(1, 2:n_y_local, 8, active_buffer)

                f(0, 1:n_y_local-1, 6, active_buffer)[domain_info%right_image_id] = &
                    f(n_x_local, 1:n_y_local-1, 6, active_buffer)
                f(0, 2:n_y_local, 9, active_buffer)[domain_info%right_image_id] = &
                    f(n_x_local, 2:n_y_local, 9, active_buffer)
            end if
        end if

        ! write bottom/top halos on neighboring images directly from owned borders
        f(1:n_x_local, n_y_local+1, 5, active_buffer)[domain_info%bottom_image_id] = &
            f(1:n_x_local, 1, 5, active_buffer)
        f(1:n_x_local, 0, 3, active_buffer)[domain_info%top_image_id] = &
            f(1:n_x_local, n_y_local, 3, active_buffer)

        if (n_x_local > 1) then
            f(2:n_x_local, n_y_local+1, 8, active_buffer)[domain_info%bottom_image_id] = &
                f(2:n_x_local, 1, 8, active_buffer)
            f(1:n_x_local-1, n_y_local+1, 9, active_buffer)[domain_info%bottom_image_id] = &
                f(1:n_x_local-1, 1, 9, active_buffer)

            f(1:n_x_local-1, 0, 6, active_buffer)[domain_info%top_image_id] = &
                f(1:n_x_local-1, n_y_local, 6, active_buffer)
            f(2:n_x_local, 0, 7, active_buffer)[domain_info%top_image_id] = &
                f(2:n_x_local, n_y_local, 7, active_buffer)
        end if

        ! write corner halos directly to diagonal neighbors
        f(0, 0, 6, active_buffer)[domain_info%top_right_image_id] = &
            f(n_x_local, n_y_local, 6, active_buffer)
        f(n_x_local+1, 0, 7, active_buffer)[domain_info%top_left_image_id] = &
            f(1, n_y_local, 7, active_buffer)
        f(n_x_local+1, n_y_local+1, 8, active_buffer)[domain_info%bottom_left_image_id] = &
            f(1, 1, 8, active_buffer)
        f(0, n_y_local+1, 9, active_buffer)[domain_info%bottom_right_image_id] = &
            f(n_x_local, 1, 9, active_buffer)
    end subroutine exchange_halos_direct_put


    subroutine exchange_x_halos_direct_scalar_put( &
        domain_info, n_x_local, n_y_local, active_buffer, f &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        integer(int32), intent(in) :: active_buffer

        ! read/write inputs
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS, 2)[*]

        ! locals
        integer(int32) :: y

        ! write left halos on neighboring images directly from owned borders
        do y = 1, n_y_local
            f(n_x_local+1, y, 4, active_buffer)[domain_info%left_image_id] = &
                f(1, y, 4, active_buffer)
        end do

        do y = 1, n_y_local - 1
            f(n_x_local+1, y, 7, active_buffer)[domain_info%left_image_id] = &
                f(1, y, 7, active_buffer)
        end do

        do y = 2, n_y_local
            f(n_x_local+1, y, 8, active_buffer)[domain_info%left_image_id] = &
                f(1, y, 8, active_buffer)
        end do

        ! write right halos on neighboring images directly from owned borders
        do y = 1, n_y_local
            f(0, y, 2, active_buffer)[domain_info%right_image_id] = &
                f(n_x_local, y, 2, active_buffer)
        end do

        do y = 1, n_y_local - 1
            f(0, y, 6, active_buffer)[domain_info%right_image_id] = &
                f(n_x_local, y, 6, active_buffer)
        end do

        do y = 2, n_y_local
            f(0, y, 9, active_buffer)[domain_info%right_image_id] = &
                f(n_x_local, y, 9, active_buffer)
        end do
    end subroutine exchange_x_halos_direct_scalar_put


    subroutine finish_halos_direct_put( &
        domain_info &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info

        ! locals
        integer(int32) :: n_neighbor_images
        integer(int32) :: neighbor_images(8)

        if (domain_info%n_images == 1) then
            return
        end if

        n_neighbor_images = 0
        call append_unique_image(domain_info%left_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%right_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%bottom_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%top_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%bottom_left_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%bottom_right_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%top_left_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%top_right_image_id, neighbor_images, n_neighbor_images)

        if (n_neighbor_images == 1) then
            sync images(neighbor_images(1))
        else
            sync images(neighbor_images(1:n_neighbor_images))
        end if
    end subroutine finish_halos_direct_put


    subroutine exchange_halos_hybrid_start( &
        domain_info, halo_buffers, n_x_local, n_y_local, active_buffer, f &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        integer(int32), intent(in) :: active_buffer

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS, 2)[*]

        ! locals
        integer(int32) :: y

        if (domain_info%n_images == 1) then
            call exchange_halos_direct_local(n_x_local, n_y_local, f(:, :, :, active_buffer))
            return
        end if

        ! pack strided left/right columns into compact x-direction buffers
        do y = 1, n_y_local
            halo_buffers%send_left_hybrid(y, 1, active_buffer) = f(1, y, 4, active_buffer)
            halo_buffers%send_left_hybrid(y, 2, active_buffer) = f(1, y, 7, active_buffer)
            halo_buffers%send_left_hybrid(y, 3, active_buffer) = f(1, y, 8, active_buffer)

            halo_buffers%send_right_hybrid(y, 1, active_buffer) = f(n_x_local, y, 2, active_buffer)
            halo_buffers%send_right_hybrid(y, 2, active_buffer) = f(n_x_local, y, 6, active_buffer)
            halo_buffers%send_right_hybrid(y, 3, active_buffer) = f(n_x_local, y, 9, active_buffer)
        end do

        ! write contiguous bottom/top rows on neighboring images directly from owned borders
        f(1:n_x_local, n_y_local+1, 5, active_buffer)[domain_info%bottom_image_id] = &
            f(1:n_x_local, 1, 5, active_buffer)
        f(1:n_x_local, 0, 3, active_buffer)[domain_info%top_image_id] = &
            f(1:n_x_local, n_y_local, 3, active_buffer)

        if (n_x_local > 1) then
            f(2:n_x_local, n_y_local+1, 8, active_buffer)[domain_info%bottom_image_id] = &
                f(2:n_x_local, 1, 8, active_buffer)
            f(1:n_x_local-1, n_y_local+1, 9, active_buffer)[domain_info%bottom_image_id] = &
                f(1:n_x_local-1, 1, 9, active_buffer)

            f(1:n_x_local-1, 0, 6, active_buffer)[domain_info%top_image_id] = &
                f(1:n_x_local-1, n_y_local, 6, active_buffer)
            f(2:n_x_local, 0, 7, active_buffer)[domain_info%top_image_id] = &
                f(2:n_x_local, n_y_local, 7, active_buffer)
        end if

        ! write corner halos directly to diagonal neighbors
        f(0, 0, 6, active_buffer)[domain_info%top_right_image_id] = &
            f(n_x_local, n_y_local, 6, active_buffer)
        f(n_x_local+1, 0, 7, active_buffer)[domain_info%top_left_image_id] = &
            f(1, n_y_local, 7, active_buffer)
        f(n_x_local+1, n_y_local+1, 8, active_buffer)[domain_info%bottom_left_image_id] = &
            f(1, 1, 8, active_buffer)
        f(0, n_y_local+1, 9, active_buffer)[domain_info%bottom_right_image_id] = &
            f(n_x_local, 1, 9, active_buffer)
    end subroutine exchange_halos_hybrid_start


    subroutine exchange_halos_hybrid_finish( &
        domain_info, halo_buffers, n_x_local, n_y_local, active_buffer, f &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        integer(int32), intent(in) :: active_buffer

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS, 2)[*]

        ! locals
        integer(int32) :: y
        integer(int32) :: n_neighbor_images
        integer(int32) :: neighbor_images(8)

        if (domain_info%n_images == 1) then
            return
        end if

        n_neighbor_images = 0
        call append_unique_image(domain_info%left_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%right_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%bottom_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%top_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%bottom_left_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%bottom_right_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%top_left_image_id, neighbor_images, n_neighbor_images)
        call append_unique_image(domain_info%top_right_image_id, neighbor_images, n_neighbor_images)

        if (n_neighbor_images == 1) then
            sync images(neighbor_images(1))
        else
            sync images(neighbor_images(1:n_neighbor_images))
        end if

        ! unpack left/right halos from compact neighboring x-direction buffers
        halo_buffers%recv_left(:, :) = halo_buffers%send_right_hybrid(:, :, active_buffer)[domain_info%left_image_id]
        halo_buffers%recv_right(:, :) = halo_buffers%send_left_hybrid(:, :, active_buffer)[domain_info%right_image_id]

        do y = 1, n_y_local
            f(0, y, 2, active_buffer) = halo_buffers%recv_left(y, 1)
            f(0, y, 6, active_buffer) = halo_buffers%recv_left(y, 2)
            f(0, y, 9, active_buffer) = halo_buffers%recv_left(y, 3)

            f(n_x_local+1, y, 4, active_buffer) = halo_buffers%recv_right(y, 1)
            f(n_x_local+1, y, 7, active_buffer) = halo_buffers%recv_right(y, 2)
            f(n_x_local+1, y, 8, active_buffer) = halo_buffers%recv_right(y, 3)
        end do
    end subroutine exchange_halos_hybrid_finish


    subroutine append_unique_image( &
        image_id, neighbor_images, n_neighbor_images &
        )
        ! inputs
        integer(int32), intent(in) :: image_id

        ! read/write inputs
        integer(int32), intent(inout) :: neighbor_images(8)
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


    subroutine exchange_halos_direct_local( &
        n_x_local, n_y_local, f &
        )
        ! inputs
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local

        ! read/write inputs
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)

        ! left/right periodic halos
        f(n_x_local+1, 1:n_y_local, 4) = f(1, 1:n_y_local, 4)
        f(0, 1:n_y_local, 2) = f(n_x_local, 1:n_y_local, 2)

        if (n_y_local > 1) then
            f(n_x_local+1, 1:n_y_local-1, 7) = f(1, 1:n_y_local-1, 7)
            f(n_x_local+1, 2:n_y_local, 8) = f(1, 2:n_y_local, 8)

            f(0, 1:n_y_local-1, 6) = f(n_x_local, 1:n_y_local-1, 6)
            f(0, 2:n_y_local, 9) = f(n_x_local, 2:n_y_local, 9)
        end if

        ! bottom/top periodic halos
        f(1:n_x_local, n_y_local+1, 5) = f(1:n_x_local, 1, 5)
        f(1:n_x_local, 0, 3) = f(1:n_x_local, n_y_local, 3)

        if (n_x_local > 1) then
            f(2:n_x_local, n_y_local+1, 8) = f(2:n_x_local, 1, 8)
            f(1:n_x_local-1, n_y_local+1, 9) = f(1:n_x_local-1, 1, 9)

            f(1:n_x_local-1, 0, 6) = f(1:n_x_local-1, n_y_local, 6)
            f(2:n_x_local, 0, 7) = f(2:n_x_local, n_y_local, 7)
        end if

        ! corner periodic halos
        f(0, 0, 6) = f(n_x_local, n_y_local, 6)
        f(n_x_local+1, 0, 7) = f(1, n_y_local, 7)
        f(n_x_local+1, n_y_local+1, 8) = f(1, 1, 8)
        f(0, n_y_local+1, 9) = f(n_x_local, 1, 9)
    end subroutine exchange_halos_direct_local

end module exchange
