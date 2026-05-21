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

        ! x-direction macro field send buffers
        real(FP), allocatable :: send_macro_left(:,:)[:]
        real(FP), allocatable :: send_macro_right(:,:)[:]

        ! x-direction macro field receive buffers
        real(FP), allocatable :: recv_macro_left(:,:)
        real(FP), allocatable :: recv_macro_right(:,:)

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
        allocate(halo_buffers%send_left(domain_info%n_y_local, 3)[*])
        allocate(halo_buffers%send_right(domain_info%n_y_local, 3)[*])
        allocate(halo_buffers%send_bottom(0:domain_info%n_x_local+1, 3)[*])
        allocate(halo_buffers%send_top(0:domain_info%n_x_local+1, 3)[*])

        allocate(halo_buffers%send_macro_left(domain_info%n_y_local, 3)[*])
        allocate(halo_buffers%send_macro_right(domain_info%n_y_local, 3)[*])
        allocate(halo_buffers%recv_macro_left(domain_info%n_y_local, 3))
        allocate(halo_buffers%recv_macro_right(domain_info%n_y_local, 3))
    end subroutine allocate_halo_buffers


    subroutine exchange_halos( &
        domain_info, halo_buffers, n_x_local, n_y_local, f, exchange_pressure_macros &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        logical, intent(in), optional :: exchange_pressure_macros

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
        logical :: do_exchange_pressure_macros

        ! check if left/right images are the same
        x_neighbor_images(1) = domain_info%left_image_id
        if (domain_info%right_image_id == domain_info%left_image_id) then
            n_x_neighbor_images = 1
        else
            x_neighbor_images(2) = domain_info%right_image_id
            n_x_neighbor_images = 2
        end if

        ! check if bottom/top images are the same
        y_neighbor_images(1) = domain_info%bottom_image_id
        if (domain_info%top_image_id == domain_info%bottom_image_id) then
            n_y_neighbor_images = 1
        else
            y_neighbor_images(2) = domain_info%top_image_id
            n_y_neighbor_images = 2
        end if

        do_exchange_pressure_macros = .false.
        if (present(exchange_pressure_macros)) then
            do_exchange_pressure_macros = exchange_pressure_macros
        end if

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
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
        f(0, 1:n_y_local, 2) = halo_buffers%send_right(:, 1)[domain_info%left_image_id]
        f(0, 1:n_y_local, 6) = halo_buffers%send_right(:, 2)[domain_info%left_image_id]
        f(0, 1:n_y_local, 9) = halo_buffers%send_right(:, 3)[domain_info%left_image_id]

        f(n_x_local+1, 1:n_y_local, 4) = halo_buffers%send_left(:, 1)[domain_info%right_image_id]
        f(n_x_local+1, 1:n_y_local, 7) = halo_buffers%send_left(:, 2)[domain_info%right_image_id]
        f(n_x_local+1, 1:n_y_local, 8) = halo_buffers%send_left(:, 3)[domain_info%right_image_id]

        ! pressure-periodic macro strips are maintained by the poiseuille kernels
        if (do_exchange_pressure_macros) then
            if (domain_info%at_left_boundary) then
                halo_buffers%recv_macro_left(:, :) = halo_buffers%send_macro_right(:, :)[domain_info%left_image_id]
            end if

            if (domain_info%at_right_boundary) then
                halo_buffers%recv_macro_right(:, :) = halo_buffers%send_macro_left(:, :)[domain_info%right_image_id]
            end if
        end if

        call sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
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
        f(0:n_x_local+1, 0, 3) = halo_buffers%send_top(0:n_x_local+1, 1)[domain_info%bottom_image_id]
        f(0:n_x_local+1, 0, 6) = halo_buffers%send_top(0:n_x_local+1, 2)[domain_info%bottom_image_id]
        f(0:n_x_local+1, 0, 7) = halo_buffers%send_top(0:n_x_local+1, 3)[domain_info%bottom_image_id]

        f(0:n_x_local+1, n_y_local+1, 5) = halo_buffers%send_bottom(0:n_x_local+1, 1)[domain_info%top_image_id]
        f(0:n_x_local+1, n_y_local+1, 8) = halo_buffers%send_bottom(0:n_x_local+1, 2)[domain_info%top_image_id]
        f(0:n_x_local+1, n_y_local+1, 9) = halo_buffers%send_bottom(0:n_x_local+1, 3)[domain_info%top_image_id]

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
