module exchange
    ! imports
    use iso_fortran_env, only: int32
    use domain, only: domain_t
    use settings, only: N_DIRS, FP, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID
    implicit none
    private

    public :: halo_buffers_t
    public :: exchange_plan_t
    public :: allocate_halo_buffers
    public :: BUF_SEND_LEFT
    public :: BUF_SEND_RIGHT
    public :: BUF_MACRO_LEFT
    public :: BUF_MACRO_RIGHT

    ! halo buffer indices in the exchange buffer window
    integer(int32), parameter :: BUF_SEND_LEFT = 1
    integer(int32), parameter :: BUF_SEND_RIGHT = 2
    integer(int32), parameter :: BUF_MACRO_LEFT = 3
    integer(int32), parameter :: BUF_MACRO_RIGHT = 4
    public :: build_exchange_plan
    public :: exchange_halos

    type :: halo_buffers_t

        ! distribution function send buffers for staged left/right halo exchange
        real(FP), allocatable :: window(:,:,:)[:]

        ! distribution function receive buffers for staged left/right halo exchange
        real(FP), allocatable :: recv_left(:,:)
        real(FP), allocatable :: recv_right(:,:)
        
        ! left/right macro field receive buffers
        real(FP), allocatable :: recv_macro_left(:,:)
        real(FP), allocatable :: recv_macro_right(:,:)

    end type halo_buffers_t

    type :: exchange_plan_t

        ! distribution function halo exchange
        logical :: left
        logical :: right
        logical :: bottom
        logical :: top

        ! macro exchange for poiseuille flow
        logical :: macro_left
        logical :: macro_right

    end type exchange_plan_t

contains

    subroutine build_exchange_plan( &
        domain_info, sim_mode, exchange_plan &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: sim_mode

        ! output
        type(exchange_plan_t), intent(out) :: exchange_plan

        exchange_plan%left = .false.
        exchange_plan%right = .false.
        exchange_plan%bottom = .false.
        exchange_plan%top = .false.
        exchange_plan%macro_left = .false.
        exchange_plan%macro_right = .false.

        select case (sim_mode)
        case (SIM_SHEAR_WAVE)
            exchange_plan%left = .true.
            exchange_plan%right = .true.
            exchange_plan%bottom = .true.
            exchange_plan%top = .true.

        case (SIM_COUETTE_FLOW)
            exchange_plan%left = .true.
            exchange_plan%right = .true.
            exchange_plan%bottom = .not. domain_info%at_bottom_boundary
            exchange_plan%top = .not. domain_info%at_top_boundary

        case (SIM_POISEUILLE_FLOW)
            exchange_plan%left = .not. domain_info%at_left_boundary
            exchange_plan%right = .not. domain_info%at_right_boundary
            exchange_plan%bottom = .not. domain_info%at_bottom_boundary
            exchange_plan%top = .not. domain_info%at_top_boundary
            exchange_plan%macro_left = domain_info%at_left_boundary
            exchange_plan%macro_right = domain_info%at_right_boundary

        case (SIM_SLIDING_LID)
            exchange_plan%left = .not. domain_info%at_left_boundary
            exchange_plan%right = .not. domain_info%at_right_boundary
            exchange_plan%bottom = .not. domain_info%at_bottom_boundary
            exchange_plan%top = .not. domain_info%at_top_boundary

        case default
            error stop "error: unknown sim mode in build_exchange_plan()"
        end select
    end subroutine build_exchange_plan


    subroutine allocate_halo_buffers( &
        domain_info, halo_buffers &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers

        if (allocated(halo_buffers%recv_macro_left)) then
            error stop "error: halo buffers are already allocated"
        end if

        allocate(halo_buffers%window(domain_info%n_y, 3, 4)[*])
        allocate(halo_buffers%recv_left(domain_info%n_y, 3))
        allocate(halo_buffers%recv_right(domain_info%n_y, 3))
        allocate(halo_buffers%recv_macro_left(domain_info%n_y, 3))
        allocate(halo_buffers%recv_macro_right(domain_info%n_y, 3))
    end subroutine allocate_halo_buffers


    subroutine exchange_halos( &
        domain_info, halo_buffers, n_x_local, n_y_local, f, exchange_plan &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        type(exchange_plan_t), intent(in) :: exchange_plan

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f(0:n_x_local+1, 0:n_y_local+1, N_DIRS)[*]

        ! locals
        integer(int32) :: n_x_neighbor_images
        integer(int32) :: n_y_neighbor_images
        integer(int32) :: x_neighbor_images(4)
        integer(int32) :: y_neighbor_images(2)

        n_x_neighbor_images = 0
        n_y_neighbor_images = 0

        if (exchange_plan%left) then
            call add_neighbor_image(x_neighbor_images, n_x_neighbor_images, domain_info%left_image_id)
        end if
        if (exchange_plan%right) then
            call add_neighbor_image(x_neighbor_images, n_x_neighbor_images, domain_info%right_image_id)
        end if
        if (exchange_plan%macro_left) then
            call add_neighbor_image(x_neighbor_images, n_x_neighbor_images, domain_info%left_image_id)
        end if
        if (exchange_plan%macro_right) then
            call add_neighbor_image(x_neighbor_images, n_x_neighbor_images, domain_info%right_image_id)
        end if
        if (exchange_plan%bottom) then
            call add_neighbor_image(y_neighbor_images, n_y_neighbor_images, domain_info%bottom_image_id)
        end if
        if (exchange_plan%top) then
            call add_neighbor_image(y_neighbor_images, n_y_neighbor_images, domain_info%top_image_id)
        end if

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! pack strided left/right boundaries into contiguous send buffers locally
        if (exchange_plan%left) then
            halo_buffers%window(:, 1, BUF_SEND_LEFT) = f(1, 1:n_y_local, 4)
            halo_buffers%window(:, 2, BUF_SEND_LEFT) = f(1, 1:n_y_local, 7)
            halo_buffers%window(:, 3, BUF_SEND_LEFT) = f(1, 1:n_y_local, 8)
        end if

        if (exchange_plan%right) then
            halo_buffers%window(:, 1, BUF_SEND_RIGHT) = f(n_x_local, 1:n_y_local, 2)
            halo_buffers%window(:, 2, BUF_SEND_RIGHT) = f(n_x_local, 1:n_y_local, 6)
            halo_buffers%window(:, 3, BUF_SEND_RIGHT) = f(n_x_local, 1:n_y_local, 9)
        end if

        call sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! exchange packed left/right halos through staging buffers
        if (exchange_plan%left) then
            halo_buffers%recv_left(:, :) = halo_buffers%window(:, :, BUF_SEND_RIGHT)[domain_info%left_image_id]
            f(0, 1:n_y_local, 2) = halo_buffers%recv_left(:, 1)
            f(0, 1:n_y_local, 6) = halo_buffers%recv_left(:, 2)
            f(0, 1:n_y_local, 9) = halo_buffers%recv_left(:, 3)
        end if

        if (exchange_plan%right) then
            halo_buffers%recv_right(:, :) = halo_buffers%window(:, :, BUF_SEND_LEFT)[domain_info%right_image_id]
            f(n_x_local+1, 1:n_y_local, 4) = halo_buffers%recv_right(:, 1)
            f(n_x_local+1, 1:n_y_local, 7) = halo_buffers%recv_right(:, 2)
            f(n_x_local+1, 1:n_y_local, 8) = halo_buffers%recv_right(:, 3)
        end if

        ! pressure-periodic macros for poiseuille flow
        if (exchange_plan%macro_left) then
            halo_buffers%recv_macro_left(:, :) = halo_buffers%window(:, :, BUF_MACRO_RIGHT)[domain_info%left_image_id]
        end if

        if (exchange_plan%macro_right) then
            halo_buffers%recv_macro_right(:, :) = halo_buffers%window(:, :, BUF_MACRO_LEFT)[domain_info%right_image_id]
        end if

        call sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)
        call sync_neighbor_images(y_neighbor_images, n_y_neighbor_images)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! exchange contiguous bottom/top halos directly, including corners
        if (exchange_plan%bottom) then
            f(0:n_x_local+1, 0, 3) = f(0:n_x_local+1, n_y_local, 3)[domain_info%bottom_image_id]
            f(0:n_x_local+1, 0, 6) = f(0:n_x_local+1, n_y_local, 6)[domain_info%bottom_image_id]
            f(0:n_x_local+1, 0, 7) = f(0:n_x_local+1, n_y_local, 7)[domain_info%bottom_image_id]
        end if

        if (exchange_plan%top) then
            f(0:n_x_local+1, n_y_local+1, 5) = f(0:n_x_local+1, 1, 5)[domain_info%top_image_id]
            f(0:n_x_local+1, n_y_local+1, 8) = f(0:n_x_local+1, 1, 8)[domain_info%top_image_id]
            f(0:n_x_local+1, n_y_local+1, 9) = f(0:n_x_local+1, 1, 9)[domain_info%top_image_id]
        end if

        call sync_neighbor_images(y_neighbor_images, n_y_neighbor_images)

    contains

        subroutine add_neighbor_image( &
            neighbor_images, n_neighbor_images, neighbor_image &
            )
            ! inputs
            integer(int32), intent(in) :: neighbor_image

            ! read/write inputs
            integer(int32), intent(inout) :: neighbor_images(:)
            integer(int32), intent(inout) :: n_neighbor_images

            ! locals
            integer(int32) :: i

            do i = 1, n_neighbor_images
                if (neighbor_images(i) == neighbor_image) then
                    return
                end if
            end do

            n_neighbor_images = n_neighbor_images + 1
            if (n_neighbor_images > size(neighbor_images)) then
                error stop "error: exchange neighbor list is too small"
            end if
            neighbor_images(n_neighbor_images) = neighbor_image
        end subroutine add_neighbor_image


        subroutine sync_neighbor_images( &
            neighbor_images, n_neighbor_images &
            )
            ! inputs
            integer(int32), intent(in) :: neighbor_images(:)
            integer(int32), intent(in) :: n_neighbor_images

            if (n_neighbor_images == 0) then
                return
            else if (n_neighbor_images == 1) then
                sync images(neighbor_images(1))
            else
                sync images(neighbor_images(1:n_neighbor_images))
            end if
        end subroutine sync_neighbor_images
    end subroutine exchange_halos


end module exchange
