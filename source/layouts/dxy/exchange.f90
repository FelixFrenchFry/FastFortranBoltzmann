module exchange
    ! imports
    use iso_fortran_env, only: int32, int64, real64
    use domain, only: domain_t
    use settings, only: N_DIRS, FP, &
        SIM_SHEAR_WAVE, SIM_COUETTE_FLOW, SIM_POISEUILLE_FLOW, SIM_SLIDING_LID
    implicit none
    private

    public :: halo_buffers_t
    public :: exchange_plan_t
    public :: exchange_timing_t
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
    integer(int32), parameter :: BUF_SEND_BOTTOM = 1
    integer(int32), parameter :: BUF_SEND_TOP = 2
    public :: build_exchange_plan
    public :: exchange_halos
    public :: exchange_poiseuille_macro_halos

    type :: halo_buffers_t

        ! distribution function send buffers for staged left/right halo exchange
        real(FP), allocatable :: window(:,:,:)[:]

        ! distribution function send buffers for staged bottom/top halo exchange
        real(FP), allocatable :: window_y(:,:,:)[:]

        ! distribution function receive buffers for staged left/right halo exchange
        real(FP), allocatable :: recv_left(:,:)
        real(FP), allocatable :: recv_right(:,:)
        real(FP), allocatable :: recv_bottom(:,:)
        real(FP), allocatable :: recv_top(:,:)


        real(FP), allocatable :: recv_macro_left(:,:)
        real(FP), allocatable :: recv_macro_right(:,:)

    end type halo_buffers_t

    type :: exchange_plan_t

        ! distribution function halo exchange
        logical :: left
        logical :: right
        logical :: bottom
        logical :: top

        ! pressure-periodic macro halos for poiseuille flow
        logical :: macro_left
        logical :: macro_right

    end type exchange_plan_t

    type :: exchange_timing_t

        real(real64) :: halo_sync_seconds
        real(real64) :: halo_transfer_seconds

    end type exchange_timing_t

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
            exchange_plan%left = .true.
            exchange_plan%right = .true.
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
        allocate(halo_buffers%window_y(0:domain_info%n_x+1, 3, 2)[*])
        allocate(halo_buffers%recv_left(domain_info%n_y, 3))
        allocate(halo_buffers%recv_right(domain_info%n_y, 3))
        allocate(halo_buffers%recv_bottom(0:domain_info%n_x+1, 3))
        allocate(halo_buffers%recv_top(0:domain_info%n_x+1, 3))
        allocate(halo_buffers%recv_macro_left(domain_info%n_y, 3))
        allocate(halo_buffers%recv_macro_right(domain_info%n_y, 3))
    end subroutine allocate_halo_buffers


    subroutine exchange_halos( &
        domain_info, halo_buffers, n_x_local, n_y_local, f, exchange_plan, clock_rate, exchange_timing &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        integer(int32), intent(in) :: n_x_local
        integer(int32), intent(in) :: n_y_local
        type(exchange_plan_t), intent(in) :: exchange_plan
        integer(int64), intent(in) :: clock_rate

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers
        real(FP), intent(inout) :: f(N_DIRS, 0:n_x_local+1, 0:n_y_local+1)[*]

        ! output
        type(exchange_timing_t), intent(out) :: exchange_timing

        ! locals
        integer(int64) :: clock_section_start
        integer(int64) :: clock_section_end
        integer(int32) :: n_x_neighbor_images
        integer(int32) :: n_y_neighbor_images
        integer(int32) :: x_neighbor_images(4)
        integer(int32) :: y_neighbor_images(2)

        exchange_timing%halo_sync_seconds = 0.0_real64
        exchange_timing%halo_transfer_seconds = 0.0_real64

        call system_clock(clock_section_start)

        n_x_neighbor_images = 0
        n_y_neighbor_images = 0

        if (exchange_plan%left) then
            call add_neighbor_image(x_neighbor_images, n_x_neighbor_images, domain_info%left_image_id)
        end if
        if (exchange_plan%right) then
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
            halo_buffers%window(:, 1, BUF_SEND_LEFT) = f(4, 1, 1:n_y_local)
            halo_buffers%window(:, 2, BUF_SEND_LEFT) = f(7, 1, 1:n_y_local)
            halo_buffers%window(:, 3, BUF_SEND_LEFT) = f(8, 1, 1:n_y_local)
        end if

        if (exchange_plan%right) then
            halo_buffers%window(:, 1, BUF_SEND_RIGHT) = f(2, n_x_local, 1:n_y_local)
            halo_buffers%window(:, 2, BUF_SEND_RIGHT) = f(6, n_x_local, 1:n_y_local)
            halo_buffers%window(:, 3, BUF_SEND_RIGHT) = f(9, n_x_local, 1:n_y_local)
        end if

        call system_clock(clock_section_end)
        call add_elapsed_seconds(exchange_timing%halo_transfer_seconds, clock_section_start, clock_section_end)

        call timed_sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! exchange packed left/right halos through staging buffers
        call system_clock(clock_section_start)

        if (exchange_plan%left) then
            if (domain_info%left_image_id == domain_info%image_id) then ! read from own image when I_X=1
                halo_buffers%recv_left(:, :) = halo_buffers%window(:, :, BUF_SEND_RIGHT)
            else
                halo_buffers%recv_left(:, :) = halo_buffers%window(:, :, BUF_SEND_RIGHT)[domain_info%left_image_id]
            end if
            f(2, 0, 1:n_y_local) = halo_buffers%recv_left(:, 1)
            f(6, 0, 1:n_y_local) = halo_buffers%recv_left(:, 2)
            f(9, 0, 1:n_y_local) = halo_buffers%recv_left(:, 3)
        end if

        if (exchange_plan%right) then
            if (domain_info%right_image_id == domain_info%image_id) then ! read from own image when I_X=1
                halo_buffers%recv_right(:, :) = halo_buffers%window(:, :, BUF_SEND_LEFT)
            else
                halo_buffers%recv_right(:, :) = halo_buffers%window(:, :, BUF_SEND_LEFT)[domain_info%right_image_id]
            end if
            f(4, n_x_local+1, 1:n_y_local) = halo_buffers%recv_right(:, 1)
            f(7, n_x_local+1, 1:n_y_local) = halo_buffers%recv_right(:, 2)
            f(8, n_x_local+1, 1:n_y_local) = halo_buffers%recv_right(:, 3)
        end if

        call system_clock(clock_section_end)
        call add_elapsed_seconds(exchange_timing%halo_transfer_seconds, clock_section_start, clock_section_end)

        call timed_sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! pack strided bottom/top boundaries into contiguous send buffers locally
        call system_clock(clock_section_start)

        if (exchange_plan%bottom) then
            halo_buffers%window_y(:, 1, BUF_SEND_BOTTOM) = f(5, 0:n_x_local+1, 1)
            halo_buffers%window_y(:, 2, BUF_SEND_BOTTOM) = f(8, 0:n_x_local+1, 1)
            halo_buffers%window_y(:, 3, BUF_SEND_BOTTOM) = f(9, 0:n_x_local+1, 1)
        end if

        if (exchange_plan%top) then
            halo_buffers%window_y(:, 1, BUF_SEND_TOP) = f(3, 0:n_x_local+1, n_y_local)
            halo_buffers%window_y(:, 2, BUF_SEND_TOP) = f(6, 0:n_x_local+1, n_y_local)
            halo_buffers%window_y(:, 3, BUF_SEND_TOP) = f(7, 0:n_x_local+1, n_y_local)
        end if

        call system_clock(clock_section_end)
        call add_elapsed_seconds(exchange_timing%halo_transfer_seconds, clock_section_start, clock_section_end)

        call timed_sync_neighbor_images(y_neighbor_images, n_y_neighbor_images)

        ! exchange packed bottom/top halos through staging buffers, including corners
        call system_clock(clock_section_start)

        if (exchange_plan%bottom) then
            if (domain_info%bottom_image_id == domain_info%image_id) then ! read from own image when I_Y=1
                halo_buffers%recv_bottom(:, :) = halo_buffers%window_y(:, :, BUF_SEND_TOP)
            else
                halo_buffers%recv_bottom(:, :) = halo_buffers%window_y(:, :, BUF_SEND_TOP)[domain_info%bottom_image_id]
            end if
            f(3, 0:n_x_local+1, 0) = halo_buffers%recv_bottom(:, 1)
            f(6, 0:n_x_local+1, 0) = halo_buffers%recv_bottom(:, 2)
            f(7, 0:n_x_local+1, 0) = halo_buffers%recv_bottom(:, 3)
        end if

        if (exchange_plan%top) then
            if (domain_info%top_image_id == domain_info%image_id) then ! read from own image when I_Y=1
                halo_buffers%recv_top(:, :) = halo_buffers%window_y(:, :, BUF_SEND_BOTTOM)
            else
                halo_buffers%recv_top(:, :) = halo_buffers%window_y(:, :, BUF_SEND_BOTTOM)[domain_info%top_image_id]
            end if
            f(5, 0:n_x_local+1, n_y_local+1) = halo_buffers%recv_top(:, 1)
            f(8, 0:n_x_local+1, n_y_local+1) = halo_buffers%recv_top(:, 2)
            f(9, 0:n_x_local+1, n_y_local+1) = halo_buffers%recv_top(:, 3)
        end if

        call system_clock(clock_section_end)
        call add_elapsed_seconds(exchange_timing%halo_transfer_seconds, clock_section_start, clock_section_end)

        call timed_sync_neighbor_images(y_neighbor_images, n_y_neighbor_images)

    contains

        subroutine add_elapsed_seconds( &
            total_seconds, section_start, section_end &
            )
            ! inputs
            integer(int64), intent(in) :: section_start
            integer(int64), intent(in) :: section_end

            ! read/write inputs
            real(real64), intent(inout) :: total_seconds

            total_seconds = total_seconds + real(section_end - section_start, real64) / real(clock_rate, real64)
        end subroutine add_elapsed_seconds


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

            if (neighbor_image == int(this_image(), int32)) then
                return
            end if

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


        subroutine timed_sync_neighbor_images( &
            neighbor_images, n_neighbor_images &
            )
            ! inputs
            integer(int32), intent(in) :: neighbor_images(:)
            integer(int32), intent(in) :: n_neighbor_images

            if (n_neighbor_images == 0) then
                return
            end if

            call system_clock(clock_section_start)
            call sync_neighbor_images(neighbor_images, n_neighbor_images)
            call system_clock(clock_section_end)
            call add_elapsed_seconds(exchange_timing%halo_sync_seconds, clock_section_start, clock_section_end)
        end subroutine timed_sync_neighbor_images
    end subroutine exchange_halos


    subroutine exchange_poiseuille_macro_halos( &
        domain_info, halo_buffers, exchange_plan, clock_rate, exchange_timing &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info
        type(exchange_plan_t), intent(in) :: exchange_plan
        integer(int64), intent(in) :: clock_rate

        ! read/write inputs
        type(halo_buffers_t), intent(inout) :: halo_buffers

        ! output
        type(exchange_timing_t), intent(out) :: exchange_timing

        ! locals
        integer(int64) :: clock_section_start
        integer(int64) :: clock_section_end
        integer(int32) :: n_x_neighbor_images
        integer(int32) :: x_neighbor_images(2)

        exchange_timing%halo_sync_seconds = 0.0_real64
        exchange_timing%halo_transfer_seconds = 0.0_real64

        n_x_neighbor_images = 0

        if (exchange_plan%macro_left) then
            call add_neighbor_image(x_neighbor_images, n_x_neighbor_images, domain_info%left_image_id)
        end if
        if (exchange_plan%macro_right) then
            call add_neighbor_image(x_neighbor_images, n_x_neighbor_images, domain_info%right_image_id)
        end if

        call timed_sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

        ! pressure-periodic macros for poiseuille flow
        call system_clock(clock_section_start)

        if (exchange_plan%macro_left) then
            if (domain_info%left_image_id == domain_info%image_id) then ! read from own image when I_X=1
                halo_buffers%recv_macro_left(:, :) = halo_buffers%window(:, :, BUF_MACRO_RIGHT)
            else
                halo_buffers%recv_macro_left(:, :) = halo_buffers%window(:, :, BUF_MACRO_RIGHT)[domain_info%left_image_id]
            end if
        end if

        if (exchange_plan%macro_right) then
            if (domain_info%right_image_id == domain_info%image_id) then ! read from own image when I_X=1
                halo_buffers%recv_macro_right(:, :) = halo_buffers%window(:, :, BUF_MACRO_LEFT)
            else
                halo_buffers%recv_macro_right(:, :) = halo_buffers%window(:, :, BUF_MACRO_LEFT)[domain_info%right_image_id]
            end if
        end if

        call system_clock(clock_section_end)
        call add_elapsed_seconds(exchange_timing%halo_transfer_seconds, clock_section_start, clock_section_end)

        call timed_sync_neighbor_images(x_neighbor_images, n_x_neighbor_images)

    contains

        subroutine add_elapsed_seconds( &
            total_seconds, section_start, section_end &
            )
            ! inputs
            integer(int64), intent(in) :: section_start
            integer(int64), intent(in) :: section_end

            ! read/write inputs
            real(real64), intent(inout) :: total_seconds

            total_seconds = total_seconds + real(section_end - section_start, real64) / real(clock_rate, real64)
        end subroutine add_elapsed_seconds


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

            if (neighbor_image == int(this_image(), int32)) then
                return
            end if

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


        subroutine timed_sync_neighbor_images( &
            neighbor_images, n_neighbor_images &
            )
            ! inputs
            integer(int32), intent(in) :: neighbor_images(:)
            integer(int32), intent(in) :: n_neighbor_images

            if (n_neighbor_images == 0) then
                return
            end if

            call system_clock(clock_section_start)
            call sync_neighbor_images(neighbor_images, n_neighbor_images)
            call system_clock(clock_section_end)
            call add_elapsed_seconds(exchange_timing%halo_sync_seconds, clock_section_start, clock_section_end)
        end subroutine timed_sync_neighbor_images
    end subroutine exchange_poiseuille_macro_halos


end module exchange
