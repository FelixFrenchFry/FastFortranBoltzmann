module domain
    ! imports
    use iso_fortran_env, only: int32, real64
    use settings, only: N_X, N_Y
    implicit none
    private

    public :: domain_t, initialize_domain, print_domain_summary

    type :: domain_t

        ! global domain info
        integer(int32) :: n_images
        integer(int32) :: n_images_x
        integer(int32) :: n_images_y
        integer(int32) :: x_global_start
        integer(int32) :: x_global_end
        integer(int32) :: y_global_start
        integer(int32) :: y_global_end

        ! local domain coordinates and size
        integer(int32) :: image_id
        integer(int32) :: image_x
        integer(int32) :: image_y
        integer(int32) :: n_x_local
        integer(int32) :: n_y_local

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! neighbor info
        integer(int32) :: bottom_image_id ! 5
        integer(int32) :: top_image_id ! 3
        integer(int32) :: left_image_id ! 4
        integer(int32) :: right_image_id ! 2
        integer(int32) :: bottom_left_image_id ! 8
        integer(int32) :: bottom_right_image_id ! 9
        integer(int32) :: top_left_image_id ! 7
        integer(int32) :: top_right_image_id ! 6

        logical :: at_left_boundary
        logical :: at_right_boundary
        logical :: at_bottom_boundary
        logical :: at_top_boundary
        
    end type domain_t

contains

    subroutine initialize_domain( &
        domain_info &
        )
        ! output
        type(domain_t), intent(out) :: domain_info

        ! locals
        integer(int32) :: n_images_sqrt
        integer(int32) :: bottom_y
        integer(int32) :: top_y
        integer(int32) :: left_x
        integer(int32) :: right_x

        domain_info%n_images = int(num_images(), int32)
        domain_info%image_id = int(this_image(), int32)

        n_images_sqrt = int(sqrt(real(domain_info%n_images, real64)), int32)

        ! check image count and domain size constraints
        if (n_images_sqrt * n_images_sqrt /= domain_info%n_images) then
            error stop "error: number of coarray images must be a square number"
        else if (mod(N_X, n_images_sqrt) /= 0) then
            error stop "error: N_X must be divisible by sqrt(num_images)"
        else if (mod(N_Y, n_images_sqrt) /= 0) then
            error stop "error: N_Y must be divisible by sqrt(num_images)"
        end if

        ! derive all domain decomposition infos
        domain_info%n_images_x = n_images_sqrt
        domain_info%n_images_y = n_images_sqrt

        domain_info%image_x = modulo(domain_info%image_id - 1, domain_info%n_images_x) + 1
        domain_info%image_y = (domain_info%image_id - 1) / domain_info%n_images_x + 1

        domain_info%n_x_local = N_X / domain_info%n_images_x
        domain_info%n_y_local = N_Y / domain_info%n_images_y

        domain_info%x_global_start = (domain_info%image_x - 1) * domain_info%n_x_local + 1
        domain_info%x_global_end = domain_info%image_x * domain_info%n_x_local
        domain_info%y_global_start = (domain_info%image_y - 1) * domain_info%n_y_local + 1
        domain_info%y_global_end = domain_info%image_y * domain_info%n_y_local

        domain_info%at_left_boundary = domain_info%image_x == 1
        domain_info%at_right_boundary = domain_info%image_x == domain_info%n_images_x
        domain_info%at_bottom_boundary = domain_info%image_y == 1
        domain_info%at_top_boundary = domain_info%image_y == domain_info%n_images_y

        left_x = wrap_image_coordinate(domain_info%image_x - 1, domain_info%n_images_x)
        right_x = wrap_image_coordinate(domain_info%image_x + 1, domain_info%n_images_x)
        bottom_y = wrap_image_coordinate(domain_info%image_y - 1, domain_info%n_images_y)
        top_y = wrap_image_coordinate(domain_info%image_y + 1, domain_info%n_images_y)

        domain_info%left_image_id = image_id_from_coordinates(left_x, domain_info%image_y, domain_info%n_images_x)
        domain_info%right_image_id = image_id_from_coordinates(right_x, domain_info%image_y, domain_info%n_images_x)
        domain_info%bottom_image_id = image_id_from_coordinates(domain_info%image_x, bottom_y, domain_info%n_images_x)
        domain_info%top_image_id = image_id_from_coordinates(domain_info%image_x, top_y, domain_info%n_images_x)

        domain_info%bottom_left_image_id = image_id_from_coordinates(left_x, bottom_y, domain_info%n_images_x)
        domain_info%bottom_right_image_id = image_id_from_coordinates(right_x, bottom_y, domain_info%n_images_x)
        domain_info%top_left_image_id = image_id_from_coordinates(left_x, top_y, domain_info%n_images_x)
        domain_info%top_right_image_id = image_id_from_coordinates(right_x, top_y, domain_info%n_images_x)
    end subroutine initialize_domain


    subroutine print_domain_summary( &
        domain_info &
        )
        ! inputs
        type(domain_t), intent(in) :: domain_info

        ! locals
        real(real64) :: halo_cell_percent

        if (this_image() /= 1) then
            return
        end if

        halo_cell_percent = 100.0_real64 * &
            real((domain_info%n_x_local + 2) * (domain_info%n_y_local + 2) - &
            domain_info%n_x_local * domain_info%n_y_local, real64) / &
            real(domain_info%n_x_local * domain_info%n_y_local, real64)

        print '(A)', ""
        print '(A)', "--- [ domain decomposition ] ----------------------------------------------"
        print '(A,T27,I0)',       "coarray images:", domain_info%n_images
        print '(A,T27,A,I0,A,I0,A)', "image grid [X/Y]:", "[ ", &
            domain_info%n_images_x, " / ", domain_info%n_images_y, " ]"
        print '(A,T27,A,I0,A,I0,A)', "local domain [X/Y]:", "[ ", &
            domain_info%n_x_local, " / ", domain_info%n_y_local, " ]"
        print '(A,T24,F8.3,A)',   "halo cells:", halo_cell_percent, " %"
    end subroutine print_domain_summary


    pure function wrap_image_coordinate( &
        image_coordinate, n_images_in_dir &
        ) result(wrapped_coordinate)
        ! inputs
        integer(int32), intent(in) :: image_coordinate
        integer(int32), intent(in) :: n_images_in_dir

        ! output
        integer(int32) :: wrapped_coordinate

        wrapped_coordinate = modulo(image_coordinate - 1, n_images_in_dir) + 1
    end function wrap_image_coordinate


    pure function image_id_from_coordinates( &
        image_x, image_y, n_images_x &
        ) result(image_id)
        ! inputs
        integer(int32), intent(in) :: image_x
        integer(int32), intent(in) :: image_y
        integer(int32), intent(in) :: n_images_x

        ! output
        integer(int32) :: image_id

        image_id = image_x + (image_y - 1) * n_images_x
    end function image_id_from_coordinates

end module domain
