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
        integer(int32) :: n_x
        integer(int32) :: n_y

        ! ---------
        ! | 7 3 6 |
        ! | 4 1 2 |
        ! | 8 5 9 |
        ! ---------
        ! neighbor info
        integer(int32) :: left_image_id ! 4
        integer(int32) :: right_image_id ! 2
        integer(int32) :: bottom_image_id ! 5
        integer(int32) :: top_image_id ! 3
        integer(int32) :: top_left_image_id ! 7
        integer(int32) :: top_right_image_id ! 6
        integer(int32) :: bottom_left_image_id ! 8
        integer(int32) :: bottom_right_image_id ! 9

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
        integer(int32) :: n_images_x_env
        integer(int32) :: n_images_y_env
        integer(int32) :: bottom_y
        integer(int32) :: top_y
        integer(int32) :: left_x
        integer(int32) :: right_x
        logical :: has_n_images_x_env
        logical :: has_n_images_y_env

        call validate_coarray_image_count_source()

        domain_info%n_images = int(num_images(), int32)
        domain_info%image_id = int(this_image(), int32)

        call read_image_grid_override( &
            "I_X", n_images_x_env, has_n_images_x_env)
        call read_image_grid_override( &
            "I_Y", n_images_y_env, has_n_images_y_env)

        ! derive image grid
        if (.not. has_n_images_x_env .or. .not. has_n_images_y_env) then
            error stop "error: I_X and I_Y must be set"
        end if

        call validate_image_grid(domain_info%n_images, n_images_x_env, n_images_y_env)
        domain_info%n_images_x = n_images_x_env
        domain_info%n_images_y = n_images_y_env

        domain_info%image_x = modulo(domain_info%image_id - 1, domain_info%n_images_x) + 1
        domain_info%image_y = (domain_info%image_id - 1) / domain_info%n_images_x + 1

        domain_info%n_x = N_X / domain_info%n_images_x
        domain_info%n_y = N_Y / domain_info%n_images_y

        domain_info%x_global_start = (domain_info%image_x - 1) * domain_info%n_x + 1
        domain_info%x_global_end = domain_info%image_x * domain_info%n_x
        domain_info%y_global_start = (domain_info%image_y - 1) * domain_info%n_y + 1
        domain_info%y_global_end = domain_info%image_y * domain_info%n_y

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
            real((domain_info%n_x + 2) * (domain_info%n_y + 2) - &
            domain_info%n_x * domain_info%n_y, real64) / &
            real(domain_info%n_x * domain_info%n_y, real64)

        print '(A)', ""
        print '(A)', "--- [ domain decomposition ] ----------------------------------------------"
        print '(A,T27,I0)',       "coarray images:", domain_info%n_images
        print '(A,T27,A,I0,A,I0,A)', "image grid [X/Y]:", "[ ", &
            domain_info%n_images_x, " / ", domain_info%n_images_y, " ]"
        print '(A,T27,A,I0,A,I0,A)', "local domain [X/Y]:", "[ ", &
            domain_info%n_x, " / ", domain_info%n_y, " ]"
        print '(A,T24,F8.3,A)',   "halo cells:", halo_cell_percent, " %"
    end subroutine print_domain_summary


    subroutine validate_coarray_image_count_source()
        ! locals
        integer :: n_images_env_length
        integer :: n_images_env_status
        integer :: config_env_length
        integer :: config_env_status
        logical :: has_n_images_env
        logical :: has_config_env

        call get_environment_variable( &
            "FOR_COARRAY_NUM_IMAGES", length=n_images_env_length, status=n_images_env_status)
        call get_environment_variable( &
            "FOR_COARRAY_CONFIG_FILE", length=config_env_length, status=config_env_status)

        if (n_images_env_status < 0 .or. config_env_status < 0) then
            error stop "error: invalid coarray image count configuration"
        end if

        has_n_images_env = n_images_env_status == 0 .and. n_images_env_length > 0
        has_config_env = config_env_status == 0 .and. config_env_length > 0

        if (.not. has_n_images_env .and. .not. has_config_env) then
            error stop "error: FOR_COARRAY_NUM_IMAGES or FOR_COARRAY_CONFIG_FILE must be set"
        end if
    end subroutine validate_coarray_image_count_source


    subroutine validate_image_grid( &
        n_images, n_images_x, n_images_y &
        )
        ! inputs
        integer(int32), intent(in) :: n_images
        integer(int32), intent(in) :: n_images_x
        integer(int32), intent(in) :: n_images_y

        if (n_images_x <= 0 .or. n_images_y <= 0) then
            error stop "error: I_X and I_Y must be positive"
        else if (n_images_x * n_images_y /= n_images) then
            error stop "error: I_X * I_Y must match coarray images"
        else if (mod(N_X, n_images_x) /= 0) then
            error stop "error: N_X must be divisible by I_X"
        else if (mod(N_Y, n_images_y) /= 0) then
            error stop "error: N_Y must be divisible by I_Y"
        end if
    end subroutine validate_image_grid


    subroutine read_image_grid_override( &
        env_name, n_images_override, has_override &
        )
        ! inputs
        character(len=*), intent(in) :: env_name

        ! output
        integer(int32), intent(out) :: n_images_override
        logical, intent(out) :: has_override

        ! locals
        character(len=32) :: env_value
        integer :: env_status
        integer :: read_status

        call get_environment_variable(env_name, env_value, status=env_status)

        n_images_override = 0

        if (env_status > 0) then
            has_override = .false.
            return
        else if (env_status < 0) then
            error stop "error: invalid image grid override"
        end if

        has_override = .true.

        read(env_value, *, iostat=read_status) n_images_override
        if (read_status /= 0) then
            error stop "error: invalid image grid override"
        end if
    end subroutine read_image_grid_override


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
