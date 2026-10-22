program test_fields
    ! imports
    use iso_fortran_env, only: int32, int64, real64, file_storage_size
    use ieee_arithmetic, only: ieee_is_finite
    use settings, only: FP, FP_DTYPE, N_X, N_Y, N_STEPS, U_MAX, N_SIN, OMEGA, &
        SIM_SHEAR_WAVE, SIM_SLIDING_LID, sim_mode_to_string
    use assertions, only: check
    implicit none

    ! exported fields
    integer(int32), parameter :: FIELD_DENSITY = 1
    integer(int32), parameter :: FIELD_VELOCITY_X = 2
    integer(int32), parameter :: FIELD_VELOCITY_Y = 3
    integer(int32), parameter :: FIELD_VELOCITY_MAG = 4
    integer(int32), parameter :: N_FIELDS = 4
    character(len=12), parameter :: FIELD_NAMES(N_FIELDS) = &
        [ character(len=12) :: "density", "velocity_x", "velocity_y", "velocity_mag" ]

    ! run configuration
    integer(int32) :: sim_mode
    integer(int32) :: n_images_x
    integer(int32) :: n_images_y
    character(len=:), allocatable :: reference_dir
    integer :: argument_length
    integer :: argument_status

    ! field data
    real(FP) :: fields(N_X, N_Y, N_FIELDS)
    real(FP) :: initial_fields(N_X, N_Y, N_FIELDS)
    real(FP) :: reference_fields(N_X, N_Y, N_FIELDS)

    ! comparison tolerances
    real(real64) :: relative_tolerance
    real(real64) :: absolute_tolerance
    real(real64) :: mass_tolerance

    ! temp
    integer(int32) :: endpoint, step, field_id

    call check(command_argument_count() == 3 .or. command_argument_count() == 4, &
        "usage: test_fields mode images_x images_y [reference_dir]")
    sim_mode = read_integer_argument(1)
    n_images_x = read_integer_argument(2)
    n_images_y = read_integer_argument(3)
    call check(sim_mode >= SIM_SHEAR_WAVE .and. sim_mode <= SIM_SLIDING_LID, "invalid flow mode")

    reference_dir = ""
    if (command_argument_count() == 4) then
        call get_command_argument(4, length=argument_length, status=argument_status)
        call check(argument_status == 0 .and. argument_length > 0, "invalid reference directory")
        deallocate(reference_dir)
        allocate(character(len=argument_length) :: reference_dir)
        call get_command_argument(4, reference_dir, status=argument_status)
        call check(argument_status == 0, "could not read reference directory")
    end if

    if (FP_DTYPE == "real32") then
        relative_tolerance = 2.0e-5_real64
        absolute_tolerance = 5.0e-7_real64
        mass_tolerance = 2.0e-5_real64
    else
        relative_tolerance = 2.0e-11_real64
        absolute_tolerance = 2.0e-12_real64
        mass_tolerance = 2.0e-12_real64
    end if

    ! verify metadata against the compiled test settings and requested decomposition
    call check_metadata(".", n_images_x, n_images_y)
    if (len(reference_dir) > 0) then
        call check_metadata(reference_dir, 1_int32, 1_int32)
    end if

    ! validate complete initial and final exports
    do endpoint = 0, 1
        step = endpoint * N_STEPS
        call read_fields(".", step, fields)
        call check(all(fields(:, :, FIELD_DENSITY) > 0.0_FP), "nonpositive exported density")
        call compare_fields( &
            "velocity magnitude", fields(:, :, FIELD_VELOCITY_MAG), &
            sqrt(fields(:, :, FIELD_VELOCITY_X)**2 + fields(:, :, FIELD_VELOCITY_Y)**2))

        if (endpoint == 0) then
            initial_fields = fields
        end if

        if (len(reference_dir) > 0) then
            call read_fields(reference_dir, step, reference_fields)
            do field_id = 1, N_FIELDS
                call compare_fields(trim(FIELD_NAMES(field_id)), fields(:, :, field_id), reference_fields(:, :, field_id))
            end do
        end if
    end do

    if (sim_mode == SIM_SHEAR_WAVE) then
        call check_shear_wave()
    end if
    write(*, '(A)') "passed: " // sim_mode_to_string(sim_mode) // " exported field checks"

contains

    function read_integer_argument( &
        position &
        ) result(value)
        ! inputs
        integer, intent(in) :: position

        ! output
        integer(int32) :: value

        ! temp
        character(len=32) :: argument
        integer :: io_stat

        call get_command_argument(position, argument, status=io_stat)
        call check(io_stat == 0, "could not read integer argument")
        read(argument, *, iostat=io_stat) value
        call check(io_stat == 0, "invalid integer argument")
    end function read_integer_argument


    subroutine read_fields( &
        run_dir, step, field_data &
        )
        ! inputs
        character(len=*), intent(in) :: run_dir
        integer(int32), intent(in) :: step

        ! write destinations
        real(FP), intent(out) :: field_data(N_X, N_Y, N_FIELDS)

        ! temp
        character(len=:), allocatable :: file_path
        character(len=9) :: suffix
        integer(int32) :: field_id
        integer(int64) :: file_size
        integer(int64) :: expected_size
        integer :: unit
        integer :: io_stat
        logical :: exists

        expected_size = int(N_X, int64) * int(N_Y, int64) * storage_size(0.0_FP) / file_storage_size
        write(suffix, '(I9.9)') step

        do field_id = 1, N_FIELDS
            file_path = run_dir // "/output/test/" // trim(FIELD_NAMES(field_id)) // "_" // suffix // ".bin"
            inquire(file=file_path, exist=exists, size=file_size, iostat=io_stat)
            call check(io_stat == 0 .and. exists, "missing field: " // file_path)
            call check(file_size == expected_size, "incorrect field size: " // file_path)

            open(newunit=unit, file=file_path, access="stream", form="unformatted", &
                status="old", action="read", iostat=io_stat)
            call check(io_stat == 0, "could not open field: " // file_path)
            read(unit, iostat=io_stat) field_data(:, :, field_id)
            close(unit)
            call check(io_stat == 0, "could not read field: " // file_path)
            call check(all(ieee_is_finite(field_data(:, :, field_id))), "nonfinite field: " // file_path)
        end do
    end subroutine read_fields


    subroutine compare_fields( &
        name, actual, expected &
        )
        ! inputs
        character(len=*), intent(in) :: name
        real(FP), intent(in) :: actual(N_X, N_Y)
        real(FP), intent(in) :: expected(N_X, N_Y)

        ! temp
        real(real64) :: difference(N_X, N_Y)
        real(real64) :: allowed(N_X, N_Y)
        integer :: location(2)

        call check(all(ieee_is_finite(actual)) .and. all(ieee_is_finite(expected)), "nonfinite comparison: " // name)
        difference = abs(real(actual, real64) - real(expected, real64))
        allowed = absolute_tolerance + relative_tolerance * abs(real(expected, real64))

        if (any(difference > allowed)) then
            location = maxloc(difference - allowed)
            write(*, *) "mismatch in ", name, " at x/y: ", location
            write(*, *) "actual, expected: ", actual(location(1), location(2)), expected(location(1), location(2))
            write(*, *) "error, allowed: ", difference(location(1), location(2)), allowed(location(1), location(2))
            call check(.false., "field mismatch: " // name)
        end if
    end subroutine compare_fields


    function read_metadata_value( &
        run_dir, key &
        ) result(value)
        ! inputs
        character(len=*), intent(in) :: run_dir
        character(len=*), intent(in) :: key

        ! output
        character(len=:), allocatable :: value

        ! temp
        character(len=1024) :: line
        integer :: unit
        integer :: io_stat
        integer :: colon
        integer :: value_length

        open(newunit=unit, file=run_dir // "/output/test/config.json", &
            status="old", action="read", iostat=io_stat)
        call check(io_stat == 0, "could not open run metadata: " // run_dir)

        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) then
                exit
            end if
            colon = index(line, ":")
            if (colon == 0) then
                cycle
            end if
            if (trim(adjustl(line(:colon-1))) /= '"' // key // '"') then
                cycle
            end if

            value = trim(adjustl(line(colon+1:)))
            value_length = len(value)
            if (value_length > 0) then
                if (value(value_length:value_length) == ",") then
                    value = trim(value(:value_length-1))
                end if
            end if
            value_length = len(value)
            if (value_length >= 2) then
                if (value(1:1) == '"' .and. value(value_length:value_length) == '"') then
                    value = value(2:value_length-1)
                end if
            end if
            close(unit)
            return
        end do

        close(unit)
        call check(.false., "missing or unreadable metadata key: " // key)
    end function read_metadata_value


    subroutine check_metadata( &
        run_dir, n_images_x, n_images_y &
        )
        ! inputs
        character(len=*), intent(in) :: run_dir
        integer(int32), intent(in) :: n_images_x
        integer(int32), intent(in) :: n_images_y

        ! expected metadata
        character(len=14), parameter :: KEYS(6) = [ character(len=14) :: &
            "N_X", "N_Y", "N_STEPS", "coarray_images", "image_grid_x", "image_grid_y" ]

        ! temp
        integer(int32) :: expected(6)
        integer(int32) :: actual
        integer :: key_id
        integer :: io_stat
        character(len=:), allocatable :: value

        call check(read_metadata_value(run_dir, "SIM_MODE") == sim_mode_to_string(sim_mode), "incorrect flow metadata")
        call check(read_metadata_value(run_dir, "file_dtype") == FP_DTYPE, "incorrect precision metadata")
        expected = [ N_X, N_Y, N_STEPS, n_images_x * n_images_y, n_images_x, n_images_y ]
        do key_id = 1, size(KEYS)
            value = read_metadata_value(run_dir, trim(KEYS(key_id)))
            read(value, *, iostat=io_stat) actual
            call check(io_stat == 0, "invalid integer metadata: " // trim(KEYS(key_id)))
            call check(actual == expected(key_id), "incorrect metadata: " // trim(KEYS(key_id)))
        end do
    end subroutine check_metadata


    subroutine check_shear_wave()
        ! temp
        real(real64) :: wave_number
        real(real64) :: viscosity
        real(real64) :: mode(N_Y)
        real(real64) :: profile(N_Y)
        real(real64) :: measured_amplitude
        real(real64) :: analytical_amplitude
        real(real64) :: relative_error
        real(real64) :: mass_error
        real(FP) :: expected_field(N_X, N_Y)
        integer(int32) :: y

        wave_number = 2.0_real64 * acos(-1.0_real64) * real(N_SIN, real64) / real(N_Y, real64)
        do y = 1, N_Y
            mode(y) = sin(wave_number * real(y - 1, real64))
            expected_field(:, y) = real(real(U_MAX, real64) * mode(y), FP)
        end do
        call compare_fields("initial shear profile", initial_fields(:, :, FIELD_VELOCITY_X), expected_field)

        ! project the measured velocity onto the prescribed sinusoidal mode
        profile = sum(real(fields(:, :, FIELD_VELOCITY_X), real64), dim=1) / real(N_X, real64)
        measured_amplitude = sum(profile * mode) / sum(mode * mode)
        viscosity = (1.0_real64 / real(OMEGA, real64) - 0.5_real64) / 3.0_real64
        analytical_amplitude = real(U_MAX, real64) * exp(-viscosity * wave_number**2 * real(N_STEPS, real64))
        relative_error = abs(measured_amplitude / analytical_amplitude - 1.0_real64)
        write(*, *) "shear wave relative amplitude error: ", relative_error
        call check(relative_error < 0.01_real64, "incorrect analytical shear wave decay")

        ! only the periodic case is required to conserve total mass
        mass_error = abs(sum(real(fields(:, :, FIELD_DENSITY), real64)) / &
            sum(real(initial_fields(:, :, FIELD_DENSITY), real64)) - 1.0_real64)
        write(*, *) "shear wave relative mass error: ", mass_error
        call check(mass_error < mass_tolerance, "shear wave mass conservation")
        call check(all(abs(real(fields(:, :, FIELD_VELOCITY_Y), real64)) <= absolute_tolerance), &
            "nonzero transverse shear velocity")
    end subroutine check_shear_wave

end program test_fields
