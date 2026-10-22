program test_domain
    ! imports
    use domain, only: domain_t, initialize_domain
    use settings, only: N_X, N_Y
    use assertions, only: check
    implicit none

    ! domain info
    type(domain_t) :: domain_info

    call initialize_domain(domain_info)
    call check(domain_info%n_images == 1, "serial domain image count")
    call check(domain_info%n_images_x == 1 .and. domain_info%n_images_y == 1, "serial domain grid")
    call check(domain_info%n_x == N_X .and. domain_info%n_y == N_Y, "serial domain size")
end program test_domain
