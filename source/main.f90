program main
    implicit none

    integer :: x[*]

    x = this_image()
    sync all

    write (*,*) 'Hello world, this is image ', this_image(), ' of ', num_images()

end program main
