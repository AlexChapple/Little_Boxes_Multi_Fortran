! Little Boxes Multi with 4th order Runge-Kutta method adapted into Fortran code
!
!
!
!

program main

    implicit none

    complex :: z 
    real :: c

    z = cmplx(1,2)

    c = modulo_func(z)

    print *, c 

    contains

    function modulo_func(z) result(c)

        implicit none

        ! Takes in a complex number z and returns its modulus

        ! Declare var types
        real :: a, b, c 
        complex , intent(in) :: z

        a = real(z)
        b = aimag(z)

        c = a**2 + b**2

    end function 

end program main 



