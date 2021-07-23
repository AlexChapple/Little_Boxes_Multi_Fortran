
program test

    implicit none

    ! Type declaration 
    complex, dimension(2) :: a, b, c

    a(1) = cmplx(1,0)
    a(2) = cmplx(0,1)

    b(1) = cmplx(1,0)
    b(2) = cmplx(0,1)

    c = a + b 

    print *, c    

end program test