
program test

    implicit none

    real (kind = 4) :: a 
    real (kind = 8) :: b 

    a = 0.0000000001
    b = 0.000000000000001

    print *, a, b

end program test