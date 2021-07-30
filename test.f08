
program test

    implicit none

    ! Type declaration 
    integer :: i, j 
    real :: a, b, c

    do i = 1,100

        print *, "start:", a 

        ! a = 0
        
        do j = 1, 100

            b = real(1)
            a = a + b 

        end do 

        print *, "end:" ,a

    end do 

end program test