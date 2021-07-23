
program test

    implicit none

    ! Type declaration 
    real, dimension(4) :: a,b 
    integer :: i 

    do i = 1,4 
        a(i) = i 
        b(i) = i ** 2
    end do 

    open(1, file="data.txt", status="new")
    do i = 1,4
        write(1,*) a(i), b(i)
    end do 


end program test