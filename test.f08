
program test

    implicit none

    ! Type declaration 
    real, dimension(4) :: a 
    integer :: i 

    do i = 1,4 
        write(*, '(1x,i0)', advance="no"), i
    end do 

    

end program test