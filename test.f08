
program test

    implicit none

    ! Type declaration 
    integer :: beginning, end, i
    integer, dimension(10000) :: a 

    call system_clock(beginning)

    do i = 1,10000
        a(i) = exp(real(i))
    end do 

    call system_clock(end)

    print *, real(end - beginning)


end program test