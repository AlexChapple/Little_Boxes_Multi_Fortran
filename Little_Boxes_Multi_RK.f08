! Little Boxes Multi with 4th order Runge-Kutta method adapted into Fortran code
!
!
!
!

program main

    implicit none

    integer, parameter :: N = 20
    






    ! Defines all the functions that are used in this program below
    contains

    function return_total (N, g_0, u_0, g_1, u_1, g_2, u_2) result(total)

        ! Takes in all the coefficient arrays and returns the total magnitude

        implicit none
        
        real :: total
        integer :: N, j, k
        complex :: g_0, u_0
        complex, dimension(N) :: g_1, u_1
        complex, dimension(N,N) :: g_2, u_2

        total = modulo_func(g_0) + modulo_func(u_0)

        do j = 1,N
            total = total + modulo_func(g_1(j)) + modulo_func(u_1(j))
        end do 

        do j = 1,N 
            do k = 1,N
                total = total + modulo_func(g_2(j,k)) + modulo_func(u_2(j,k))
            end do 
        end do

    end function 

    function modulo_func(z) result(c)

        ! Takes in a complex number z and returns its modulus

        implicit none

        ! Declare var types
        real :: a, b, c 
        complex , intent(in) :: z

        a = real(z)
        b = aimag(z)

        c = a**2 + b**2

    end function 


end program main 



