! Little Boxes Multi with 4th order Runge-Kutta method adapted into Fortran code
! This calculates the waiting time distribution. 

program main

    ! ----------------------------------------------------------------------------------
    ! 
    ! Main program that simulates the open quantum system. 
    ! Also calculates the g2 distribution  
    ! and updates the waiting time list accordingly. 
    !
    ! ----------------------------------------------------------------------------------

    implicit none

    ! Declare variables and parameters
    integer, parameter :: N = 20
    integer, parameter :: time_steps = 10000
    integer, parameter :: end_time = 8
    integer, parameter :: num_of_simulations = 1500 
    integer, parameter :: bin_width = 10 * time_steps ! Creates 800 bins for waiting time distribution 
    real, parameter :: pi = 3.1415927
    real, parameter :: phase = 0  
    real, parameter :: gammaL = 0.5 
    real, parameter :: gammaR = 0.5 
    complex, parameter :: Omega = cmplx(10 * pi, 0)
    real, parameter :: dt = real(end_time) / real(time_steps) 
    real, parameter :: tau = 10.0 * dt * real(N) 
    real :: total, waiting_total, current_time ! Last time the photon was found, total for normalisation purposes 
    integer :: sim, index, j, k, beginning, end, rate
    real, dimension(time_steps) :: time_list, rand_list
    real, dimension(bin_width) :: waiting_time_list, reduced_time_list 

    ! The coefficients (g for ground, u for up)
    complex :: g_0, g_0_k1, g_0_k2, g_0_k3, g_0_k4, g_0_new, e_0 ,e_0_k1, e_0_k2, e_0_k3, e_0_k4, e_0_new
    complex, dimension(N) :: g_1, g_1_k1, g_1_k2, g_1_k3, g_1_k4, g_1_new, e_1 ,e_1_k1, e_1_k2, e_1_k3, e_1_k4, e_1_new 
    complex, dimension(N,N) :: g_2, g_2_k1, g_2_k2, g_2_k3, g_2_k4, g_2_new, e_2 ,e_2_k1, e_2_k2, e_2_k3, e_2_k4, e_2_new 
    complex :: lambdaL, lambdaR
    real :: psi_0, psi_1, prob, rand_num 
    
    ! Program execution time tracking 
    call system_clock(beginning, rate)

    ! Construct time_list
    call linspace(start=0.0, end=end_time, time_list=time_list) ! This makes the time list 
    call linspace(start=0.0, end=end_time, time_list=waiting_time_list) ! Creates the waiting time list 

    ! Initialise lambdaL and lambdaR
    lambdaL = exp(cmplx(0, phase / 2)) * sqrt(gammaL) * sqrt(N/tau)
    lambdaR = exp(cmplx(0, -phase / 2)) * sqrt(gammaR) * sqrt(N/tau)

    ! A do loop will go through and do the simulations here
    do sim = 1, num_of_simulations

        ! Initialise arrays
        call initialise_arrays(N, g_0, g_0_k1, g_0_k2, g_0_k3, g_0_k4, g_0_new, e_0 , & 
                            e_0_k1, e_0_k2, e_0_k3, e_0_k4, e_0_new, g_1, g_1_k1, & 
                            g_1_k2, g_1_k3, g_1_k4, g_1_new, e_1 ,e_1_k1, e_1_k2, & 
                            e_1_k3, e_1_k4, e_1_new, g_2, g_2_k1, g_2_k2, g_2_k3, &
                            g_2_k4, g_2_new, e_2 ,e_2_k1, e_2_k2, e_2_k3, e_2_k4, e_2_new)

        ! Construct random number list 
        call random_number(rand_list)

        do index = 1, size(time_list)

            psi_0 = 0; psi_1 = 0; prob = 0; rand_num = 0; total = 0

            ! Calculates k1 values 
            g_0_k1 = (-Omega/2) * e_0
            e_0_k1 = (lambdaL * g_1(1)) + (lambdaR * g_1(N)) + ((Omega/2) * g_0) 

            g_1_k1(1) = (lambdaL * e_0) - (Omega/2)*e_1(1) 
            g_1_k1(N) = (lambdaR * e_0) - (Omega/2)*e_1(N) 

            e_1_k1(1) = (lambdaR * g_2(1,1)) + (Omega/2)*g_1(1)
            e_1_k1(N) = (lambdaL * g_2(1,N)) + (Omega/2)*g_1(N)

            do j = 2, (N-1)
                g_1_k1(j) = (-Omega/2) * e_1(j)
                e_1_k1(j) = (lambdaL * g_2(j,1)) + (lambdaR * g_2(j,N)) + ((Omega/2)*g_1(j))
            end do 

            do j = 2,N 
                g_2_k1(j,1) = (lambdaL * e_1(j)) - (Omega/2)*e_2(j,1)
            end do 

            do j = 1,(N-1)
                g_2_k1(j,N) = (lambdaR * e_1(j)) - (Omega/2)*e_2(j,N)
            end do 

            do j = 1,(N-2)
                do k = (j+1),(N-1) 
                    g_2_k1(j,k) = (-Omega/2) * e_2(j,k)
                end do 
            end do

            do j = 1,(N-1)
                do k = (j+1), N 
                    e_2_k1(j,k) = (Omega/2) * g_2(j,k)
                end do 
            end do 

            ! Calculates k2 values
            g_0_k2 = (-Omega/2) * (e_0 + (dt * g_0_k1/2))
            e_0_k2 = (lambdaL * (g_1(1) + (dt * g_1_k1(1)/2))) + & 
            (lambdaR * (g_1(N) + (dt * g_1_k1(N)/2))) + ((Omega/2) * (g_0 + (dt * g_0_k1/2)))

            g_1_k2(1) = (lambdaL * (e_0 + (dt * g_0_k1/2))) - (Omega/2)*(e_1(1) + (dt*e_1_k1(1)/2))
            g_1_k2(N) = (lambdaR * (e_0 + (dt * g_0_k1/2))) - (Omega/2)*(e_1(N) + (dt*e_1_k1(N)/2))

            e_1_k2(1) = (lambdaR * (g_2(1,1) + (dt * g_2_k1(1,1)/2))) + (Omega/2)*(g_1(1) + (dt * g_1_k1(1)/2))
            e_1_k2(N) = (lambdaL * (g_2(1,N) + (dt * g_2_k1(1,N)/2))) + (Omega/2)*(g_1(N) + (dt * g_1_k1(N)/2))

            do j = 2, (N-1)
                g_1_k2(j) = (-Omega/2) * (e_1(j) + (dt * e_1_k1(j)/2))
                e_1_k2(j) = (lambdaL * (g_2(j,1) + (dt * g_2_k1(j,1) / 2))) + &
                (lambdaR * (g_2(j,N) + (dt * g_2_k1(j,N) / 2))) + (Omega/2)*(g_1(j) + (dt * g_1_k1(j)/2))
            end do 

            do j = 2,N 
                g_2_k2(j,1) = (lambdaL * (e_1(j) + (dt * e_1_k1(j)/2))) - (Omega/2) * (e_2(j,1) + (dt * e_2_k1(j,1)/2))
            end do 

            do j = 1,(N-1)
                g_2_k2(j,N) = (lambdaR * (e_1(j) + (dt * e_1_k1(j)/2))) - (Omega/2) * (e_2(j,1) + (dt * e_2_k1(j,1)/2))
            end do 

            do j = 1,(N-2)
                do k = (j+1),(N-1) 
                    g_2_k2(j,k) = (-Omega/2) * (e_2(j,k) + (dt * e_2_k1(j,k)/2))
                end do 
            end do

            do j = 1,(N-1)
                do k = (j+1), N 
                    e_2_k2(j,k) = (Omega/2) * (g_2(j,k) + (dt * g_2_k1(j,k)/2))
                end do 
            end do 

            ! Calculates k3 values
            g_0_k3 = (-Omega/2) * (e_0 + (dt * g_0_k2/2))
            e_0_k3 = (lambdaL * (g_1(1) + (dt * g_1_k2(1)/2))) + &
            (lambdaR * (g_1(N) + (dt * g_1_k2(N)/2))) + ((Omega/2) * (g_0 + (dt * g_0_k2/2)))

            g_1_k3(1) = (lambdaL * (e_0 + (dt * g_0_k2/2))) - (Omega/2)*(e_1(1) + (dt*e_1_k2(1)/2))
            g_1_k3(N) = (lambdaR * (e_0 + (dt * g_0_k2/2))) - (Omega/2)*(e_1(N) + (dt*e_1_k2(N)/2))

            e_1_k3(1) = (lambdaR * (g_2(1,1) + (dt * g_2_k2(1,1)/2))) + (Omega/2)*(g_1(1) + (dt * g_1_k2(1)/2))
            e_1_k3(N) = (lambdaL * (g_2(1,N) + (dt * g_2_k2(1,N)/2))) + (Omega/2)*(g_1(N) + (dt * g_1_k2(N)/2))

            do j = 2, (N-1)
                g_1_k3(j) = (-Omega/2) * (e_1(j) + (dt * e_1_k2(j)/2))
                e_1_k3(j) = (lambdaL * (g_2(j,1) + (dt * g_2_k2(j,1) / 2))) + &
                (lambdaR * (g_2(j,N) + (dt * g_2_k2(j,N) / 2))) + (Omega/2)*(g_1(j) + (dt * g_1_k2(j)/2))
            end do 

            do j = 2,N 
                g_2_k3(j,1) = (lambdaL * (e_1(j) + (dt * e_1_k2(j)/2))) - (Omega/2) * (e_2(j,1) + (dt * e_2_k2(j,1)/2))
            end do 

            do j = 1,(N-1)
                g_2_k3(j,N) = (lambdaR * (e_1(j) + (dt * e_1_k2(j)/2))) - (Omega/2) * (e_2(j,1) + (dt * e_2_k2(j,1)/2))
            end do 

            do j = 1,(N-2)
                do k = (j+1),(N-1) 
                    g_2_k3(j,k) = (-Omega/2) * (e_2(j,k) + (dt * e_2_k2(j,k)/2))
                end do 
            end do

            do j = 1,(N-1)
                do k = (j+1), N 
                    e_2_k3(j,k) = (Omega/2) * (g_2(j,k) + (dt * g_2_k2(j,k)/2))
                end do 
            end do 

            ! Calculates k4 values 
            g_0_k4 = (-Omega/2) * (e_0 + (dt * g_0_k2))
            e_0_k4 = (lambdaL * (g_1(1) + (dt * g_1_k2(1)))) + &
            (lambdaR * (g_1(N) + (dt * g_1_k2(N)))) + ((Omega/2) * (g_0 + (dt * g_0_k2)))

            g_1_k4(1) = (lambdaL * (e_0 + (dt * g_0_k2))) - (Omega/2)*(e_1(1) + (dt*e_1_k2(1)))
            g_1_k4(N) = (lambdaR * (e_0 + (dt * g_0_k2))) - (Omega/2)*(e_1(N) + (dt*e_1_k2(N)))

            e_1_k4(1) = (lambdaR * (g_2(1,1) + (dt * g_2_k2(1,1)))) + (Omega/2)*(g_1(1) + (dt * g_1_k2(1)))
            e_1_k4(N) = (lambdaL * (g_2(1,N) + (dt * g_2_k2(1,N)))) + (Omega/2)*(g_1(N) + (dt * g_1_k2(N)))

            do j = 2, (N-1)
                g_1_k4(j) = (-Omega/2) * (e_1(j) + (dt * e_1_k2(j)))
                e_1_k4(j) = (lambdaL * (g_2(j,1) + (dt * g_2_k2(j,1)))) + &
                (lambdaR * (g_2(j,N) + (dt * g_2_k2(j,N)))) + (Omega/2)*(g_1(j) + (dt * g_1_k2(j)))
            end do 

            do j = 2,N 
                g_2_k4(j,1) = (lambdaL * (e_1(j) + (dt * e_1_k2(j)))) - (Omega/2) * (e_2(j,1) + (dt * e_2_k2(j,1)))
            end do 

            do j = 1,(N-1)
                g_2_k4(j,N) = (lambdaR * (e_1(j) + (dt * e_1_k2(j)))) - (Omega/2) * (e_2(j,1) + (dt * e_2_k2(j,1)))
            end do 

            do j = 1,(N-2)
                do k = (j+1),(N-1) 
                    g_2_k4(j,k) = (-Omega/2) * (e_2(j,k) + (dt * e_2_k2(j,k)))
                end do 
            end do

            do j = 1,(N-1)
                do k = (j+1), N 
                    e_2_k4(j,k) = (Omega/2) * (g_2(j,k) + (dt * g_2_k2(j,k)))
                end do 
            end do  

            ! Collect them all together
            g_0_new = g_0 + (dt/6)*(g_0_k1 + 2*g_0_k2 + 2*g_0_k3 + g_0_k4)
            e_0_new = e_0 + (dt/6)*(e_0_k1 + 2*e_0_k2 + 2*e_0_k3 + e_0_k4)

            g_1_new = g_1 + (dt/6)*(g_1_k1 + 2*g_1_k2 + 2*g_1_k3 + g_1_k4)
            e_1_new = e_1 + (dt/6)*(e_1_k1 + 2*e_1_k2 + 2*e_1_k3 + e_1_k4)

            g_2_new = g_2 + (dt/6)*(g_2_k1 + 2*g_2_k2 + 2*g_2_k3 + g_2_k4)
            e_2_new = e_2 + (dt/6)*(e_2_k1 + 2*e_2_k2 + 2*e_2_k3 + e_2_k4)

            ! Normalise everything here
            call normalise_new(total, N, g_0_new, e_0_new, g_1_new, e_1_new, g_2_new, e_2_new)

            ! Check photon 
            if (mod(index, 10) == 0) then ! Only checks for a photon every 10 time steps 

                ! Do statistics here 
                psi_0 = modulo_func(g_0_new)**2 + modulo_func(e_0_new)**2

                do j = 1,(N-1) 
                    psi_0 = psi_0 + modulo_func(g_1_new(j))**2 + modulo_func(e_1_new(j))**2
                end do 

                do j = 1,(N-2)
                    do k = (j+1),(N-1) 
                        psi_0 = psi_0 + modulo_func(g_2_new(j,k))**2 + modulo_func(e_2_new(j,k))**2
                    end do 
                end do 

                psi_1 = modulo_func(g_1_new(N))**2 + modulo_func(e_1_new(N))**2

                do j = 1,(N-1)
                    psi_1 = psi_1 + modulo_func(g_2_new(j,N))**2 + modulo_func(e_2_new(j,N))**2
                end do 

                prob = psi_1 / (psi_1 + psi_0)
                
                ! Grab a random number from a pre existing list of random numbers
                rand_num = rand_list(index)

                ! Check if photon is in the Nth box
                if (rand_num <= prob) then 

                    g_0 = g_1_new(N)
                    e_0 = e_1_new(N)

                    g_1(1) = 0
                    e_1(1) = 0

                    do j = 2,N 
                        g_1(j) = g_2_new(j-1,N)
                        e_1(j) = e_2_new(j-1,N)
                    end do 

                    g_2 = 0.0; e_2 = 0.0

                    call normalise(total, N, g_0, e_0, g_1, e_1, g_2, e_2)

                    ! Does all the photon waiting time distribution here 
                    current_time = time_list(index)
                    call calc_waiting_time(bin_width, time_steps, waiting_time_list, time_list, current_time)

                else ! Photon not found 
                    
                    g_0 = g_0_new
                    e_0 = e_0_new 

                    g_1(1) = 0 
                    e_1(1) = 0

                    do j = 2,N 
                        g_1(j) = g_1_new(j-1)
                        e_1(j) = e_1_new(j-1)
                    end do 

                    g_2 = 0; e_2 = 0

                    do k = 2,N 
                        g_2(1,k) = 0 
                        e_2(1,k) = 0 
                    end do  

                    do j = 2,(N-1)
                        do k = (j+1),N 
                            g_2(j,k) = g_2_new(j-1, k-1)
                            e_2(j,k) = e_2_new(j-1, k-1)
                        end do 
                    end do 
                
                    ! Call normalisation here 
                    call normalise(total, N, g_0, e_0, g_1, e_1, g_2, e_2)

                end if 

            else ! It's not a multiple of 10 time step
                ! In this case we just update it with Runge-Kutta 

                g_0 = g_0_new 
                e_0 = e_0_new
                g_1 = g_1_new 
                e_1 = e_1_new
                g_2 = g_2_new 
                e_2 = e_2_new 


            end if 

        end do 

        if (mod(sim, 100) == 0) then 
            print *, sim ,' simulations completed.'
        end if         

    end do 

    ! Normalise the waiting distribution array 
    waiting_total = sum(waiting_time_list)
    waiting_time_list = waiting_time_list / waiting_total 
    
    ! Make time array for waiting time list with correct length 
    call linspace(start=0.0, end=end_time, time_list=reduced_time_list)

    ! Write out final result to a txt file
    open(1, file="g2.txt", status="replace")
    do index = 1,size(reduced_time_list)
        write(1,*) reduced_time_list(index), waiting_time_list(index)
    end do 
    
    ! Find end time of simulation 
    call system_clock(end)

    print *, "All simulations completed. Execution time: ", real(end - beginning) / real(rate), " seconds."

    ! -------------------------------------------------------------------------------------------
    ! 
    !   Functions and Subroutines 
    !
    !-------------------------------------------------------------------------------------------
    
    contains

    subroutine linspace(start, end, time_list)

        ! Linspace is created by  

        real, intent(in) :: start
        integer, intent(in) :: end 
        real, intent(out) :: time_list(:)
        real :: range
        integer :: n, i

        n = size(time_list)
        range = end - start
        
        do i = 1,n
            time_list(i) = start + (range * (i - 1) / (n - 1))
        end do

    end subroutine

    ! function return_total (N, g_0, e_0, g_1, e_1, g_2, e_2) result(total)

    !     ! Takes in all the coefficient arrays and returns the total magnitude

    !     implicit none
        
    !     real:: total
    !     integer :: N, j, k
    !     complex :: g_0, e_0
    !     complex, dimension(N) :: g_1, e_1
    !     complex, dimension(N,N) :: g_2, e_2

    !     total = modulo_func(g_0) + modulo_func(e_0)

    !     do j = 1,N
    !         total = total + modulo_func(g_1(j)) + modulo_func(e_1(j))
    !     end do 

    !     do j = 1,N 
    !         do k = 1,N
    !             total = total + modulo_func(g_2(j,k)) + modulo_func(e_2(j,k))
    !         end do 
    !     end do

    ! end function 

    subroutine normalise (total, N, g_0, e_0, g_1, e_1, g_2, e_2)

        implicit none 

        real :: total 
        integer :: N, j, k
        complex :: g_0, e_0
        complex, dimension(N) :: g_1, e_1
        complex, dimension(N,N) :: g_2, e_2

        total = modulo_func(g_0) + modulo_func(e_0)

        do j = 1,N
            total = total + modulo_func(g_1(j)) + modulo_func(e_1(j))
        end do 

        do j = 1,N 
            do k = 1,N
                total = total + modulo_func(g_2(j,k)) + modulo_func(e_2(j,k))
            end do 
        end do

        g_0 = g_0 / total 
        e_0 = e_0 / total 
        g_1 = g_1 / total 
        e_1 = e_1 / total 
        g_2 = g_2 / total 
        e_2 = e_2 / total 

    end subroutine

    subroutine normalise_new (total, N, g_0_new, e_0_new, g_1_new, e_1_new, g_2_new, e_2_new)

        implicit none 

        real :: total 
        integer :: N, j, k
        complex :: g_0_new, e_0_new
        complex, dimension(N) :: g_1_new, e_1_new
        complex, dimension(N,N) :: g_2_new, e_2_new

        total = modulo_func(g_0_new) + modulo_func(e_0_new)

        do j = 1,N
            total = total + modulo_func(g_1_new(j)) + modulo_func(e_1_new(j))
        end do 

        do j = 1,N 
            do k = 1,N
                total = total + modulo_func(g_2_new(j,k)) + modulo_func(e_2_new(j,k))
            end do 
        end do

        g_0_new = g_0_new / total 
        e_0_new = e_0_new / total 

        g_1_new = g_1_new / total 
        e_1_new = e_1_new / total 

        g_2_new = g_2_new / total 
        e_2_new = e_2_new / total 

    end subroutine

    function modulo_func(z) result(c)

        ! Takes in a complex number z and returns its modulus

        implicit none

        ! Declare var types
        real :: a, b, c
        complex , intent(in) :: z

        a = real(z)
        b = aimag(z)

        c = sqrt(a**2 + b**2)

    end function 

    subroutine initialise_arrays(N, g_0, g_0_k1, g_0_k2, g_0_k3, g_0_k4, g_0_new, e_0 , & 
                                e_0_k1, e_0_k2, e_0_k3, e_0_k4, e_0_new, g_1, g_1_k1, & 
                                g_1_k2, g_1_k3, g_1_k4, g_1_new, e_1 ,e_1_k1, e_1_k2, & 
                                e_1_k3, e_1_k4, e_1_new, g_2, g_2_k1, g_2_k2, g_2_k3, &
                                g_2_k4, g_2_new, e_2 ,e_2_k1, e_2_k2, e_2_k3, e_2_k4, &
                                e_2_new)

        ! This stores the initial condition of the simulations

        ! Declare types 
        integer :: N
        complex :: g_0, g_0_k1, g_0_k2, g_0_k3, g_0_k4, g_0_new, e_0 ,e_0_k1, e_0_k2, e_0_k3, e_0_k4, e_0_new
        complex, dimension(N) :: g_1, g_1_k1, g_1_k2, g_1_k3, g_1_k4, g_1_new, e_1 ,e_1_k1, e_1_k2, e_1_k3, e_1_k4, e_1_new 
        complex, dimension(N,N) :: g_2, g_2_k1, g_2_k2, g_2_k3, g_2_k4, g_2_new, e_2 ,e_2_k1, e_2_k2, e_2_k3, e_2_k4, e_2_new 

        g_0 = 1
        g_0_k1 = 0; g_0_K2 = 0; g_0_K3 = 0; g_0_K4 = 0; g_0_new = 0
        e_0 = 0; e_0_k1 = 0; e_0_K2 = 0; e_0_K3 = 0; e_0_K4 = 0; e_0_new = 0

        g_1 = 0; g_1_k1 = 0; g_1_K2 = 0; g_1_K3 = 0; g_1_K4 = 0; g_1_new = 0
        e_1 = 0; e_1_k1 = 0; e_1_K2 = 0; e_1_K3 = 0; e_1_K4 = 0; e_1_new = 0

        g_2 = 0; g_2_k1 = 0; g_2_K2 = 0; g_2_K3 = 0; g_2_K4 = 0; g_2_new = 0
        e_2 = 0; e_2_k1 = 0; e_2_K2 = 0; e_2_K3 = 0; e_2_K4 = 0; e_2_new = 0

    end subroutine

    subroutine calc_waiting_time (bin_width, time_steps, waiting_time_list, time_list, current_time)

        integer :: bin_width, i
        integer :: time_steps 
        real :: current_time
        real, dimension(bin_width) :: waiting_time_list
        real, dimension(time_steps) :: time_list 

        ! Finds the right bin to add the statistic 

        if (current_time <= waiting_time_list(2)) then 

            waiting_time_list(1) = waiting_time_list(1) + 1

        else if (current_time > waiting_time_list(size(waiting_time_list) - 1) .AND. & 
            current_time <= waiting_time_list(size(waiting_time_list))) then 

            waiting_time_list(size(waiting_time_list)) = waiting_time_list(size(waiting_time_list)) + 1

        else

            ! Go into a do loop here 
            do i = 2, (size(waiting_time_list) - 2)

                if (((current_time > waiting_time_list(i) .AND. & 
                    current_time <= waiting_time_list(i+1)))) then 

                    waiting_time_list(i) = waiting_time_list(i) + 1

                    exit 

                end if 

            end do 

        end if  

    end subroutine


end program main 



