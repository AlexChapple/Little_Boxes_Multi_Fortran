! Little Boxes Multi with 4th order Runge-Kutta method adapted into Fortran code
!
!
!
!

program main

    ! The equivalent of the average simulations method in the julia code
    ! Once difference is instead of plotting the function using julia, 
    ! we're going to write out the final result and use python to plot it

    implicit none

    ! Declare variables and parameters
    integer, parameter :: N = 20
    integer, parameter :: time_steps = 10000
    integer, parameter :: end_time = 8
    integer, parameter :: num_of_simulations = 500 
    real, parameter :: pi = 3.1415927
    real, parameter :: phase = 0.0 
    real, parameter :: gammaL = 0.5 
    real, parameter :: gammaR = 0.5 
    complex, parameter :: Omega = cmplx(10 * pi, 0)
    real, parameter :: dt = real(end_time) / real(time_steps) 
    real, parameter :: tau = 10.0 * dt * real(N) 
    real :: total 
    integer :: sim, index, j, k, beginning, end, rate 
    real, dimension(time_steps) :: time_list, spin_up_list, spin_down_list, rand_list 

    ! The coefficients (g for ground, u for up)
    complex :: g_0, g_0_k1, g_0_k2, g_0_k3, g_0_k4, g_0_new, e_0 ,e_0_k1, e_0_k2, e_0_k3, e_0_k4, e_0_new
    complex, dimension(N) :: g_1, g_1_k1, g_1_k2, g_1_k3, g_1_k4, g_1_new, e_1 ,e_1_k1, e_1_k2, e_1_k3, e_1_k4, e_1_new 
    complex, dimension(N,N) :: g_2, g_2_k1, g_2_k2, g_2_k3, g_2_k4, g_2_new, e_2 ,e_2_k1, e_2_k2, e_2_k3, e_2_k4, e_2_new 
    complex :: lambdaL, lambdaR
    real :: psi_0, psi_1, prob, rand_num, spin_up_prob, spin_down_prob, spin_total ! spin_total is just the total probability of spin up and down for normalisation purposes 
    
    spin_up_list = 0 
    spin_down_list = 0

    ! Program execution time 
    call system_clock(beginning, rate)

    ! Construct time_list
    call linspace(start=0.0, end=end_time, time_list=time_list) ! This makes the time list 

    ! A do loop will go through and do the simulations here
    do sim = 1, num_of_simulations

        ! Instead of a function called each time here, instead will do everything in this program itself
        call initialise_arrays(N, g_0, g_0_k1, g_0_k2, g_0_k3, g_0_k4, g_0_new, e_0 , & 
                            e_0_k1, e_0_k2, e_0_k3, e_0_k4, e_0_new, g_1, g_1_k1, & 
                            g_1_k2, g_1_k3, g_1_k4, g_1_new, e_1 ,e_1_k1, e_1_k2, & 
                            e_1_k3, e_1_k4, e_1_new, g_2, g_2_k1, g_2_k2, g_2_k3, &
                            g_2_k4, g_2_new, e_2 ,e_2_k1, e_2_k2, e_2_k3, e_2_k4, e_2_new)

        ! Construct random number list 
        call random_number(rand_list)

        lambdaL = exp(cmplx(0, phase / 2)) * sqrt(gammaL) * sqrt(N/tau)
        lambdaR = exp(cmplx(0, -phase / 2)) * sqrt(gammaR) * sqrt(N/tau)

        psi_0 = 0; psi_1 = 0; prob = 0; rand_num = 0; spin_up_prob = 0; spin_down_prob = 0; spin_total = 0; total = 0
  
        ! Change some variables into integer 

        do index = 1, size(time_list)

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
                    g_2_k1(j,k) = (-Omega/2)*e_2(j,k)
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
                do k = (j+1),N 
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
                do k = (j+1),N 
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
                do k = (j+1),N 
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
            total = return_total(N, g_0_new, e_0_new, g_1_new, e_1_new, g_2_new, e_2_new)

            ! if (index <= 15) then 
            !     print *, total, e_2_new 
            ! end if 

            g_0_new = g_0_new / total 
            e_0_new = e_0_new / total 

            g_1_new = g_1_new / total 
            e_1_new = e_1_new / total 

            g_2_new = g_2_new / total 
            e_2_new = e_2_new / total 

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

                prob = psi_1 / (psi_0 + psi_1)
                
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

                    g_2 = 0.0
                    e_2 = 0.0

                    ! Normalisation 
                    ! total = modulo_func(g_0) + modulo_func(e_0)
                    
                    ! do j = 1,N 
                    !     total = total + modulo_func(g_1(j)) + modulo_func(e_1(j))
                    ! end do 

                    total = return_total(N, g_0, e_0, g_1, e_1, g_2, e_2)

                    g_0 = g_0 / total 
                    e_0 = e_0 / total 
                    g_1 = g_1 / total 
                    e_1 = e_1 / total 
                    g_2 = g_2 / total 
                    e_2 = e_2 / total 

                    ! Calculate spin up and down probability here 
                    spin_down_prob = modulo_func(g_0)
                    spin_up_prob = modulo_func(e_0) 
 
                    do j = 1,N 
                        spin_down_prob = spin_down_prob + modulo_func(g_1(j))**2
                        spin_up_prob = spin_up_prob + modulo_func(e_1(j))**2
                    end do 

                    do j = 1,N ! NOTE: This can possibly be optimised slightly and reduce redundancy 
                        do k = 1,N 
                            spin_down_prob = spin_down_prob + modulo_func(g_2(j,k))**2
                            spin_up_prob = spin_up_prob + modulo_func(e_2(j,k))**2
                        end do 
                    end do

                    ! Normalisation 
                    spin_total = spin_up_prob + spin_down_prob
                    spin_up_prob = spin_up_prob / spin_total 
                    spin_down_prob = spin_down_prob / spin_total

                    ! Store the information in the local spin_up/down_list 
                    spin_up_list(index) = spin_up_list(index) + spin_up_prob 
                    spin_down_list(index) = spin_down_list(index) + spin_down_prob

                else ! Photon not found 
                    
                    g_0 = g_0_new
                    e_0 = e_0_new 

                    g_1(1) = 0 
                    e_1(1) = 0

                    do j = 2,N 
                        g_1(j) = g_1_new(j-1)
                        e_1(j) = e_1_new(j-1)
                    end do 

                    do k = 1,N 
                        g_2(1,k) = 0 
                        e_2(1,k) = 0 
                    end do  

                    do j = 2,(N-1)
                        do k = (j+1),N 
                            g_2(j,k) = g_2_new(j-1, k-1)
                            e_2(j,k) = e_2_new(j-1, k-1)
                        end do 
                    end do 
                
                    ! Normalisation 
                    ! total = modulo_func(g_0) + modulo_func(e_0)

                    ! do j = 2,N 
                    !     total = total + modulo_func(g_1(j)) + modulo_func(e_1(j))
                    ! end do 

                    ! do j = 2,(N-1)
                    !     do k = (j+1),N 
                    !         total = total + modulo_func(g_2(j,k)) + modulo_func(e_2(j,k))
                    !     end do 
                    ! end do 

                    total = return_total(N, g_0, e_0, g_1, e_1, g_2, e_2)

                    ! print *, g_0, e_0, g_1, e_1, g_2, e_2 

                    g_0 = g_0 / total 
                    e_0 = e_0 / total 
                    g_1 = g_1 / total 
                    e_1 = e_1 / total 
                    g_2 = g_2 / total 
                    e_2 = e_2 / total 

                    ! Calculate up and down probability here
                    spin_down_prob = modulo_func(g_0) 
                    spin_up_prob = modulo_func(e_0) 

                    do j = 1,N 
                        spin_down_prob = spin_down_prob + modulo_func(g_1(j))**2 
                        spin_up_prob = spin_up_prob + modulo_func(e_1(j))**2
                    end do 

                    do j = 1,N ! NOTE: This can possibly be optimised slightly and reduce redundancy 
                        do k = 1,N 
                            spin_down_prob = spin_down_prob + modulo_func(g_2(j,k))**2
                            spin_up_prob = spin_up_prob + modulo_func(e_2(j,k))**2
                        end do 
                    end do

                    ! Normalisation 
                    spin_total = spin_up_prob + spin_down_prob
                    spin_up_prob = spin_up_prob / spin_total
                    spin_down_prob = spin_down_prob / spin_total

                    ! Store the information in the local spin_up/down_list 
                    spin_up_list(index) = spin_up_list(index) + spin_up_prob 
                    spin_down_list(index) = spin_down_list(index) + spin_down_prob

                end if 

            else ! It's not a multiple of 10 time step
                ! In this case we just update it with Runge-Kutta 
                ! Here I'm not renormalising it because I'm fairly certain it doesn't have to be 

                g_0 = g_0_new 
                e_0 = e_0_new
                g_1 = g_1_new 
                e_1 = e_1_new
                g_2 = g_2_new 
                e_2 = e_2_new 

                ! Find probability system is in spin up / down 
                spin_down_prob = modulo_func(g_0) 
                spin_up_prob = modulo_func(e_0) 

                do j = 1,N 
                    spin_down_prob = spin_down_prob + modulo_func(g_1(j))**2
                    spin_up_prob = spin_up_prob + modulo_func(e_1(j))**2
                end do 

                do j = 1,N 
                    do k = 1,N 
                        spin_down_prob = spin_down_prob + modulo_func(g_2(j,k))**2
                        spin_up_prob = spin_up_prob + modulo_func(e_2(j,k))**2
                    end do
                end do 

                ! Normalisation 
                spin_total = spin_up_prob + spin_down_prob
                spin_up_prob = spin_up_prob / spin_total 
                spin_down_prob = spin_down_prob / spin_total 

                ! Store the information in the local spin_up/down_list 
                spin_up_list(index) = spin_up_list(index) + spin_up_prob 
                spin_down_list(index) = spin_down_list(index) + spin_down_prob

            end if 

        end do 

        if (mod(sim, 100) == 0) then 
            print *, sim ,' simulations completed.'
        end if         

    end do 

    ! Normalise spin up and down list by number of simulations 
    spin_up_list = spin_up_list / num_of_simulations
    spin_down_list = spin_down_list / num_of_simulations

    !!! Write out final result to a txt file
    open(1, file="spin_up.txt", status="replace")
    open(2, file="spin_down.txt", status="replace")
    do index = 1,size(time_list)
        write(1,*) time_list(index), spin_up_list(index)
        write(2,*) time_list(index), spin_down_list(index)
    end do 
    
    call system_clock(end)

    print *, "All simulations completed. Execution time: ", real(end - beginning) / real(rate), " seconds."

    ! -------------------------------------------------------------------------------------------
    ! 
    !   Functions and Subroutines 
    !
    !-------------------------------------------------------------------------------------------
    
    contains

    subroutine linspace(start, end, time_list) 

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

    function return_total (N, g_0, e_0, g_1, e_1, g_2, e_2) result(total)

        ! Takes in all the coefficient arrays and returns the total magnitude

        implicit none
        
        real:: total
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

    end function 

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


end program main 



