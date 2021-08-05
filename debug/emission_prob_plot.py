import matplotlib.pyplot as plt 
import numpy as np 

tau_list = [0.02, 0.05 ,0.1, 0.15, 0.2, 0.4, 0.8]
first_emission = [4.4 * 10 ** -5, 5.7*10**-4 ,2.5*10**-3, 2.7*10**-3, 2.5*10**-4, 1.2*10**-3, 6*10**-3]
second_emission = [1.25 * 10 ** -3, 4.5*10**-3, 8.63*10**-3, 4.11*10**-3, 0.1, 8*10**-2, 7*10**-2]


plt.figure(1)
for i in range(7):

    plt.scatter(first_emission[i], second_emission[i], label="tau = " + str(tau_list[i]), s=20, marker="x")
    plt.title("Probability of first emission vs second emission \n for different $\\tau$")
    plt.xlabel("rough average of first emission probability")
    plt.ylabel("rough average emission of subsequent emission")
    plt.legend()
    plt.savefig("emission_prob_1.png")

first_emission_log = [np.log(i) for i in first_emission]
second_emission_log = [np.log(i) for i in second_emission]

plt.figure(2)
for i in range(7):

    plt.scatter(first_emission_log[i], second_emission_log[i], label="tau = " + str(tau_list[i]), s=20, marker="x")
    plt.title("Log probability of first emission vs second emission \n for different $\\tau$")
    plt.xlabel("rough average of first emission log probability")
    plt.ylabel("rough average emission of subsequent log emission")
    plt.legend()
    plt.savefig("emission_prob_2.png")
