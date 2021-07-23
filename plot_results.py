import matplotlib.pyplot as plt 
import numpy as np 

# Open data files 
spin_up_data =  np.loadtxt("spin_up.txt")
spin_down_data = np.loadtxt("spin_down.txt")

# Initialises arrays 
time_list = []
spin_up_prob_list = []
spin_down_prob_list = []

time_list = spin_down_data[:, 0]
spin_down_prob_list = spin_down_data[:, 1]
spin_up_prob_list = spin_up_data[:,1]

# Plotting function 
plt.figure(1)
plt.plot(time_list, spin_down_prob_list)
plt.xlabel("Time (s)")
plt.ylabel("Probability spin down")
plt.savefig("Figures/spin_down.png", dpi=600)

plt.figure(2)
plt.plot(time_list, spin_up_prob_list)
plt.xlabel("Time (s)")
plt.ylabel("Probability spin up")
plt.savefig("Figures/spin_up.png", dpi=600)