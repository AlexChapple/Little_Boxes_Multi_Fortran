import matplotlib.pyplot as plt 
import numpy as np 


data = np.loadtxt("photon_sequences.txt")

time_list = data[:,0]
photon_sequences = data[:,1]
within_list = data[:,2]

for i in range(len(photon_sequences)):

    if photon_sequences[i] == 1 and within_list[i] == 1:

        plt.axvline(time_list[i], c="red", lw=0.1)
    
    elif photon_sequences[i] == 1 and within_list[i] != 1:

        plt.axvline(time_list[i], c="blue", lw=0.1)



plt.xlabel("Time (s)")
plt.title("photon sequences of one simulation")
plt.savefig("photon_sequences.png", dpi=600)
