import matplotlib.pyplot as plt 
import numpy as np 


waiting_data = np.loadtxt("waiting_distribution.txt")

time_list = []
waiting_time_list = []

time_list = waiting_data[:,0]
waiting_time_list = waiting_data[:,1]

plt.plot(time_list, waiting_time_list)
plt.xlabel("Time (s)")
plt.ylabel("Waiting time")
plt.savefig("Figures/waiting_distribution.png", dpi=600)