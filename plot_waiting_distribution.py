import matplotlib.pyplot as plt 
import numpy as np 


waiting_data = np.loadtxt("waiting_distribution.txt")

time_list = []
waiting_time_list = []

time_list = waiting_data[:,0]
waiting_time_list = waiting_data[:,1]

plt.scatter(time_list, waiting_time_list, linewidth=0.5, marker="x")
plt.xlabel("Time (s)")
plt.ylabel("Waiting time")
plt.title("waiting time distribution")
plt.savefig("Figures/waiting_distribution.png", dpi=600)