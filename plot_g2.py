import matplotlib.pyplot as plt 
import numpy as np 


waiting_data = np.loadtxt("g2.txt")

time_list = []
waiting_time_list = []

time_list = waiting_data[:,0]
waiting_time_list = waiting_data[:,1]

reduced_time_list = []
reduced_waiting_time_list = []

for i in range(len(waiting_time_list)):

    if i <= 10:

        reduced_time_list.append(time_list[i])
        reduced_waiting_time_list.append(waiting_time_list[i])

    elif i % 1000 == 0:

        reduced_time_list.append(time_list[i])
        reduced_waiting_time_list.append(waiting_time_list[i])



plt.scatter(reduced_time_list, reduced_waiting_time_list, linewidth=0.5, marker=".")
plt.xlabel("Time (s)")
plt.ylabel("Waiting time")
plt.title("waiting time distribution")
plt.savefig("Figures/g2.png", dpi=600)