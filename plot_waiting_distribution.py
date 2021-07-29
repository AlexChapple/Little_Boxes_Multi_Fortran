import matplotlib.pyplot as plt 
import numpy as np 


waiting_data = np.loadtxt("waiting_distribution.txt")

time_list = []
waiting_time_list = []

time_list = waiting_data[:,0]
waiting_time_list = waiting_data[:,1]

reduced_time_list = []
reduced_waiting_time_list = []

for i in range(len(waiting_time_list)):

    if i <= 1:

        reduced_time_list.append(time_list[i])
        reduced_waiting_time_list.append(waiting_time_list[i])

    elif i % 20000 == 0:

        reduced_time_list.append(time_list[i])
        reduced_waiting_time_list.append(waiting_time_list[i])


plt.figure(1)
plt.scatter(reduced_time_list, reduced_waiting_time_list, linewidth=0.5, marker=".")
plt.xlabel("Waiting Time (s)")
plt.ylabel("Waiting time distribution")
plt.title("waiting time distribution")
plt.savefig("Figures/waiting_distribution_s.png", dpi=600)

plt.figure(2)
plt.plot(reduced_time_list, reduced_waiting_time_list, linewidth=1)
plt.scatter(reduced_time_list, reduced_waiting_time_list, linewidth=1, marker="o", s=10, facecolors="none", edgecolors="red")
plt.xlabel("Waiting Time (s)")
plt.ylabel("Waiting time distribution")
plt.title("waiting time distribution")
plt.savefig("Figures/waiting_distribution_p.png", dpi=600)