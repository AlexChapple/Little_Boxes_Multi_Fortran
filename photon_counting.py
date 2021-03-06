import matplotlib.pyplot as plt 
import numpy as np 


photon_data = np.loadtxt("photon_counting.txt")

photon_count = []
photon_count_hist = []

photon_count_data = photon_data

x_list = [i for i in range(0,8)]



plt.bar(x_list, photon_count_data)
plt.xlabel("Photon number")
plt.ylabel("frequency (a.u)")
plt.title("photon counting distribution")
plt.savefig("Figures/photon_counting.png", dpi=600)