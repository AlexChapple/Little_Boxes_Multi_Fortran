import matplotlib.pyplot as plt 

# Open data files 
spin_up_data =  open("spin_down.txt", "r") 
spin_down_data = open("spin_up.txt", "r")

# Initialises arrays 
time_list = []
spin_up_prob_list = []
spin_down_prob_list = []

# Extract data from data files
for line in spin_up_data:
    p = line.split()
    time_list.append(p[1])
    spin_up_prob_list.append(p[2])

for line in spin_down_data:
    p = line.split()
    spin_down_prob_list.append(p[2])

spin_up_data.close()
spin_down_data.close()

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