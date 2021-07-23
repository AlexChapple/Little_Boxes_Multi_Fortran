import matplotlib.pyplot as plt 





data =  open("data.txt", "r") 

x = []
y = []

for line in data:
    p = line.split()
    x.append(float(p[0]))
    y.append(float(p[1]))

data.close()


plt.plot(x,y)
plt.show()