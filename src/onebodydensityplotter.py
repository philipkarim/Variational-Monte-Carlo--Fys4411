import matplotlib.pyplot as plt
import numpy as np

fonts = {"font.size": 14}
plt.rcParams.update(fonts)

#Defining the filenames
filenames=["jastrow/N10Dim3Jold.txt", "nonjastrow/N10Dim3old.txt", "jastrow/N50Dim3new.txt", "nonjastrow/N50Dim3new.txt"]

#Making som empty array to fill up with values from the files
lists_r=[np.empty(shape=(0)) for _ in range(4)]
lists_rho=[np.empty(shape=(0)) for _ in range(4)]

#Reading the files
for files in range(0,4):
    with open("Results/onebodydensity/"+filenames[files], "r") as infile:
        lines = infile.readlines()
        for i in range(len(lines)):
            line = lines[i]
            value = line.split()
            lists_r[files]=np.append(lists_r[files],float(value[1]))
            lists_rho[files]=np.append(lists_rho[files],float(value[0]))

#Just a file to remove too high values
def finetuning(rad, rho):
    rho = rho.tolist()
    rad = rad.tolist()
    for j in range(len(rho)):
        if rho[j] > 1.0:
            del rad[j]

    for i in rho[:]:
        if i > 1.0:
            rho.remove(i)
    return rad, rho

#making lists ready to be plotted
rad10J, rho10J=finetuning(lists_r[0], lists_rho[0])
rad10, rho10=finetuning(lists_r[1], lists_rho[1])
rad50J, rho50J=finetuning(lists_r[2], lists_rho[2])
rad50, rho50=finetuning(lists_r[3], lists_rho[3])


#Plotting from here on
fig, (ax1, ax2) = plt.subplots(2)
ax1.scatter(rad10J, rho10J, s=7, label="N=10, With Jastrow",  color='#d62728')
ax1.scatter(rad10, rho10, s=7, label="N=10, Without Jastrow", color='#1f77b4')
ax2.scatter(rad50J, rho50J, s=7, label="N=50, With Jastrow",  color='#d62728')
ax2.scatter(rad50, rho50, s=7, label="N=50, Without Jastrow", color='#1f77b4')

ax1.set_xlabel("r")
ax1.set_ylabel(r'$\rho(r)$')
ax1.legend(loc="best")
ax1.grid()
ax2.set_xlabel("r")
ax2.set_ylabel(r'$\rho(r)$')
ax2.legend(loc="best")
ax2.grid()

plt.show()

plt.scatter(rad50J, rho50J, s=7, label="N=50, With Jastrow",  color='#d62728')
plt.scatter(rad50, rho50, s=7, label="N=50, Without Jastrow", color='#1f77b4')
plt.xlabel("r")
plt.ylabel(r'$\rho(r)$')
plt.legend(loc="best")
plt.tight_layout()
plt.grid()
plt.show()

