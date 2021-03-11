import matplotlib.pyplot as plt

fonts = {"font.size": 14}
plt.rcParams.update(fonts)

def read_file(filename):
    prob = []; r = []
    with open("Results/onebodydensity/"+filename, "r") as infile:
        lines = infile.readlines()
        for i in range(len(lines)):
            line = lines[i]
            vals = line.split()
            prob.append(float(vals[0]))
            r.append(float(vals[1]))

    for j in range(len(prob)):
        if prob[j] > 1.0:
            del r[j]
    for i in prob[:]:
        if i > 1.0:
            prob.remove(i)

    return prob, r


probEJ10, rEJ10 = read_file("jastrow/N10D3")
probE10, rE10 = read_file("jastrow/N10D3")


plt.figure()
plt.title("One-body density for N=10 with \n and without Jastrow factor")
plt.plot(rE10, probE10, "+", label="N=10, no-Jastrow")
plt.plot(rEJ10, probEJ10, "+", label="N=10, Jastrow")
plt.xlabel("r/$a_{ho}$")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/onebody_N_10.png")


plt.show()