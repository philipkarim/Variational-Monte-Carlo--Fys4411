import matplotlib.pyplot as plt
import numpy as np

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline, make_pipeline

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


"""
probEJ10, rEJ10 = read_file("jastrow/N10Dim3Jold.txt")
probE10, rE10 = read_file("nonjastrow/N10Dim3old.txt")
"""
probEJ10, rEJ10 = read_file("jastrow/N50Dim3new.txt")
probE10, rE10 = read_file("nonjastrow/N50Dim3new.txt")

plt.figure()

plt.plot(rE10, probE10, "ro", markersize=3,label="N=10, no-Jastrow")
plt.plot(rEJ10, probEJ10, "bo",markersize=3, label="N=10, Jastrow")


#x=rEJ10
#y=probEJ10

x=np.array(rEJ10)
y=np.array(probEJ10)

x = x.reshape(-1,1)
y = y.reshape(-1,1)

#plt.plot(rE10, probE10, label="N=10, no-Jastrow")
#plt.plot(rEJ10, probEJ10,  label="N=10, Jastrow")


plt.xlabel("r/$a_{ho}$")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/onebody_N_10.png")
plt.grid()
plt.show()

#probEJ10, rEJ10 = read_file("jastrow/N10Dim3Jold.txt")

x=rEJ10
y=probEJ10

x2=rE10
y2=probE10

x=np.array(x)
y=np.array(y)

x2=np.array(x)
y2=np.array(y)


x = x.reshape(-1,1)
y = y.reshape(-1,1)

model_1 = make_pipeline(PolynomialFeatures(degree = 6),LinearRegression())
model_1.fit(x.reshape(-1,1),y)
y=model_1.predict(x.reshape(-1,1))

model_2 = make_pipeline(PolynomialFeatures(degree = 7),LinearRegression())
model_2.fit(x2.reshape(-1,1),y2)
y2=model_2.predict(x2.reshape(-1,1))
#plt.figure(figsize = (8,5))
#plt.scatter(x,y, label = 'Data')
plt.plot(x,y, color = 'red', label = 'jastrow')
plt.plot(x2,y2, color = 'blue', label = 'no jastrow')
plt.xlabel("r/$a_{ho}$")
plt.ylabel("Probability")
plt.legend(), plt.show()