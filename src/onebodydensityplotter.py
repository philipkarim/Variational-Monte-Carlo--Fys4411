import matplotlib.pyplot as plt
import numpy as np

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


probEJ10, rEJ10 = read_file("jastrow/N10Dim3Jold.txt")
probE10, rE10 = read_file("nonjastrow/N10Dim3old.txt")


plt.figure()

plt.plot(rE10, probE10, "ro", markersize=3,label="N=10, no-Jastrow")
plt.plot(rEJ10, probEJ10, "bo",markersize=3, label="N=10, Jastrow")


#plt.plot(rE10, probE10, label="N=10, no-Jastrow")
#plt.plot(rEJ10, probEJ10,  label="N=10, Jastrow")
plt.xlabel("r/$a_{ho}$")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/onebody_N_10.png")
plt.show()

"""
def plot_one_body_density():
    
    Plot the radial one-body density based on discrete bins with recordings
    from long simulations of 2**20 samples. The positions were recorded to
    these bins. The files are loaded from the directory ./Data/onebodydensity
    Returns: None, but plots and shows the result.
    
    max = 4.0
    DIR = "../Data/onebodydensity/"
    fns = ["Jastrow.csv", "NoJastrow.csv"]
    labs = ["Correlation", "Without correlation"]
    fmts = ["--", ":"]
    for fn, label, fmt in zip(fns, labs, fmts):
        df = pd.read_csv(DIR + fn)
        px, py, pz = df["x"].values, df["y"].values, df["z"].values
        bins = px.shape[0]
        num = int(bins / 2)
        pr = np.zeros(num)
        r = np.linspace(0, max, num)
        total = 0
        for i in range(num):
            # We have to do this "strange" couting method because we defined
            # the bins to start from -4, hence the middle indices will actually
            # be where the radius is the least.
            bin = (px[i] + px[-1 - i])**2  # Add from both sides of the origin
            bin += (py[i] + py[-1 - i])**2
            bin += (pz[i] + pz[-1 - i])**2
            rval = np.sqrt(bin)
            total += rval  # Scale factor so that integral sums up to 1
            pr[-i - 1] = rval
        pr /= total
        plt.plot(r, pr, fmt, label=label)
    plt.legend()
    plt.xlabel(r"$|r|$")
    plt.ylabel(r"$p(r)$", rotation=0)
    plt.show()
    return None


import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import numpy as np


#create data
#x = np.array([5.00E-07, 1.40E-06, 4.10E-06, 1.25E-05, 3.70E-05, 1.11E-04, 3.33E-04, 1.00E-03])
#y= np.array([494.55, 333.4666667, 333.3333333, 333.1, 303.4966667, 197.7533333, 66.43333333, 67.715])

x, y = read_file("jastrow/N10Dim3Jold.txt")
print(y)
print(x)

#define x as 200 equally spaced values between the min and max of original x 
xnew = np.linspace(0, 4, 400) 

#xnew=x

#define spline
spl = make_interp_spline(x, y, k=3)
y_smooth = spl(xnew)


#create smooth line chart 
plt.plot(xnew, y_smooth)
plt.show()
"""

"""
x, y = read_file("jastrow/N10Dim3Jold.txt")

from sklearn.svm import SVC
#... load the data into X,y
model = SVC(kernel='poly')
model.fit(x,y)


from scipy import interpolate
rE10, probEJ10 = read_file("jastrow/N10Dim3Jold.txt")

#x = np.arange(0, 10)
#y = np.exp(-x/3.0)
f = interpolate.interp1d(x, y)

xnew = np.arange(0,4, 0.1)
ynew = f(xnew)   # use interpolation function returned by `interp1d`
plt.plot(x, y, 'o', xnew, ynew, '-')
plt.show()
"""