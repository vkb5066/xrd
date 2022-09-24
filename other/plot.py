#Read
n, x, y = -1, [], []
with open("graph.xrd", 'r') as infile:
    res = infile.readline().split()
    n = int(res[0])
    for i in range(1, n + 1):
        x.append(float(res[i]))
        y.append(float(res[n + i]))

#Transform
if(0):
    from numpy import arcsin as asin
    pi = 3.14159
    wavelength = 1.5406 ##angstroms
    for i in range(0, n):
        arg = x[i]*wavelength/4./pi ## Q -> sin(theta)
        if(-1. < arg < +1.):
            x[i] = 2. * asin(arg) * 180./pi ## sin(theta) -> twoTheta in degrees
        else:
            print("invalid value at i =", i, "for arcsin()")
            y[i] = 0.


#Plot
from matplotlib import pyplot as plt
plt.plot(x, y, 'k')
plt.show()
