import matplotlib.pyplot as plt
import numpy as np
import OutRead as OR

kp,bs = OR.read_EIGENVAL('Test/EIGENVAL')
plt.plot(bs)
#plt.show()
plt.figure()
DOS = np.loadtxt('Test/DOSCAR',skiprows = 6, max_rows=1000)
plt.plot(DOS[:1000,0],DOS[:1000,1:])
plt.show()
print(DOS.shape)
