import matplotlib.pyplot as plt
import OutRead as OR

kp,bs = OR.read_EIGENVAL('Test/EIGENVAL')
plt.plot(bs)
plt.show()
