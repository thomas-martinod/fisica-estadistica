import numpy as np

hbar = 1.05457182 * 10**(-34)
kb = 1.380649 * 10**(-23)
T = 298.15
m = 1.673*10**(-27)
L = 1.8605*10**(-10)
sigma = 24

I = 4*m*L**2
throt = hbar**2 / (2 * I * kb)
qrot = np.sqrt(np.pi) / sigma * (T / throt)**1.5

print(I)
print(throt)
print(qrot)