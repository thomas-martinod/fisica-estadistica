import numpy as np

hbar = 1.05457182 * 10**(-34)
kb = 1.380649 * 10**(-23)
T = 298.15

th = np.deg2rad(32.049)
phi = np.deg2rad([90, 210, 330, 90, 210, 330])
m = 1.673*10**(-27)
L = 1.8605*10**(-10)

sum_A = 0
for i in phi:
    sum_A += np.sin(th)**2 * np.sin(i)**2 + np.cos(th)**2

I_A = m * L**2 * sum_A
print(I_A)


sum_B = 0
for i in phi:
    sum_B += np.sin(th)**2 * np.cos(i)**2 + np.cos(th)**2

I_B = m * L**2 * sum_B
print(I_B)

I_C = 6 * m * L**2 * np.sin(th)**2
print(I_C)


throt_A = hbar**2 / (2 * I_A * kb)
throt_B = hbar**2 / (2 * I_B * kb)
throt_C = hbar**2 / (2 * I_C * kb)

print(throt_A)
print(throt_B)
print(throt_C)

q_rot = np.sqrt(np.pi)/6 * (T**3 / (throt_A*throt_B*throt_C))**(1/2)
print(q_rot)