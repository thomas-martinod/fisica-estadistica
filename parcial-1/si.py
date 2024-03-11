import numpy as np

pi = np.pi
e = np.exp(1)
h = 6.62607e-34
kb = 1.38065e-23

T = 298
V = 4.11e-26

m = np.array([6.64e-27, 3.35e-26, 4.49e-26, 2.28e-25, 6.63e-26, 1.49e-26, 1.32e-25, 1.99e-26])
Qe = np.array([1.00, 1.00, 2.00, 1.00, 1.00, 1.00, 2.00, 1.00])

N = 6.023e23

Gt = np.zeros(8)
Ut = np.zeros(8)
St = np.zeros(8)

Ge = np.zeros(8)
Ue = np.zeros(8)
Se = np.zeros(8)

Gte = np.zeros(8)
Ute = np.zeros(8)
Ste = np.zeros(8)

# Calculations with translational partition function
with open('Funcion_Traslacional.dat', 'w') as f:
    f.write("          G                          U                         S\n")
    for i in range(8):
        G = 8.43e-5 * T * (-np.log((((2 * pi * m[i] * kb * T) / (h ** 2)) ** (1.5)) * V) + 1)
        U = N * kb * T * 1.5 / 1000
        S = N * ((1.5 * kb) + (kb * np.log((((2 * pi * m[i] * kb * T) / (h ** 2)) ** (1.5)) * V)))
        f.write(f"{G} {U} {S}\n")
        Gt[i] = G
        Ut[i] = U
        St[i] = S

# Calculations with electronic partition function
with open('Funcion_Electronica.dat', 'w') as f:
    f.write("          G                          U                         S\n")
    for i in range(8):
        G = ((-kb * T * np.log(Qe[i])) + (kb * T))
        U = kb * (T ** 2) * (1 / Qe[i]) * 10e9
        S = N / 1000 * kb * T * ((1 / Qe[i]) + (kb * np.log(Qe[i])))
        f.write(f"{G} {U} {S}\n")
        Ge[i] = G
        Ue[i] = U
        Se[i] = S

# Calculations with both translational and electronic partition functions
with open('Funcion_Traslacional+Electronica.dat', 'w') as f:
    f.write("          G                          U                         S\n")
    for i in range(8):
        G = 8.31e-3 * T * np.log((((((2 * pi * m[i] * kb * T) / (h ** 2)) ** (1.5)) * V) * Qe[i]) / N)
        U = Ue[i] + Ut[i]
        S = Se[i] + St[i]
        f.write(f"{G} {U} {S}\n")
        Gte[i] = G
        Ute[i] = U
        Ste[i] = S
