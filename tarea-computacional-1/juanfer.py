import numpy as np
N_A = 6.02214076e+23        #s Avogadro's Number (1 / mol)
k_B = 1.380649e-23          # Boltzmann constant (J / K)
h = 6.62607015e-34          # Plank's constant (J s)

print(N_A*  k_B)



m = 4.0026 / 1000 / N_A

T = 298.15

P = 10**5

V = N_A * k_B * T / P
print(V)

q_tras = (2* np.pi * m* k_B * T / h**2)**(1.5) * V
print(q_tras)

S_th = N_A * k_B * np.log(q_tras * np.exp(5/2) * 1 / N_A)
print(S_th)

S_jf = 2.5 * N_A * k_B + N_A * k_B* np.log(q_tras)
print(S_jf)

S_total_jf =  S_th * 1
print(S_total_jf)