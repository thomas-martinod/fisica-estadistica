import numpy as np

# Constantes
pi = np.pi                  # Pi
h = 6.62607015e-34          # Plank's constant (J s)
k_B = 1.380649e-23          # Boltzmann constant (J / K)
N_A = 6.02214076e+23        #s Avogadro's Number (1 / mol)


degeneracy_table = {'s^2': 1, 'p^6': 1, 'd^10': 1, 'p^2': 1, 'd^4': 1, 'p^1': 2, 'p^5': 4, 'd^1': 4,
                    'd^3': 4, 'p^4': 5, 'd^2': 5, 'p^3': 6, 'd^9': 6, 'd^8': 9, 'd^6': 9, 'd^7': 10, 'd^5': 10}


class Tras:
    def __init__(self, m, N, T, P) -> None:
        self.m = m / 1000 / N_A  # En kg
        self.N = N
        self.T = T
        self.P = P
        self.V = self.N * k_B * self.T / self.P
        self.q = (2 * pi * self.m * k_B * self.T / h**2)**(3/2) * self.V
        self.U = 1.5 * self.N * k_B * self.T
        self.S = self.N * k_B * (2.5 + np.log(self.q / self.N))
        self.A = -self.N * k_B * self.T * (np.log(self.q / self.N) + 1)
        self.G = -self.N * k_B * self.T * np.log(self.q / self.N)
        self.H = 2.5 * self.N * k_B * self.T
        self.mu = - k_B * self.T * np.log(self.q / self.N)

    def get_q(self):
        return self.q

    def get_U(self):
        return self.U

    def get_S(self):
        return self.S

    def get_A(self):
        return self.A

    def get_G(self):
        return self.G

    def get_H(self):
        return self.H

    def get_mu(self):
        return self.mu

    def print_all(self, num_decimals):
        print("------------------------------------------------------------------------")
        print(f"q =\t{self.q:.{num_decimals}f} [unit]")
        print(f"U =\t{self.U:.{num_decimals}f} [unit]")
        print(f"S =\t{self.S:.{num_decimals}f} [unit]")
        print(f"A =\t{self.A:.{num_decimals}f} [unit]")
        print(f"G =\t{self.G:.{num_decimals}f} [unit]")
        print(f"H =\t{self.H:.{num_decimals}f} [unit]")
        print(f"mu =\t{self.mu:.{num_decimals}f} [unit]")


class Ele:
    def __init__(self, ge1, N, T) -> None:
        self.ge1 = ge1
        self.N = N
        self.T = T
        self.U = 0
        self.S = self.N * k_B * self.T * np.log(self.ge1 * np.exp(1) / self.N)
        self.A = 0
        self.G = 0
        self.H = 0
        self.mu = 0

    def get_q(self):
        return self.ge1

    def get_U(self):
        return self.U

    def get_S(self):
        return self.S

    def get_A(self):
        return self.A

    def get_G(self):
        return self.G

    def get_H(self):
        return self.H

    def get_mu(self):
        return self.mu


class Tras_and_ele:
    def __init__(self, m, N, T, P, ge1) -> None:
        self.m = m / 1000 / N_A  # En kg
        self.N = N
        self.T = T
        self.P = P
        self.V = self.N * k_B * self.T / self.P
        self.q_tras = (2 * pi * self.m * k_B * self.T / h**2)**(3/2) * self.V
        self.q_ele = ge1
        self.U = 1.5 * self.N * k_B * self.T
        self.S = self.N * k_B * np.log(self.q_tras * self.q_ele * np.exp(2.5) / self.N)
        self.A = -self.N * k_B * self.T * (np.log(self.q_tras * self.q_ele / self.N) + 1)
        self.G = -self.N * k_B * self.T * np.log(self.q_tras * self.q_ele / self.N)
        self.H = 2.5 * self.N * k_B * self.T
        self.mu = - k_B * self.T * np.log(self.q_tras * self.q_ele / self.N)

    def get_q_tras(self):
        return self.q_tras

    def get_q_ele(self):
        return self.q_ele

    def get_U(self):
        return self.U

    def get_S(self):
        return self.S

    def get_A(self):
        return self.A

    def get_G(self):
        return self.G

    def get_H(self):
        return self.H

    def get_mu(self):
        return self.mu

def main():
    elem = input('Ingrese el elemento del gas monoatómico a considerar: ')

    aux_P = input('Ingrese la presión del sistema en [bar] (al ingresar \'st\' se asumen condiciones estándar): ')
    P_sys = 10e5 if aux_P == 'st' else float(aux_P) * 10e5

    aux_T = input('Ingrese la temperatura del sistema en [K] (al ingresar \'st\' se asumen condiciones estándar): ')
    T_sys = 298.15 if aux_T == 'st' else float(aux_T)

    aux_N = input('Ingrese el número de moles del sistema (al ingresar \'st\' se asume un mol): ')
    N_sys = N_A if aux_N == 'st' else float(aux_N)

    mass = float(input('Ingrese la masa atómica del ' + elem + ' en [g/mol]: '))

    electron_config = str(input('Ingrese la terminación electrónica del elemento (ej: para el azufre S \'3p^4\'): '))[1:]
    try:
        degeneracy = degeneracy_table[electron_config]
    except:
        print('Pon una terminación electrónica válida la próxima :)')
        quit()

    sys_tras = Tras(m = mass, N = N_sys, T = T_sys, P = P_sys)
    sys_ele = Ele(ge1 = degeneracy, N = N_sys, T = T_sys)
    sys_both = Tras_and_ele(m = mass, N = N_sys, T = T_sys, P = P_sys, ge1 = degeneracy)




if __name__ == "__main__":
    main()

