import numpy as np
eV = 1.60218e-19  # 1 eV in joules
kB = 1.380649e-23  # Boltzmann constant in J/K


def fixed_point_iteration(f, x0, tolerance=1e-6, max_iterations=100):
    """
    Solves x = f(x) using the fixed-point iteration method.

    Parameters:
    f (function): The function for which we are trying to find a fixed point.
    x0 (float): Initial guess for the fixed point.
    tolerance (float): The tolerance level for convergence.
    max_iterations (int): Maximum number of iterations.

    Returns:
    float: Approximate solution to x = f(x).
    int: Number of iterations performed.
    """
    x_current = x0
    for iteration in range(max_iterations):
        x_next = f(x_current)
        if abs(x_next - x_current) < tolerance:
            return x_next, iteration + 1
        x_current = x_next
    raise ValueError(f"Did not converge within {max_iterations} iterations")

# Constants
mu_F = 8.852369811566442e-19  # Fermi energy, you need to set the actual value
T = 298.15  # Temperature in K, you need to set the actual value

# Define the function f(mu)
def f(mu):
    return mu_F * (1 - (np.pi**2 * kB**2 * T**2) / (12 * mu**2))

# Initial guess
mu0 = mu_F  # A reasonable initial guess might be the Fermi energy

# Solve the equation mu = f(mu)
solution, iterations = fixed_point_iteration(f, mu0)

print(f"Solution: mu = {solution}")
print(f"Number of iterations: {iterations}")
print(f"In electron-volts: {solution/eV}")
print(f"Initial: {mu0/eV}")