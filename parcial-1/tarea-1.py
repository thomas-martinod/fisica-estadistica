# Python code to extract the orbital term from an electronic configuration

# Read the electronic configuration.
electron_config = input("Enter the electronic configuration: ")

# Find the length of the string.
len_config = len(electron_config)

# Find the term referring to the orbital.
orbital_start = len_config
for i in range(len_config-1, -1, -1):
    if electron_config[i] == ' ':
        break
    if electron_config[i] == '^':
        orbital_start = i - 1
        break

# Extract the last term (orbital).
last_term = electron_config[orbital_start:len_config]

# Print the last term (orbital).
print("The orbital is:", last_term)

# Python code to determine and print the degeneracy of the last orbital term

# Define the degeneracy for each orbital configuration
degeneracy = {
    's^2': 1,
    'p^6': 1,
    'd^10': 1,
    'p^2': 1,
    'd^4': 1,
    'p^1': 2,
    'p^5': 4,
    'd^1': 4,
    'd^3': 4,
    'p^4': 5,
    'd^2': 5,
    'p^3': 6,
    'd^9': 6,
    'd^8': 9,
    'd^6': 9,
    'd^7': 10,
    'd^5': 10
}

# Find the degeneracy of the last orbital term
if last_term in degeneracy:
    print("The degeneracy in this configuration is", degeneracy[last_term])
else:
    print("Degeneracy for this configuration is not defined.")

