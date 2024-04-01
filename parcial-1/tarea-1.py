degeneracy_table = {'s^2': 1, 'p^6': 1, 'd^10': 1, 'p^2': 1, 'd^4': 1, 'p^1': 2, 'p^5': 4, 'd^1': 4,
                    'd^3': 4, 'p^4': 5, 'd^2': 5, 'p^3': 6, 'd^9': 6, 'd^8': 9, 'd^6': 9, 'd^7': 10, 'd^5': 10}

electron_config = str(input('Ingrese la terminación electrónica del elemento (ej: para el azufre S \'3p^4\'): '))[1:]
try:
    degeneracy = degeneracy_table[electron_config]
    print(degeneracy)
except:
    print('Pon una terminación electrónica válida la próxima :)')


