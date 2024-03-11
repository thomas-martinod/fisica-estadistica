PROGRAM PropiedadesTermodinamicas
    implicit none

    ! Constantes
    real, parameter :: pi = acos(-1.0)
    real, parameter :: R = 8.314462618 ! J/(mol*K) - Constante de los gases ideales
    real, parameter :: T = 298.0       ! K - Temperatura
    real, parameter :: P = 1.0         ! atm - Presión
    real, parameter :: h = 6.62607015e-34
    real, parameter :: kB = 1.380649e-23


    ! Función de partición translacional
    real function ParticionTranslacional(m)
        real, intent(in) :: m
        ParticionTranslacional = (2.0 * pi * m * kB * T / h**2) ** (3.0 / 2.0) * (kB * T / P)
    end function ParticionTranslacional

    ! Variable para la masa
    real :: masa

    ! Asigna un valor a la masa
    masa = 28.97e-3 ! Masa molar del nitrógeno en kg/mol

    ! Calcula la función de partición translacional para la masa dada
    print *, "Masa:", masa, "kg/mol"
    print *, "Función de partición translacional:", ParticionTranslacional(masa)

END PROGRAM PropiedadesTermodinamicas
