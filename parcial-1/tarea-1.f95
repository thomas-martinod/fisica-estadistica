program Thermodynamics
  implicit none
  real :: pi, T, P, h, k_B, N_A, R, V, U
  real, dimension(8) :: molar_masses
  real, dimension(8) :: degeneracies ! List of degeneracies (ge1)
  real :: atom_mass, partition_function, element_entropy, element_gibbs
  integer :: i
  character(len=26) :: filename
 
  ! Constants
  pi = 3.14159265359                 ! Pi
  T = 298.0                           ! Standard Temperature (K)
  P = 1.0E+5                          ! Standard pressure (Pa)
  h = 6.62607015E-34                 ! Plank's constant (J s)
  k_B = 1.380649E-23                 ! Boltzmann constant (J / K)
  N_A = 6.02214076E+23                ! Avogadro's Number (1 / mol)
 
  R = k_B * N_A                       ! Ideal gas constant (J / mol K)
  V = k_B * T / P                     ! Standard volume (m^3)
  U = (3.0/2.0) * R * T / 1000.0      ! Standard internal energy U_298 (kJ)
 
  print*, 'Standard Volume: ', V, 'm^3'
  print*, 'Standard internal energy:', U, 'kJ'
 
  molar_masses = [4.00260128, 20.18004638, 26.98153841, 137.32667172, &
                   39.94779856, 9.01218306, 79.90432616, 12.01063556] ! List of molar masses of the elements listed (g / mol)
 
  degeneracies = [1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 1.0]             ! List of degeneracies (ge1)
 
  ! Open file for writing results
  filename = 'thermodynamics-results.csv'
  open(unit=10, file=filename, status='replace', action='write')
 
  ! Write header to the CSV file
  write(10, '(A, A, A)') 'Element', 'S_298 (J/mol K)', 'G_298 (J/mol)'
   
  ! Calculate entropy and Gibbs free energy for each element
  do i = 1, size(molar_masses)
       ! Calculate mass of atom
       atom_mass = mass_of_an_atom(molar_masses(i))
       print*, 'Standard Volume: ', atom_mass
       
       ! Calculate translational partition function
       partition_function = translational_partition(molar_masses(i))
       print*, 'Standard Volume: ', partition_function

       
       ! Calculate entropy
       element_entropy = entropy(partition_function)
       print*, 'Standard Volume: ', element_entropy

       
       ! Calculate Gibbs free energy
       element_gibbs = gibbs(partition_function)
       print*, 'Standard Volume: ', element_gibbs


       ! Write results to the CSV file
       write(10, '(A, F20.10, F20.10)') 'Element ' // trim(adjustl(int2str(i))), element_entropy, element_gibbs
  end do
   
  ! Close the file
  close(unit=10)
 
 contains
 
  ! Function to calculate mass of an atom (kg / atom)
  real function mass_of_an_atom(m)
       real, intent(in) :: m
       mass_of_an_atom = m / N_A / 1000.0
  end function mass_of_an_atom
 
  ! Function to calculate translational partition function
  real function translational_partition(m)
       real, intent(in) :: m
       translational_partition = ((2.0 * pi * m * k_B * T) / (h**2))**(3.0/2.0) * V
  end function translational_partition
 
  ! Function to calculate entropy
  real function entropy(q)
       real, intent(in) :: q
       entropy = N_A * ((3.0/2.0) * k_B + k_B * log(q))
  end function entropy
 
  ! Function to calculate Gibbs free energy
  real function gibbs(q)
       real, intent(in) :: q
       gibbs = k_B * T * log(exp(1.0) / q)
  end function gibbs
 
  ! Function to convert integer to string
  function int2str(i) result(str)
       integer, intent(in) :: i
       character(len=:), allocatable :: str
       character(len=10) :: tmp
       write(tmp, '(I0)') i
       str = trim(adjustl(tmp))
  end function int2str
 
 end program Thermodynamics
 