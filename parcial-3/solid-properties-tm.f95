program main
    implicit none
    ! Declare constants
    real(8) :: hbar, kB, me, eV
  
    ! Declare variables
    real(8) :: V, rho, T, Ef_J, Vf, mu_J, E_J, Ef_eV, mu_eV, E_eV
    character(len=20) :: bravais_cell
    real(8) :: valence_electrons
    integer :: it1, it2
  
    ! Constants
    hbar = 1.054571817d-34  ! Planck's constant over 2 pi in J·s
    kB = 1.380649d-23       ! Boltzmann constant in J/K
    me = 9.1093837d-31      ! Mass of electron in kg
    eV = 1.60218d-19        ! 1 eV in joules
  
    ! Read inputs
    print *, "Enter the type of Bravais lattice cell:"
    read *, bravais_cell
    print *, "Enter the number of valence electrons Z:"
    read *, valence_electrons
  
    ! Call the function with user inputs
    call calculate_volume_and_density(bravais_cell, valence_electrons, V, rho)
    print *, "Enter the temperature in K:"
    read *, T
    print *, ''
  
    ! Call functions to calculate Fermi energy and velocity
    Ef_J = calculate_fermi_energy(rho)  ! Fermi energy in joules
    Vf = calculate_fermi_velocity(rho)  ! Fermi velocity in m/s
    ! Call fixed-point iteration to calculate chemical potential and internal energy
    call calculate_energy(Ef_J, T, E_J, it1)
    call calculate_chemical_potential(Ef_J, T, mu_J, it2)
  
    ! Convert energies from joules to eV
    Ef_eV = Ef_J / eV
    mu_eV = mu_J / eV
    E_eV = E_J / eV
  
    ! Print results
    print *, '-------------------------------------------------------------------------'
    print *, "Volume:", V, "m^3"
    print *, "Electron Density:", rho, "electrons/m^3"
    print *, "Fermi Energy:", Ef_eV, "eV"
    print *, "Fermi Velocity:", Vf, "m/s"
    print *, "Chemical Potential at T =", T, "K:", mu_eV, "eV"
    print *, "Energy at T =", T, "K:", E_eV, "eV"
    print *, '-------------------------------------------------------------------------'
  
  contains
  
  subroutine calculate_volume_and_density(bravais, Z, V, rho)
    implicit none
    character(len=20), intent(in) :: bravais
    real(8), intent(in) :: Z
    real(8), intent(out) :: V, rho
    real(8) :: a, b, c, alpha, beta, gamma
    real(8) :: a_meters, b_meters, c_meters, alpha_radians, beta_radians, gamma_radians
    integer :: N
  
    print *, 'Enter the lattice parameter a (in Ångströms):'
    read *, a
    a_meters = a * 1.0d-10  ! Convert a from Ångströms to meters
  
    select case (trim(bravais))
      case ('cubic', 'fcc', 'bcc')
        if (bravais == 'cubic') then
          N = 1
        elseif (bravais == 'fcc') then
          N = 4
        elseif (bravais == 'bcc') then
          N = 2
        endif
        V = a_meters**3
      case ('tetra', 'bc-tetra')
        print *, 'Enter the lattice parameter c (in Ångströms):'
        read *, c
        c_meters = c * 1.0d-10
        if (bravais == 'tetra') then
          N = 1
          V = a_meters**2 * c_meters
        else if (bravais == 'bc-tetra') then
          N = 2
          V = a_meters**2 * c_meters / 2.0
        endif
      case ('ortho', 'bc-ortho', 'fc-ortho')
        print *, 'Enter the lattice parameter b (in Ångströms):'
        read *, b
        print *, 'Enter the lattice parameter c (in Ångströms):'
        read *, c
        b_meters = b * 1.0d-10
        c_meters = c * 1.0d-10
        if (bravais == 'ortho') then
          N = 1
          V = a_meters * b_meters * c_meters
        else if (bravais == 'bc-ortho') then
          N = 2
          V = a_meters * b_meters * c_meters / 2.0
        else if (bravais == 'fc-ortho') then
          N = 4
          V = a_meters * b_meters * c_meters / 4.0
        endif
      case ('hex')
        print *, 'Enter the lattice parameter c (in Ångströms):'
        read *, c
        c_meters = c * 1.0d-10
        N = 1
        V = (3.0 * sqrt(3.0) / 2.0) * a_meters**2 * c_meters
      case ('rhombo')
        print *, 'Enter the angle alpha (in degrees):'
        read *, alpha
        alpha_radians = alpha * acos(-1.0d0) / 180.0d0
        N = 1
        V = a_meters**3 * sqrt(1.0 - 3.0 * cos(alpha_radians)**2 + 2.0 * cos(alpha_radians)**3)
      case ('mono', 'bc-mono')
        print *, 'Enter the lattice parameter b (in Ångströms):'
        read *, b
        print *, 'Enter the lattice parameter c (in Ångströms):'
        read *, c
        print *, 'Enter the angle beta (in degrees):'
        read *, beta
        b_meters = b * 1.0d-10
        c_meters = c * 1.0d-10
        beta_radians = beta * acos(-1.0d0) / 180.0d0
        if (bravais == 'mono') then
          N = 1
          V = a_meters * b_meters * c_meters * sin(beta_radians)
        else if (bravais == 'bc-mono') then
          N = 2
          V = a_meters * b_meters * c_meters * sin(beta_radians) / 2.0
        endif
      case ('tri')
        print *, 'Enter the lattice parameter b (in Ångströms):'
        read *, b
        print *, 'Enter the lattice parameter c (in Ångströms):'
        read *, c
        print *, 'Enter the angle alpha (in degrees):'
        read *, alpha
        print *, 'Enter the angle beta (in degrees):'
        read *, beta
        print *, 'Enter the angle gamma (in degrees):'
        read *, gamma
        b_meters = b * 1.0d-10
        c_meters = c * 1.0d-10
        alpha_radians = alpha * acos(-1.0d0) / 180.0d0
        beta_radians = beta * acos(-1.0d0) / 180.0d0
        gamma_radians = gamma * acos(-1.0d0) / 180.0d0
        N = 1
        V = a_meters * b_meters * c_meters * sqrt(1.0 - cos(alpha_radians)**2 - cos(beta_radians)**2 - cos(gamma_radians)**2 + 2.0 &
             * cos(alpha_radians) * cos(beta_radians) * cos(gamma_radians))
      case default
        print *, "Unsupported Bravais lattice type"
        stop
    end select
  
    rho = N * Z / V  ! Electron density in electrons/m^3
  end subroutine calculate_volume_and_density
  
    ! Function to calculate Fermi energy
    real(8) function calculate_fermi_energy(rho)
      implicit none
      real(8), intent(in) :: rho
      ! Declare constants
      real(8) :: hbar, me
      ! Define constants
      hbar = 1.054571817d-34
      me = 9.1093837d-31
      ! Calculate Fermi energy
      calculate_fermi_energy = hbar**2 / (2.0 * me) * (3.0 * acos(-1.0d0)**2 * rho)**(2.0/3.0)
    end function calculate_fermi_energy
  
    ! Function to calculate Fermi velocity
    real(8) function calculate_fermi_velocity(rho)
      implicit none
      real(8), intent(in) :: rho
      ! Declare constants
      real(8) :: hbar, me
      ! Define constants
      hbar = 1.054571817d-34
      me = 9.1093837d-31
      ! Calculate Fermi velocity
      calculate_fermi_velocity = hbar / me * (3.0 * acos(-1.0d0)**2 * rho)**(1.0/3.0)
    end function calculate_fermi_velocity
  
    ! Subroutine to calculate chemical potential
    subroutine calculate_chemical_potential(x0, T, mu, iterations)
      implicit none
      ! Declare variables
      real(8), intent(in) :: x0, T
      real(8), intent(out) :: mu
      integer, intent(out) :: iterations
      ! Declare local variables
      real(8) :: x_current, tolerance, mu_F
  
      ! Constants
      real(8), parameter :: pi = acos(-1.0d0)
      real(8), parameter :: kB = 1.380649d-23  ! Boltzmann constant in J/K
  
      ! Initialize variables
      x_current = x0
      tolerance = 1.0d-8
      iterations = 0
      mu_F = x0  ! Assuming initial guess is Fermi energy
  
      ! Perform fixed-point iteration
      do
        mu = mu_F * (1.0d0 - (pi**2 * kB**2 * T**2) / (12.0d0 * x_current**2))
        iterations = iterations + 1
        ! Check convergence
        if (abs(mu - x_current) < tolerance .or. iterations > 100) exit
        x_current = mu
      end do
  
      ! Handle non-convergence
      if (iterations > 100) then
        print *, "Did not converge within 100 iterations"
        stop
      end if
    end subroutine calculate_chemical_potential
  
    ! Subroutine to calculate energy
    subroutine calculate_energy(x0, T, E, iterations)
      implicit none
      ! Declare variables
      real(8), intent(in) :: x0, T
      real(8), intent(out) :: E
      integer, intent(out) :: iterations
      ! Declare local variables
      real(8) :: x_current, tolerance, E_F
  
      ! Constants
      real(8), parameter :: pi = acos(-1.0d0)
      real(8), parameter :: kB = 1.380649d-23  ! Boltzmann constant in J/K
  
      ! Initialize variables
      x_current = x0
      tolerance = 1.0d-8
      iterations = 0
      E_F = x0  ! Assuming initial guess is Fermi energy
  
      ! Perform fixed-point iteration
      do
        E = E_F * (1.0d0 + (5.0d0 * pi**2 * kB**2 * T**2) / (12.0d0 * x_current**2))
        iterations = iterations + 1
        ! Check convergence
        if (abs(E - x_current) < tolerance .or. iterations > 100) exit
        x_current = E
      end do
  
      ! Handle non-convergence
      if (iterations > 100) then
        print *, "Did not converge within 100 iterations"
        stop
      end if
    end subroutine calculate_energy

  end program main
