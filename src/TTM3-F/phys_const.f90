module phys_const
      !module for physical constants

      double precision, parameter :: hbar = 1.d0
      double precision, parameter :: pi = 4.d0*datan(1.d0)	! because 4*atan(1) = pi !
      double precision,parameter :: auang=0.5291772083d0	! convert length from Bohr to Angstrom; 1 Bohr = 0.5291772083 Angstrom
      double precision,parameter :: aucm=219474.6313705d0	! convert energy from Hartree to cm^{-1}
      double precision,parameter :: aukcal=627.509469d0		! convert energy from Hartree to kcal/mol; 1 Hartree = 627.509469 kcal/mol
      double precision,parameter :: auamu=1822.88848d0		! convert mass from amu to me; 1 amu = 1822.88848 me

      !Atomic masses
      double precision,parameter :: Hamu = 1.00782503207d0	! Mass of hydrogen in amu (IUPAC Green Book Third Edition 2007 p.122)
      double precision,parameter :: Damu = 2.0141017778d0	! Mass of deuterium in amu (IUPAC Green Book Third Edition 2007 p.122)
      double precision,parameter :: Oamu = 15.99491461956d0 	! Mass of oxygen-16 in amu (IUPAC Green Book Third Edition 2007 p.122)
end module
