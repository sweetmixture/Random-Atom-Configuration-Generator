# Random-Atom-Configuration-Generator

A python code to generate randomised atom configuration of binary materials (e.g., MgO, Al2O3 ... ) output will be written in standard 'xyz' format.

## how to run the code?

use the following command
```bash
python /path/to/(MAIN_RACG.py) 
```
Beware that the directory where "MAIN_RACG.py" is going to be run must include input file "RACG.in"


- Contents in the inputfile: "RACG.in"

Example "RACG.in" for (Bi2O3)n (n=8) randomised nanocluster generator
(lines start with '#' are comments)

```bash
  1	#################################################
  2	# GENERAL INFO
  3	#################################################
  4	BOX_SIZE:		12.0	6.0	6.0
  5	BASE_ELEMENT:		Bi	O
  6	UNIT_STOCHIO:		2	3	
  7	ATOM_CHARGE:		+3.	-2.
  8	UNIT_SIZE:		8
  9	#################################################
 10	# TOLERACES FOR GENERATING CONFIGURATIONS
 11	#################################################
 12	INTER_ATOMIC_DISTANCE_TOL(Angs):	2.0
 13	INTER_ATOMIC_FORCE_TOL(eV/Angs):	20.
 14	#################################################
 15	# OPTION WORDS IF USE PRE-OPTIMISATION
 16	#
 17	# BORN-MAYER POTENTIAL BETWEEN CATION - ANION
 18	# IS USED TO CIRCUMVENT ATOM CLASHES
 19	# 
 20	# BY TWEAKING "GNORM_TOL" USER CAN KEEP A
 21	# RANDOM CONFIG IN LOOSELY RELAXED STRUCTURE
 22	#################################################
 23	BM_A(eV):		12435.
 24	BM_RHO(/Angs):		0.257
 25	LJ_(Angs^6):		32.00
 26	GNORM_TOL(eV/Angs):	0.000050
 27	
 28	
 29	# Dummy Comment
```

## Essential Arguments 

 1. BOX_SIZE: X Y Z
 
 A size of container or a box (x/y/z) where atoms randomly go in (lengths are in Angstrom unit)

 2. BASE_ELEMENT: A B

 Species A and B, (only supports binary system in this version but feel free to edit the source code).

 3. UNIT_STOCHIO: a b

 Stochiometry of a building unit, e.g., for the case of (Bi2O3)n, a = 2 and b = 3.

 4. ATOM_CHARGE: qa qb

 User can feed formal charges of atoms in a building unit (in the sample above, Bi = +3 and O = -2, used). 

 5. UNIT_SIZE: n

 Size of an atom configuration or cluster, e.g., for (Bi2O3)8, n = 8.

 - Optional arguements

 6. INTER_ATOMIC_DISTANCE_TOL(Angs):

 * During the runtime of the program (while try to randomly put an atom in the container),
   it will meausre interatomic distances between the newly put atom and the existing others in the system, 
   if one of the distances is within tolerance the trial move will be rejected (else accepted).
 
   (There is a default tolerance for the distance (1.5 Angs), but user can feed it - see line 12 in the sample above)

 7. INTER_ATOMIC_FORCE_TOL(eV/Angs):

 * Here, especially feeding the charges of atoms in (4) would be useful. 
   During runtime, program will measure interatomic force (based on the coulomb law) acting on the newply put atom, 
   if the force is within tolerance the trial move will be rejected (else accetped).    
   The situation may happen when a collection of anions or cations are in a pack local region in the cotainer.

   (There is a default tolerance for the force (10.0 eV/Angs), but user can feed it - see line 13 in the sample above)

 8. Parameters for short-range interaction & pre-optimiser:

    The parameters are generally in effect if user wants to pre-optimised the randomised structure. The pre-optimisation step would be useful when a user wants to 
    use the radomised configurations for a post-process with high level of theory. 
    (i.e., if the interatomic distances of atoms in the random configuration are too far each other, the geometry optimistion step in the post-process will take a bit of hours)
    
    (1) BM_A(eV):

    Born-Mayer 'A' parameter (default value with 3000).
 
    (2) BM_RHO(/Angs):

    Born-Mayer 'Rho' parameter (default value with 0.30).

    (3) LJ_(Angs^6):

    Attractive term in Lennard-Jones (default value with 36, if you put negative sign ... the interaction will turns intp repulsive).

    (4) GNORM_TOL(eV/Angs):

    Gnorm tolerance for terminating a pre-geometric-optimisation process (default value with 0.035).

 * Born-Mayer potential is used between cation and anion, the parameters (BM_A & BM_RHO) are used to describe short-range repulsive interactions 
   to circumvent collision of atoms during pre-geometric-optimisation (e.g., anion - cation collision by attractive coulomb interaction).
 
 * Also attractive term in Lennard-Jones potential is used to avoid a randomised cluster in a sparsed configuration after pre-geometric-optimisation.
   
 * Note that a user can custom all those values. For intance, if you want a pure randomised configuration, then can specify a large value in "GNORM_TOL(eV/Angs):" (e.g., 5.0),
   then the initial random configuration will be written as a result (since the gnorm value of the initial configration is already less then the tolerance).
