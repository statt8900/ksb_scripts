# traj -> x -> y -> z -> adsorber coordinates ->  feff with no absorber defined

header = """EDGE      L3
S02       1

* SETEDGE

*i         pot    xsph  fms   paths genfmt ff2chi
CONTROL   1      1     1     1     1      1
PRINT     1      0     0     0     0      3

EXCHANGE  0 0 1
* CORRECTIONS   1.0

SCF       5.5 0 25

XANES     6 0.05 0.3

FMS       7.5 0

* LDOS      -30   20     0.1

* POLARIZATION  0   0   0

COREHOLE none

* The following lines are needed to account for the many-pole self-energy
OPCONS
MPSE 2
SFCONV

* This card accounts for dynamic screening of the x-ray and core-hole fields
TDLDA 1

* This card calculates the Debye-Waller factor
* The parameters are [temperature] [Debye temperature] and [idwopt -- i.e., model]
DEBYE  190 315


 POTENTIALS
 *   ipot   Z  element        l_scmt  l_fms   stoichiometry
     0     79     Au           3       3       0.001
     1     79     Au           3       3       2
     2     46     Pd           3       3       2
     3     1      H            1       1       1


 ATOMS                  * this list contains 293 atoms
 *   x          y          z     ipot tag           distance"""