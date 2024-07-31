import sys
import numpy as np

from   FELiCS.Parameters.parameters import parameters


from calcBaseFlow                      import calculateBaseFlow

from calcModesAndSensitivity           import calculateModesAndSensitivity
from calcNSResiduum                    import calculateNavierStokesResiduum


#-----------------------------------------------------------------------
## READ PARAMETER FILES AND PARAMETER VALUES FOR OPTIMIZER
#-----------------------------------------------------------------------
# Replace spurious double quotation marks which occur in MobaXterm when executing shell scripts
for i_argument,argument in enumerate(sys.argv):
        sys.argv[i_argument] = argument.replace(chr(8221), '')

# The first parameter file is the one for the base flow calculations, the second parameter file is for the eigensystem.
settingsFileName_baseFlow    = sys.argv[2]
settingsFileName_eigenSystem = sys.argv[3]


# Read the given parameters. They have to fit the case in the GeometryDeformer
a_i  = np.load("params.npy")

#-----------------------------------------------------------------------
# MAIN PART 
#-----------------------------------------------------------------------
print('~~~~~ base flow ~~~~~~~')
baseFlow_array                                  = calculateBaseFlow            (settingsFileName_baseFlow)
print('~~~~~ base flow sensitivity ~~~~~~~')
residuum                                        = calculateNavierStokesResiduum (settingsFileName_baseFlow, baseFlow_array, a_i) 


# write out solutions s.t. the shapedescender can read them
np.save("F.npy",  residuum)

print("residuum: ", residuum)




