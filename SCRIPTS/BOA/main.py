#-----------------------------------------------------------------------
## READ PARAMETER FILES
#-----------------------------------------------------------------------
import sys
from   FELiCS.Parameters.parameters import parameters

# Replace spurious double quotation marks which occur in MobaXterm when executing shell scripts
for i_argument,argument in enumerate(sys.argv):
        sys.argv[i_argument] = argument.replace(chr(8221), '')

# The first parameter file is the one for the base flow calculations, the second parameter file is for the eigensystem.
settingsFileName_baseFlow    = sys.argv[2]
settingsFileName_eigenSystem = sys.argv[3]


#-----------------------------------------------------------------------
## MAIN PART 
#-----------------------------------------------------------------------
from calcBaseFlow                import calculateBaseFlow
from calcSensitivityFromModes    import calculateSensitivityFromModes
from calcSensitivityFromBaseFlow import calculateSensitivityFromBaseFlow


baseFlow                                        = calculateBaseFlow                (settingsFileName_baseFlow)

baseFlow_array = baseFlow.getCoefficientArray()


sensitivity1, eigenValue, rhs_adjoint_baseFlow  = calculateSensitivityFromModes    (settingsFileName_eigenSystem, baseFlow_array)

sensitivity2                                    = calculateSensitivityFromBaseFlow (settingsFileName_baseFlow,    baseFlow_array, rhs_adjoint_baseFlow) 


# write out solutions s.t. the shapedescender can read them
np.save("f.npy",  eigenValue)
np.save("df.npy", sensitivity1 + sensitivity2)

