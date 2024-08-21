import sys
import numpy as np

from   FELiCS.Parameters.parameters import parameters


from calcBaseFlow                import calculateBaseFlow
from calcMeanFlow                import calculateMeanFlow

from calcModesAndSensitivity     import calculateModesAndSensitivity
from calcBaseFlowSensitivity     import calculateBaseFlowSensitivity


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
#baseFlow_array                                  = calculateMeanFlow            (settingsFileName_baseFlow)
print('~~~~~ base flow sensitivity ~~~~~~~')
baseFlowSensitivity_array                       = calculateBaseFlowSensitivity (settingsFileName_baseFlow, baseFlow_array, a_i) 
print('~~~~~ modes and sensitivity ~~~~~~~')
sensitivity1, sensitivity2, eigenValue          = calculateModesAndSensitivity (settingsFileName_eigenSystem, baseFlow_array, baseFlowSensitivity_array, a_i)



# write out solutions s.t. the shapedescender can read them
np.save("f.npy",  eigenValue)
np.save("df.npy", (sensitivity1 + sensitivity2))

print("##################", eigenValue)
print("##################", (sensitivity1)) 
print("##################", (sensitivity2) )
print("##################", (sensitivity1 + sensitivity2))



##-----------------------------------------------------------------------
### MAIN PART - OLD VERSION WITH SOLVING THE BASE FLOW EQUATION 
##-----------------------------------------------------------------------
#from calcSensitivityFromModes    import calculateSensitivityFromModes
#from calcSensitivityFromBaseFlow import calculateSensitivityFromBaseFlow
#baseFlow                                        = calculateBaseFlow                (settingsFileName_baseFlow)
#
#baseFlow_array = baseFlow.getCoefficientArray()
#
#
#sensitivity1, eigenValue, rhs_adjoint_baseFlow  = calculateSensitivityFromModes    (settingsFileName_eigenSystem, baseFlow_array)
#
#sensitivity2                                    = calculateSensitivityFromBaseFlow (settingsFileName_baseFlow,    baseFlow_array, rhs_adjoint_baseFlow) 
# 
#
## write out solutions s.t. the shapedescender can read them
#np.save("f.npy",  eigenValue)
#np.save("df.npy", (sensitivity1 + sensitivity2))
#
#
#print("##################", eigenValue)
#print("##################", (sensitivity1 + sensitivity2))

