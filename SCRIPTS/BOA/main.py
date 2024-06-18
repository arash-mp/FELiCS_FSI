import sys
import numpy as np

from   FELiCS.Parameters.parameters import parameters


from calcBaseFlow                import calculateBaseFlow
from calcSensitivityFromModes    import calculateSensitivityFromModes
from calcSensitivityFromBaseFlow import calculateSensitivityFromBaseFlow

from calcModesAndSensitivity     import calculateModesAndSensitivity
from calcBaseFlowSensitivity     import calculateBaseFlowSensitivity


#-----------------------------------------------------------------------
## READ PARAMETER FILES
#-----------------------------------------------------------------------
# Replace spurious double quotation marks which occur in MobaXterm when executing shell scripts
for i_argument,argument in enumerate(sys.argv):
        sys.argv[i_argument] = argument.replace(chr(8221), '')

# The first parameter file is the one for the base flow calculations, the second parameter file is for the eigensystem.
settingsFileName_baseFlow    = sys.argv[2]
settingsFileName_eigenSystem = sys.argv[3]


##-----------------------------------------------------------------------
### MAIN PART 
##-----------------------------------------------------------------------
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

#-----------------------------------------------------------------------
# MAIN PART VERSION 2 
#-----------------------------------------------------------------------
baseFlow                                        = calculateBaseFlow            (settingsFileName_baseFlow)
baseFlow_array = baseFlow.getCoefficientArray()

baseFlowSensitivity                             = calculateBaseFlowSensitivity (settingsFileName_baseFlow,    baseFlow_array) 

baseFlowSensitivity_array=[]
for sens in baseFlowSensitivity:
    baseFlowSensitivity_array.append(sens.getCoefficientArray())

sensitivity1, sensitivity2, eigenValue          = calculateModesAndSensitivity (settingsFileName_eigenSystem, baseFlow_array, baseFlowSensitivity_array)



# write out solutions s.t. the shapedescender can read them
np.save("f.npy",  eigenValue)
np.save("df.npy", (sensitivity1 + sensitivity2))

print("##################", eigenValue)
print("##################", (sensitivity1)) 
print("##################", (sensitivity2) )
print("##################", (sensitivity1 + sensitivity2))

