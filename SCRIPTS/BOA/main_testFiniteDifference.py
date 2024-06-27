import sys
import numpy as np

from   FELiCS.Parameters.parameters import parameters


from calcBaseFlow                import calculateBaseFlow

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

## Read the given parameters. They have to fit the case in the GeometryDeformer
#a_i  = np.load("params.npy")

#-----------------------------------------------------------------------
# MAIN PART - TEST FINITE DIFFERENCE FOR PROFILE CASE 
#-----------------------------------------------------------------------

# Set control points
a_0 = [ -0.00242300886715957411, \
        -0.00042898338465717500, \
         0.00114762923408253878, \
         0.00239516041194041962, \
         0.00336637871645898909, \
         0.00409472049133493068, \
         0.00460176647276058781, \
         0.00490103870112774259]

baseFlow_array                                  = calculateBaseFlow            (settingsFileName_baseFlow)
baseFlowSensitivity_array                       = calculateBaseFlowSensitivity (settingsFileName_baseFlow, baseFlow_array, a_0) 
sensitivity1, sensitivity2, eigenValue_0        = calculateModesAndSensitivity (settingsFileName_eigenSystem, baseFlow_array, baseFlowSensitivity_array, a_0)

# write out solutions s.t. the shapedescender can read them
np.save("f.npy",  eigenValue)
np.save("df.npy", (sensitivity1 + sensitivity2))
epsilon = 1.e-4

a_i = a_0
FD  = np.empty(len(a_0),dtype=complex)
for i in range(len(a_0)):
    a_i[i] = a_i[i] + epsilon 
    baseFlow_array                            = calculateBaseFlow            (settingsFileName_baseFlow)
    baseFlowSensitivity_array                 = calculateBaseFlowSensitivity (settingsFileName_baseFlow, baseFlow_array, a_i) 
    sensitivity1, sensitivity2, eigenValue_i  = calculateModesAndSensitivity (settingsFileName_eigenSystem, baseFlow_array, baseFlowSensitivity_array, a_i)
    FD[i] = (eigenValue_i - eigenValue_0)/epsilon
    a_i[i] = a_i[i] - epsilon 

np.save("FD_-4.npy", FD)





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

