import sys
import numpy as np

from   FELiCS.Parameters.parameters    import parameters


from calcBaseFlow                      import calculateBaseFlow

from calcResolventModesAndSensitivity  import calculateResolventModesAndSensitivity
from calcBaseFlowSensitivity           import calculateBaseFlowSensitivity


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

baseFlow_array                          = calculateBaseFlow            (settingsFileName_baseFlow)
baseFlowSensitivity_array               = calculateBaseFlowSensitivity (settingsFileName_baseFlow, baseFlow_array, a_0) 
sensitivity1, sensitivity2, gain_0      = calculateResolventModesAndSensitivity (settingsFileName_eigenSystem, baseFlow_array, baseFlowSensitivity_array, a_0)

print("################## Resolvent gain:", gain_0)

np.save("OutputTest/f.npy",  gain_0)
np.save("OutputTest/df.npy", (sensitivity1 + sensitivity2))


for j in range(5,12):

    epsilon = 10**(-j)
    a_i = a_0
    FD  = np.empty(len(a_0),dtype=complex)
    for i in range(len(a_0)):
        a_i[i] = a_i[i] + epsilon 
        baseFlow_array                          = calculateBaseFlow            (settingsFileName_baseFlow, a_i, deformed = True)
        baseFlowSensitivity_array               = calculateBaseFlowSensitivity (settingsFileName_baseFlow, baseFlow_array, a_i, deformed = True) 
        sensitivity1_i, sensitivity2_i, gain_i  = calculateResolventModesAndSensitivity (settingsFileName_eigenSystem, baseFlow_array, baseFlowSensitivity_array, a_i, deformed = True)
        FD[i] = (gain_i - gain_0)/epsilon
        a_i[i] = a_i[i] - epsilon  
        print("################## Resolvent gain:", gain_0)
        print("################## Sensitivities:")
        print(i, np.real(sensitivity1[i]+sensitivity2[i]), np.imag(sensitivity1[i]+sensitivity2[i]))
        print('################## Finite Difference, epsilon= ', epsilon)
        print(i, np.real(FD[i]), np.imag(FD[i]))

    
    np.save("OutputTest/FD_-"+str(j)+".npy", FD)
    
    
    print("################## Resolvent Gain:", gain_0)
    print("################## Sensitivities:")
    for i in range(len(a_0)):
        print(i, np.real(sensitivity1[i]+sensitivity2[i]), np.imag(sensitivity1[i]+sensitivity2[i]))
    print('################## Finite Difference, epsilon= ', epsilon)
    for i in range(len(a_0)):
        print(i, np.real(FD[i]), np.imag(FD[i]))


