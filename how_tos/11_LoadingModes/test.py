#-----------------------------------------------------------------------
from 	FELiCS.Parameters.config	        import 	config
from 	FELiCS.Misc.logging			        import  Logger
from    FELiCS.IO.Reader                    import  Reader
from    FELiCS.Fields.Field                 import  Field
from    FELiCS.IO.Writer                    import  Writer
from    FELiCS.SpaceDisc.FEMSpaces          import  FEMSpaces
from    FELiCS.Fields.ModeCollection        import  ModeCollection
import  ufl
from    FELiCS.Misc.tensorUtils             import Tensor, iDot, iConj

# Get the logger
logger      = Logger(True, False, "felics")
logger      = Logger.get_logger("felics")

# Get the config of the case
param       = config()
param.importFromFile("Input/config.json")

#-----------------------------------------------------------------------
## INITIALIZATION
#-----------------------------------------------------------------------
# Get the same FELiCS objects as the ones used in the analysis

# Mesh
mesh        = param.getMesh()

# FEMSpaces
FEMSpaces   = FEMSpaces(param, mesh)

# Create Reader object
reader      = Reader(
    sourceDir="Input",
    needInterpolation=False,
)

# for the Writer object, the default export directory is "Output"
writer      = Writer()

#-----------------------------------------------------------------------
## Loading modes from file
#-----------------------------------------------------------------------
# Create empty mode collection to import into
importSolution  = ModeCollection(
    FEMSpaces.VMixed, 
    mesh
)

# Using the import this way loads all modes in the directory
importSolution.importData(
    reader,
    "Input",
)

# Print the contents of the imported solution
importSolution.describe()

#-----------------------------------------------------------------------
## Plot the streamwise velocity of each mode
#-----------------------------------------------------------------------
# Separate direct and adjoint modes
directMode      = [mode for mode in importSolution.modeList if mode.modeType.name == "DIRECT"][0]
adjointMode     = [mode for mode in importSolution.modeList if mode.modeType.name == "ADJOINT"][0]
# We can check the mode type with directMode.modeType.name

# Keep only the streamwise velocity component 
# The velocity vector is assumed to be the first subfield of the state vector
# and the streamwise component is assumed to be the first subfield of the velocity vector
uxFieldDirect   = directMode.getListOfSubFields()[0].getListOfSubFields()[0]
uxFieldAdjoint  = adjointMode.getListOfSubFields()[0].getListOfSubFields()[0]

# Plot the streamwise velocity of the direct mode
uxFieldDirect.plot(xlim=(-5, 20), ylim=(0, 5))

# Plot the streamwise velocity of the adjoint mode
uxFieldAdjoint.plot(xlim=(-20, 5), ylim=(0, 5))


#-----------------------------------------------------------------------
## Computing the structural sensitivity
#-----------------------------------------------------------------------
# Get the velocity fluctuations
uFieldDirect            = directMode.getListOfSubFields()[0]
uFieldAdjoint           = adjointMode.getListOfSubFields()[0]

# Test space for the structural sensitivity (scalar field) and a test function
V_scalar                = FEMSpaces.P2
v                       = ufl.TestFunction(V_scalar)

# Create empty field for the structural sensitivity
sS                      = Field(V_scalar, mesh)
sS.name                 = "StructuralSensitivity"

# Define the expression for the structural sensitivity
J_hat                   = sS.mesh.coordinateSystem.J_hat
v_tens                  = Tensor(v, CoordSys=sS.mesh.coordinateSystem, m = sS.m, mayHaveSpectralDimension=sS.hasSpectralDimension)
exprSS                  = (
        iDot(uFieldDirect.getTensor(), iConj(uFieldDirect.getTensor()))**0.5
        *
        iDot(uFieldAdjoint.getTensor(), iConj(uFieldAdjoint.getTensor()))**0.5
        * iConj(v_tens)
        ).ufl_tens * J_hat* ufl.dx
sS.evaluateUflTensorExpression(exprSS)

# Plot the structural sensitivity
sS.plot(xlim=(-2, 6), ylim=(0, 4))
