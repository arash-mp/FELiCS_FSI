#-----------------------------------------------------------------------
from 	FELiCS.Parameters.config	        import 	config
from 	FELiCS.Misc.logging			        import  Logger
from    FELiCS.IO.Reader                    import  Reader
from    FELiCS.Fields.Mode                  import  Mode
from    FELiCS.IO.Writer                    import  Writer
from    FELiCS.SpaceDisc.FEMSpaces          import  FEMSpaces
from    FELiCS.Fields.ModeCollection        import  ModeCollection

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

