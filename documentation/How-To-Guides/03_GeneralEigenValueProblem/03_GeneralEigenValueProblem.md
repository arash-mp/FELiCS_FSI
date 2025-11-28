# FELiCS scripts explained: Modal analysis

**Step 1. Import of Python packages.**


```python
import numpy as np
```

**Step 2. Import of FELiCS packages.**


```python
from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
from    FELiCS.Fields.meanFlowClass         import meanFlowClass
from    FELiCS.Fields.Field                 import Field
from    FELiCS.Parameters.config            import config 
from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
from    FELiCS.Solvers.LinearSolver         import LinearSolver
from    FELiCS.Fields.ModeCollection        import ModeCollection
from    FELiCS.IO.Writer                    import Writer
```

**Step 3: Read the settings from the setting file**

In this tutorial the Modal Analysis has been done for base flow with $Re = 50$. Hence, the flow field is imported through the file <code style="color : Darkorange">base_flow_for_FELiCS.fel</code>, with its path also defined in the settings file <code style="color : Darkorange">modal.json</code>. 


```python
# Read parameters
settingsFileName = "modal.json"
param            = config()
param.importFromFile(settingsFileName)

mesh             = param.getMesh()
```

    Info     | config.py              | importFromFile             (line 205 ) : Loading configuration from modal.json
    Warning  | config.py              | importFromFile             (line 224 ) : Field "needInterpolation" missing from file, setting default: True
    Warning  | config.py              | importFromFile             (line 224 ) : Field "MolViscPerturbModel" missing from file, setting default: {'type': 'Constant', 'Constants': {'Viscosity': 1.0}}
    Warning  | config.py              | importFromFile             (line 224 ) : Field "PrandtlNumber" missing from file, setting default: 0.72
    Warning  | config.py              | importFromFile             (line 224 ) : Field "ForcingCoeff" missing from file, setting default: []
    Warning  | config.py              | importFromFile             (line 224 ) : Field "ForcingNorm" missing from file, setting default: TKE
    Warning  | config.py              | importFromFile             (line 224 ) : Field "ResponseNorm" missing from file, setting default: TKE
    Info     | logging.py             | change_log_location        (line 286 ) : Log files moved to: output_dir/log
    Info     | MixtureClass.py        | __init__                   (line 75  ) : No mixture file Mixture.mix, using defaults.
    Info     | FELiCSMesh.py          | __init__                   (line 94  ) : Opening mesh file: cylinder_wake.msh
    Info     | FELiCSMesh.py          | __init__                   (line 102 ) : Mesh contains 3613 nodes and 7224 elements
    Info     | config.py              | importFromFile             (line 267 ) : Configuration loaded successfully


**Step 4: Define FEM spaces and read in the baseflow**

We provide a baseflow file, that already contains all the necessary information of the flow field. This baseflow file gets imported and gets stored in the FELiCS meanFlowClass, which contains all the variables of the flow field. However, to import it we need to define the necassary finite element function space via the <code style="color : Cyan">FEMSpaces</code> class.


```python
# Get FEMSpaces
spaces   = FEMSpaces(
    param,
    mesh,
)

# Initialize the writer
writer = Writer(mesh, param.Export.ExportFolder)

# Initialize mean flow class & import from file
meanFlow    = meanFlowClass(
    param, 
    spaces, 
    mesh
)
meanFlow.importDataFromFileAndExportToH5(writer)
```

    Info     | FEMSpaces.py           | __init__                   (line 132 ) : Defining FEM-spaces.
    Info     | meanFlowClass.py       | importDataFromFileAndExportToH5 (line 155 ) : Reading input flow from: 'base_flow_for_FELiCS.fel'
    Warning  | Reader.py              | _check_variable_availability_and_type (line 818 ) : No variables for field 'rho' found in file. Set to default values.
    Warning  | Reader.py              | _set_arrays_to_field       (line 769 ) : Variable 'rho' not found in loaded arrays for scalar field. Set to default values.
    Warning  | Reader.py              | _check_variable_availability_and_type (line 818 ) : No variables for field 'spg' found in file. Set to default values.
    Warning  | Reader.py              | _set_arrays_to_field       (line 769 ) : Variable 'spg' not found in loaded arrays for scalar field. Set to default values.


**Step 5: Define the Equations for the linear problem**

Setting up an equation is straightforward in FELiCS. The EquationCollectionClass contains a range of predefined equations like the Navier-Stokes-Equations, that we will solve today. You can find a detailed list and information about the equations [here](./)


```python
# Set up quations
equation    = EquationCollectionClass(
    param,
    spaces,
    meanFlow,
    mesh
)
```

    Info     | EquationCollection.py  | __init__                   (line 137 ) : Initializing the equation collection class.


**Step 7: Define and solve the general Eigenproblem**

The general Eigenproblem is defined by getting the linear operator and the weight matrix from the <code style="color : Cyan">EquationCollectionClass</code>. Mathematically, splitting up the temporal and spatial parts of the solution as $\mathbf{u}' = \hat{\mathbf{u}} e^{-i \omega t}$ and $p' = \hat{p} e^{-i\omega t}$ which yeilds the following generalized Eigenvalue problem (GEVP)

$$
    \mathrm{Re} (-i \omega \hat{\mathbf{u}} + \underbrace{(\mathbf{u}_b \cdot \nabla ) \hat{\mathbf{u}}  + (\hat{\mathbf{u}} \cdot \nabla ) \mathbf{u}_b}_{C \hat{\mathbf{u}}} ) = \underbrace{- \nabla \hat{p}}_{G \hat{p}} + \underbrace{\nabla^2 \hat{\mathbf{u}}}_{D \hat{\mathbf{u}}}
$$

In the matrix form the GEVP can be formulated as

$$ 
\omega\underbrace{\begin{bmatrix}
        -i*\mathrm{Re} & 0 \\
        0     & 0
    \end{bmatrix}}_{B}
    \begin{bmatrix}
        \hat{\mathbf{u}} \\
        \hat{p}
    \end{bmatrix}
    =
    \underbrace{\begin{bmatrix}
        \mathrm{Re}*C+D & 0 \\
        0 & G
    \end{bmatrix}}_{A}
    \begin{bmatrix}
        \hat{\mathbf{u}} \\
        \hat{p}
    \end{bmatrix}
$$

After initializing the <code style="color : Cyan">ModeCollection</code> class, which will store the solution of our modal analysis, we are ready to solve the general Eigenproblem.The general Eigenproblem is solved using the <code style="color : Cyan">LinearSolver</code> class. The solution of the general Eigenproblem may take some time.


```python
# Get matrices for eigenproblem
A       = equation.getLinearOperator(meanFlow)
B       = equation.getWeightMatrix  (meanFlow)

guesses  = param.Numerics.EigenValueGuess
nSol    = param.Numerics.nSolut
# Solve direct eigenproblem for each guess
solution    = ModeCollection(
    spaces.VMixed, 
    mesh
)
for guess in guesses:
    tmp = LinearSolver.solveGeneralEigenproblem(
        A,
        B,
        guess,
        nSol,
    )
    solution.appendSolutionOfEigenProblem(tmp, guess)

# Get leading eigenvalue
eigenValue  = solution.getLeadingMode().eigenValue
print(f"Leading eigenvalue: {str(eigenValue)}")
```

    Leading eigenvalue: (0.744756877836706+0.013262932365236848j)


After solving the problem we can get the leading mode and its respective eigenvalue from the <code style="color : Cyan">ModeCollection</code> class. The eigenvalue printed should be approximatley 0.74+0.013j.
