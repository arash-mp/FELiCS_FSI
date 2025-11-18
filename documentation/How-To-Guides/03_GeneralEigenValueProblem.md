# How to solve a general eigenvalue problem

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

**Step 4: Define FEM spaces and read in the baseflow**

We provide a baseflow file, that already contains all the necessary information of the flow field. This baseflow file gets imported and gets stored in the FELiCS meanFlowClass, which contains all the variables of the flow field. However, to import it we need to define the necassary finite element function space via the <code style="color : Cyan">FEMSpaces</code> class.


```python
# Get FEMSpaces
spaces   = FEMSpaces(
    param,
    mesh,
)

# Initialize mean flow class & import from file
meanFlow    = meanFlowClass(
    param, 
    spaces, 
    mesh
)
meanFlow.importDataFromFileAndExportToH5()
```

**Step 5. Export of flow fields and mesh.**

The meanflow and the mesh need to be saved as <code style="color : Darkorange">.hdf5</code> so that the mesh can be utilized again while saving the eigenmodes in <code style="color : Darkorange">.xmf</code> and <code style="color : Darkorange">.hdf5</code> formats.


```python
baseFlow_array = meanFlow._fieldDict['u'].function.x.array[:]

# Load base flow into meanFlow object
if np.linalg.norm(meanFlow._fieldDict['u'].function.x.array[:]) < 1.e-8:
    baseFlow    = Field(spaces.VMixed, mesh)
    baseFlow.setCoefficientArray(baseFlow_array)
    [u,p]       = baseFlow.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u.function

# Export mean flow in "h5" file
meanflowFilename = 'meanflow.h5'
meanFlow.mapToExportMeshAndExport(spaces, meanflowFilename)
```

**Step 6: Define the Equations for the linear problem**

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
eigenValue  = solution.getLeadingMode().getEigenValue()
print(f"Leading eigenvalue: {str(eigenValue)}")
```

After solving the problem we can get the leading mode and its respective eigenvalue from the <code style="color : Cyan">ModeCollection</code> class. The eigenvalue printed should be approximatley 0.74+0.013j.
