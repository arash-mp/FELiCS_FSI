# FELiCS scripts explained: Modal analysis

**Step 1. Import of Python packages.**


```python
import numpy as np
```

**Step 2. Import of FELiCS packages.**


```python
from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
from    FELiCS.Fields.MeanFlowClass         import MeanFlowClass
from    FELiCS.Fields.Field                 import Field
from    FELiCS.Parameters.Config            import Config 
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
param            = Config()
param.import_from_file(settingsFileName)

mesh             = param.get_mesh()
```

    Info     | Config.py              | import_from_file           (line 245 ) : Loading configuration from modal.json
    Warning  | Config.py              | import_from_file           (line 264 ) : Field "needInterpolation" missing from file, setting default: True
    Warning  | Config.py              | import_from_file           (line 264 ) : Field "MolViscPerturbModel" missing from file, setting default: {'type': 'Constant', 'Constants': {'Viscosity': 1.0}}
    Warning  | Config.py              | import_from_file           (line 264 ) : Field "PrandtlNumber" missing from file, setting default: 0.72
    Warning  | Config.py              | import_from_file           (line 264 ) : Field "ForcingCoeff" missing from file, setting default: []
    Warning  | Config.py              | import_from_file           (line 264 ) : Field "ForcingNorm" missing from file, setting default: TKE
    Warning  | Config.py              | import_from_file           (line 264 ) : Field "ResponseNorm" missing from file, setting default: TKE
    Info     | logging.py             | change_log_location        (line 364 ) : Log files moved to: output_dir/log
    Info     | MixtureClass.py        | __init__                   (line 101 ) : No mixture file Mixture.mix, using defaults.
    Info     | FELiCSMesh.py          | __init__                   (line 124 ) : Opening mesh file: cylinder_wake.msh
    Info     | FELiCSMesh.py          | __init__                   (line 132 ) : Mesh contains 3613 nodes and 7224 elements
    Info     | Config.py              | import_from_file           (line 314 ) : Configuration loaded successfully


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
meanFlow    = MeanFlowClass(
    param, 
    spaces, 
    mesh
)
meanFlow.import_data_from_file_and_export_to_h5(writer)
```

    Info     | FEMSpaces.py           | __init__                   (line 174 ) : Defining FEM-spaces.
    Info     | MeanFlowClass.py       | import_data_from_file_and_export_to_h5 (line 192 ) : Reading input flow from: 'base_flow_for_FELiCS.fel'
    Info     | Reader.py              | _interpolate_to_calc_mesh  (line 706 ) : Linear interpolation took 0.8 seconds. (Mesh size: (14058, 2))
    Warning  | Reader.py              | _check_variable_availability_and_type (line 938 ) : No variables for field 'rho' found in file. Set to default values.
    Warning  | Reader.py              | _set_arrays_to_field       (line 889 ) : Variable 'rho' not found in loaded arrays for scalar field. Set to default values.
    Warning  | Reader.py              | _check_variable_availability_and_type (line 938 ) : No variables for field 'spg' found in file. Set to default values.
    Warning  | Reader.py              | _set_arrays_to_field       (line 889 ) : Variable 'spg' not found in loaded arrays for scalar field. Set to default values.


**Step 5: Define the Equations for the linear problem**

Setting up an equation is straightforward in FELiCS. The EquationCollectionClass contains a range of predefined equations like the Navier-Stokes-Equations, that we will solve today. You can find a detailed list and information about the equations [here](https://felics-d43476.gitlab.io/GoverningEquations/index.html)


```python
# Set up quations
equation    = EquationCollectionClass(
    param,
    spaces,
    meanFlow,
    mesh
)
```

    Info     | EquationCollection.py  | __init__                   (line 165 ) : Initializing the equation collection class.


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
A       = equation.get_linear_operator(meanFlow)
B       = equation.get_weight_matrix  (meanFlow)

guesses  = param.Numerics.EigenValueGuess
nSol    = param.Numerics.nSolut
# Solve direct eigenproblem for each guess
solution    = ModeCollection(
    spaces.VMixed, 
    mesh
)
for guess in guesses:
    tmp = LinearSolver.solve_general_eigenproblem(
        A,
        B,
        guess,
        nSol,
    )
    solution.append_solution_of_eigen_problem(tmp, guess)

# Get leading eigenvalue
eigenValue  = solution.get_leading_mode().eigen_value
print(f"Leading eigenvalue: {str(eigenValue)}")
```

    Leading eigenvalue: (0.744756877836705+0.013262932365236025j)


After solving the problem we can get the leading mode and its respective eigenvalue from the <code style="color : Cyan">ModeCollection</code> class. The eigenvalue printed should be approximatley 0.74+0.013j.
