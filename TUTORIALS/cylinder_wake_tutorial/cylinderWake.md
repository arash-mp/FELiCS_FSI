0. ***Please make sure to switch to a FELiCS branch that has a base-flow solver***

    This tutorial can be found inside the jupyter notebook ```cylinderWake.ipynb```. We highly encourage readers to try it out.

1. Importing all the necassary Python packages:

```
import os
import time
import numpy as np
import matplotlib.pyplot as plt
```

2. Import FEniCSx packages

```
import  dolfinx
from    mpi4py import MPI
from    dolfinx.fem import Function, FunctionSpace
from    ufl import VectorElement, SpatialCoordinate, exp, FiniteElement

```

3. Importing all the necassary FELiCS packages

```
from FELiCS.Parameters.config import config
from FELiCS.SpaceDisc import DefineFEMSpaces
from FELiCS.Fields.meanFlowClass import meanFlowClass
from FELiCS.Equation.EquationCollection import EquationCollectionClass
from FELiCS.Misc.functions import printDebug
from FELiCS.Solvers.LinearSolver import LinearSolver 
from FELiCS.Fields.Field import Field
```

4. Importing all the associated Python scripts

```
import generate_cylinder_mesh
```

5. Generating the mesh file. 

    Here, the idea is create 2D mesh for a flow around a cylinder with different levels of refinement directly via a python script. The flow around the cylinder is described by a large rectangular domain spanning from (-100, 0) to (200, 25), consisting of two smaller rectangular regions, such that the mesh in those regions is refined, which looks something like :  

![Cylinder Wake Mesh](./CylinderWakeMesh.png)

The region of the mesh outlined by a red marker represents the sphere, which upon magnification looks like:

![Region outlined by the red marker, magnified](./magCylinderWake.png)

The length and height of the refined regions can be adjusted by changing the magnitudes values of <code style="color : Cyan">RefineMinX</code>, <code style="color : Cyan">RefineMaxX</code>, <code style="color : Cyan">RefineMinY</code> and <code style="color : Cyan">RefineMaxY</code>, in the script <code style="color : Darkorange">create_mesh_cylinder_wake.py</code>, which define the boundaries of the outer refinement region, and parameters <code style="color : Cyan">Refine2MinX</code>, <code style="color : Cyan">Refine2MaxX</code>, <code style="color : Cyan">Refine2MinY</code> and <code style="color : Cyan">Refine2MaxY</code>, define the boundaries of the inner refinement region as given below:

```py
# Geometry parameters for refinment regions

RefineMinX                  = -7
RefineMaxX                  = 20
RefineMaxY                  = 5
Refine2MinX                 = -1.5
Refine2MaxX                 = 15
Refine2MaxY                 = 1.5
```

The magnitudes of <code style="color : Darkorange">RefinementFactorCoarse</code>, <code style="color : Darkorange">RefinementFactorFine</code> and <code style="color : Darkorange">RefinementFactorFinest</code>, define the element sizes inside of the three regions of the mesh. Furthermore, all of them these have been scaled via the parameter <code style="color : Darkorange">globalRefinementFactor</code>, such that $\texttt{element size}= \texttt{globalRefinementFactor} * \texttt{RefinementFactor}$ (Coarse, Fine or Finest). 

```py
# Mesh sizes

globalRefinementFactor      = 1.5
RefinementFactorCoarse      = 2
RefinementFactorFine        = 0.1
RefinementFactorFinest      = .05
```

We generate the mesh file : <code style="color : Darkorange">CylinderWakeMesh.msh</code> directly using the python file <code style="color : Darkorange">generate_cylinder_mesh.py</code>. However, it is extremely important to make sure that <code style="color : Cyan">RefineMinX</code> < <code style="color : Cyan">Refine2MinX</code> and <code style="color : Cyan">Refine2MaxX</code> < <code style="color : Cyan">RefineMaxX</code> as well as <code style="color : Cyan">Refine2MaxY</code> < <code style="color : Cyan">RefineMaxY</code>, otherwise it wouldn't work. 

```py
# Make sure that RefMinX < Refine2MinX and Refine2MaxX < RefineMaxX as well as Refine2MaxY < RefineMaxY 


mesh_definitions = generate_cylinder_mesh.create_custom_mesh(RefineMinX, RefineMaxX, RefineMaxY, Refine2MinX, Refine2MaxX, Refine2MaxY,globalRefinementFactor, RefinementFactorCoarse, RefinementFactorFine, RefinementFactorFinest)

generate_cylinder_mesh.create_custom_mesh.gen_cylinder_mesh(mesh_definitions)
```

The following piece of code imports all the parameters in the script <code style="color : Darkorange">test.json</code>. The file <code style="color : Darkorange">CylinderWakeMesh.msh</code>, is a spatial discretization of a **symmetric** rectangular domain around a cylinder with a **unit** diameter. The viscosity of the fluid around the cylinder is given by <code style="color : Darkorange">MolVisc</code>, which is modelled as <code style="color : Darkorange">Constant</code>. We only solve the steady-state Navier-Stokes equations, without any source terms that in-general describe a base-flow. These non-dimensionalised equations are given by, wherein the pressure has been non-dimensionalized via $p = \frac{\mu \bold{u}}{d} \tilde{p}$

$$
    \nabla \cdot \bold{u} = 0 ,
$$

$$
    \mathrm{Re}  \left(\bold{u} \cdot \nabla \right) \bold{u} = -\nabla p +  \nabla^2 \bold{u} .
$$



Where $\bold{u} = \{u_x, u_y\}$ and the Reynolds number, $\mathrm{Re} = 50$ since $\nu = 0.02$ and $U_\infty = 1$. 

```py
# reading parameters from the json file.  

param = config()
param.importFromFile('test.json')
```

The flow-field throughout the domain is initialised as $U_\infty = 1$ inside the domain, but at the boundaries: 

1. Inlet &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; :  &nbsp; &nbsp;  $u_x = 1$, $u_y = 0$  &  $\eta \cdot \nabla p = 0$,

2. Symmetry  &nbsp; : &nbsp; &nbsp;  $\eta \cdot \nabla u_x = 0$, $u_y = 0$ & $\eta \cdot \nabla p = 0$,

3. Outlet &nbsp; &nbsp; &ensp; &nbsp; : &nbsp; &nbsp; $\eta \cdot \nabla u_x = 0$, $u_y = 0$ & $\eta \cdot \nabla p = 0$,

4. Top &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &ensp; : &nbsp; &nbsp; $\eta \cdot \nabla u_x = 0$, $u_y = 0$ & $\eta \cdot \nabla p = 0$,

5. Wall &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; : &nbsp; &nbsp; $u_x, u_y = 0$ & $\eta \cdot \nabla p = 0$. 

However, in the file <code style="color : Darkorange">bc.json</code>, the specified boundary conditions may seem different than what they should ideally be. Especially, the Dirichlet conditions for velocity at the **Inlet**, **Outlet** and **Top** of the domain. This is because, these values are superimposed on the ones present in <code style="color : Darkorange">meanflow.fel</code>    

A mixed finite element function space is defined via $W_h = \bold{V}_h^2 \times V_h^1$, where $\bold{V}_h^k, V_h^k$ are vector and scalar spaces respectively for a triangulation $\mathcal{T}_h$ that decomposes the domain $ \Omega = \bigcup_{T \in \mathcal{T}_h } T $ consisting of $k^{th}$ order compactly supported polynomials $P_k (T)$.

However, if one needs to visualize meanflow field it has been exported to the h5 file : <code style="color : Darkorange">meanflow.h5</code>  with a higher resolution, via mapping the $P_2$ finite elements to $P_1$ finite elements.


```py
# Reading the mesh based on the filename stated in the json file.

mesh=param.__mesh__

# Constructs a mixed finite element function space based on P2 and P1 scalar finite elements.

FEMSpaces = DefineFEMSpaces.FEMSpacesClass(param,mesh,)

# Initialize the solution via importing data from meanflow.fel

meanFlow = meanFlowClass(param, FEMSpaces, mesh)
meanFlow.importDataFromFile()

# Exporting the meanflow as a h5 file and then mapping a P2 FE space to a P1 FE space

meanflowFilename = 'meanflow.h5'
meanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)
```

With $\bold{v} \in \bold{V}_h^2$ and $q \in V_h^1$ we operate the strong form and propose the following variational formulation, 

$$
Re \int_\Omega  \bold{u} \cdot \nabla \bold{u} \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x - \int_\Omega p \nabla \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x + \int_\Omega \nabla \bold{u} : \nabla \bold{v} \hspace{0.2cm} \mathrm{d}x  + \int_\Omega \nabla \cdot \bold{u} q  \hspace{0.2cm} \mathrm{d}x = 0 .
$$

It is important to note that we shift all the terms to the left hand side, and perform integration by parts on the pressure term. Furthermore, we will discretize this variational problem using the corresponding basis functions from the mixed finite element space, and solve the discretized problem using the Newton's method. Hence,

$$ 
    f(x) = Re \int_\Omega  \bold{u}^{(k)} \cdot \nabla \bold{u}^{(k)} \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x - \int_\Omega p^{(k)} \nabla \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x + \int_\Omega \nabla \bold{u}^{(k)} : \nabla \bold{v} \hspace{0.2cm} \mathrm{d}x  + \int_\Omega \nabla \cdot \bold{u^{(k)}} q  \hspace{0.2cm} \mathrm{d}x .
$$

$$
    f'(x) \delta x = \underbrace{\int_\Omega  \bold{u}^{(k)} \cdot \nabla \delta \bold{u} \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x}_{C_1\delta \bold{u}} + \underbrace{\int_\Omega \delta \bold{u} \cdot \nabla \bold{u}^{(k)} \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x}_{C_2\delta \bold{u}} - \underbrace{\int_\Omega \delta p \nabla \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x}_{B^T\delta p} + \underbrace{\int_\Omega \nabla \delta \bold{u} : \nabla \bold{v} \hspace{0.2cm} \mathrm{d}x}_{D\delta \bold{u}}  + \underbrace{\int_\Omega (\nabla \cdot \delta \bold{u}) q  \hspace{0.2cm} \mathrm{d}x}_{B\delta \bold{u}} .
$$

This can be assembled into a system matrix given by,

<!-- $$    
    \begin{bmatrix}
        C_1 + C_2 + D & B^T \\ 
        B   & 0       
    \end{bmatrix}

    \begin{bmatrix}
        \delta \bold{u} \\ 
        \delta p       
    \end{bmatrix}
      = - 
    \begin{bmatrix}
        r_\bold{u} \\
        r_p
    \end{bmatrix}
$$ -->
```math
\left[
\begin{array}{cc}
C_1 + C_2 + D & B^T \\
B & 0
\end{array}
\right]
\left[
\begin{array}{c}
\delta \bold{u} \\
\delta p
\end{array}
\right]
= -
\left[
\begin{array}{c}
r_{\bold{u}} \\
r_p
\end{array}
\right]
```

This is of the form of the Newton's method which is formulated as $x_{new} = x_{old} + \delta x$, where $\delta x = - \frac{f(x)}{f'(x)}$. Here, $r_\bold{u}$ and $r_p$ together represent $f(x^{(k)})$, where $^k$ is the current iterate. Therefore,

$$
    r_\bold{u} = Re \int_\Omega  \bold{u}^{(k)} \cdot \nabla \bold{u}^{(k)} \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x - \int_\Omega p^{(k)} \nabla \cdot \bold{v} \hspace{0.2cm} \mathrm{d}x + \int_\Omega \nabla \bold{u}^{(k)} : \nabla \bold{v} \hspace{0.2cm} \mathrm{d}x.
$$
$$
    r_p = \int_\Omega \nabla \cdot \bold{u}^{(k)} q  \hspace{0.2cm} \mathrm{d}x 
$$

```py
# Getting the variational formulation for the Navier-Stokes continuity and the momentum equation.

equation = EquationCollectionClass(param, FEMSpaces, meanFlow, mesh)   

# Creation of an instance of the function space called baseFlow for saving the solution

baseFlow = Field(FEMSpaces.VMixed, mesh)

# Initialization of baseFlow using the  meanFlow

baseFlow.interpolateFieldsFromMeanfield(meanFlow)

# ---------------- START LOOP ---------------------------------------------- 

## start Newton solver
target_residuum = 1.e-11
# track time
start= time.time()
# calculate the base flow
i=0
residuum=1.
[u,p] = baseFlow.getListOfSingleFields()
meanFlow._fieldDict['u'] = u.function
meanFlow._fieldDict['p'] = p.function
N = equation.getNonlinearExpression(meanFlow)
while(residuum > target_residuum and i<10):
    i+=1

    # solve equation system
    L = equation.getLinearOperator(meanFlow)
    newtonSummand_array = LinearSolver.solveEquationSystem(L,N)
    L.destroy()
    N.destroy()

    # update baseFlow
    baseFlow_array = baseFlow.getCoefficientArray() + newtonSummand_array
    baseFlow.setCoefficientArray(baseFlow_array)
    
    # calculate nonlinear expression & residuum
    [u,p] = baseFlow.getListOfSingleFields()
    meanFlow._fieldDict['u'] = u.function
    meanFlow._fieldDict['p'] = p.function
    N = equation.getNonlinearExpression(meanFlow)
    residuum = np.linalg.norm(N.getArray())
    printDebug(True, "-------------------------------------------------------------" )
    printDebug(True, "-- Base flow iteration: "+str(i)+"; Residuum: %4g " % residuum)
    printDebug(True, "-------------------------------------------------------------" )



# end tracking time
end = time.time() - start
printDebug(True, '-- Solving the base flow problem took %4g s' % end)
printDebug(True, '-- Residuum:  %12g' % (residuum))

```

6. Saving the solution in <code style="color : Darkorange">.xdmf</code> format, and the velocity field in Paraview can be visualized as:

![Base-flow](./BaseFlow.png)

The same tutorial can be found inside the jupyter notebook ```cylinderWake.ipynb```. We highly encourage readers to try it out.