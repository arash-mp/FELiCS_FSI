# 1. Goals of the tutorial.
In this first tutorial, an introduction to setting up a FELiCS case will be given via solving for the base flow around a cylinder for a Reynolds number, $\mathrm{Re} = 50$.

# 2. Requirements.

* FELiCS should be installed as per the [FELiCS installation guide](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/development/DOCUMENTATION/installation_guide.md?ref_type=heads) with an active conda environment.

* Choose a branch of the FELiCS that contains the base flow solver.

# 3. Case setup.

In general it is recommended to create a folder for every case, where the source code, input and output files will be located. Here, the case folder can be found in ```felics2.0/TUTORIALS/cylinder_wake_tutorial```. Therefore, simply copy this tutorial to your working directory ```workDir```
```sh
cp -r felics2.0/TUTORIALS/cylinder_wake_tutorial workDir/

cd workDir/cylinder_wake_tutorial/ 
```

# 4. Mesh generation.

Since the code is using the Finite Element Continuous Galerkin
approach, a computational grid needs to be created, which spatially discretizes the domain. For 2D computations like the one in this tutorial, the code uses triangular elements only. In order to create the mesh, the package ```python-gmsh``` is used. An input file for this program called ```cylinder_wake_mesh.py``` is already prepared in the folder. Upon exceuting it in ```python``` it generates the mesh file, ```cylinder_wake.msh``` in Version 2 *ASCII* format. The documentation for using the ```python-gmsh``` package can be found [here](). 

```sh
python cylinder_wake_mesh.py
```
After generating the mesh it will look like:

![Image1](../../TUTORIALS/cylinder_wake_tutorial/CylinderWakeMesh.png) <a id="fig:Mesh"></a>
Figure 1. Mesh.


Upon magnifying in the vicinity of the cylinder it looks like

![Image2](../../TUTORIALS/cylinder_wake_tutorial/magCylinderWake.png) <a id="fig:MagMesh"></a>
Figure 2. Magnified version of [Figure 1](#fig:Mesh)

The file ```cylinder_wake.msh``` consists of the following domains.

```sh
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
6
1 1 "Inlet"
1 2 "Symmetry"
1 3 "Outlet"
1 4 "Top"
1 5 "Wall"
2 6 "all"
$EndPhysicalNames
```
In [Figure 1](#fig:Mesh), <code style="color : Darkorange">Inlet</code> is the left vertical boundary of the domain, <code style="color : Darkorange">Outlet</code> is the right vertical boundary of the domain, <code style="color : Darkorange">Top</code> is the top horizontal boundary of the domain. <code style="color : Darkorange">Symmetry</code> is the bottom horizontal boundary of the domain, except the half-cylinder in [Figure 2](#fig:MagMesh), which is the <code style="color : Darkorange">Wall</code>.


# 5. Obtaining the base flow. 
The base flow for FELiCS can be obtained by various sources depending on the case: Numerical simulations, such as RANS, URANS, LES, DNS, as well as experimental results or analytical models. In this case, we will calculate our base flow ourselves. But don't worry, everything is prepared: a finite element Newton solver called ```solveBaseFlow.py``` can be found in the working directory. The newton flow solver doesn't need to be adapted, nevertheless, the interested reader will find that it is not hard at all to solve the flow equations for different Reynolds numbers. Before the solver can be started, the conda environment FELiCS created during the installation needs to be activated:

```sh
conda activate felics2025_dolfin9
```

The settings and the boundary conditions can be found in ```Re50.json``` and ```bc.json```. For $\mathrm{Re} = 50$,  the diameter of the cylinder, $d = 1$ and the bulk velocity $U_\infty = 1$, gives the viscostiy, $\nu = 0.02$. Therefore in the file ```Re50.json``` we have

```json
"Molvisc":0.02
```
Furthermore, since we only solve the non-dimensionalized Navier-stokes continuity and the momentum equation we set, 

```json
"Reaction": false
```
and in the <code style="color : Cyan">"SetofEquations"</code>, we set everything else to <code style="color : Darkorange">None</code>, except for <code style="color : Cyan">"Momentum"</code> and <code style="color : Cyan">"Mass"</code> equations. 

Additionally we import the boundary conditions from the ```bc.json``` file. Overall the boundary conditions can be formulated as : 

| Boundary | $u_x$     | $u_y$      | $p$       |
|:----------|:-----------|:-----------|:-----------|
| <code style="color : Darkorange">Inlet</code>    | Dirichlet | Dirichlet | Neumann   |
| <code style="color : Darkorange">Symmetry</code> | Neumann   | Dirichlet | Neumann   |
| <code style="color : Darkorange">Outlet</code>   | Neumann   | Neumann   | Dirichlet |
| <code style="color : Darkorange">Top</code>      | Neumann   | Dirichlet | Neumann   |
| <code style="color : Darkorange">Wall</code>     | Dirichlet | Dirichlet | Neumann   |


Finally, the mean-flow field is imported via
```json
"MeanFlowFilePath": "meanflow.fel"
```

Lastly, we obtain the base flow via solving the python script ```solveBaseFlow.py``` as

```sh
python solveBaseFlow.py
```
This will output the file ```base_flow.xdmf```, which in paraview can be visualised as 

![Image3](../../TUTORIALS/cylinder_wake_tutorial/BaseFlow.png)
Figure 3. Magnitude of the flow-field. 


<!-- # References
1.   Thomas Ludwig Kaiser, Thierry Poinsot, and Kilian Oberleithner. “Stability and Sensitivity Analysis of
Hydrodynamic Instabilities in Industrial Swirled Injection Systems”. In: Journal of Engineering for Gas
Turbines and Power 140.5 (Jan. 2018), p. 051506. ISSN: 0742-4795. DOI: 10.1115/1.4038283. 
URL: http://gasturbinespower.asmedigitalcollection.asme.org/article.aspx?DOI=10.1115/1.4038283.

 2. Thomas Ludwig Kaiser et al. “Examining the Effect of Geometry Changes in Industrial Fuel Injection Systems On Hydrodynamic Structures with Biglobal Linear Stability Analysis”. In: Journal of
Engineering for Gas Turbines and Power (Sept. 2019). ISSN: 0742-4795. DOI: 10.1115/1.4045018. URL: https://asmedigitalcollection.asme.org/gasturbinespower/article/doi/10.1115/1.4045018/1047002/Examining-the-Effect-of-Geometry-Changes-in.

3. Thomas Ludwig Kaiser et al. “Impact of symmetry breaking on the Flame Transfer Function of a laminar premixed flame”. In: Proceedings of the Combustion Institute 37.2 (2019), pp. 1953–1960. ISSN: 1540-7489. DOI: https://doi.org/10.1016/j.proci.2018.06.047. URL: http://www.sciencedirect.com/science/article/pii/S154074891830230X. -->