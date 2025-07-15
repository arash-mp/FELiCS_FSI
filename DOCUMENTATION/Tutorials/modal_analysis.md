# Tutorial 2: Modal Analysis
## Goals of the tutorial
In this tutorial, we will do an eigenvalue decomposition of the base flow obtained in the tutorial of [cylinder wake](./cylinder_wake.md). By the end of this tutorial, you will be able to:

- Run a modal analysis case.
- Postprocess modal analysis results with paraview.
## Requirements
Before you begin this tutorial make sure to 
* have completed the [base flow tutorial](./cylinder_wake.md).
* access the case folder ```felics2.0/TUTORIALS/modal_analysis_tutorial``` and copy it into your working directory.

## Modal analysis settings
### Boundary conditions
We define the boundary conditions (BCs) for the modal analysis in [```bc_Modal.json```](./../../TUTORIALS/modal_analysis_tutorial/bc_modal.json). The file structure is detailed in [Setting files](hhttps://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/development/DOCUMENTATION/Running_FELiCS/FELiCS_settings.md?ref_type=heads).
The boundary conditions for the perturbations are: 

| Boundary | $u'_x$     | $u'_y$      | $p'$       |
|:----------|:-----------|:-----------|:-----------|
| <code style="color : Darkorange">Inlet</code>    | Dirichlet | Dirichlet | Dirichlet   |
| <code style="color : Darkorange">Symmetry</code> | Dirichlet   | Neumann | Dirichlet   |
| <code style="color : Darkorange">Outlet</code>   | Dirichlet   | Dirichlet   | Dirichlet |
| <code style="color : Darkorange">Top</code>      | Dirichlet   | Dirichlet | Dirichlet   |
| <code style="color : Darkorange">Wall</code>     | Dirichlet | Dirichlet | Neumann   |

**Note:** The BCs for the base flow in [base flow tutorial](./cylinder_wake.md) and for modal analysis are different. 

### Settings
The setting file, [```modal.json```](./../../TUTORIALS/modal_analysis_tutorial/modal.json) contains all the information relevant for modal analysis. The complete architecture of [```modal.json```](./../../TUTORIALS/modal_analysis_tutorial/bc_modal.json) can be found inside [FELiCS settings](./../Running_FELiCS/FELiCS_settings.md). 

Since the modal analysis is done for the base flow around the cylinder, we use the same mesh [```cylinder_wake.msh```](./../../TUTORIALS/cylinder_wake_tutorial/cylinder_wake.msh). Hence, the setting file [```modal.json```](./../../TUTORIALS/modal_analysis_tutorial/modal.json) mentions:
```json
"MeshFilePath":"cylinder_wake.msh"
```

Furthermore, we use the base flow file [```base_flow_for_FELiCS```](./../../TUTORIALS/modal_analysis_tutorial/base_flow_for_FELiCS.fel), which was get generated in the [base flow tutorial](./cylinder_wake.md). It is defined:
```json
"MeanFlowFilePath": "base_flow_for_FELiCS.fel"
```
The BCs path as:
```json
"BCsFilePath": "bc_modal.json"
```
Lastly, we initialize the solver using the eigenvalue guess, and the number of eigenvalues we want to compute:
```json
"EigenValueGuess": ["0.7"], 
"nSolut": 100
```
## Running the analysis
You can run the modal analysis via 
```sh
FELiCS -f modal.json
```
## Postprocessing
After running the modal analysis, the working directory should look like:

```bash
.
└── logs
└── output_dir
├ base_flow_for_FELiCS.fel
├ bc_modal.json
├ cylinder_wake.msh
├ modal.json
├ modal.json
├ PlotScatter.py
└ ...
```

Inside the ```output_dir``` directory, the all the eigenmodes in ```.h5``` and ```.xmf``` format can be found, along with eigenspectrum in ```spectrum.csv``` file. 
Run the pyhton script [PlotScatter.py](./../../TUTORIALS/modal_analysis_tutorial/PlotScatter.py) to plot the complete eigenvalues spectrum:
![](./../../TUTORIALS/modal_analysis_tutorial/eigenspectrum.png) <a id="fig:EigSpec"></a>
Figure 1. Eigenspectrum

In [Figure 1](#EigSpec), an eigenvalue with a positive imaginary part stands out. We use paraview to visualize the corresponding real part of the eigenmode stored in `ModalSolution_Omega_Direct_(0.745+0.013j).xmf`. 

![](./../../TUTORIALS/modal_analysis_tutorial/ux_real.png) <a id="fig:RealUx"></a>
Figure 2. Real eigenmode, $u'_x$

![](./../../TUTORIALS/modal_analysis_tutorial/uy_real.png) <a id="fig:RealUy"></a>
Figure 3. Real eigenmode, $u'_y$

![](./../../TUTORIALS/modal_analysis_tutorial/p_real.png) <a id="fig:Realp"></a>
Figure 4. Real eigenmode, $p'$
