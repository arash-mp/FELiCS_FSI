# Overview of FELiCS workflow

## 0. Disclaimer

⚠️ The present file represents a *work in progress* version of a guide for the FELiCS workflow. 

🚧 Current limitations and ToDos:
* Update the name of the methods/variables to fit the recent refactoring of the FELiCS code.
* Update and complete the links to the documentation pages
* Complete the description of the workflow (some parts are not yet covered)
* Figure out an easier way to visualize the guide (collapsible list?, ...)

## 1. Initialisation & Meshing

```python
param = Config()
```
*Initialise a FELiCS [`Config()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#module-FELiCS.Parameters.Config) class instance*

1. Fills in defaults parameters from [`get_all_settings_dict()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_all_settings_dict)


```python
param.import_from_file(filename)
```
*Reads in parameters from json file and compute needed extra parameters*

1. Reads in data from json file
2. Overwrite default parameters with file values
    
    > 💡 Only default parameters are read in! To add new parameters, make sure they are also added to defaults stored in [`Config.get_all_settings_dict()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_all_settings_dict)

3. Check presence of FELiCS [mandatory files](https://felics-d43476.gitlab.io/Running_FELiCS/index.html#overview) in running directory
4. Initialize a FELiCS [`MixtureClass`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/MixtureClass/index.html#module-FELiCS.Equation.MixtureClass) class  instance
    * Sets default physical parameters
    * Reads in physical parameters from a FELiCS [mixture file](https://felics-d43476.gitlab.io/Running_FELiCS/FELiCS_settings.html#structure-of-the-mixture-json-file) if present in the running directory.

5. Read the domain data (*not conceptually clear*) with [`read_domain_data()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.read_domain_data)
    *  Initialize a [`FELiCSMesh`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FELiCSMesh/index.html#module-FELiCS.SpaceDisc.FELiCSMesh) class instance.
        * Reads in the [.msh](https://felics-d43476.gitlab.io/Running_FELiCS/index.html#mesh-file) mesh file and setup DolfinX mesh objects
        * Initialize [`CoordinateSystem`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Misc/tensorUtils/index.html#FELiCS.Misc.tensorUtils.CoordinateSystem) class instance, which defines geometric and spectral dimensions of the problem.

    * Extracts equation list from parameters, [`get_equation_list()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_equation_list), and gets corresponding transported quantities, [`get_state_vector_variables()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_state_vector_variables). (SHOULD MOVE TO [`calculate_parameters()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.calculate_parameters)??)

6. Compute additional parameters: [`calculate_parameters()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.calculate_parameters)

    * Number of velocity components: [`get_velocity_components()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_velocity_components) (spectral and resolved)
    * Get variables of state vector [`get_transported_quantity_list()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_transported_quantity_list) (**DUPLICATED from** [`get_state_vector_variables()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_state_vector_variables)**??**)
    * Sets total number of dimensions in [`FELiCSMesh`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FELiCSMesh/index.html#module-FELiCS.SpaceDisc.FELiCSMesh) and the [`CoordinateSystem`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Misc/tensorUtils/index.html#FELiCS.Misc.tensorUtils.CoordinateSystem) objects (**WEIRD TO REDO HERE**)

```python
mesh = param.get_mesh()
```
*For easy access to the [`FELiCSMesh`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FELiCSMesh/index.html#module-FELiCS.SpaceDisc.FELiCSMesh) instance in the rest of the script.*

## 2. Defining FEM Spaces

This step requires the [`param`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#module-FELiCS.Parameters.Config) and [`FELiCSMesh`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FELiCSMesh/index.html#module-FELiCS.SpaceDisc.FELiCSMesh) objects from the previous step and prepares the FEM spaces and objects needed to setup the problem. 

```python
FEMSpaces = FEMSpaces(param, mesh)
```
*Instantiate the FELiCS [`FEMSpaces`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FEMSpaces/index.html#module-FELiCS.SpaceDisc.FEMSpaces) class [REF API]*
1. Creates default FEM spaces (with [`create_function_space()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FEMSpaces/index.html#FELiCS.SpaceDisc.FEMSpaces.create_function_space)) for scalar and vector quantities on the computation mesh with specs from [`param`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#module-FELiCS.Parameters.Config)
2. Prepares similar spaces on the FELiCS export mesh (order 1 only)
    * If not yet instantiated, creates an [`ExportMesh()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FELiCSMesh/index.html#FELiCS.SpaceDisc.FELiCSMesh.ExportMesh) instance
3. Creates the mixed function space corresponding to state vector [`VMixed`]. One on the computation mesh and one on the export mesh.

> 💡 At this stage, it is useful to introduce the different meshes involved in FELiCS runs:
>    1. The **computation mesh**, defined by the `.msh` file given in `config.Case["MeshFilePath"]`. It is used to solve the problem. Currently, both `P1` and `P2` spaces can be used on that mesh.
>    2. The **mean flow mesh**, given in the `.fel` file given in `config.Case["MeshFilePath"]`
>    3. The **export mesh**, saved as `mesh.h5` in the `config.Case[ExportFolder]`. Currently, FELiCS only exports variables on `P1` meshes.  

### Step 3: Setting up mean flow functions
```python
meanFlow = MeanFlowClass(param, FEMSpaces, mesh)
```
*Initialize a [`MeanFlowClass`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/MeanFlowClass/index.html#module-FELiCS.Fields.MeanFlowClass) instance, containing all MeanFlow variables, which are stored in FELiCS [`Fields`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/Field/index.html#module-FELiCS.Fields.Field).*

1. Sets default attributes
2. Creates helper fields (zero field, one field)

```python
meanFlow.import_data_from_file_and_export_to_h5(writer)
```
*Read the required mean fields for the problem from a file and interpolate on the FELiCS mesh (if necessary). Also exports the meanflow fields to a .h5 file.*

1. Collects the names of meanFields needed for the problem: [`meanFlow._getMeanFieldsToBeRead()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/MeanFlowClass/index.html#FELiCS.Fields.MeanFlowClass.MeanFlowClass.add_custom_mean_flow_quantity)

2. Instantiate a [`Reader`](https://felics-d43476.gitlab.io/autoapi/FELiCS/IO/Reader/index.html#FELiCS.IO.Reader.Reader) class object [REF API] linked to the meanFlow path in [`param`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#module-FELiCS.Parameters.Config).

3. For each mean flow variable, creates an empty [`Field`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/Field/index.html#module-FELiCS.Fields.Field) and reads in the values from the meanflow file [`Field.import_data()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/Field/index.html#FELiCS.Fields.Field.Field.import_data)
    * In the [`Field`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/Field/index.html#module-FELiCS.Fields.Field), call the reader instance with [`reader.import_in_field()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/IO/Reader/index.html#FELiCS.IO.Reader.Reader.import_in_field)
        * For more details, see [`Reader`](https://felics-d43476.gitlab.io/autoapi/FELiCS/IO/Reader/index.html#module-FELiCS.IO.Reader) [REF API]

4. Compute additional mean flow variables based on models selected by user with [`init_lam_diff()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/MeanFlowClass/index.html#FELiCS.Fields.MeanFlowClass.MeanFlowClass.init_lam_diff) and [`init_thermodynamic_quantities()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/MeanFlowClass/index.html#FELiCS.Fields.MeanFlowClass.MeanFlowClass.init_thermodynamic_quantities)

5. Using the [`Writer`](https://felics-d43476.gitlab.io/autoapi/FELiCS/IO/Writer/index.html#FELiCS.IO.Writer.Writer) instance, map the meanFlow [`Fields`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/Field/index.html) onto the [`ExportMesh`](https://felics-d43476.gitlab.io/autoapi/FELiCS/SpaceDisc/FELiCSMesh/index.html#FELiCS.SpaceDisc.FELiCSMesh.ExportMesh) and save them to `meanflow.h5` file (+xdmf). 
    * For more details see [`Writer`](https://felics-d43476.gitlab.io/autoapi/FELiCS/IO/Writer/index.html#FELiCS.IO.Writer.Writer)


### Step 4: Setting up the weak form and test/trial functions
```python
equation = EquationCollectionClass(param, FEMSpaces, meanFlow, mesh)
```
*Instantiate the FELiCS [`EquationCollectionClass`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/EquationCollection/index.html#module-FELiCS.Equation.EquationCollection).*

1. Store needed FELiCS objects and ufl geometric objects (`SpatialCoordinate`)
2. Instantiate the FELiCS [`BoundaryHandler`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/Boundary/index.html#FELiCS.Equation.Boundary.BoundaryHandler)
    * Read the FELiCS [`boundary file`](https://felics-d43476.gitlab.io/Running_FELiCS/FELiCS_settings.html#structure-of-the-boundaries-json-file)
    * Set list of BC at each mesh boundary
        * Instantiate a FELiCS [`BoundaryCondition()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/Boundary/index.html#FELiCS.Equation.Boundary.BoundaryCondition) object
        * Set the specific type of boundary
3. Instantiate the FELiCS [`FluctuationClass`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/FluctuationClass/index.html#module-FELiCS.Fields.FluctuationClass), containing the trial functions corresponding to unknowns of the problem.
    * Inherits the FELiCS [`FieldProperties`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/FieldProperties/index.html#module-FELiCS.Fields.FieldProperties) class. **IMPORTANT:** all variables (mean or fluctuations) used in FELiCS equations must appear here (?)
    * Sets up useful FELiCS [`Fields`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/Field/index.html#module-FELiCS.Fields.Field) (zeros, ones, ...)
    * Instantiate a ufl `TrialFunctions` over the `FEMSpaces.VMixed`
    * For each transported variable (defined in [`param.get_transported_quantity_list()`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#FELiCS.Parameters.Config.Config.get_transported_quantity_list)), instantiate a FELiCS [`Tensor`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Misc/tensorUtils/index.html#FELiCS.Misc.tensorUtils.Tensor) with the corresponding trial function. Stores all in `fluc._fieldDict`.
    * Gather list of all fluctuation variables from each of the equations (e.g. [`MomentumHandler._getNeededFieldsForLinearMomentum`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/dependentVariables/MomentumHandler/index.html#module-FELiCS.Equation.dependentVariables.MomentumHandler))
    * Go over all equations handlers (e.g. see [`MomentumHandler`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/dependentVariables/MomentumHandler/index.html)) to connect/substitute all fluctuation variables to state vector fluctuations, using models selected in [`param`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Parameters/Config/index.html#module-FELiCS.Parameters.Config). Repeats process untill all variables are substituted/defined.
    * If needed, instantiate the FELiCS [`ReactionHandler`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/dependentVariables/ReactionHandler/index.html)
    * NOTE: This part of the code will be overhauled/streamlined soon.
4. Instantiate a ufl `TrialFunctions` over the `FEMSpaces.VMixed` (NOTE: redone from fluctuation class)
5. Instantiate a ufl `TestFunctions` also over the `FEMSpaces.VMixed`
6. Fill the EquationCollectionClass.equationList with equations corresponding to param
7. Add sponge term to define damping
    > 💡 Sponge term is only added if spg field was given in mean flow file!


## Step 5: Matrix Assembly
```python
A = equation.getLinearOperator(meanFlow)
B = equation.getWeightMatrix  (meanFlow)
```
*Construct the linear operator matrix for the equation system. Then construct the weight matrix. The workflow is the same for both cases and only presented once*
1. Instantiate empty [`UflDecorator`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/UflDecorator/index.html#module-FELiCS.Equation.UflDecorator)
2. Iterate through the `equationList` and add the linear expression for each equation
3. [`get_assembled_matrix`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Equation/UflDecorator/index.html#FELiCS.Equation.UflDecorator.UflDecorator.get_assembled_matrix) from all linear expressions


## Step 6: Solution & Export
```python
solution = ModeCollection(FEMSpaces.VMixed, mesh, analysisType = "modal")
```
*Instantiate a [`ModeCollection`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/ModeCollection/index.html#module-FELiCS.Fields.ModeCollection) object. This object will serve as a container for the solutions*

```python
for guess in guesses:
    tmp = LinearSolver.solve_general_eigenproblem(A, B, guess, nSol)
    solution.append_solution_of_eigen_problem(tmp, guess)
    n_calculated = nSol

    if adjoint:
        tmp = LinearSolver.solve_general_eigenproblem(A, B, guess, nSol, adjoint=True)
        solution.append_solution_of_eigen_problem(tmp, guess, adjoint=True)
        n_calculated += nSol

    solution.export_spectrum_to_csv(writer)
    solution.export_modes(writer, onlyNewN = n_calculated)
```
*Cycle through eigenvalue guesses and compute N closest eigenvalues for each. Then export the solutions as spectrum and modes.*
1. Solve the eigenvalue problem using SLEPc and PETSc [`LinearSolver.solve_general_eigenproblem`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Solvers/LinearSolver/index.html#FELiCS.Solvers.LinearSolver.LinearSolver.solve_general_eigenproblem)
2. Add the modes to the [`ModeCollection`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/ModeCollection/index.html#module-FELiCS.Fields.ModeCollection) with [`append_solution_of_eigen_problem`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/ModeCollection/index.html#FELiCS.Fields.ModeCollection.ModeCollection.append_solution_of_eigen_problem)
4. Repeat the process for adjoint calculation if stated in the configuration file
5. Update and export the spectrum with all solutions so far
6. Export the Mode for respective guess to [`ModeCollection`](https://felics-d43476.gitlab.io/autoapi/FELiCS/Fields/ModeCollection/index.html#module-FELiCS.Fields.ModeCollection)
    > 💡 This is handled by the [`Writer`](https://felics-d43476.gitlab.io/autoapi/FELiCS/IO/Writer/index.html#module-FELiCS.IO.Writer) Class
