#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
import  numpy as np
from ufl import (
    Dx,
    TestFunctions,
    TrialFunctions,
    dx,
    SpatialCoordinate,
    FacetNormal,
    Measure,
    conj,
    lhs,
    rhs,
    # as_tensor, i, j,
)
from dolfinx.fem import (
    Function,
    dirichletbc,
    Constant,
    form,
    locate_dofs_topological,
    assemble_scalar,
)
from dolfinx.fem.petsc import (
    assemble_matrix,
    assemble_vector,
    set_bc

)
from FELiCS.Misc.tensorUtils import (
    Tensor,
    as_vector,
    iDot,
    iConj,
    # iInner,
)
from    petsc4py.PETSc               import ScalarType
from    FELiCS.Equation.UflDecorator import UflDecorator
from    FELiCS.Equation.Boundary     import BoundaryHandler
from 	FELiCS.Misc.logging          import Logger


# Get the logger
logger = Logger.get_logger("felics")

class EquationCollectionClass():
    """
    Constructs variational formulations for equations and manages boundary conditions.

    This class generates variational formulations for various matrices, including:
    - A (real and imaginary parts of the linear operator)
    - B (real and imaginary parts of the weight/time-derivative matrix)
    - B_forcing (resolvent forcing norm)
    - B_response (resolvent response norm)

    It also manages test and trial functions, boundary conditions, and mesh-based
    coordinate transformations. The design supports multiple physical models
    (momentum, mass, energy, species) and enables configuration through an external
    parameter object.

    **Initialize the EquationCollectionClass object**

    Parameters
    ----------
    param : object
        Parameters for the simulation, including boundary conditions, equations, and I/O settings.
    FEMSpaces : object
        Collection of Finite Element Method spaces for the mixed formulation.
    mean : object
        Mean flow fields used in variational formulations.
    mesh : object
        Mesh object defining the spatial domain and coordinate system.

    Attributes
    ----------
    x : ufl.SpatialCoordinate
        Spatial coordinates of the computational domain.
    ds : ufl.Measure
        Measure for integration over boundary subdomains.
    all_ds : ufl.Measure
        Combined measure over all boundary subdomains.
    n_BC : ufl.FacetNormal
        Outward-pointing unit normal on the domain boundary.
    n : Tensor
        Tensor form of the boundary normal vector in the chosen coordinate system.
    BCs : list
        List of Dirichlet boundary condition objects.
    trialFunctionsFEM : ufl.TrialFunctions
        Mixed trial functions used in the variational problem.
    testFunctionsFEM : ufl.TestFunctions
        Mixed test functions used in the variational problem.
    X : list of Tensor
        List of structured test functions with tensor-awareness.
    R : ufl.Coefficient
        Radial coordinate value for cylindrical systems or unity for Cartesian.
    equationList : list
        List of equation objects that define the weak formulations.
    resolventResponseIndices : numpy.ndarray
        Indices selecting solution components used in the resolvent response norm.
    resolventForcingIndices : numpy.ndarray
        Indices selecting solution components used in the resolvent forcing norm.

    Notes
    -----
    The class dynamically loads and constructs equation objects depending on
    the configuration provided in `param.Case.Equations`. It supports tensor-based
    formulations and multiple coordinate systems.

    Examples
    --------
    >>> eq_collection = EquationCollectionClass(param, FEMSpaces, mean, mesh)
    >>> A = eq_collection.getLinearOperator(mean)
    >>> B = eq_collection.getWeightMatrix(mean)
    """


    def __init__(self, param, FEMSpaces, mean, mesh):
        """
        Initialize the EquationCollectionClass with simulation parameters and spaces.

        This method sets up the computational domain, coordinate systems, 
        boundary conditions, and prepares the equation list based on the 
        specified set of equations.

        Parameters
        ----------
        param : object
            Parameters for the simulation, including boundary conditions and equations.
        FEMSpaces : object
            Finite Element Method spaces for the simulation.
        mean : object
            Mean flow fields used in variational formulations.
        mesh : object
            Mesh defining the domain.
        """
        
        from FELiCS.Fields.fluctuationClass import fluctuationClass
        logger.info('Initializing the equation collection class.')

        # Store input parameters as private attributes
        self._param            = param
        self._FEMSpaces        = FEMSpaces
        self._mean             = mean
        self._mesh             = mesh

        # Get spatial coordinates
        self.x                  = SpatialCoordinate(mesh.dolfinxMesh)
        self._coordinateSystem  = mesh.coordinateSystem
        # self._coordinateSystem.setTrueDimension(len(self._param.getVelocityComponents())) # NOTE: moved to config.py

        ## BOUNDARIES
        # Initialize boundary handler
        self.variables          = param.Case.StateVectorVariables
        self.boundaryHandler    = BoundaryHandler(self.variables, mesh, param.BoundaryCondition.BCsFilePath)
        # Initialize ds: Get all boundaries (So far hard coded)
        self.ds                 = Measure("ds", subdomain_data=mesh.facet_tags)
        self.all_ds             = self.ds(self.boundaryHandler.IDs[0])
        for i in range(1,len(self.boundaryHandler.IDs)):
            self.all_ds += self.ds(self.boundaryHandler.IDs[i])
        # Get boundary normals
        self.n_BC               = FacetNormal(self._FEMSpaces.P2.mesh)
        self.n                  = Tensor(
            as_vector((self.n_BC[0], self.n_BC[1])),
            self._coordinateSystem
        )
        # Get Dirichlet boundary conditions
        self.BCs                = self.boundaryHandler.getListOfDirichletBCsForDolfinx(FEMSpaces.VMixed)


        ## TEST AND TRIAL FUNCTIONS
        # Define test and trial functions
        fluctuationC = fluctuationClass(
            param,
            mean,
            FEMSpaces,
            self._coordinateSystem,
        )
        self.trialFunctionsFEM  = TrialFunctions(FEMSpaces.VMixed)
        self.hat                = fluctuationC.fluc
        self.fluctuationC       = fluctuationC

        fluc = {}
        for sol in self._param.Case.SolutionList:
            fluc[sol]           = self.hat[self._param.Case.SolutionList.index(sol)]
        
        XTemp                   = TestFunctions(self._FEMSpaces.VMixed)
        self.testFunctionsFEM   = XTemp
        
        # Test function
        X = []
        for i in list(XTemp):
            X.append(Tensor(
                i,
                self._coordinateSystem,
                mayHaveSpectralDimension=True,
                ))
        self.X = X
            
        # Get radial coordinate
        if self._param.Case.CoordinateSystem in ['Cylindrical']:
            self.R = self.x[1]
        else:
            from petsc4py import PETSc
            self.R = Constant(self._FEMSpaces.P2.mesh, PETSc.ScalarType(1.0))

        logger.debug('State vector: %s.' % param.Case.SolutionList)
 
        # Create equation list from parameters
        # TODO: the following could be in a method for readability
        logger.debug("Creating equation list.")
        self.equations = param.Case.Equations
        self.equationList = []

        for equation in self.equations:
            index = self.equations.index(equation) 
            if   equation[0]  == "Momentum" and equation[1]["Equation"] == "NSPrimitive":  
                from FELiCS.Equation.Equations.MomentumEquation import MomentumEquation
                logger.debug('Adding momentum equation for u-fluc -> X[%d].' % index)     # Hardcoded u' for mom eq.
                eqObject    = MomentumEquation(index,self,fluctuationC,X[index],param)
                self.equationList.append(eqObject)

            elif equation[0]  == "Mass" and equation[1]["Equation"] == "Continuity":  
                from FELiCS.Equation.Equations.MassEquation import MassEquation
                logger.debug('Adding mass equation for %s-fluc -> X[%d].' % (index, index))
                eqObject    = MassEquation(index,self,fluctuationC,X[index],self._param)
                self.equationList.append(eqObject)

            elif equation[0]  == "Energy" and equation[1]["Equation"] == "Enthalpy":  
                from FELiCS.Equation.Equations.EnthalpyEquation import EnthalpyEquation
                logger.debug('Adding enthalpy-energy equation for %s-fluc -> X[%d].' % (index, index))
                eqObject    = EnthalpyEquation(index,self,fluctuationC,X[index],self._param)
                self.equationList.append(eqObject)

            elif equation[0]  == "Energy" and equation[1]["Equation"] == "primitive-p":  
                from FELiCS.Equation.Equations.EnergyPressureEquation import EnergyPressureEquation
                logger.debug('Adding pressure-energy equation for %s-fluc -> X[%d].' % (index, index))
                eqObject     = EnergyPressureEquation(index,self,fluctuationC,X[index],self._param)
                self.equationList.append(eqObject)

            # TODO: Jens: put species equations back in FELiCS 
            # Notes from Sophie:
            # - Now for every species, the equation "Species" has to given in the main json file with the appropriate variable name.
            # - Info on the species can still be found in the "Mixture.json" file.  
            elif equation[0]  == "Species" and equation[1]["Equation"] == "Non-conservative":
                from FELiCS.Equation.Equations.SpeciesEquation import SpeciesEquation
                species       = equation[1]["Variable"]
                logger.debug(f"Adding non-conservative species equation for %s-fluc -> X[%d], species name: {species}." % (index, index))
                eqObject      = SpeciesEquation(index,self,fluctuationC,X[index],species,self._param)
                self.equationList.append(eqObject)
                            
                #        #elif self._param.Case.SetOfEquations['Species']['Equation'] == 'Conservative':
                #        #    # This eq has not been derived in tensor framework yet.
                #        #    from FELiCS.Equation.speciesConservative.addSpeciesConservativeEq import addSpeciesConservativeEq
                #        #    print('-- Adding equation for species '+specie +' in conservative form')
                #        #    addSpeciesConservativeEq(self,fluctuationC,X[i_eqn],self.mean,specie,self._param)
                #        else:
                #            logger.error('Species transport equation type ' + self._param.Case.SetOfEquations['Species']['Equation'] + ' unknown.' )
                #            raise Exception('Species transport equation type ' + self._param.Case.SetOfEquations['Species']['Equation'] + ' unknown.' )

            elif equation[0] not in  ["EquationOfState", "ProgressVariableLinear"] :
                logger.error('Equation type ' + str(equation)  + ' unknown.' )
                raise Exception('Equation type ' + str(equation)  + ' unknown.' )

            
        # Add sponge region to equation list only if the field was given in the mean flow file
        if 'spg' not in mean._notInFileList:
            from FELiCS.Equation.Equations.SpongeTerm import SpongeTerm
            logger.debug('Adding sponge damping.')
            eqObject      = SpongeTerm(self.equations,self,fluctuationC,X,self._param) # give equationsList-Dictionary as "index"
            self.equationList.append(eqObject)


        if self._param.Case.AnalysisMode in ['Resolvent']:
            logger.debug('Defining specific forms for resolvent analysis.')
            self.computeResolventNorms     (X,self._param,mean,fluctuationC)
            self.computeResolventFEMWeights(X,self._param,mean,fluctuationC)

            # Get indices for forcing and response, depending on used norm, 
            # to use when creating the shrinker matrices
            index_u = param.Case.SolutionList.index('u')
            
            # Get indices for response norm
            if param.IOResolvent.ResponseNorm == 'Chu':
                # Variables used: all
                self.resolventResponseIndices = np.arange(self._FEMSpaces.VMixed.dofmap.index_map.local_range[1]) # whole size of VMixed
            elif param.IOResolvent.ResponseNorm == 'TKE':
                # Variables used: u
                self.resolventResponseIndices = self._FEMSpaces.VMixed.sub(index_u).collapse()[1]
            else:
                logger.error(f"Response norm type '{param.IOResolvent.ResponseNorm}' not implemented.")
                raise Exception(f"Response norm type '{param.IOResolvent.ResponseNorm}' not implemented.")
                
            # Same for the forcing norm
            if param.IOResolvent.ForcingNorm == 'Chu':
                # Variables used: all
                self.resolventForcingIndices = np.arange(self._FEMSpaces.VMixed.dofmap.index_map.local_range[1]) # whole size of VMixed
            elif param.IOResolvent.ForcingNorm == 'TKE':
                # u
                self.resolventForcingIndices = self._FEMSpaces.VMixed.sub(index_u).collapse()[1]
            else:
                logger.error(f"Forcing norm type '{param.IOResolvent.ForcingNorm}' not implemented.")
                raise Exception(f"Forcing norm type '{param.IOResolvent.ForcingNorm}' not implemented.")

    def getLinearOperator(self, meanFlow):
        """
        Construct the linear operator matrix for the equation system.

        This method assembles the linear part of the weak formulation 
        by iterating through all equations in the equation list and 
        adding their linear expressions.

        Parameters
        ----------
        meanFlow : object
            Mean flow properties used in constructing the linear operator

        Returns
        -------
        petsc4py.PETSc.Mat
            Assembled linear operator matrix with boundary conditions applied
        """
        # create ufl object with the linear equation system 
        A_ufl = UflDecorator()
        for equation in self.equationList:
            equation.addLinearExpression(A_ufl, meanFlow)

        # assemble petsc matrix
        return A_ufl.getAssembledMatrix(self._mesh, bcs=self.BCs)


    def getWeightMatrix(self, meanFlow):
        """
        Construct the weight matrix representing the time derivative term.  

        Assembles the matrix by summing contributions from all equations.
        Boundary conditions are included for resolvent analysis mode only.  

        Parameters
        ----------
        meanFlow : object
            Mean flow properties used in constructing the weight matrix.    

        Returns
        -------
        petsc4py.PETSc.Mat
            Assembled weight matrix.
       
        Notes
        -----
        - For Resolvent analysis, boundary conditions are applied
        - For other analysis modes, no boundary conditions are applied to 
          avoid computational issues during eigenvalue problem solving
        """
        # create ufl object with the weight matrix expression ("time derivative")
        B_ufl = UflDecorator()
        for equation in self.equationList:
            equation.addWeightMatrixExpression(B_ufl, meanFlow)

        # assemble petsc matrix
        if self._param.Case.AnalysisMode in ['Resolvent']:
            return B_ufl.getAssembledMatrix(self._mesh, bcs=self.BCs)
        else:
            return B_ufl.getAssembledMatrix(self._mesh, bcs=[]) # no boundaries applied, else there is a bug when computing the eigenvalue problem


    def getFEMWeightMatrix(self):
        """
        Construct the full Finite Element Method (FEM) weight matrix.

        This method creates a weight matrix by performing inner products 
        of test and trial functions across all function spaces in the 
        mixed function space.

        Returns
        -------
        petsc4py.PETSc.Mat
            Assembled FEM weight matrix without boundary conditions
        
        Notes
        -----
        - Handles both scalar and vector function spaces
        - Uses conjugate of test functions for matrix construction
        - Includes a workaround for mesh object compatibility with different versions of DOLFINx
        """
        # create ufl object with the full FEM weight matrix expression
        W_ufl     = UflDecorator()
        test_FEM  = self.testFunctionsFEM
        trial_FEM = self.trialFunctionsFEM

        # go through all (scalar) function spaces in the mixed function space 
        i=0
        for test in test_FEM:
            try:    # try if the function space is a "VectorFunctionSpace"
                j=0
                for subTest in test:
                    W_ufl.add(conj(subTest)*trial_FEM[i][j]*dx)
                    j+=1

            except: # function space is scalar
                W_ufl.add(conj(test)*trial_FEM[i]*dx)
            i+=1

        # assemble petsc matrix
        return W_ufl.getAssembledMatrix(self._mesh)

  
    # TODO Sophie: This is a very quick implementation. Tensor framework needed!
    def getFEMDiffusionMatrix(self, diffusionFactor, sponge=None):
        """
        Construct the Finite Element Method (FEM) diffusion matrix.

        This method creates a diffusion matrix by computing the weak form 
        of the diffusion operator across all function spaces in the mixed 
        function space.

        Parameters
        ----------
        diffusionFactor : float
            Coefficient multiplying the diffusion terms
        sponge : float or None, optional
            Additional damping term to be applied to the matrix (default is None)

        Returns
        -------
        petsc4py.PETSc.Mat
            Assembled FEM diffusion matrix with boundary conditions applied
        
        Notes
        -----
        - Handles both scalar and vector function spaces
        - Computes second-order derivatives in x and y directions
        - Currently a quick implementation; a more comprehensive tensor framework is needed for future improvements

        """
        # TODO: This implementation is considered a temporary solution and 
        # requires a more robust tensor framework in future iterations.
        # Rest of the existing implementation remains unchanged
        D_ufl = UflDecorator()
        test_FEM = self.testFunctionsFEM
        trial_FEM = self.trialFunctionsFEM

        i = 0
        for test in test_FEM:
            try:
                j = 0
                for subTest in test:
                    D_ufl.add(conj(subTest)*trial_FEM[i][j]*dx)
                    D_ufl.add(diffusionFactor*(Dx(conj(subTest),0)*Dx(trial_FEM[i][j],0)+Dx(conj(subTest),1)*Dx(trial_FEM[i][j],1))*dx)
                    if sponge is not None:
                        D_ufl.add(sponge*conj(subTest)*trial_FEM[i][j]*dx)
                    j += 1
            except: # function space is scalar
                D_ufl.add(conj(test)*trial_FEM[i]*dx)
                D_ufl.add(diffusionFactor*(Dx(conj(test),0)*Dx(trial_FEM[i],0)+Dx(conj(test),1)*Dx(trial_FEM[i],1))*dx)
                if sponge is not None:
                    D_ufl.add(sponge*conj(test)*trial_FEM[i]*dx)
            i += 1

        # Assemble petsc matrix
        return D_ufl.getAssembledMatrix(self._mesh, self.BCs)


    def getFullRHS(self, func):
        """
        Construct the full FEM right-hand side (RHS) vector.    

        Forms the RHS by integrating each test function against the 
        corresponding part of the provided function.    

        Parameters
        ----------
        func : dolfinx.fem.Function
            Mixed function representing the source term or solution.    

        Returns
        -------
        petsc4py.PETSc.Vec
            Assembled RHS vector.
        """

        # create ufl object with the full FEM weight matrix expression
        rhs_ufl   = UflDecorator()
        test_FEM  = self.testFunctionsFEM
        func_i    = func.split()

        # go through all (scalar) function spaces in the mixed function space 
        i=0
        for test in test_FEM:
            try:    # try if the function space is a "VectorFunctionSpace"
                j=0
                for subTest in test:
                    rhs_ufl.add(conj(subTest)*func_i[i][j]*dx)
                    j+=1

            except: # function space is scalar
                rhs_ufl.add(conj(test)*func_i[i]*dx)
            i+=1

        # assemble petsc matrix
        return rhs_ufl.getAssembledVector(self._mesh)

    def getNonlinearExpression(self, meanFlow, setBC=True):
        """
        Construct the nonlinear vector expression for the system.

        Assembles the nonlinear contributions from all equations.

        Parameters
        ----------
        meanFlow : object
            Mean flow field used in the nonlinear formulation.
        setBC : bool, optional
            Whether to apply boundary conditions (default is True).

        Returns
        -------
        petsc4py.PETSc.Vec
            Assembled nonlinear vector.
        """

        N_ufl = UflDecorator()
        for equation in self.equationList:
            equation.addNonlinearExpression(N_ufl, meanFlow)
        if setBC:
            return N_ufl.getAssembledVector(self._mesh, self.BCs)
        else:
            return N_ufl.getAssembledVector(self._mesh)

    def getBilinearOperator(self, meanFlow):
        """
        Construct the bilinear operator matrix.

        Assembles all bilinear terms across the system equations.

        Parameters
        ----------
        meanFlow : object
            Mean flow field used in the bilinear formulation.

        Returns
        -------
        petsc4py.PETSc.Mat
            Assembled bilinear operator matrix with boundary conditions applied.
        """

        # create ufl object with the linear equation system 
        BL_ufl = UflDecorator()
        for equation in self.equationList:
            equation.addBilinearExpression(BL_ufl, meanFlow)
        # assemble petsc matrix
        return BL_ufl.getAssembledMatrix(self._mesh, bcs=self.BCs)

    def getForcingForInputOutput(self,meanFlow):
        """
        Construct the forcing vector for input-output analysis.

        Creates ufl object with the linear equation system. 

        Parameters
        ----------
        meanFlow : object
            Mean flow properties used to construct the forcing vector.

        Returns
        -------
        petsc4py.PETSc.Vec
            Assembled complex forcing vector.
        """

        f_ufl = UflDecorator()
        for equation in self.equationList:
            equation.addLinearExpression(f_ufl, meanFlow)
        # assemble forcing vector
        forcing = f_ufl.getAssembledVector(self._mesh, self.BCs)
        forcing.scale(-1j)
        return forcing


    def getResolventNorm_response(self, meanFlow):
        """
        Construct the resolvent response norm matrix.

        Parameters
        ----------
        meanFlow : object
            Mean flow field. Currently not used directly here, but kept for
            interface consistency with other assembly routines.

        Returns
        -------
        petsc4py.PETSc.Mat
            PETSc matrix representing the response norm, assembled from the
            UFL form stored in ``self.response_vf``.
        """
        # assemble petsc matrix
        return self.response_vf.getAssembledMatrix(self._mesh)


    def getResolventNorm_forcing(self, meanFlow):
        """
        Construct the resolvent forcing norm matrix.

        Parameters
        ----------
        meanFlow : object
            Mean flow field. Currently not used directly here, but kept for
            interface consistency with other assembly routines.

        Returns
        -------
        petsc4py.PETSc.Mat
            PETSc matrix representing the forcing norm, assembled from the
            UFL form stored in ``self.forcing_vf``.
        """

        # assemble petsc matrix
        return self.forcing_vf.getAssembledMatrix(self._mesh)


    def getResolventWeighting_FEM(self, meanFlow):
        """
        Assemble the FEM weighting matrix for resolvent analysis.

        The underlying UFL form for the weighting is constructed beforehand
        by :meth:`computeResolventFEMWeights` and stored in
        :attr:`self.fem_weighting`. This method only assembles the
        corresponding PETSc matrix.

        Parameters
        ----------
        meanFlow : object
            Mean flow field. Currently not used directly here, but included
            for API compatibility with other assembly routines.

        Returns
        -------
        petsc4py.PETSc.Mat
            PETSc matrix with FEM weights based on the previously defined
            variable projections.
        """

        # create ufl object with the linear equation system 
        W_FEM_ufl = UflDecorator()

        # assemble petsc matrix
        return self.fem_weighting.getAssembledMatrix(self._mesh)


    ########################### Resolvent Norm  ############################
    def computeResolventNorms(self, X, param, mean, fluc):
        """
        Compute energy norms for resolvent forcing and response.

        Parameters
        ----------
        X : list of Tensor
            Test functions with tensor awareness.
        param : object
            Simulation configuration parameters.
        mean : object
            Mean flow fields.
        fluc : object
            Fluctuation fields object.

        Notes
        -----
        Only 'TKE' and 'Chu' norm types are currently implemented.
        """

        #Initialize forcing and response
        self.forcing_vf     = UflDecorator()
        self.response_vf    = UflDecorator()

        # If the density field is inhomogeneous, the mean density
        # field must be taken into account, if not it is set to 1
        barrho = mean.rho

        #In body forcing, forcing is allowed in the entire domain (later restricted by P matrix)
        if param.IOResolvent.ForcingMode == 'Body':
            # Loop through forcing coefficients (The coefficients that are chosen by the user,
            # corresponding to the respective equations)

            if param.IOResolvent.ForcingNorm == 'Chu':
                logger.debug("Using Chu's disturbance energy (rho-T) for forcing norm.")
                idu   = param.Case.SolutionList.index('u')
                idrho = param.Case.SolutionList.index('rho')
                idT   = param.Case.SolutionList.index('T')
                self.forcing_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx     # TKE term
                self.forcing_vf += (mean.R_spe*mean.T/mean.rho * fluc.rho*iConj(X[idrho])).ufl_tens*self._coordinateSystem.J_hat*dx     # density term
                self.forcing_vf += (mean.rho*mean.cp/(mean.T*mean.gamma) * fluc.T*iConj(X[idT])).ufl_tens*self._coordinateSystem.J_hat*dx       # Temperature term
            elif param.IOResolvent.ForcingNorm == 'TKE':
                logger.debug("Using TKE energy for forcing norm.")
                idu   = param.Case.SolutionList.index('u')
                self.forcing_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx
                self.__forcing_coeff = self._param.IOResolvent.ForcingCoeff
            else:
                logger.error(f"Forcing norm type '{param.IOResolvent.ForcingNorm}' not implemented.")
                raise Exception(f"Forcing norm type '{param.IOResolvent.ForcingNorm}' not implemented.")
                
        # In boundary forcing, forcing is allowed only on the specific boundaries
        elif param.IOResolvent.ForcingMode=='Boundary':
            raise Exception("Boundary forcing not implemented for Resolvent analysis in Tensor notation")

        # Same for the response
        if param.IOResolvent.ResponseNorm == 'Chu':
            logger.debug("Using Chu's disturbance energy (rho-T) for response norm.")
            idu =  param.Case.SolutionList.index('u')
            idrho = param.Case.SolutionList.index('rho')
            idT = param.Case.SolutionList.index('T')
            self.response_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx     # TKE term
            self.response_vf += (mean.R_spe*mean.T/mean.rho * fluc.rho*iConj(X[idrho])).ufl_tens*self._coordinateSystem.J_hat*dx     # density term
            self.response_vf += (mean.rho*mean.cp/(mean.T*mean.gamma) * fluc.T*iConj(X[idT])).ufl_tens*self._coordinateSystem.J_hat*dx       # Temperature term
        elif param.IOResolvent.ResponseNorm == 'TKE':
            logger.debug("Using TKE energy for response norm.")
            idu =  param.Case.SolutionList.index('u')
            self.response_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx
        else:
            logger.error(f"Response norm type '{param.IOResolvent.ResponseNorm}' not implemented.")
            raise Exception(f"Response norm type '{param.IOResolvent.ResponseNorm}' not implemented.")

    def computeResolventFEMWeights(self, X, param, mean, fluc):
        """
        Compute FEM weighting matrix for resolvent input/output scaling.

        Parameters
        ----------
        X : list of Tensor
            Test functions.
        param : object
            Parameter configuration.
        mean : object
            Mean flow fields.
        fluc : object
            Fluctuation fields.
        """
        
        # Initialize the matrix
        self.fem_weighting = UflDecorator()
        J_hat = self._coordinateSystem.J_hat
        
        # Loop over all eqs, and multiply fluctuation
        # with corresponding test function
        for eqID in param.Case.SetOfEquations :
            if not (param.Case.SetOfEquations[eqID]['Equation'] == 'None' or \
                param.Case.SetOfEquations[eqID]['Variable'] == 'None'):
                
                #  Get the corresponding variable and its index in X
                varName = param.Case.SetOfEquations[eqID]['Variable']
                varIndex = param.Case.SolutionList.index(varName)
                
                # Dynamically get the corresponding fluctuation field
                fluc_var = getattr(fluc, '%s' % varName)
                
                # Multiply by corresponding X*
                if varName == 'u': # For u we need the dot product with X
                    self.fem_weighting += ( iDot(fluc_var, iConj(X[varIndex])) ).ufl_tens*J_hat*dx
                else:
                    self.fem_weighting += ( fluc_var*iConj(X[varIndex]) ).ufl_tens*J_hat*dx
        


    def getRestrictorMatResponse(self):
        """
        Construct the response restrictor matrix for the resolvent analysis.

        This matrix is diagonal and applies spatial restriction (if defined) to the 
        response vector based on the `responseDomain` field in the mean flow.

        Returns
        -------
        petsc4py.PETSc.Mat
            PETSc matrix with diagonal entries corresponding to the spatial restriction mask.

        Notes
        -----
        - Provides a quadratic matrix, with the size of the solution space (VMixed).
        - Has the response restrictor values, given with the mean field, on the diagonal.
        """

        from petsc4py import PETSc
        
        # First we check if response Domain is zero everywhere = no spatial limiter
        if max(self._mean.getVertexValues().responseDomain, key=abs) == 0:
            responseRestrictor_scalarP2             = Function(self._FEMSpaces.P2)
            responseRestrictor_scalarP1             = Function(self._FEMSpaces.P1)
            responseRestrictor_scalarP2.x.array[:]  = 1. # Setting 1 to everywhere
            responseRestrictor_scalarP1.x.array[:]  = 1. # Setting 1 to everywhere
            logger.debug('No spatial restriction of Resolvent response.')
        else:
            responseRestrictor_scalarP2             = self._mean.responseDomain.function    # using actual values
            responseRestrictor_scalarP1             = Function(self._FEMSpaces.P1)
            responseRestrictor_scalarP1.interpolate(responseRestrictor_scalarP2)
            logger.debug('Imposing spatial restriction of Resolvent response.')

        # Crude way to go over all scalar spaces and get their indices
        # (some of them are in the vector space for the velocity) 
        notFinished     = True
        i               = 0
        j               = 0
        indices         = []
        
        while notFinished:
            try:
                indices.append(self._FEMSpaces.VMixed.sub(i).sub(j).collapse()[1])
                j +=1
            except:
                j = 0
                i += 1 
                try:
                    indices.append(self._FEMSpaces.VMixed.sub(i).collapse()[1]) 
                except:
                    notFinished = False

        # create quadratic petsc matrix and fill it with the restrictor values
        range_all = self._FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed
        m         = len(np.arange(*range_all))
        array     = np.empty(m,dtype=complex)
        for index in indices:
            if len(index) == len(responseRestrictor_scalarP1.x.array[:]):
                array[index] = responseRestrictor_scalarP1.x.array[:]
            elif len(index) == len(responseRestrictor_scalarP2.x.array[:]):
                array[index] = responseRestrictor_scalarP2.x.array[:]

        vec_diag = PETSc.Vec().createSeq(m)
        vec_diag.setValues(np.arange(m,dtype=np.int32),array[:])
        vec_diag.assemble()

        P_petsc = PETSc.Mat().createAIJ([m,m])
        P_petsc.setUp()
        P_petsc.setDiagonal(vec_diag)
        P_petsc.assemble()

        return P_petsc

    def getRestrictorMatForcing(self):
        """
        Construct the forcing restrictor matrix for the resolvent analysis.

        This matrix is diagonal and applies spatial restriction (if defined)
        to the forcing vector based on the ``forcingDomain`` field in the
        mean flow. If the forcing domain is zero everywhere, the matrix
        reduces to the identity.

        Returns
        -------
        petsc4py.PETSc.Mat
            PETSc matrix with diagonal entries corresponding to the spatial
            restriction mask.

        Notes
        -----
        - The matrix is quadratic with the size of the solution space
          (``VMixed``).
        - The diagonal entries are given by the forcing restrictor values
          provided with the mean field (or ones everywhere if no spatial
          restriction is prescribed).
        """

        from petsc4py import PETSc
        
        # First we check if forcingDomain is zero everywhere = no spatial limiter
        if max(self._mean.getVertexValues().forcingDomain, key=abs) == 0:
            forcingRestrictor_scalarP1              = Function(self._FEMSpaces.P1)
            forcingRestrictor_scalarP2              = Function(self._FEMSpaces.P2)
            forcingRestrictor_scalarP1.x.array[:]   = 1. # Setting 1 to everywhere
            forcingRestrictor_scalarP2.x.array[:]   = 1. # Setting 1 to everywhere
            logger.debug('No spatial restriction of Resolvent forcing.')
        else:
            # invert values, if non-zero
            forcingRestrictor_scalarP2              = self._mean.forcingDomain.function    # using actual values
            forcingRestrictor_scalarP1              = Function(self._FEMSpaces.P1)
            forcingRestrictor_scalarP1.interpolate(forcingRestrictor_scalarP2)
            logger.debug('Imposing spatial restriction of Resolvent forcing.')

        # crude way to go over all scalar spaces and get their indices (some of them are in the vector space for the velocity) 
        notFinished = True
        i           = 0
        j           = 0
        indices     = []
        while notFinished:
            try:
                indices.append(self._FEMSpaces.VMixed.sub(i).sub(j).collapse()[1])
                j +=1
            except:
                j = 0
                i += 1 
                try:
                    indices.append(self._FEMSpaces.VMixed.sub(i).collapse()[1]) 
                except:
                    notFinished = False

        # create quadratic petsc matrix and fill it with the restrictor values
        range_all = self._FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed
        m         = len(np.arange(*range_all))
        array     = np.empty(m,dtype=complex)
        for index in indices:
            if len(index) == len(forcingRestrictor_scalarP1.x.array[:]):
                array[index] = forcingRestrictor_scalarP1.x.array[:]
            elif len(index) == len(forcingRestrictor_scalarP2.x.array[:]):
                array[index] = forcingRestrictor_scalarP2.x.array[:]

        vec_diag = PETSc.Vec().createSeq(m)
        vec_diag.setValues(np.arange(m,dtype=np.int32),array[:])
        vec_diag.assemble()

        P_petsc   = PETSc.Mat().createAIJ([m,m])
        P_petsc.setUp()
        P_petsc.setDiagonal(vec_diag)
        P_petsc.assemble()

        return P_petsc

    def getShrinkerMatResponse(self):
        """
        Construct a shrinker matrix for the resolvent response vector.

        Returns
        -------
        petsc4py.PETSc.Mat
            PETSc matrix projecting full response to a reduced subspace.
        """
    
        from petsc4py import PETSc

        indices_response    = self.resolventResponseIndices  # depending on the response norm
        range_all           = self._FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed

        n                   = len(indices_response)
        m                   = len(np.arange(*range_all))
        array               = np.empty(n)
        array[:]            = 1.

        P_petsc             = PETSc.Mat().createAIJ([m,n])
        P_petsc.setUp()
        #P_petsc.setValues(indices_response, np.arange(n,dtype=np.int32), array)
        j                   = 0
        for i in indices_response:
            P_petsc.setValue(i,j,1.)
            j               += 1
        P_petsc.assemble()

        return P_petsc

    def getShrinkerMatForcing(self):
        """
        Creates a (possibly rectangular) matrix which serves the purpose to "shrink" the forcing vector to the requested size 
        (e.g. only consider the velocity components when using the forcing norm "TKE").

        Returns
        -------
        petsc4py.PETSc.Mat
            PETSc matrix projecting full forcing vector to reduced subspace.
        """
    
        from petsc4py import PETSc

        indices_forcing     = self.resolventForcingIndices  # depending on the forcing norm
        range_all           = self._FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed

        n                   = len(indices_forcing)
        m                   = len(np.arange(*range_all))

        P_petsc             = PETSc.Mat().createAIJ([m,n])
        P_petsc.setUp()
        j                   = 0
        for i in indices_forcing:
            P_petsc.setValue(i,j,1.)
            j               += 1
        P_petsc.assemble()

        return P_petsc



