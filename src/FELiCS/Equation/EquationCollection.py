from dolfinx import __version__


from ufl import (
                i,
                j,
                Dx,
                TestFunctions,
                TrialFunctions,
                dx,
                SpatialCoordinate,
                FacetNormal,
                as_tensor,
                Measure,
                rhs,
                lhs,
                conj,
                )

from dolfinx.fem import (
                Function,
                dirichletbc,
                Constant,
                form,
                locate_dofs_topological,
)
from dolfinx.fem.petsc import (
                assemble_matrix,
                assemble_vector,

)
from dolfinx.cpp.la import (
                SparsityPattern,

)
from dolfinx.cpp.la.petsc import (
                create_matrix,
)

# SLEPc modification
from petsc4py.PETSc import ScalarType
import pdb

import FELiCS.Solvers.LinearSystem as LinearSystem
from   FELiCS.Misc.functions import *
from   FELiCS.Misc.tensorUtils import (
    Tensor,
    as_vector,
    iInner,
    iDot,
    iConj,
)
from .WeakForm import WeakForm

class EquationCollectionClass():
    """
    A class to build variational formulations for relevant matrices in a computational fluid dynamics context.

    This class constructs weak formulations for various matrices including:
    - A matrix (real and imaginary parts)
    - B matrix (real and imaginary parts)
    - B_forcing (Resolvent forcing norm)
    - B_response (Resolvent response norm)

    The B matrix (time derivative) is always positive and real.

    Attributes
    ----------
    m : int
        Crossstreamwise wave number
    x : ufl.SpatialCoordinate
        Spatial coordinates of the mesh
    boundaries : object
        Facet tags for defining boundary conditions
    ds : ufl.Measure
        Measure for surface integration with boundary data
    n_BC : ufl.FacetNormal
        Boundary normals
    n : Tensor
        Boundary normals as a tensor
    BCs : list
        List of Dirichlet boundary conditions
    trialFunctionsFEM : tuple
        Trial functions for finite element method
    testFunctionsFEM : list
        Test functions for finite element method
    R : float or Constant
        Radial coordinate (1.0 for non-cylindrical coordinates)
    A_vf : WeakForm
        Weak form for the A matrix
    B_vf : WeakForm
        Weak form for the B matrix
    equationList : list
        List of equations to be solved

    Parameters
    ----------
    param : object
        Parameter object containing simulation configuration
    FEMSpaces : object
        Finite element method spaces
    mean : object
        Mean flow properties
    mesh : object
        Computational mesh
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
            Configuration parameters for the simulation
        FEMSpaces : object
            Finite element method spaces
        mean : object
            Mean flow properties
        mesh : object
            Computational mesh
        """
        from itertools import compress
        from FELiCS.Fields.fluctuationClass import fluctuationClass
        from FELiCS.Misc.tensorUtils import (
            CoordinateSystem,
            )

        printDebug(True,'--------------------------------')
        printDebug(True,'-- Initializing the equations...')

        # Store input parameters as private attributes
        self.__param = param
        self.__FEMSpaces = FEMSpaces
        self.__mean = mean
        self.__mesh = mesh

        # Get crossstreamwise wave number
        self.m = self.__param.Case.m
               
        # Get spatial coordinates
        self.x = SpatialCoordinate(mesh)
        self._coordinateSystem = mesh.coordinateSystem

        ## BOUNDARIES
        # Get class for integrating along boundaries
        self.boundaries = mesh.facet_tags
        self.ds = Measure("ds", subdomain_data=self.boundaries)

        # Get all boundaries (So far hard coded)
        first_BC_flag=True
        for Boundary in self.__param.BCs.getBCsDict()[list(self.__param.BCs.getBCsDict().keys())[0]]:
            if first_BC_flag:
                self.all_ds  = self.ds(Boundary['ID'])
                first_BC_flag = False
            else:
                self.all_ds += self.ds(Boundary['ID'])
 
        # Get boundary normals
        self.n_BC=FacetNormal(self.__FEMSpaces.P2.mesh)
        self.n = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), self._coordinateSystem)

        # initialize Dirichlet boundary conditions
        self.BCs = self.__getListOfDirichletBCs()

        ## TEST AND TRIAL FUNCTIONS
        # Define test and trial functions
        fluctuationC = fluctuationClass(
                                   param,
                                   mean,
                                   FEMSpaces,
                                   self._coordinateSystem,
                                   )
        self.trialFunctionsFEM = TrialFunctions(FEMSpaces.VMixed)
        self.hat=fluctuationC.fluc
        self.fluctuationC = fluctuationC

        fluc={}
        for sol in self.__param.SolutionList:
            fluc[sol]=self.hat[self.__param.SolutionList.index(sol)]
        
        XTemp = TestFunctions(self.__FEMSpaces.VMixed)
        self.testFunctionsFEM = XTemp
        
        X=[]
        for i in list(XTemp):
            X.append(Tensor(
                i,
                self._coordinateSystem,
                containsTestFunction=True,
                ))
        self.X = X
            
        # Get radial coordinate
        if self.__param.Case.CoordinateSystem in ['Cylindrical']:
            self.R = self.x[1]
        else:
            from petsc4py import PETSc
            self.R=Constant(self.__FEMSpaces.P2.mesh, PETSc.ScalarType(1.0))

        ## Initialize variational formulations
        self.A_vf = WeakForm()
        self.B_vf = WeakForm()
        printDebug(True, '-- Primary fluctuations: %s.' % param.SolutionList)
 
        ## create equation list from parameters
        self.equationList = []

        if self.__param.Case.SetOfEquations['Momentum']['Equation'] == 'NSPrimitive':
            from FELiCS.Equation.Equations.MomentumEquation import MomentumEquation
            printDebug(True, '-- Adding momentum equation for u-fluc -> X[0].')     # Hardcoded u' for mom eq.

            momentum = MomentumEquation(self,fluctuationC,X[0],param)
            self.equationList.append(momentum)

        if self.__param.Case.SetOfEquations['Mass']['Equation'] == 'Continuity':
            from FELiCS.Equation.Equations.MassEquation import MassEquation
            varEq = self.__param.Case.SetOfEquations['Mass']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding mass-balance equation for %s-fluc -> X[%d].' % (varEq,idVar))

            mass = MassEquation(self,fluctuationC,X[idVar],self.__param)
            self.equationList.append(mass)

        if self.__param.Case.SetOfEquations['Energy']['Equation'] == 'Enthalpy':
            from FELiCS.Equation.Equations.EnthalpyEquation import EnthalpyEquation
            varEq = self.__param.Case.SetOfEquations['Energy']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding enthalpy-energy equation for %s-fluc -> X[%d].' % (varEq,idVar))

            enthalpy = EnthalpyEquation(self,fluctuationC,X[idVar],self.__param)
            self.equationList.append(enthalpy)
        ### CODE BLOCK 1
        if self.__param.Case.SetOfEquations['Energy']['Equation'] == 'primitive-p':
            # printError('Energy equation in primitive form is not ready to use!!! Ask Simon Demange for updates.')
            from FELiCS.Equation.Equations.EnergyPressureEquation import EnergyPressureEquation
            varEq = self.__param.Case.SetOfEquations['Energy']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding pressure-energy equation for %s-fluc -> X[%d].' % (varEq,idVar))

            energyP = EnergyPressureEquation(self,fluctuationC,X[idVar],self.__param)
            self.equationList.append(energyP)

            #energyP.addLinearExpression(self.A_vf,mean)
            #energyP.addWeightMatrixExpression(self.B_vf,mean)
            
        # Add sponge region only if the field was given in the mean flow file
        if not('spg' in mean._meanFlowClass__notInFileList):
            from FELiCS.Equation.Equations.SpongeTerm import SpongeTerm
            printDebug(True, '-- Adding sponge damping.')

            sponge = SpongeTerm(self,fluctuationC,X,self.__param)
            self.equationList.append(sponge)

            #sponge.addLinearExpression(self.A_vf,mean)
            #sponge.addWeightMatrixExpression(self.B_vf,mean)

        # Add species transport equation for all transported species
        transportedSpecies=self.__param.Case.Mixture.getSpeciesList('transported')
        for specie in transportedSpecies:
            i_eqn=self.__param.SolutionList.index(specie)
            if self.__param.Case.SetOfEquations['Species']['Equation'] == 'Non-conservative':
                from FELiCS.Equation.Equations.SpeciesEquation import SpeciesEquation
                print('-- Adding equation for species '+specie + ' in non-conservative form')

                species = SpeciesEquation(self,fluctuationC,X[i_eqn],specie,self.__param)
                self.equationList.append(species)

                #species.addLinearExpression(self.A_vf,mean)
                #species.addWeightMatrixExpression(self.B_vf,mean)
                
            #elif self.__param.Case.SetOfEquations['Species']['Equation'] == 'Conservative':
            #    # This eq has not been derived in tensor framework yet.
            #    from FELiCS.Equation.speciesConservative.addSpeciesConservativeEq import addSpeciesConservativeEq
            #    print('-- Adding equation for species '+specie +' in conservative form')
            #    addSpeciesConservativeEq(self,fluctuationC,X[i_eqn],self.mean,specie,self.__param)
            else:
                raise Exception('Species transport equation type ' + self.__param.Case.SetOfEquations['Species']['Equation'] + ' unknown' )

        if self.__param.Case.AnalysisMode in ['Resolvent']:
            self.computeResolventNorms     (X,self.__param,mean,fluctuationC)
            self.computeResolventFEMWeights(X,self.__param,mean,fluctuationC)

            # get indices for forcing and response, depending on used norm, to use when creating the shrinker matrices
            index_u =  param.SolutionList.index('u')
            if param.IOResolvent.ResponseNorm == 'Chu':
                # TODO: (next step) initialize the indices with the name of the variables! Here: all are used, hard-coded, as a quick fix for Simon
                #index_rho = param.SolutionList.index('rho')
                #index_T   = param.SolutionList.index('T')
                ## add up index lists 
                #size = 0
                #for i in indexList:
                #    size += len(i)
                #index = np.zeros(size, dtype=int)
                #scalarSize = len(indexList[0])
                #for i in range(len(indexList)):
                #    index[i*scalarSize:(i+1)*scalarSize] = indexList[i][:]

                # u, rho, p
                self.resolventResponseIndices = np.arange(self.__FEMSpaces.VMixed.dofmap.index_map.local_range[1]) # whole size of VMixed
            elif param.IOResolvent.ResponseNorm == 'TKE':
                # u
                self.resolventResponseIndices = self.__FEMSpaces.VMixed.sub(index_u).collapse()[1]
            else:
                print("ERROR")
                exit()
            if param.IOResolvent.ForcingNorm == 'Chu':
                # u, rho, p
                self.resolventForcingIndices = np.arange(self.__FEMSpaces.VMixed.dofmap.index_map.local_range[1]) # whole size of VMixed
            elif param.IOResolvent.ForcingNorm == 'TKE':
                # u
                self.resolventForcingIndices = self.__FEMSpaces.VMixed.sub(index_u).collapse()[1]
            else:
                print("ERROR")
                exit()

            # TODO: raise Error if the norm is set with a wrong value!!!

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
            
            Notes
            -----
            Includes a workaround for mesh object compatibility with different 
            versions of DOLFINx.
            """
            # create ufl object with the linear equation system 
            A_ufl = WeakForm()
            for equation in self.equationList:
                equation.addLinearExpression(A_ufl, meanFlow)

            # Mesh compatibility workaround
            try:
                A_ufl.setCorrectMeshObject(self.__mesh)
            except:
                printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")

            # assemble petsc matrix
            A = assemble_matrix(form(A_ufl.lhs), bcs=self.BCs)
            A.assemble()

            return A


    def getWeightMatrix(self, meanFlow):
        """
        Construct the weight matrix representing the time derivative term.

        This method assembles the weight matrix by iterating through all 
        equations in the equation list and adding their weight matrix 
        expressions.

        Parameters
        ----------
        meanFlow : object
            Mean flow properties used in constructing the weight matrix

        Returns
        -------
        petsc4py.PETSc.Mat
            Assembled weight matrix, with or without boundary conditions 
            depending on the analysis mode
        
        Notes
        -----
        - For Resolvent analysis, boundary conditions are applied
        - For other analysis modes, no boundary conditions are applied to 
          avoid computational issues during eigenvalue problem solving
        - Includes a workaround for mesh object compatibility with different 
        versions of DOLFINx
        """
        # create ufl object with the weight matrix expression ("time derivative")
        B_ufl = WeakForm()
        for equation in self.equationList:
            equation.addWeightMatrixExpression(B_ufl, meanFlow)

        # Mesh compatibility workaround
        try:
            B_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")

        # assemble petsc matrix
        if self.__param.Case.AnalysisMode in ['Resolvent']:
            B = assemble_matrix(form(B_ufl.lhs), bcs=self.BCs)
        else:
            B = assemble_matrix(form(B_ufl.lhs), bcs=[]) # no boundaries applied, else there is a but when computing the eigenvalue problem
        B.assemble()

        return B


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
        - Includes a workaround for mesh object compatibility with different 
        versions of DOLFINx
        """
        # create ufl object with the full FEM weight matrix expression
        W_ufl = WeakForm()
        test_FEM = self.testFunctionsFEM
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

        # Mesh compatibility workaround
        try:
            W_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")

        # assemble petsc matrix
        W = assemble_matrix(form(W_ufl.lhs), [])#self.BCs) #without BCs
        W.assemble()

        return W

  
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
        - Includes a workaround for mesh object compatibility with different 
        versions of DOLFINx
        - Currently a quick implementation; a more comprehensive tensor 
        framework is needed for future improvements

        Warnings
        --------
        TODO: This implementation is considered a temporary solution and 
        requires a more robust tensor framework in future iterations.
        """
        # Rest of the existing implementation remains unchanged
        D_ufl = WeakForm()
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

        # Mesh compatibility workaround
        try:
            D_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")

        # Assemble petsc matrix
        D = assemble_matrix(form(D_ufl.lhs), self.BCs)
        D.assemble()
        return D


    def getFullRHS(self,func):

        # create ufl object with the full FEM weight matrix expression
        rhs_ufl   = WeakForm()
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

        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            rhs_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        # assemble petsc matrix
        rhs = assemble_vector(form(-rhs_ufl.rhs))
        rhs.assemble()

        return rhs



    def getNonlinearExpression(self, meanFlow, setBC=True):
        from dolfinx.fem.petsc import set_bc
        N_ufl = WeakForm()
        for equation in self.equationList:
            equation.addNonlinearExpression(N_ufl, meanFlow)

        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            N_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        N = assemble_vector(form(N_ufl.rhs))
        if setBC:
            set_bc(N, self.BCs)

        return N

    def getBilinearOperator(self, meanFlow):

        # create ufl object with the linear equation system 
        BL_ufl = WeakForm()
        for equation in self.equationList:
            try:
                equation.addBilinearExpression(BL_ufl, meanFlow)
            except:
                pass

        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            BL_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        # assemble petsc matrix
        BL = assemble_matrix(form(BL_ufl.lhs), bcs=self.BCs)
        BL.assemble()

        return BL

    def getForcingForInputOutput(self,meanFlow):
        from dolfinx.fem.petsc import set_bc

        # create ufl object with the linear equation system 
        A_ufl = WeakForm()
        for equation in self.equationList:
            equation.addLinearExpression(A_ufl, meanFlow)
        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            A_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        # assemble forcing vector
        forcing = assemble_vector(form(A_ufl.rhs))
        forcing.assemble()
        set_bc(forcing, self.BCs)
        forcing.scale(1j)

        return forcing


        B_forcing = assemble_matrix(form(self.forcing_vf))
        B_forcing.assemble()
        self.__matrix_dict_petsc['B_forcing'] = B_forcing 
        del B_forcing

        B_response = assemble_matrix(form(self.response_vf))
        B_response.assemble()
        self.__matrix_dict_petsc['B_response'] = B_response 
        del B_response
 


    def getResolventNorm_response(self, meanFlow):

        # create ufl object with the linear equation system 
        W_r_ufl = WeakForm()
        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            W_r_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        # assemble petsc matrix
        W_response = assemble_matrix(form(self.response_vf))
        W_response.assemble()

        return W_response


    def getResolventNorm_forcing(self, meanFlow):

        # create ufl object with the linear equation system 
        W_f_ufl = WeakForm()
        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            W_f_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        # assemble petsc matrix
        W_forcing = assemble_matrix(form(self.forcing_vf))
        W_forcing.assemble()

        return W_forcing


    def getResolventWeighting_FEM(self, meanFlow):

        # create ufl object with the linear equation system 
        W_FEM_ufl = WeakForm()
        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            W_FEM_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        # assemble petsc matrix
        W_FEM = assemble_matrix(form(self.fem_weighting))
        W_FEM.assemble()

        return W_FEM


########################### Resolvent Norm  ############################
    def computeResolventNorms(self,X,param,mean,fluc):
        ''' This function yields the norms for the resolvent analysis.
        Note that both the forcing and response norm must be real!'''
        #Initialize forcing and response
        self.forcing_vf = 0
        self.response_vf = 0

        # If the density field is inhomogeneous, the mean density
        # field must be taken into account, if not it is set to 1
        barrho=mean.rho

        #In body forcing, forcing is allowed in the entire domain (later restricted by P matrix)
        if param.IOResolvent.ForcingMode=='Body':
            # Loop through forcing coefficients (The coefficients that are chosen by the user,
            # corresponding to the respective equations)
            printWarning("  -- Currently only the L2 norm is implemented for both forcing and response in a resolvent analysis. Here, ALL velocity components are taken into account, no matter the choices in the settings file.")
            if param.IOResolvent.ForcingNorm == 'Chu':
                printDebug(True, "-- Using Chu's disturbance energy (rho-T) for forcing norm.")
                idu   = param.SolutionList.index('u')
                idrho = param.SolutionList.index('rho')
                idT   = param.SolutionList.index('T')
                self.forcing_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx     # TKE term
                self.forcing_vf += (mean.R_spe*mean.T/mean.rho * fluc.rho*iConj(X[idrho])).ufl_tens*self._coordinateSystem.J_hat*dx     # density term
                self.forcing_vf += (mean.rho*mean.cp/(mean.T*mean.gamma) * fluc.T*iConj(X[idT])).ufl_tens*self._coordinateSystem.J_hat*dx       # Temperature term
            elif param.IOResolvent.ForcingNorm == 'TKE':
                printDebug(True, "-- Using TKE energy for forcing norm.")
                idu   = param.SolutionList.index('u')
                self.forcing_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx
                self.__forcing_coeff = self.__param.IOResolvent.ForcingCoeff
            #TODO: raise Error!!!
                
        # In boundary forcing, forcing is allowed only on the specific boundaries
        elif param.IOResolvent.ForcingMode=='Boundary':
            raise Exception("Boundary forcing not implemented for Resolvent analysis in Tensor notation")

          
        if param.IOResolvent.ResponseNorm == 'Chu':
            printDebug(True, "-- Using Chu's disturbance energy (rho-T) for response norm.")
            idu =  param.SolutionList.index('u')
            idrho = param.SolutionList.index('rho')
            idT = param.SolutionList.index('T')
            self.response_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx     # TKE term
            self.response_vf += (mean.R_spe*mean.T/mean.rho * fluc.rho*iConj(X[idrho])).ufl_tens*self._coordinateSystem.J_hat*dx     # density term
            self.response_vf += (mean.rho*mean.cp/(mean.T*mean.gamma) * fluc.T*iConj(X[idT])).ufl_tens*self._coordinateSystem.J_hat*dx       # Temperature term
        elif param.IOResolvent.ResponseNorm == 'TKE':
            printDebug(True, "-- Using TKE energy for response norm.")
            idu =  param.SolutionList.index('u')
            self.response_vf += (barrho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx
        #TODO: raise Error!!
         

        ## Prompt variational formulations in debug mode
        #printDebug(param.debug,'-- Resolvent forcing norm is '+ str(self.forcing_vf))
        #printDebug(param.debug,'-- Resolvent response norm is '+ str(self.response_vf))

    def computeResolventFEMWeights(self,X,param,mean,fluc):
        ''' 
        This function yields a matrix containing the FEM weights
        corresponding to a diagonal unit matrix.
        (Need to chat with Sophie to better define this)
        '''
        
        # Initialize the matrix
        self.fem_weighting = 0
        J_hat = self._coordinateSystem.J_hat
        
        # Loop over all eqs, and multiply fluctuation
        # with corresponding test function
        for eqID in param.Case.SetOfEquations :
            if not (param.Case.SetOfEquations[eqID]['Equation'] == 'None' or \
                param.Case.SetOfEquations[eqID]['Variable'] == 'None'):
                
                #  Get the corresponding variable and its index in X
                varName = param.Case.SetOfEquations[eqID]['Variable']
                varIndex = param.SolutionList.index(varName)
                
                # Dynamically get the corresponding fluctuation field
                fluc_var = getattr(fluc, '%s' % varName)
                
                # Multiply by corresponding X*
                if varName == 'u': # For u we need the dot product with X
                    self.fem_weighting += ( iDot(fluc_var, iConj(X[varIndex])) ).ufl_tens*J_hat*dx
                else:
                    self.fem_weighting += ( fluc_var*iConj(X[varIndex]) ).ufl_tens*J_hat*dx

        

    def __getListOfDirichletBCs(self, listOfComponents=None):
        BClist = []
        mesh=self.__FEMSpaces.P2.mesh

        printDebug(True, '--------------------------------')
        printDebug(True, '-- Setting boundary conditions...')
        self.__boundaries = self.__param.BCs.getBoundaries()
        self.__bcDict = self.__param.BCs.getBCsDict()
        VelocityComponents=self.__param.Case.getVelocityComponents()
        SolutionList=self.__param.SolutionList#Case.getTransportedQuantityList()
        for k,m in zip(list(self.__bcDict.keys()),range(0,len(self.__bcDict.keys()))):
            # Get index of equation/variable i_eqn and if needed the index of the velocity component
            if k[0]=='u' and k[1] in VelocityComponents:
                i_eqn=0
                i_component=VelocityComponents.index(k[1])
            elif k in SolutionList:
                i_eqn = SolutionList.index(k)
            else:
                continue

            for mm in range(0,len(self.__bcDict[k])) :
                if self.__bcDict[k][mm]['type']=='Dirichlet':


                    printDebug(self.__param.debug,"-- Adding Dirichlet BC for "+str(k)+ " in equation "+str(i_eqn)+" with value "+str(self.__bcDict[k][mm]['value'])+" on boundary with index "+str(self.__bcDict[k][mm]['ID']))
                    if k in ['u'+ component for component in VelocityComponents] :

                        #if __version__.find('0.4.1') >= 0:
                        BClist.append( dirichletbc(ScalarType(self.__bcDict[k][mm]['value']), locate_dofs_topological(self.__FEMSpaces.VMixed.sub(i_eqn).sub(i_component), 1, self.boundaries.indices[self.boundaries.values==self.__bcDict[k][mm]['ID']]), self.__FEMSpaces.VMixed.sub(i_eqn).sub(i_component)) )

                        #else:
                            #BClist.append( dirichletbc(ScalarType(self.__bcDict[k][mm]['value']), locate_dofs_topological(self.__FEMSpaces.VMixed.sub(i_eqn).sub(i_component), 1, self.boundaries.find(self.__bcDict[k][mm]['ID'])), self.__FEMSpaces.VMixed.sub(i_eqn).sub(i_component)) )

                    else:
                        if len(self.__param.Case.getTransportedQuantityList()) == 1:

                            #if __version__.find('0.4.1') >= 0:
                            BClist.append( dirichletbc(ScalarType(self.__bcDict[k][mm]['value']), locate_dofs_topological(self.__FEMSpaces.FunctionSpaceList[i_eqn], 1, self.boundaries.indices[self.boundaries.values==self.__bcDict[k][mm]['ID']]), self.__FEMSpaces.FunctionSpaceList[i_eqn]) )
                            #else:
                            #    BClist.append( dirichletbc(ScalarType(self.__bcDict[k][mm]['value']), locate_dofs_topological(self.__FEMSpaces.FunctionSpaceList[i_eqn], 1, self.boundaries.find(self.__bcDict[k][mm]['ID'])), self.__FEMSpaces.FunctionSpaceList[i_eqn]) )
                        else:
                            #if __version__.find('0.4.1') >= 0:
                            BClist.append( dirichletbc(ScalarType(self.__bcDict[k][mm]['value']), locate_dofs_topological(self.__FEMSpaces.VMixed.sub(i_eqn), 1, self.boundaries.indices[self.boundaries.values==self.__bcDict[k][mm]['ID']]), self.__FEMSpaces.VMixed.sub(i_eqn)) )
                            #else:
                            #    BClist.append( dirichletbc(ScalarType(self.__bcDict[k][mm]['value']), locate_dofs_topological(self.__FEMSpaces.VMixed.sub(i_eqn), 1, self.boundaries.find(self.__bcDict[k][mm]['ID'])), self.__FEMSpaces.VMixed.sub(i_eqn)) )
        return BClist



    def getRestrictorMatResponse(self):
        # provides a quadratic matrix, with the size of the solution space (VMixed)
        # has the response restrictor values, given with the mean field, on the diagonal
        from petsc4py import PETSc
        
        
        # First we check if forcingDom is zero everywhere = no spatial limiter
        if max(self.__mean.getVertexValues().responseDomain, key=abs) == 0:
            responseRestrictor_scalarP2 = Function(self.__FEMSpaces.P2)
            responseRestrictor_scalarP1 = Function(self.__FEMSpaces.P1)
            responseRestrictor_scalarP2.x.array[:] = 1. # Setting 1 to everywhere
            responseRestrictor_scalarP1.x.array[:] = 1. # Setting 1 to everywhere
            printDebug(True, '-- No spatial restriction of response.')
        else:
            responseRestrictor_scalarP2 = self.__mean.responseDomain    # using actual values
            responseRestrictor_scalarP1 = Function(self.__FEMSpaces.P1)
            responseRestrictor_scalarP1.interpolate(responseRestrictor_scalarP2)
            printDebug(True, '-- Spatial restriction of response.')

        # crude way to go over all scalar spaces and get their indices (some of them are in the vector space for the velocity) 
        notFinished = True
        i = 0; j = 0
        indices = []
        while notFinished:
            try:
                indices.append(self.__FEMSpaces.VMixed.sub(i).sub(j).collapse()[1])
                j +=1
            except:
                j = 0
                i += 1 
                try:
                    indices.append(self.__FEMSpaces.VMixed.sub(i).collapse()[1]) 
                except:
                    notFinished = False


        # create quadratic petsc matrix and fill it with the restrictor values
        range_all = self.__FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed
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

        P_petsc   = PETSc.Mat().createAIJ([m,m])
        P_petsc.setUp()
        P_petsc.setDiagonal(vec_diag)
        P_petsc.assemble()

        return P_petsc


    def getRestrictorMatForcing(self):
        # provides a quadratic matrix, with the size of the solution space (VMixed)
        # has the inverse of the forcing restrictor values, given with the mean field, on the diagonal
        from petsc4py import PETSc
        
        # First we check if forcingDom is zero everywhere = no spatial limiter
        if max(self.__mean.getVertexValues().forcingDomain, key=abs) == 0:
            forcingRestrictor_scalarP1 = Function(self.__FEMSpaces.P1)
            forcingRestrictor_scalarP2 = Function(self.__FEMSpaces.P2)
            forcingRestrictor_scalarP1.x.array[:] = 1. # Setting 1 to everywhere
            forcingRestrictor_scalarP2.x.array[:] = 1. # Setting 1 to everywhere
            printDebug(True, '-- No spatial restriction of forcing.')
        else:
            # invert values, if non-zero
            array = self.__mean.forcingDomain.x.array
            np.where(array[:] != 0., 1./ array[:], 0.)
            self.__mean.forcingDomain.x.array[:] = array[:]
            forcingRestrictor_scalarP2 = self.__mean.forcingDomain    # using actual values
            forcingRestrictor_scalarP1 = Function(self.__FEMSpaces.P1)
            forcingRestrictor_scalarP1.interpolate(forcingRestrictor_scalarP2)
            printDebug(True, '-- Spatial restriction of forcing.')


        # crude way to go over all scalar spaces and get their indices (some of them are in the vector space for the velocity) 
        notFinished = True
        i = 0; j = 0
        indices = []
        while notFinished:
            try:
                indices.append(self.__FEMSpaces.VMixed.sub(i).sub(j).collapse()[1])
                j +=1
            except:
                j = 0
                i += 1 
                try:
                    indices.append(self.__FEMSpaces.VMixed.sub(i).collapse()[1]) 
                except:
                    notFinished = False


        # create quadratic petsc matrix and fill it with the restrictor values
        range_all = self.__FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed
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
        # creates a (possibly rectangular) matrix which serves the purpose to "shrink" the solution (response) vector to the requested size 
        # (e.g. only consider the velocity components when using the response norm "TKE")
    
        from petsc4py import PETSc

        indices_response = self.resolventResponseIndices  # depending on the response norm
        range_all        = self.__FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed

        n        = len(indices_response)
        m        = len(np.arange(*range_all))
        array    = np.empty(n)
        array[:] = 1.

        P_petsc = PETSc.Mat().createAIJ([m,n])
        P_petsc.setUp()
        #P_petsc.setValues(indices_response, np.arange(n,dtype=np.int32), array)
        j = 0
        for i in indices_response:
            P_petsc.setValue(i,j,1.)
            j+=1
        P_petsc.assemble()

        return P_petsc


    def getShrinkerMatForcing(self):
        # creates a (possibly rectangular) matrix which serves the purpose to "shrink" the forcing vector to the requested size 
        # (e.g. only consider the velocity components when using the forcing norm "TKE")
    
        from petsc4py import PETSc

        indices_forcing  = self.resolventForcingIndices  # depending on the forcing norm
        range_all        = self.__FEMSpaces.VMixed.dofmap.index_map.local_range # whole size of VMixed

        n = len(indices_forcing)
        m = len(np.arange(*range_all))

        P_petsc = PETSc.Mat().createAIJ([m,n])
        P_petsc.setUp()
        j = 0
        for i in indices_forcing:
            P_petsc.setValue(i,j,1.)
            j+=1
        P_petsc.assemble()

        return P_petsc


