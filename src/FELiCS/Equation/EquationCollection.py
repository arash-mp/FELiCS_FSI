#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * opermission of the copyright owner, the FLOW group at TU Berlin
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de
# */
'''
# **********************************************************************
# * This file contains the Weak formulation class, which builds all
# * necessary weak formulations for the later use
# *
# * This file created by Thomas L. Kaiser. Significant contributions
# * were made by
# * Mario Casel, transcribing the Low Mach Equations from FreeFEM
# * Chuhan Wang, including adding the enthalpy equation
# *
# ********************
'''

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
    '''This class build the variational formulations for all relevant matrices
    Currently these are:
    -A (imag and real)
    -B (imag and real)
    -B_forcing (Resolvent forcing norm)
    -B_response (Resolvent response norm)
    The convention is such that the B matrix (time derivative) is always positive and real
    '''
    def __init__(self,param,FEMSpaces,mean,mesh):
        #from fenics import Function
        from itertools import compress
        from FELiCS.Fields.fluctuationClass import fluctuationClass
        from FELiCS.Misc.tensorUtils import (
            CoordinateSystem,
            )

        printDebug(True,'--------------------------------')
        printDebug(True,'-- Initializing the equations...')

        # add the parameters of the constructor as attributs of the class to use them in DiscretizeFlow-method:
        self.__param = param
        self.__FEMSpaces = FEMSpaces
        self.__mean = mean

        self.__mesh = mesh

        # Get crossstreamwise wave number
        self.m = self.__param.Case.m
               
        # Get spatial coordinates
        self.x = SpatialCoordinate(mesh)
        self._coordinateSystem = mesh.coordinateSystem
        ## Define tensor coordinate system, we always assume the third dimension to be homogenous
        #if param.Case.CoordinateSystem =='Cartesian':
        #    self._coordinateSystem = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        #elif param.Case.CoordinateSystem =='Cylindrical':
        #    self._coordinateSystem = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        #else:
        #    printError('Coord. syst not yet implemented in tensor framework.')
       # self.coordinateSystem = self.coord_sys
 

        ## BOUNDARIES
        # Get class for integrating along boundaries
        self.boundaries = mesh.facet_tags
        self.ds = Measure("ds", subdomain_data=self.boundaries)

        # Get all boundaries (So far hard coded)
        first_BC_flag=True
        for Boundary in self.__param.BCs.getBCsDict()[list(self.__param.BCs.getBCsDict().keys())[0]]:
            if first_BC_flag:
                self.all_ds = self.ds(Boundary['ID'])
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
        #self.hat=TrialFunctions(FEMSpaces.VMixed)
        self.hat=fluctuationC.fluc
        self.fluctuationC = fluctuationC

        fluc={}
        for sol in self.__param.SolutionList:
            fluc[sol]=self.hat[self.__param.SolutionList.index(sol)]
        
        XTemp = TestFunctions(self.__FEMSpaces.VMixed)
        
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
            #self.ThirdVelCompIndex = self.__param.SolutionList.index('ut')
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

            #momentum.addLinearExpression(self.A_vf,mean)
            #momentum.addWeightMatrixExpression(self.B_vf,mean)
            
        if self.__param.Case.SetOfEquations['Mass']['Equation'] == 'Continuity':
            from FELiCS.Equation.Equations.MassEquation import MassEquation
            varEq = self.__param.Case.SetOfEquations['Mass']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding mass-balance equation for %s-fluc -> X[%d].' % (varEq,idVar))

            mass = MassEquation(self,fluctuationC,X[idVar],self.__param)
            self.equationList.append(mass)

            #mass.addLinearExpression(self.A_vf,mean)
            #mass.addWeightMatrixExpression(self.B_vf,mean)

        if self.__param.Case.SetOfEquations['Energy']['Equation'] == 'Enthalpy':
            from FELiCS.Equation.Equations.EnthalpyEquation import EnthalpyEquation
            varEq = self.__param.Case.SetOfEquations['Energy']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding enthalpy-energy equation for %s-fluc -> X[%d].' % (varEq,idVar))

            enthalpy = EnthalpyEquation(self,fluctuationC,X[idVar],self.__param)
            self.equationList.append(enthalpy)

            #enthalpy.addLinearExpression(self.A_vf,mean)
            #enthalpy.addWeightMatrixExpression(self.B_vf,mean)
        
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

        ## Add reactions
        ## Reaction eqs not derived in tensor framework yet
        #if self.__param.Case.Reaction:
        #    if self.__param.Case.Mixture.ReactionMechanism['type']=='WestbrookDryer_Max':
        #        from FELiCS.Equation.Reactions.GlobalReaction import GlobalReaction
        #        ReactionModelName="WestbrookDryer_Max" #to be put in param
        #        Reaction=GlobalReaction(self.__param.Case.Mixture.ReactionMechanism)
        #        reactionRateMean=Reaction.computeMeanField(self.mean,self.__FEMSpaces.P2)
        #        reactionForm=Reaction.addReaction(self.mean, X, fluctuationC, self.__param.SolutionList)
        #        self.A_vf.add(1j * reactionForm)
        #    elif self.__param.Case.Mixture.ReactionMechanism['type']=='TwoStep':
        #        from FELiCS.Equation.Reactions.TwoStepReaction import TwoStepReaction
        #        ReactionModelName="BFER" #to be put in param
        #        Reaction=TwoStepReaction(ReactionModelName)
        #        Reaction.computeMeanField(MF,self.__FEMSpaces.P2)
        #        Reaction.testM()
        #        reactionForm=Reaction.addReaction(MF, X, fluc, self.__param.SolutionList,self.__FEMSpaces.P2)
        #        self.A_vf.add(1j * reactionForm)
        #    elif self.__param.Case.Mixture.ReactionMechanism['type']=='2S-SM2':
        #        from FELiCS.Equation.Reactions.c2sm2 import C2SM2
        #        ReactionModelName="2S-SM2" #to be put in param
        #        YCH4_lim=0.043*1e-4
        #        #c2=C2SM2(YCH4_lim,2)
        #        c2=self.mean.reaction
        #        self.TR=fluctuationC.T
        #        self.rhoR=fluctuationC.rho
        #        self.YCH4R=fluctuationC.Y('CH4')
        #        self.YO2R=fluctuationC.Y('O2')
        #        self.YCOR=fluctuationC.Y('CO')
        #        self.YCO2R=fluctuationC.Y('CO2')
        #        self.v_eneR=X[self.__param.Case.getTransportedQuantityList().index('rho')]
        #        self.v_YCH4R=X[self.__param.Case.getTransportedQuantityList().index('CH4')]
        #        self.v_YO2R=X[self.__param.Case.getTransportedQuantityList().index('O2')]
        #        self.v_YH2OR=X[self.__param.Case.getTransportedQuantityList().index('H2O')]
        #        self.v_YCOR=X[self.__param.Case.getTransportedQuantityList().index('CO')]
        #        self.v_YCO2R=X[self.__param.Case.getTransportedQuantityList().index('CO2')]
        #        self.dQMean=self.mean.dQ
        #        self.order=2
        #        self.dx=dx

        #        #c2.computeSensitivities(MeanFlow.T,
        #        #                        MeanFlow.rho,
        #        #                        MeanFlow.Y('CH4'),
        #        #                        MeanFlow.Y('CO'),
        #        #                        MeanFlow.Y('O2'),
        #        #                        MeanFlow.Y('CO2'))

        #        self.A_vf.add(1j * -c2.add_source_to_weak_form(self))
        #    elif self.__param.Case.Mixture.ReactionMechanism['type']=='NOx':
        #        reaction=self.mean.reaction
        #        self.v_NO  = X[self.__param.Case.getTransportedQuantityList().index('NO')]
        #        self.v_NO2 = X[self.__param.Case.getTransportedQuantityList().index('NO2')]
        #        self.T = self.mean.T
        #        self.phi = self.mean.phi
        #        self.A_vf.add(1j * -reaction.add_source_to_weak_form(self))

        if self.__param.Case.AnalysisMode in ['Resolvent']:
            self.getResolventNorms(X,self.__param,mean,fluctuationC)
            self.getResolventFEMWeights(X,self.__param,mean,fluctuationC)


#################################################################################

    def getLinearOperator(self, meanFlow):

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

        # assemble petsc matrix
        A = assemble_matrix(form(A_ufl.lhs), bcs=self.BCs)
        A.assemble()

        return A


    def getWeightMatrix(self, meanFlow):

        # create ufl object with the weight matrix expression ("time derivative")
        B_ufl   = WeakForm()
        for equation in self.equationList:
            equation.addWeightMatrixExpression(B_ufl, meanFlow)

        #####################################################################################################
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            B_ufl.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
        #####################################################################################################

        # assemble petsc matrix
        B = assemble_matrix(form(B_ufl.lhs), bcs=self.BCs)
        B.assemble()

        return B


    def getNonlinearExpression(self, meanFlow):
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


########################### Resolvent Norm  ############################
    def getResolventNorms(self,X,param,mean,fluc):
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
            u_f = fluc.u
            self.forcing_vf += (barrho*iDot(u_f,iConj(X[0]))).ufl_tens*self._coordinateSystem.J_hat*dx
            
            # printDebug(True, "-- Using Chu's disturbance energy for forcing norm!")
            # idu = 0
            # idrho = param.SolutionList.index('rho')
            # idp = param.SolutionList.index('p')
            # self.forcing_vf += (mean.rho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx
            # self.forcing_vf += (mean.p/(mean.rho*mean.rho)*mean.gamma/(mean.gamma-1)*\
            #     fluc.rho*iConj(X[idrho])).ufl_tens*self._coordinateSystem.J_hat*dx
            # self.forcing_vf += (fluc.p*iConj(X[idp])/(mean.p*(mean.gamma-1))).ufl_tens*self._coordinateSystem.J_hat*dx
            # self.forcing_vf += (-1*fluc.p*iConj(X[idrho])/(mean.rho*(mean.gamma-1))).ufl_tens*self._coordinateSystem.J_hat*dx
            # self.forcing_vf += (-1*fluc.rho*iConj(X[idp])/(mean.rho*(mean.gamma-1))).ufl_tens*self._coordinateSystem.J_hat*dx
            
            # Below are arbitrary weights used for debugging resolvent considering other norms
            # than the TKE one:
            # idrhoF = param.SolutionList.index('rho')
            # self.forcing_vf += (fluc.rho*iConj(X[idrhoF])).ufl_tens*self._coordinateSystem.J_hat*dx
            # idpF = param.SolutionList.index('p')
            # self.forcing_vf += (fluc.p*iConj(X[idpF])).ufl_tens*self._coordinateSystem.J_hat*dx
            
            #velocityForcingList = [0,0,0]
            #for i in param.IOResolvent.ForcingCoeff:
            #    # In case the coefficient correspionds to a velocity, i.e. i is smaller
            #    # than the number of velocity components, the coefficient must be applied
            #    # to the corresponding (second level) subspace of u, which correspionds to the right
            #    # velocity component. If not, it is applied directly to the first level subspace,
            #    # and the index is corrected by param.nVelocityComponents+1
            #    
            #    if i<param.Case.getNVelocityComponents():
            #        velocityForcingList[i] = self.hat[0][i]
            #        #input(velocityForcingList)
            #        #self.forcing_vf += conj(self.X[0][i])*barrho*self.hat[0][i]*self.R*dx
            #        #self.forcing_vf += conj(self.X[0][i])*self.hat[0][i]*self.R*dx
            #        
            #    else:
            #        self.forcing_vf += conj(self.X[i-param.nVelocityComponents+1])*\
            #            barrho*self.hat[i-param.nVelocityComponents+1]*self.R*dx
            #velocityComponents = Tensor(as_vector(velocityForcingList),self.coord_sys)
            #u_f = Tensor(self.hat[0], self.coord_sys)
            #self.forcing_vf += conj(self.X[0][i])*barrho*self.hat[0][i]*self.R*dx
            ##self.A_vf.add(( 1j*iDot(iGrad(iConj(X), self.m),fluc_rhou) ).ufl_tens*coord.J_hat*dx)
            #velocityComponents = Tensor(as_vector(velocityForcingList),self.coord_sys)
            #self.forcing_vf += barrho *iDot iConj(X)*self.hat[0][i]*self.R*dx
            #self.forcing_vf += temporalVF
        # In boundary forcing, forcing is allowed only on the specific boundaries
        elif param.IOResolvent.ForcingMode=='Boundary':
            raise Exception("Boundary forcing not implemented for Resolvent analysis in Tensor notation")
            # Create integrator for the respective boundaries
#            Ds = ds(subdomain_data=self.boundaries)
#
#            # Loop through forcing coefficients (The coefficients that are chosen by the user,
#            # corresponding to the respective equations)
#            for i in param.IOResolvent.ForcingCoeff:
#                # Loop through the forcing boundaries specified by the user
#                for k in param.IOResolvent.ForcingBoundaryIndices:
#                    # In case the coefficient correspionds to a velocity, i.e. i is smaller
#                    # than the number of velocity components, the coefficient must be applied
#                    # to the corresponding (second level) subspace of u, which correspionds to the right
#                    # velocity component. If not, it is applied directly to the first level subspace,
#                    # and the index is corrected by param.nVelocityComponents+1
#                    if i<param.Case.getNVelocityComponents():
#                        self.forcing_vf += self.X[0][i]*barrho*\
#                            self.hat[0][i]*self.R*Ds(int(k))
#                    else:
#                        self.forcing_vf += self.X[i-param.nVelocityComponents+1]*\
#                            barrho*self.hat[i-param.nVelocityComponents+1]\
#                            *self.R*Ds(int(k))

        # Loop through response coefficients (The coefficients that are chosen by the user,
        # corresponding to the respective solutions to be maximized)
        #velocityResponseList = [0,0,0]
        #for i in param.IOResolvent.ResponseCoeff:
        #    # In case the coefficient correspionds to a velocity, i.e. i is smaller
        #    # than the number of velocity components, the coefficient must be applied
        #    # to the corresponding (second level) subspace of u, which correspionds to the right
        #    # velocity component. If not, it is applied directly to the first level subspace,
        #    # and the index is corrected by param.nVelocityComponents+1
        #    if i<param.Case.getNVelocityComponents():
        #        velocityResponseList[i] = self.hat[0][i]
        #   #     self.response_vf += conj(self.X[0][i])*barrho*self.hat[0][i]*self.R*dx
        #    else:
        #        self.response_vf += conj(self.X[i-param.nVelocityComponents+1])*barrho*\
        #            self.hat[i-param.nVelocityComponents+1]*self.R*dx
        #velocityComponents = Tensor(as_vector(velocityResponseList),self.coord_sys)
        #temporalVF = (barrho*iDot(velocityComponents,iConj(X[0]))).ufl_tens*self.coord_sys.J_hat*dx
        #self.response_vf += temporalVF
        self.response_vf += (barrho*iDot(u_f,iConj(X[0]))).ufl_tens*self._coordinateSystem.J_hat*dx
        
        # printDebug(True, "-- Using Chu's disturbance energy for response norm!")
        # idu = 0
        # idrho = param.SolutionList.index('rho')
        # idp = param.SolutionList.index('p')
        # self.response_vf += (mean.rho*iDot(fluc.u,iConj(X[idu]))).ufl_tens*self._coordinateSystem.J_hat*dx
        # self.response_vf += (mean.p/(mean.rho*mean.rho)*mean.gamma/(mean.gamma-1)*\
        #     fluc.rho*iConj(X[idrho])).ufl_tens*self._coordinateSystem.J_hat*dx
        # self.response_vf += (fluc.p*iConj(X[idp])/(mean.p*(mean.gamma-1))).ufl_tens*self._coordinateSystem.J_hat*dx
        # self.response_vf += (-1*fluc.p*iConj(X[idrho])/(mean.rho*(mean.gamma-1))).ufl_tens*self._coordinateSystem.J_hat*dx
        # self.response_vf += (-1*fluc.rho*iConj(X[idp])/(mean.rho*(mean.gamma-1))).ufl_tens*self._coordinateSystem.J_hat*dx
        
        # Below are arbitrary weights used for debugging resolvent considering other norms
        # than the TKE one:
        # idrhoF = param.SolutionList.index('rho')
        # self.response_vf += (fluc.rho*iConj(X[idrhoF])).ufl_tens*self._coordinateSystem.J_hat*dx
        # idpF = param.SolutionList.index('p')
        # self.response_vf += (fluc.p*iConj(X[idpF])).ufl_tens*self._coordinateSystem.J_hat*dx

        # Prompt variational formulations in debug mode
        printDebug(param.debug,'-- Resolvent forcing norm is '+ str(self.forcing_vf))
        printDebug(param.debug,'-- Resolvent response norm is '+ str(self.response_vf))
        
    def getResolventFEMWeights(self,X,param,mean,fluc):
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

        
    def DiscretizeFlow(self):

        #self.__DiscretizeAndSolve = DiscretizeAndSolve
        #self.__WeakForm = WeakForm
        from copy import deepcopy
        # Check for type of case
        
        # Used for debugging
        # from matspy import spy
        # import matplotlib
        # matplotlib.use('TkAGG')

        AnalysisMode = self.__param.Case.AnalysisMode
        FEMSpaces = self.__FEMSpaces
        #Get the Dirichlet BCs set by the user in a list...
        # BClist= self.__getListOfDirichletBCs()
        bcFunction = Function(self.__FEMSpaces.VMixed)
        bcFunction.x.array[:] = 0.0
        # for bc in BClist:
        #   dofs = bc.dof_indices()[0]
        #   bc_vals = 1.0
        #   bcFunction.x.array[dofs] = bc_vals
        mesh = self.__FEMSpaces.P2.mesh
        pattern = SparsityPattern(mesh.comm, [self.__FEMSpaces.VMixed.dofmap.index_map, self.__FEMSpaces.VMixed.dofmap.index_map],
                                                            [self.__FEMSpaces.VMixed.dofmap.index_map_bs, self.__FEMSpaces.VMixed.dofmap.index_map_bs])
        pattern.insert_diagonal(np.arange(len(bcFunction.vector.array), dtype=np.int32))
        try:
            pattern.finalize()
        except:
            printDeprecatedWarning("dolfinx version is <0.7.0")
            pattern.assemble()
        BC_Diriclet = create_matrix(mesh.comm, pattern)
        BC_Diriclet.setDiagonal(bcFunction.vector)
        BC_Diriclet.assemble()

        bcs= self.BCs
        n_dof=BC_Diriclet.size[0]
		
        # Sophie: This is a weird work-around, because somehow the wrong mesh object is given to the UFL-form 
        # when using a newer version of dolfinx (version >= 0.6.*).
        # I will try and understand why that is (probably has something to do with the class FelicsMesh?), 
        # but for now this works fine. 
        try:
            self.A_vf.setCorrectMeshObject(self.__mesh)
            self.B_vf.setCorrectMeshObject(self.__mesh)
        except:
            printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")


        if not self.A_vf.lhsIsZero():
            A = assemble_matrix(form(self.A_vf.lhs), bcs=bcs)
            A.assemble()
        else:
            A = BC_Diriclet.scale(0.)

        if not self.B_vf.lhsIsZero():
            B = assemble_matrix(form(self.B_vf.lhs), bcs=bcs)
            B.assemble()
        else:
            B = BC_Diriclet.scale(0.)

        self.__matrix_dict_petsc={}  

        if AnalysisMode in ['Input-Output']:
            if not self.A_vf.rhsIsZero():
                forcing_vec_petsc = assemble_vector(form(self.A_vf.rhs))
                forcing_vec_petsc.assemble()
            else:
                forcing_vec_petsc = PETScVector()
                forcing_vec_petsc.init(n_dof)

            from dolfinx.fem.petsc import set_bc
            set_bc(forcing_vec_petsc, bcs)

            #b_forcing = 1j*forcing_vec
            self.__matrix_dict_petsc['b_forcing'] = forcing_vec_petsc.copy()  
            self.__matrix_dict_petsc['b_forcing'].scale(1j)
            del forcing_vec_petsc

        # Get the BCs provided by the user
        # Aplly the BCs
        ## Boundary conditions are applied via penalisation method manually. This is in order to keep the rhs matrix invertible.
        else:

            bcFunction.x.array[:] = 0.0
            for bc in bcs:
                try:
                    dofs = bc._cpp_object.dof_indices()[0]
                except:
                    printDeprecatedWarning("Mesh module from dolfinx version <0.7.0 is used.")
                    dofs = bc.dof_indices()[0]
                bc_vals = 1.0
                bcFunction.x.array[dofs] = bc_vals
            BC_Diriclet.setDiagonal(bcFunction.vector)
            BC_Diriclet.assemble()
        from petsc4py import PETSc
        self.__matrix_dict_petsc['B'] = B  
        self.__matrix_dict_petsc['A'] = A
        del A,B

        if AnalysisMode in ['Resolvent']:
            #forcing_norm,response_norm=getResolventNorms(param,MF,FEMSpaces)
            
            B_forcing = assemble_matrix(form(self.forcing_vf))
            B_forcing.assemble()
            self.__matrix_dict_petsc['B_forcing'] = B_forcing 
            del B_forcing

            B_response = assemble_matrix(form(self.response_vf))
            B_response.assemble()
            self.__matrix_dict_petsc['B_response'] = B_response 
            del B_response
            
            # FEM weighting matrix
            B_femWeight = assemble_matrix(form(self.fem_weighting))
            B_femWeight.assemble()
            self.__matrix_dict_petsc['B_femWeight'] = B_femWeight 
            del B_femWeight

        self.__matrix_dict_petsc['bcs'] = bcs

        return self.__buildSolutionObj()

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

    def __buildSolutionObj(self):
        return LinearSystem.linearSystem(
                                        self,
                                        self.__matrix_dict_petsc,
                                        self.__FEMSpaces,
                                        self.__param,
                                        self.__mean
                                        )        

    def getPMat(self):
        '''
        This function provides the P matrix, which restricts the forcing
        in terms of variables and spatial region
        
        TODO This does not work properly if a P1-fluctuations is part
        of the forcing/response norm! Indices of the DOFs will be wrong.
        '''
    
        from petsc4py import PETSc

        self.__forcing_coeff = self.__param.IOResolvent.ForcingCoeff
        self.__nVelocityComponents = self.__param.Case.getNVelocityComponents()
        printDebug(True, '-- Building Pu matrix...')

        # Get the matrix that restricts the forcing in space
        forcingDom = self.__mean.forcingDomain

        # By default, the forcing is applied everywhere, but the corresponding
        # matrix is only zeros, so we check and convert to ones in the default setting
        if max(forcingDom.x.array[:], key=abs) == 0:
            forcingDom.x.array[:] =  1
            flagdom = False
            printDebug(True, '-- No spatial restriction of forcing.')
        else:
            forcingDom.x.array[:] = np.rint(forcingDom.x.array[:])
            flagdom = True
            printDebug(True, '-- Applying spatial restriction of forcing from MeanFlow file.')

            nfluctvar = len(self.__bcDict)      # counting velocity components as 
            nDim = self.__param.Case.nDim
            forcingDomainVMixed = self.__FEMSpaces._projectField2allFEMSpaces(forcingDom, nfluctvar, nDim) # the last two inputs don't matter

        index = np.empty(shape=(0, 0))

        for i in self.__forcing_coeff:
            if i < self.__nVelocityComponents:
                dofsIterator = self.__FEMSpaces.VMixed.sub(0).sub(i).collapse()[1]
            else:
                dofsIterator = self.__FEMSpaces.VMixed.sub(i - self.__nVelocityComponents + 1).collapse()[1]

            # Only goes through the index loop if spatial limiter given
            if flagdom:
                for index_local in dofsIterator:
                    #if round(forcingDomainVMixed.vector()[index_local]) == 1:
                    if forcingDomainVMixed.x.array[index_local] == 1:
                        index = np.append(index, index_local)
                        # Rounding added there since mesh interpolation can result
                        # in non-integer values of the limiter domain flag
            else:
                index = np.append(index, dofsIterator)

        index.sort()
        # see: https://fenicsproject.discourse.group/t/fenicsx-method-which-is-equaivalent-to-dofmap-dofs-in-fenics/9152/3
        local_range = self.__FEMSpaces.VMixed.dofmap.index_map.local_range
        m = len(np.arange(*local_range))
        n = len(index)
        row_ind = index
        col_ind = np.arange(n)

        P_petsc = PETSc.Mat().createAIJ([m,n])
        P_petsc.setUp()
        for i in range(n):
            # P_petsc.setValue(row_ind[i],col_ind[i],1.,1)
            P_petsc.setValue(row_ind[i],col_ind[i],1.)
        P_petsc.assemble()

        printDebug(True, '-- Done.')

        return P_petsc

    def getCrMat(self):
        ''' This function provides the Cr matrix, which restricts the response '''
        from petsc4py import PETSc

        self.__response_coeff = self.__param.IOResolvent.ResponseCoeff
        self.__nVelocityComponents = self.__param.Case.getNVelocityComponents()
        printDebug(True, '-- Building Cr matrix...')

        # Get the matrix that restricts the forcing in space
        responseDom = self.__mean.responseDomain

        # By default, the forcing is applied everywhere, but the corresponding
        # matrix is only zeros, so we check and convert to ones in the default setting
        nfluctvar = len(self.__bcDict)
        nDim = self.__param.Case.nDim
        if max(responseDom.x.array[:], key=abs) == 0:
            responseDom.x.array[:] = 1
            flagdom = False
            printDebug(True, '-- No spatial restriction of response.')
        else:
            responseDom.x.array[:] = np.rint(responseDom.x.array[:])
            flagdom = True
            printDebug(True, '-- Applying spatial restriction of response from MeanFlow file.')

            responseDomainVMixed = self.__FEMSpaces._projectField2allFEMSpaces(responseDom, nfluctvar, nDim)

        index = np.empty(shape=(0,0))

        # Loop over fluctuations which are part of response coeff
        for i in self.__response_coeff:
            if i < self.__nVelocityComponents:
                dofsIterator = self.__FEMSpaces.VMixed.sub(0).sub(i).collapse()[1]
            else:
                dofsIterator = self.__FEMSpaces.VMixed.sub(i - self.__nVelocityComponents + 1).collapse()[1]

            # Only goes through the index loop if spatial limiter given
            if flagdom:
                for index_local in dofsIterator:
                    if round(responseDomainVMixed.x.array[index_local]) == 1:
                        index = np.append(index, index_local)
            else:
                index = np.append(index, dofsIterator)

        index.sort()
        local_range = self.__FEMSpaces.VMixed.dofmap.index_map.local_range
        m = len(np.arange(*local_range))
        n = len(index)
        row_ind = index
        col_ind = np.arange(n)

        Cr_petsc = PETSc.Mat().createAIJ([m,m])
        Cr_petsc.setUp()
        for row in row_ind:
            Cr_petsc.setValue(row,row,1.)
        Cr_petsc.assemble()

        ## Cr = csr_matrix((np.ones(n),(col_ind,row_ind)),(n,m))     # in theory this should be the size of Cr

        printDebug(True, '-- Done.')
        
        return Cr_petsc
