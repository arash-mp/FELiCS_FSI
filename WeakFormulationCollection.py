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
                form,
                Constant,
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

#from fenics import  PETScMatrix, PETScVector,DirichletBC, as_backend_type
from scipy.sparse import (
	csr_matrix,
	csc_matrix
	)
from functions import *
#from fenics import FunctionAssigner,Constant
import pdb
import LinearSystem

from tensorUtils import (
    Tensor,
    as_vector,
    iInner,
    iDot,
    iConj,
)

class WeakFormulationCollectionClass():
    '''This class build the variational formulations for all relevant matrices
    Currently these are:
    -A (imag and real)
    -B (imag and real)
    -B_forcing (Resolvent forcing norm)
    -B_response (Resolvent response norm)
    The convention is such that the B matrix (time derivative) is always positive and real
    '''
    def __init__(self,param,FEMSpaces,mean):
        #from fenics import Function
        from itertools import compress
        from fluctuationClass import fluctuationClass
        from WeakForm import WeakForm
        from tensorUtils import (
            CoordinateSystem,
            )

        # add the parameters of the constructor as attributs of the class to use them in DiscretizeFlow-method:
        self.__param = param
        self.__FEMSpaces = FEMSpaces
        self.__mean = mean
        mesh = self.__FEMSpaces.P2.mesh

        # Get class for integrating along boundaries
        self.boundaries = mesh.facet_tags
        self.ds = Measure("ds", subdomain_data=self.boundaries)

        # Get crossstreamwise wave number
        self.m = self.__param.Case.m

        # Get all boundaries (So far hard coded)
        first_BC_flag=True
        for Boundary in self.__param.BCs.getBCsDict()[list(self.__param.BCs.getBCsDict().keys())[0]]:
            if first_BC_flag:
                self.all_ds = self.ds(Boundary['ID'])
                first_BC_flag = False
            else:
                self.all_ds += self.ds(Boundary['ID'])
                
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
        
        # Define test and trial functions
        fluctuationC = fluctuationClass(
                                   param,
                                   mean,
                                   FEMSpaces,
                                   self._coordinateSystem,
                                   )
        #self.hat=TrialFunctions(FEMSpaces.VMixed)
        self.hat=fluctuationC.fluc

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
            

        # Get radial coordinate
        if self.__param.Case.CoordinateSystem in ['Cylindrical']:
            self.R = self.x[1]
            #self.ThirdVelCompIndex = self.__param.SolutionList.index('ut')
        else:
            from petsc4py import PETSc
            self.R=Constant(self.__FEMSpaces.P2.mesh, PETSc.ScalarType(1.0))

        # Get boundary normals
        self.n_BC=FacetNormal(self.__FEMSpaces.P2.mesh)
        self.n = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), self._coordinateSystem)
        ## Initialize variatial formulations
        self.A_vf = WeakForm()
        self.B_vf = WeakForm()
        printDebug(True, '-- Primary fluctuations: %s.' % param.SolutionList)
        if self.__param.Case.SetOfEquations['Momentum']['Equation'] == 'NSPrimitive':
            from Equations.Momentum.addMomentumEq_tensorial import addMomentumEq
            printDebug(True, '-- Adding momentum equation for u-fluc -> X[0].')     # Hardcoded u' for mom eq.
            addMomentumEq(self,fluctuationC,X[0],mean,param)
            print('-- Adding momentum equation.')
            
        if self.__param.Case.SetOfEquations['Mass']['Equation'] == 'Continuity':
            from Equations.Mass.addMassEq_tensorial import addMassEq
            varEq = self.__param.Case.SetOfEquations['Mass']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding mass-balance equation for %s-fluc -> X[%d].' % (varEq,idVar))
            addMassEq(self,fluctuationC,X[idVar],mean,self.__param)

        if self.__param.Case.SetOfEquations['Energy']['Equation'] == 'Enthalpy':
            from Equations.Enthalpy.addEnthalpyEq_tensorial import addEnthalpyEq
            varEq = self.__param.Case.SetOfEquations['Energy']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding enthalpy-energy equation for %s-fluc -> X[%d].' % (varEq,idVar))
            addEnthalpyEq(self,fluctuationC,X[idVar],mean,self.__param)
        
        if self.__param.Case.SetOfEquations['Energy']['Equation'] == 'primitive-p':
            printError('Energy equation in primitive form is not ready to use!!! Ask Simon Demange for updates.')
            from Equations.Energy_Pressure.addEnergyPEq_tensorial import addEnergyPEq
            varEq = self.__param.Case.SetOfEquations['Energy']['Variable']
            idVar = param.SolutionList.index(varEq)
            printDebug(True, '-- Adding pressure-energy equation for %s-fluc -> X[%d].' % (varEq,idVar))
            addEnergyPEq(self,fluctuationC,X[idVar],mean,self.__param)
            
        # Add sponge region only if the field was gieven in the mean flow file
        if not('spg' in mean._meanFlowClass__notInFileList):
            from Equations.addSpongeEq_tensorial import addSpongeEq
            printDebug(True, '-- Adding sponge damping.')
            addSpongeEq(self,fluctuationC,X,mean,self.__param)

        if self.__param.Case.AnalysisMode in ['Resolvent']:
            self.getResolventNorms(X,self.__param,mean,fluctuationC)

        # Add species transport equation for all transported species
        transportedSpecies=self.__param.Case.Mixture.getSpeciesList('transported')
        for specie in transportedSpecies:
            i_eqn=self.__param.SolutionList.index(specie)
            if self.__param.Case.SetOfEquations['Species']['Equation'] == 'Non-conservative':
                from Equations.Species.addSpeciesEq_tensorial import addSpeciesEq
                print('-- Adding equation for species '+specie + ' in non-conservative form')
                addSpeciesEq(self,fluctuationC,X[i_eqn],mean,specie,self.__param)
                
            elif self.__param.Case.SetOfEquations['Species']['Equation'] == 'Conservative':
                # This eq has not been derived in tensor framework yet.
                from Equations.speciesConservative.addSpeciesConservativeEq import addSpeciesConservativeEq
                print('-- Adding equation for species '+specie +' in conservative form')
                addSpeciesConservativeEq(self,fluctuationC,X[i_eqn],self.mean,specie,self.__param)
            else:
                raise Exception('Species transport equation type ' + self.__param.Case.SetOfEquations['Species']['Equation'] + ' unknown' )

        # Add reactions
        # Reaction eqs not derived in tensor framework yet
        if self.__param.Case.Reaction:
            if self.__param.Case.Mixture.ReactionMechanism['type']=='WestbrookDryer_Max':
                from Reactions.GlobalReaction import GlobalReaction
                ReactionModelName="WestbrookDryer_Max" #to be put in param
                Reaction=GlobalReaction(self.__param.Case.Mixture.ReactionMechanism)
                reactionRateMean=Reaction.computeMeanField(self.mean,self.__FEMSpaces.P2)
                reactionForm=Reaction.addReaction(self.mean, X, fluctuationC, self.__param.SolutionList)
                self.A_vf.add(1j * reactionForm)
            elif self.__param.Case.Mixture.ReactionMechanism['type']=='TwoStep':
                from Reactions.TwoStepReaction import TwoStepReaction
                ReactionModelName="BFER" #to be put in param
                Reaction=TwoStepReaction(ReactionModelName)
                Reaction.computeMeanField(MF,self.__FEMSpaces.P2)
                Reaction.testM()
                reactionForm=Reaction.addReaction(MF, X, fluc, self.__param.SolutionList,self.__FEMSpaces.P2)
                self.A_vf.add(1j * reactionForm)
            elif self.__param.Case.Mixture.ReactionMechanism['type']=='2S-SM2':
                from Reactions.c2sm2 import C2SM2
                ReactionModelName="2S-SM2" #to be put in param
                YCH4_lim=0.043*1e-4
                #c2=C2SM2(YCH4_lim,2)
                c2=self.mean.reaction
                self.TR=fluctuationC.T
                self.rhoR=fluctuationC.rho
                self.YCH4R=fluctuationC.Y('CH4')
                self.YO2R=fluctuationC.Y('O2')
                self.YCOR=fluctuationC.Y('CO')
                self.YCO2R=fluctuationC.Y('CO2')
                self.v_eneR=X[self.__param.Case.getTransportedQuantityList().index('rho')]
                self.v_YCH4R=X[self.__param.Case.getTransportedQuantityList().index('CH4')]
                self.v_YO2R=X[self.__param.Case.getTransportedQuantityList().index('O2')]
                self.v_YH2OR=X[self.__param.Case.getTransportedQuantityList().index('H2O')]
                self.v_YCOR=X[self.__param.Case.getTransportedQuantityList().index('CO')]
                self.v_YCO2R=X[self.__param.Case.getTransportedQuantityList().index('CO2')]
                self.dQMean=self.mean.dQ
                self.order=2
                self.dx=dx

                #c2.computeSensitivities(MeanFlow.T,
                #                        MeanFlow.rho,
                #                        MeanFlow.Y('CH4'),
                #                        MeanFlow.Y('CO'),
                #                        MeanFlow.Y('O2'),
                #                        MeanFlow.Y('CO2'))

                self.A_vf.add(1j * -c2.add_source_to_weak_form(self))
            elif self.__param.Case.Mixture.ReactionMechanism['type']=='NOx':
                reaction=self.mean.reaction
                self.v_NO  = X[self.__param.Case.getTransportedQuantityList().index('NO')]
                self.v_NO2 = X[self.__param.Case.getTransportedQuantityList().index('NO2')]
                self.T = self.mean.T
                self.phi = self.mean.phi
                self.A_vf.add(1j * -reaction.add_source_to_weak_form(self))

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
            printWarning("Currently only the L2 norm is implemented for both forcing and response in a resolvent analysis. Here, ALL velocity components are taken into account, no matter the choices in the settings file.")
            u_f = fluc.u
            self.forcing_vf += (barrho*iDot(u_f,iConj(X[0]))).ufl_tens*self._coordinateSystem.J_hat*dx
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

        # Prompt variational formulations in debug mode
        printDebug(param.debug,'Resolvent forcing norm is '+ str(self.forcing_vf))
        printDebug(param.debug,'Resolvent response norm is '+ str(self.response_vf))

    def DiscretizeFlow(self):
        #self.__DiscretizeAndSolve = DiscretizeAndSolve
        #self.__WeakForm = WeakForm
        from copy import deepcopy
        # Check for type of case

        AnalysisMode = self.__param.Case.AnalysisMode
        FEMSpaces = self.__FEMSpaces
        #Get the Dirichlet BCs set by the user in a list...
        # BClist= self.__getListOfDirichletBCs()
        bcFunction = Function(self.__FEMSpaces.VMixed)
        bcFunction.x.array[:] = 0.0
        # for bc in BClist:
        # 	dofs = bc.dof_indices()[0]
        # 	bc_vals = 1.0
        # 	bcFunction.x.array[dofs] = bc_vals
        mesh = self.__FEMSpaces.P2.mesh
        pattern = SparsityPattern(mesh.comm, [self.__FEMSpaces.VMixed.dofmap.index_map, self.__FEMSpaces.VMixed.dofmap.index_map],
                                                            [self.__FEMSpaces.VMixed.dofmap.index_map_bs, self.__FEMSpaces.VMixed.dofmap.index_map_bs])
        pattern.insert_diagonal(np.arange(len(bcFunction.vector.array), dtype=np.int32))
        pattern.assemble()
        BC_Diriclet = create_matrix(mesh.comm, pattern)
        BC_Diriclet.setDiagonal(bcFunction.vector)
        BC_Diriclet.assemble()

        bcs= self.__getListOfDirichletBCs()
        n_dof=BC_Diriclet.size[0]
        if not self.A_vf.lhsIsZero():
            if AnalysisMode in ['Input-Output']:
                A = assemble_matrix(form(self.A_vf.lhs), bcs=bcs)
                A.assemble()
            else:
                A = assemble_matrix(form(self.A_vf.lhs))
                A.assemble()
        else:
            A = 0*BC_Diriclet

        if not self.B_vf.lhsIsZero():
            B = assemble_matrix(form(self.B_vf.lhs))
            B.assemble()
        else:
            B = 0*BC_Diriclet

        self.__matrix_dict={}


        if AnalysisMode in ['Input-Output']:
            if not self.A_vf.rhsIsZero():
                forcing_vec_petsc = assemble_vector(form(self.A_vf.rhs))
                forcing_vec_petsc.assemble()
            else:
                forcing_vec_petsc = PETScVector()
                forcing_vec_petsc.init(n_dof)

            from dolfinx.fem.petsc import set_bc
            set_bc(forcing_vec_petsc, bcs)

            b_forcing = 1j * forcing_vec_petsc.array
            self.__matrix_dict['b_forcing'] = b_forcing
            del b_forcing, forcing_vec_petsc

        # Get the BCs provided by the user
        # Aplly the BCs
        ## Boundary conditions are applied via penalisation method manually. This is in order to keep the rhs matrix invertible.
        else:

            bcFunction.x.array[:] = 0.0
            for bc in bcs:
                dofs = bc.dof_indices()[0]
                bc_vals = 1.0
                bcFunction.x.array[dofs] = bc_vals
            BC_Diriclet.setDiagonal(bcFunction.vector)
            BC_Diriclet.assemble()
        # In the next 12 lines the Imaginary and real parts of both the lhs and rhs matrix are combined to the
        # sparse matrix A and B, respectively. Not needed matrices are deleted
        self.__matrix_dict['A'] = csr_matrix(A.getValuesCSR()[::-1], shape = A.size,dtype=complex)
        del A
        if not AnalysisMode in ['Input-Output']:
            #BC_Diriclet_mat = as_backend_type(BC_Diriclet).mat()
            self.__matrix_dict['A'] = self.__matrix_dict['A'] + 10**30*(1+1j) * csr_matrix(BC_Diriclet.getValuesCSR()[::-1], shape = BC_Diriclet.size,dtype=complex)
            del BC_Diriclet
        self.__matrix_dict['B'] = csr_matrix(B.getValuesCSR()[::-1], shape = B.size,dtype=complex)
        del B
        #tempMat=1*matrix_dict['B'].transpose()
        #tempMat[1,1]=100
        #input((matrix_dict['B'] != tempMat).nnz==0)
        #input(matrix_dict['B'])
        # must be constructed. So far only the L2 norm is implemented (To be extended!)
        if AnalysisMode in ['Resolvent']:
            #forcing_norm,response_norm=getResolventNorms(param,MF,FEMSpaces)
            B_forcing = assemble_matrix(form(self.forcing_vf))
            B_forcing.assemble()
            self.__matrix_dict['B_forcing'] = csr_matrix(B_forcing.getValuesCSR()[::-1], shape = B_forcing.size,dtype=complex)
            del B_forcing
            B_response = assemble_matrix(form(self.response_vf))
            B_response.assemble()
            self.__matrix_dict['B_response'] = csr_matrix(B_response.getValuesCSR()[::-1], shape = B_response.size,dtype=complex)
            del B_response

        return self.__buildSolutionObj()

    def __getListOfDirichletBCs(self):
        BClist = []
        mesh=self.__FEMSpaces.P2.mesh

        self.__boundaries = self.__param.BCs.getBoundaries()
        self.__bcDict = self.__param.BCs.getBCsDict()
        VelocityComponents=self.__param.Case.getVelocityComponents()
        SolutionList=self.__param.Case.getTransportedQuantityList()
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
                    if k in ['u'+ component for component in VelocityComponents]:

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
                                        self.__matrix_dict,
                                        self.__FEMSpaces,
                                        self.__param,
                                        self.__mean,
                                        )

    def getPMat(self):
        '''
        This function provides the P matrix, which restricts the forcing
        '''
    
        self.__forcing_coeff = self.__param.IOResolvent.ForcingCoeff
        #self.__nVelocityComponents = self.__param.nVelocityComponents
        self.__nVelocityComponents = self.__param.Case.getNVelocityComponents()
        print('-- Building Pu matrix...')

        # Get the matrix that restricts the forcing in space
        forcingDom = self.__mean.forcingDomain

        # By default, the forcing is applied everywhere, but the corresponding
        # matrix is only zeros, so we check and convert to ones in the default setting
        if max(forcingDom.x.array[:], key=abs) == 0:
            forcingDom.x.array[:] =  1
            flagdom = False
            print('-- No spatial restriction of forcing.')
        else:
            forcingDom.x.array[:] = np.rint(forcingDom.x.array[:])
            flagdom = True
            print('-- Applying spatial restriction of forcing from MeanFlow file.')

            nfluctvar = len(self.__bcDict)
            nDim = self.__param.Case.nDim
            forcingDomainVMixed = self.__FEMSpaces._projectField2allFEMSpaces(forcingDom, nfluctvar, nDim)

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
        P = csr_matrix((np.ones(n), (row_ind,col_ind)), (m, n))
        print('-- Done.')
        return P

    def getCrMat(self):
        ''' This function provides the Cr matrix, which restricts the response '''
        self.__forcing_coeff = self.__param.IOResolvent.ForcingCoeff
        #self.__nVelocityComponents = self.__param.nVelocityComponents
        self.__nVelocityComponents = self.__param.Case.getNVelocityComponents()
        print('-- Building Cr matrix...')

        # Get the matrix that restricts the forcing in space
        responseDom = self.__mean.responseDomain

        # By default, the forcing is applied everywhere, but the corresponding
        # matrix is only zeros, so we check and convert to ones in the default setting
        nfluctvar = len(self.__bcDict)
        nDim = self.__param.Case.nDim
        if max(responseDom.x.array[:], key=abs) == 0:
            responseDom.x.array[:] = 1
            flagdom = False
            print('-- No spatial restriction of response.')
        else:
            responseDom.x.array[:] = np.rint(responseDom.x.array[:])
            flagdom = True
            print('-- Applying spatial restriction of response from MeanFlow file.')

            responseDomainVMixed = self.__FEMSpaces._projectField2allFEMSpaces(responseDom, nfluctvar, nDim)

        index = np.empty(shape=(0,0))

        for i in range(nfluctvar): # HARDCODED FOR U, V, P: incompressible 2D
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
        Cr = csr_matrix((np.ones(n),(row_ind,row_ind)),(m,m))
        print('-- Done.')
        return Cr
