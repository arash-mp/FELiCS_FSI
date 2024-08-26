import numpy as np

import dolfinx

from ufl import (
                i,
                j,
                Dx,
                TestFunctions,
                TrialFunctions,
                dx,
                conj,
                inner,
                grad,
                div,
                )

from dolfinx.fem import (
                Function,
                FunctionSpace,
                dirichletbc,
                Constant,
                form,
)

from dolfinx.fem.petsc import (
                assemble_matrix,
                assemble_vector,
                set_bc,
)


from   FELiCS.Solvers.LinearSolver import LinearSolver 

from   FELiCS.Fields.Field import Field

from .TimeStepperTemplate import TimeStepperTemplate

class CrankNicolson_inc(TimeStepperTemplate):

    def __init__(self, FEMSpace, mesh, meanFlow, equation):
        super().__init__(FEMSpace, mesh, meanFlow, equation)

        
    def doTimeStepping(self, q_init, t_start, dt, t_end):
     
        q         = Field(self.FEMSpace, self.mesh)
        dq        = Field(self.FEMSpace, self.mesh)
        q_new     = Field(self.FEMSpace, self.mesh)
        q_old     = Field(self.FEMSpace, self.mesh)
        q_true    = Field(self.FEMSpace, self.mesh)
        #q_smooth = Field(FEMSpaces.VMixed, mesh)
        q.setCoefficientArray(    q_init.getCoefficientArray()) 
        q_new.setCoefficientArray(q_init.getCoefficientArray()) 
        q_old.setCoefficientArray(q_init.getCoefficientArray()) 
    
    
        t = t_start
        i = 0
        numberOfNewtonIterations = 1

        while t < t_end:

            dq.setConstantValue(0.)  
    
            for j in range(numberOfNewtonIterations):
                q_new.setCoefficientArray(q.getCoefficientArray() + dq.getCoefficientArray()) #q_new = q_new - dq_Newton
    
                A  = self._getOperator(q, q_new, dt)
                b  = self._getRHS(q, q_new, dt)  
               
                dq_Newton = LinearSolver.solveEquationSystem(A,b,destroy=True)
    
                dq.setCoefficientArray(dq.getCoefficientArray() - dq_Newton[:]) #q_new = q_new - dq_Newton
    
            q_old.setCoefficientArray(q.getCoefficientArray())
            q.setCoefficientArray(q_old.getCoefficientArray() + dq.getCoefficientArray())
    
            q_true.setCoefficientArray((q.getCoefficientArray() + q_old.getCoefficientArray())/2.)
    
            t = t + dt 
            i = i + 1
    
            ############################################
            ####### treat solution #####################
            ############################################
    
            # handle oscillating pressure solution (due to a non-divergence-free start solution) 
            if np.mod(i,10)==0:
                q.setCoefficientArray(q_true.getCoefficientArray())
                t = t - 0.5*dt
   
            if self.monitor != None:
                self.monitor(i, q, t)

    
        return q 
    
    
    
    def _getOperator(self, q, q_new,dt):
        u,p = q_new.function.split()
        
        gamma = self.gamma 
        
        op  = 1./dt * self.uTest[i]*self.uFluc[i]*dx
        
        op -= 1./2.*(3.-gamma)*u[0]*self.uFluc[0]*Dx(self.uTest[0],0)*dx
        op -= 1./2.*(1.-gamma)*u[1]*self.uFluc[1]*Dx(self.uTest[0],0)*dx
        op -= 1./2.*(u[0]*self.uFluc[1]+u[1]*self.uFluc[0])*Dx(self.uTest[0],1)*dx
        op -= 1./2.*(u[0]*self.uFluc[1]+u[1]*self.uFluc[0])*Dx(self.uTest[1],0)*dx
        op -= 1./2.*(3.-gamma)*u[1]*self.uFluc[1]*Dx(self.uTest[1],1)*dx
        op -= 1./2.*(1.-gamma)*u[0]*self.uFluc[0]*Dx(self.uTest[1],1)*dx
        
        op += 1./2. * self.visc * Dx(self.uFluc[i],j)*Dx(self.uTest[i],j)*dx
        op -= 1./2. * (gamma-1.) * inner(div(self.uTest),conj(self.pFluc))*dx  
        # no 1./2. here: implicit Euler for divergence-free equation  
        op -= inner(grad(self.pTest),conj(self.uFluc))*dx 
        op += 1./2. * self.sponge*self.uFluc[i]*self.uTest[i]*dx
        op += 1./2. * self.sponge*self.pFluc*self.pTest*dx
    
        # set correct mesh object
        sd = op.subdomain_data()
        domain, = list(sd.keys())  
        domain._ufl_cargo = self.mesh._cpp_object._cpp_object
    
        # assemble petsc matrix
        A = assemble_matrix(form(op), bcs=self.BCs)
        A.assemble()
    
        return A
    
    
    def _getRHS(self, q, q_new, dt):
        u,p         = q.function.split() 
        u_new,p_new = q_new.function.split()
        
        gamma = self.gamma
        
        rhs  = 1./dt * self.uTest[i] * (u_new[i] - u[i]) * dx
        
        rhs -= 1./2.* (3.-gamma)/2.*u[0]*u[0]*Dx(self.uTest[0],0)*dx
        rhs -= 1./2.* (1.-gamma)/2.*u[1]*u[1]*Dx(self.uTest[0],0)*dx
        rhs -= 1./2.* (u[0]*u[1])*Dx(self.uTest[0],1)*dx
        rhs -= 1./2.* (u[0]*u[1])*Dx(self.uTest[1],0)*dx
        rhs -= 1./2.* (3.-gamma)/2.*u[1]*u[1]*Dx(self.uTest[1],1)*dx
        rhs -= 1./2.* (1.-gamma)/2.*u[0]*u[0]*Dx(self.uTest[1],1)*dx
        
        rhs -= 1./2.* (3.-gamma)/2.*u_new[0]*u_new[0]*Dx(self.uTest[0],0)*dx
        rhs -= 1./2.* (1.-gamma)/2.*u_new[1]*u_new[1]*Dx(self.uTest[0],0)*dx
        rhs -= 1./2.* (u_new[0]*u_new[1])*Dx(self.uTest[0],1)*dx
        rhs -= 1./2.* (u_new[0]*u_new[1])*Dx(self.uTest[1],0)*dx
        rhs -= 1./2.* (3.-gamma)/2.*u_new[1]*u_new[1]*Dx(self.uTest[1],1)*dx
        rhs -= 1./2.* (1.-gamma)/2.*u_new[0]*u_new[0]*Dx(self.uTest[1],1)*dx
        
        rhs += 1./2.* self.visc * Dx (u[i],j)*Dx(self.uTest[i],j)*dx
        rhs += 1./2.* self.visc * Dx (u_new[i],j)*Dx(self.uTest[i],j)*dx
        
        rhs -= 1./2.* (gamma-1.)*inner(div(self.uTest),p)*dx  
        rhs -= 1./2.* (gamma-1.)*inner(div(self.uTest),p_new)*dx  
        
        rhs += 1./2.* self.sponge*((u[0]-self.u_target[0])*self.uTest[0]+(u[1]-self.u_target[1])*self.uTest[1])*dx
        rhs += 1./2.* self.sponge*((u_new[0]-self.u_target[0])*self.uTest[0]+(u_new[1]-self.u_target[1])*self.uTest[1])*dx
        
        rhs += 1./2.* self.sponge*(p-self.p_target)*self.pTest*dx
        rhs += 1./2.* self.sponge*(p_new-self.p_target)*self.pTest*dx
        
        # no 1./2. here: implicit Euler for divergence-free equation  
        rhs -= inner(grad(self.pTest),u_new)*dx 
       
        # set correct mesh object
        sd = rhs.subdomain_data()
        domain, = list(sd.keys())  
        domain._ufl_cargo = self.mesh._cpp_object._cpp_object
    
        # assemble petsc vector 
        b = assemble_vector(form(rhs))
        set_bc(b, self.BCs)
        b.assemble()
    
        return b
    
    
#def getOperatorCrankNicolson(q_new,dt,B,equation,meanFlow):
#    # Calculates the operator of the Crank Nicolson method, calculated with Newton-Steps:
#    # op = B + 0.5*dt*LinearOperator_new
#
#    # linear operator with q_new
#    [u_new,p_new] = q_new.getListOfSingleFields()
#    meanFlow._fieldDict['u'] = u_new.function
#    meanFlow._fieldDict['p'] = p_new.function
#    L = equation.getLinearOperator(meanFlow)
#    L.scale(-0.5)
#
#    # L = B + L*dt/2
#    #L.aypx(-dt/2., B) # L is negative, due to the weak form
#    #L.aypx(-dt, B) # L is negative, due to the weak form # implicit Euler
#    L.axpy(1./dt, B) # L is negative, due to the weak form 
#     
#    return L 
#                           
#
#
#def getRightHandSideCrankNicolson(q,q_new,dt,B,equation, meanFlow):
#    # Calculates the right hand side of the Crank Nicolson method, calculated with Newton-Steps:
#    # rhs = (q_new-q) + 0.5*dt*(NavierStokes+NavierStokes_new)
#
#    [u,p] = q.getListOfSingleFields()
#    meanFlow._fieldDict['u'] = u.function
#    meanFlow._fieldDict['p'] = p.function
#    N = equation.getNonlinearExpression(meanFlow,withoutP=True)
#
#    [u_new,p_new] = q_new.getListOfSingleFields()
#    meanFlow._fieldDict['u'] = u_new.function
#    meanFlow._fieldDict['p'] = p_new.function
#    N_new = equation.getNonlinearExpression(meanFlow)
#    #N = equation.getNonlinearExpression(meanFlow) # implicit Euler
#
#    # dq = u_new - u
#    q_new  = q_new.getPetscVector()
#    dq     = q.getPetscVector()
#    q_new.axpy(-1.,dq)
#
#    B.mult(q_new, dq)
#
#    dq.scale(1./dt)
#
#    # N = N + N_new
#    N.axpy(1.,N_new)
#    N.scale(0.5)
#
#    # N = dq + N*dt/2
#    #N.aypx(dt/2., dq) # N has the opposite sign to L 
#    #N.aypx(dt, dq) # N has the opposite sign to L # implicit Euler
#    N.axpy(1., dq) # N has the opposite sign to L # implicit Euler
#   
#    #N_new.destroy()
#    q_new.destroy()
#    dq.destroy()
#
#    return N
                           

