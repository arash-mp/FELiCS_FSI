from abc import ABC, abstractmethod

import dolfinx

class TimeStepperTemplate(ABC):

    def __init__(self, FEMSpace, mesh, meanFlow, equation):
        from ufl import TestFunctions, TrialFunctions, conj

        self.FEMSpace = FEMSpace

        self.uTest, self.pTest = TestFunctions(FEMSpace)
        self.uTest, self.pTest = conj(self.uTest), conj(self.pTest) 
        self.uFluc, self.pFluc = TrialFunctions(FEMSpace) 

        self.sponge   = meanFlow._fieldDict['spg']
        self.visc     = meanFlow._fieldDict['nulam']
        self.u_target = meanFlow._fieldDict['u_target']
        self.p_target = meanFlow._fieldDict['p_target']

        self.gamma = 2.

        self.BCs   = equation.BCs
        self.mesh  = mesh


    @abstractmethod
    def doTimeStepping(self, q_init, t_start, dt, t_end):
        pass

    @abstractmethod
    def _getOperator(self, q, q_new, dt):
        pass


    @abstractmethod
    def _getRHS(self, q, q_new, dt):
        pass


    def setMonitoringMethod(self, method):
        # the monitoring method has to have 3 arguments: 
        # iter: [integer]: iteration number (counted from 0 onwards)
        # q: [Field]: time dependent solution
        # t: [float]: current time
        self.monitor = method

