#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 17:59:50 2022

@author: kai hildebrandt
Modifications: 
- Thomas L. Kaiser
- Simon Demange
"""

from typing import Union
from ufl import indices
from ufl import Identity, as_vector, as_matrix, as_tensor
from ufl import sin, cos, tan, sqrt
from ufl import det, tr
from ufl import conj


### CLASSES
class CoordinateSystem():
    def __init__(self, SpatialCoordinateObj, name: str, **kwargs):
        x = SpatialCoordinateObj
        
        # mesh_dims is a variable to handle that the ufl vector x[i] might re-
        # present for example (r,z) in cylindrical coordinates but
        # tensor_utils can only operate with (r, theta, z).
        if 'mesh_dims' in kwargs:
            self.mesh_dims = kwargs['mesh_dims']
            self.dim = len(self.mesh_dims)
            x_vec = []
            iter = 0

            # x vector is filled with 1 if i-th direction is not a dimension of
            # the mesh. Its filled with a one to avoid 1/0 in the christoffel
            # symbols etc.
            for i in self.mesh_dims:
                if i == 1:
                    x_vec.append(x[iter])
                    iter += 1
                else:
                    x_vec.append(1.0)
            self.x = as_vector(x_vec)

        # if there is no dimensionality reduction standard x can be used
        else:
            self.dim = SpatialCoordinateObj.geometric_dimension()
            self.mesh_dims = self.dim * [1]
            self.x = x

        # this is to avoid writing self.x everywhere
        x = self.x

        # Hard coded data about the coordinate system. Information for example
        # in Tensoranalysis by Schade and Neemann (2018, ISBN: 9783110404265)
        if name == "cartesian":
            self.ch = as_tensor([[[0 for col in range(self.dim)] for row in \
                                  range(self.dim)] for christoffel in \
                                 range(self.dim)])
            self.cov_metric = self.con_metric = Identity(self.dim)
            """ 
            In standard case (cartesian) the christoffel symbols are zero.
            """
        
        elif name == "polar":
            """
            See cylindrical for conventions.
            """
            self.dim = 2
            ch1 = [[0, 0],[0, -x[0]]]
            ch2 = [[0, 1/x[0]],[1/x[0], 0]]           
            self.ch = as_tensor([ch1, ch2])
            self.cov_metric = as_matrix([[1, 0], [0, x[0]*x[0]]])
            self.con_metric = as_matrix([[1, 0], [0, 1/(x[0]*x[0])]])
        
        elif name == "cylindrical":
            self.dim = 3
            """
            Given for {x,y,z} = {r\cos\phi, r\sin\phi, z}, where ordering is
            {r, \phi, z}. This results in tangent basis:
            {g_i} = {\cos\phi \vec{e}_x + \sin\phi \vec{e}_y,
                   r(-\sin\phi \vec{e}_x + \cos\phi \vec{e}_y),
                   \vec{e}_z}
            where {\vec{e}_i} denote cartesian base vectors.
            """
            ch1 = [[0, 0, 0],\
                   [0, -x[0], 0],\
                   [0, 0, 0]]
            ch2 = [[0, 1/x[0], 0],\
                   [1/x[0], 0, 0],\
                   [0, 0, 0]]
            ch3 = [[0, 0, 0],\
                   [0, 0, 0],\
                   [0, 0, 0]]
            self.ch = as_tensor([ch1, ch2, ch3])
            self.cov_metric = as_matrix([[1, 0, 0], \
                                         [0, x[0]*x[0], 0], \
                                         [0, 0, 1]]
                                        )
            self.con_metric = as_matrix([[1, 0, 0], \
                                         [0, 1/(x[0]*x[0]), 0], \
                                         [0, 0, 1]]
                                        )
            
        elif name == "cylindricalfelics":
            self.dim = 3
            """
            Given for {x,y,z} = {r\cos\phi, r\sin\phi, z}, where ordering is
            {z, r, \phi}.
            This results in tangent basis:
            {g_i} = {\vec{e}_z,
                   \cos\phi \vec{e}_x + \sin\phi \vec{e}_y,
                   r(-\sin\phi \vec{e}_x + \cos\phi \vec{e}_y)}
            where {\vec{e}_i} denote cartesian base vectors.
            """
            ch1 = [[0, 0, 0], \
                   [0, 0, 0], \
                   [0, 0, 0]]
            ch2 = [[0, 0, 0], \
                   [0, 0, 0], \
                   [0, 0, -x[1]]]
            ch3 = [[0, 0, 0], \
                   [0, 0, 1/x[1]], \
                   [0, 1/x[1], 0]]
            self.ch = as_tensor([ch1, ch2, ch3])
            self.cov_metric = as_matrix([[1, 0, 0], \
                                         [0, 1, 0], \
                                         [0, 0, x[1]*x[1]]]
                                        )
            self.con_metric = as_matrix([[1, 0, 0], \
                                         [0, 1, 0], \
                                         [0, 0, 1/(x[1]*x[1])]]
                                        )

        elif name == "spherical":
            """
            Given for {x, y, z} = {r\sin\theta\cos\phi, r\sin\theta\sin\phi,
            r\cos\theta} with order {r, \theta, \phi}
            """
            self.dim = 3
            
            ch1 = [[0,0,0], [0,-x[0],0], [0,0,-x[0]*sin(x[1])**2]]
            ch2 = [[0,1/x[0],0], [1/x[0],0,0], [0,0,-sin(x[1])*cos(x[1])]]
            ch3 = [[0,0,1/x[0]], [0,0,1/tan(x[1])], [1/x[0],1/tan(x[1]),0]]
            self.ch = as_tensor([ch1, ch2, ch3])
            self.cov_metric = as_matrix([[1, 0, 0], [0, x[0]*x[0], 0], \
                                         [0, 0, x[0]*x[0]*sin(x[1])*sin(x[1])]])
            self.con_metric = as_matrix([[1, 0, 0], [0, 1/(x[0]*x[0]), 0], \
                                         [0, 0, 1/(x[0]*x[0]*sin(x[1])**2)]])
            
        else:
            raise ValueError("The specified coordinate system isn't " \
                             "implemented.")
        
        self.g = det(self.cov_metric)
        self.J_hat = sqrt(self.g)
            
        """
        # Example for index operations on Matrices in UFL.
        ch1 = [[1, 2, 3],[4,5,6], [7,8,9]]
        ch2 = [[10, 11, 12],[13,14,15], [16,17,18]]
        ch3 = [[19, 20, 21],[22,23,24], [25,26,27]]
        ch = as_matrix([ch1, ch2, ch3])
        print(ch3)
        print(ch[2,:,0]) # first column of the last christoffel symbol
        """


class Tensor():#TestFunction):    
    def __init__(
        self, 
        ufl_tens, 
        CoordSys: CoordinateSystem, 
        containsTestFunction = False, 
        containsFluctuation = False,
        m = 0,
        **kwargs,
        ):
        self.CoordSys = CoordSys
        self.x = CoordSys.x
        self.dim = CoordSys.dim
        self.order = len(ufl_tens.ufl_shape) # scalar --> order = 0
        self.containsTestFunction = containsTestFunction
        self.containsFluctuation = containsFluctuation
        
        # if no basis is specified, ufl_tens is assumed to be a vector given
        # in  physical basis. The vector is then transformed into a natural
        # basis, namely the tangent basis. This is done because tensor-
        # analytical operators, as they are defined in tensor_utils, require
        # it.
        if 'basis' in kwargs:
            self.ufl_tens = ufl_tens
            self.basis = kwargs['basis']
        else:
            if self.order == 0:
                self.ufl_tens = ufl_tens
                self.basis = []
            elif self.order == 1:
                length = ufl_tens.ufl_shape[0]
                ##if self.dim == 2:
                if length == 2:
                    self.ufl_tens = \
                        as_vector((ufl_tens[0]/sqrt(self.CoordSys.cov_metric[0,0]), \
                                   ufl_tens[1]/sqrt(self.CoordSys.cov_metric[1,1]),
                                    0.0))
                elif length == 3:
                    self.ufl_tens = \
                        as_vector((ufl_tens[0]/sqrt(self.CoordSys.cov_metric[0,0]), \
                                   ufl_tens[1]/sqrt(self.CoordSys.cov_metric[1,1]), \
                                   ufl_tens[2]/sqrt(self.CoordSys.cov_metric[2,2])))
                #entryTuple = []
                #print (length)
                #input (self.dim)
                #if length > self.dim:
                #    for i in range(self.dim):
                #        if CoordSys.mesh_dims[i] == 1:
                #            entryTuple.append(ufl_tens[i] / sqrt(self.CoordSys.cov_metric[i,i]))
                #        else:
                #            entryTuple.append( 0.0 )
                #else: 
                #    for i in range(self.dim):
                #        entryTuple.append(ufl_tens[i] / sqrt(self.CoordSys.cov_metric[i,i]))
                #self.ufl_tens = as_vector(tuple(entryTuple))
                self.basis = [False]
            else:
                raise ValueError("Basis transformation from physical to " \
                                 "tangent of tensors with order > 1 isn't " \
                                 "implemented yet.")
            
        if 'sym' in kwargs:
            self.sym = kwargs['sym']
        else:
            # if there is a zero in ufl_tens, sym entry must be zero to avoid
            # try to differentiate the zero in iGrad. If k-th mesh_dims entry
            # is zero, no tensor element can be a function of the k-th
            # coordinate
            sym = []
            if self.order == 0:
                for i in range(self.dim):
                    if self.ufl_tens == 0 or self.CoordSys.mesh_dims[i] == 0:
                        sym.append(0)
                    else:
                        sym.append(1)
                        
            elif self.order == 1:
                for i in range(self.dim):
                    row = []
                    for j in range(self.dim):
                        if self.ufl_tens[i] == 0 or self.CoordSys.mesh_dims[j] == 0:
                            row.append(0)
                        else:
                            row.append(1)
                    sym.append(tuple(row))
            elif self.order == 2:
                for i in range(self.dim):
                    row = []
                    for j in range(self.dim):
                        column = []
                        for k in range(self.dim):
                            if self.ufl_tens[i,j] == 0 or self.CoordSys.mesh_dims[k] == 0:
                                column.append(0)
                            else:
                                column.append(1)
                        row.append(tuple(column))
                    sym.append(tuple(row))
            self.sym = tuple(sym)
            
    # addition
    def __add__(self, other):
            
        if type(other) == Tensor:
            if not (self.containsTestFunction == other.containsTestFunction):
                ValueError("In a tensor sum, both tensors must be of same order in test functions.")
            if self.basis == other.basis:
                added = self.ufl_tens + other.ufl_tens
            else:
                new_self = convertBasis(self, other.basis)
                added = new_self.ufl_tens + other.ufl_tens
        else:
            added = self.ufl_tens + other
        return Tensor(added, self.CoordSys, basis = self.basis, containsTestFunction = self.containsTestFunction)

    # division, 
    def __truediv__(self, other): # Tensor object to the left
        if type(other) == Tensor:
            if other.containsTestFunction:
                ValueError("Division by test function not possible.")
            if other.order == 0:
                return Tensor(
                            self.ufl_tens / other.ufl_tens,
                            self.CoordSys, 
                            basis = self.basis, 
                            containsTestFunction = self.containsTestFunction,
                            )
            elif self.order == 0:
                return Tensor(
                            self.ufl_tens / other.ufl_tens,
                            self.CoordSys,
                            basis = other.basis,
                            containsTestFunction = self.containsTestFunction,
                            )
            else:
                raise ValueError("Division operation between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) != Tensor:
            return Tensor(
                        self.ufl_tens * other,
                        self.CoordSys,
                        basis = self.basis,
                        containsTestFunction = self.containsTestFunction,
                        )
    
    # subtraction: A - B is the same as A.__sub__(B)
    def __sub__(self, other):
        if not (self.containsTestFunction == other.containsTestFunction):
            ValueError("In a tensor difference, both tensors must be of same order in test functions.")
        if type(other) == Tensor:
            if self.basis == other.basis:
                subtracted = self.ufl_tens - other.ufl_tens
            else:
                new_self = convertBasis(self, other.basis)
                subtracted = new_self.ufl_tens - other.ufl_tens
        else:
            subtracted = self.ufl_tens - other
        
        return Tensor(subtracted, self.CoordSys, basis = self.basis, containsTestFunction = self.containsTestFunction)

    # muliplication, both ways, because matrix mul not commutative
    def __mul__(self, other): # Tensor object to the left
        if type(other) == Tensor:
            if self.containsTestFunction and other.containsTestFunction:
                ValueError("Tensor product at least second order in test functions.")
            if other.order == 0:
                return Tensor(
                            self.ufl_tens * other.ufl_tens,
                            self.CoordSys, 
                            basis = self.basis, 
                            containsTestFunction = self.containsTestFunction or self.containsTestFunction,
                            )
            elif self.order == 0:
                return Tensor(
                            self.ufl_tens * other.ufl_tens,
                            self.CoordSys,
                            basis = other.basis,
                            containsTestFunction = self.containsTestFunction or self.containsTestFunction,
                            )
            else:
                raise ValueError("Multiplication between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) != Tensor:
            return Tensor(
                        self.ufl_tens * other,
                        self.CoordSys,
                        basis = self.basis,
                        containsTestFunction = self.containsTestFunction,
                        )
    
    def __rmul__(self, other): # Tensor object to the right
        if type(other) == Tensor:
            if self.containsTestFunction and other.containsTestFunction:
                ValueError("Tensor product at least second order in test functions.")
            if other.order == 0:
                return Tensor(
                            other.ufl_tens * self.ufl_tens, 
                            self.CoordSys,
                            basis = self.basis,
                            containsTestFunction = self.containsTestFunction or self.containsTestFunction,
                            )
            elif self.order == 0:
                return Tensor(
                            other.ufl_tens * self.ufl_tens, 
                            self.CoordSys,
                            basis = other.basis,
                            containsTestFunction = self.containsTestFunction or self.containsTestFunction,
                            )
            else:
                raise ValueError("Multiplication between Tensors only " \
                                 "defined, if one is a scalar.")
        else:
            return Tensor(
                        other * self.ufl_tens,
                        self.CoordSys,
                        basis = self.basis,
                        containsTestFunction = self.containsTestFunction,
                        )
    
    """
    TODO: dont know if and how division is implemented in ufl
    def __truediv__(self, other):
        return Tensor(self.ufl_tens / other, self.CoordSys, self.basis, self.x)
    """


### TENSOR OBJECT FUNCTIONS
def iDot(
         tensorA: Tensor, 
         tensorB: Tensor,
         ):
    if tensorA.containsTestFunction and tensorB.containsTestFunction:
        ValueError("iDot product at least second order in test functions.")
    # Metric depends on base vectors which are being contracted
    metric = getMetric(tensorA.basis[-1], tensorB.basis[0], tensorA.CoordSys)
    A = tensorA.ufl_tens
    B = tensorB.ufl_tens
    
    # Perform dot Product
    if tensorA.order == 1 and tensorB.order == 1:
        i,j = indices(2)
        dotted = A[i]*metric[i,j]*B[j]
    if tensorA.order == 1 and tensorB.order == 2:
        i,j,k = indices(3)
        dotted = as_tensor(A[i]*metric[i,j]*B[j,k], (k))
    if tensorA.order == 2 and tensorB.order == 1:
        i,j,k = indices(3)
        dotted = as_tensor(A[i,j]*metric[j,k]*B[k], (i))
    if tensorA.order == 2 and tensorB.order == 2:
        i,j,k,l = indices(4)
        dotted = as_tensor(A[i,j]*metric[j,k]*B[k,l], (i,l))
    
    return Tensor(
            dotted, 
            tensorA.CoordSys, 
            basis = tensorA.basis[:-1] + tensorB.basis[1:],
            containsTestFunction = tensorA.containsTestFunction or tensorB.containsTestFunction, 
            )

def iDotT(tensorA: Tensor, tensorB: Tensor):
    """
    
    Args:
        - A: not transposed tensor of order 2
        - B_ tensor of order 1
    
    Return:
        - Dot product of A^T * b
    """
    if tensorA.containsTestFunction and tensorB.containsTestFunction:
        ValueError("iDotT product at least second order in test functions.")
    # Metric depends on base vectors which are being contracted
    metric = getMetric(tensorA.basis[-2], tensorB.basis[0], tensorA.CoordSys)
    A = tensorA.ufl_tens
    B = tensorB.ufl_tens

    if tensorA.order == 2 and tensorB.order == 1:
        i,j,k = indices(3)
        dotted = as_tensor(A[i,j]*metric[i,k]*B[k], (j))

    else:
        raise Exception('not implemented')
    # # Perform dot Product
    # if tensorA.order == 2 and tensorB.order == 1:
    #     i,j,k = indices(3)
    #     # dotted = as_tensor(A[j,i]*metric[i,k]*B[k], (j))
    #     dotted = as_tensor(metric[i,j]*A[j,k]*B[k], (i))
    #     print("I am here")
    
    return Tensor(
            dotted, 
            tensorA.CoordSys,
            basis = [tensorA.basis[-1]],
            containsTestFunction = tensorA.containsTestFunction or tensorB.containsTestFunction,
            )


def iInner(tensorA: Tensor, tensorB: Tensor):
    if tensorA.containsTestFunction and tensorB.containsTestFunction:
        raise ValueError("iInner product at least second order in test functions.")
    if tensorA.order != tensorB.order or tensorA.order != 2:
        raise ValueError("The order of both tensors must be two.")
    
    # Metric depends on base vectors which are being contracted
    metric1 = getMetric(tensorA.basis[-2], tensorB.basis[0], tensorA.CoordSys)
    metric2 = getMetric(tensorA.basis[-1], tensorB.basis[1], tensorA.CoordSys)
    
    i,j,k,l = indices(4)
    innered = tensorA.ufl_tens[i,j] * tensorB.ufl_tens[k,l] * metric1[i,k] * \
        metric2[j,l]
    
    return Tensor(
            innered,
            tensorA.CoordSys, 
            basis = tensorA.basis[:-2] + tensorB.basis[2:],
            containsTestFunction = tensorA.containsTestFunction or tensorB.containsTestFunction,
            )


def getMetric(basisA: Union[list, bool, int], basisB: Union[list, bool, int], \
              CoordSys: CoordinateSystem):
    if type(basisA) == list:
        basisA = basisA[0]
    if type(basisB) == list:
        basisB = basisB[0]
    
    if basisA and basisB:
        metric = CoordSys.con_metric
    elif not(basisA) and not(basisB):
        metric = CoordSys.cov_metric
    else:
        metric = Identity(CoordSys.dim)
    return metric


def iGrad(
    T: Tensor, 
    m=0,
    ):
    # implementing all kind of derivatives is tideous. --> convert every tensor
    # to tangent basis and only implement grads in tangent basis.
    if sum(T.basis) != 0:
        T = convertBasis(T, T.order*[False])

    # If the expression contains a test function, the wave numbers in mean flow homogeneous directions must be multiplied by minus one for the following reason: if a gradient operator contains a test function, integration by parts has been applied. With the wave numbers, however an analytical expression for the gradient is found instead of a numerical one. The easily readable application of the tensorial framework, however will apply integration by parts 'falsely' also to the terms where derivations in mean flow homogeneous directions are applied analytically. The inversion of the sign in front of the wave number cancels this effect.

    if T.containsTestFunction:
        mLocal = -m
    else:
        mLocal = m 
    # compute gradient: first check symmetry, then add christoffel parts
    if T.order == 0:
        diffs =  []
        mesh_iter = 0
        for i in range(T.dim):
            if T.sym[i] == 1:
                diffs.append(T.ufl_tens.dx(mesh_iter))
            else:
                if mLocal==0:
                    diffs.append(0.0)
                else:
                    diffs.append(1j * mLocal * T.ufl_tens) # Added wavenumber on homogeneous direction
            if T.CoordSys.mesh_dims[i] == 1:
                mesh_iter += 1
        gradient = as_vector(diffs)

    elif T.order == 1:
        # partial derivatives part
        diffs =  []
        for i in range(T.dim):
            row = []
            mesh_iter = 0
            for j in range(T.dim):
                if T.sym[i][j] == 1:
                    row.append(T.ufl_tens[i].dx(mesh_iter))
                else:
                    if mLocal==0:
                        row.append(0.0)
                    else:
                        row.append(1j*mLocal*T.ufl_tens[i]) # Added wavenumber on homogeneous direction
                if T.CoordSys.mesh_dims[j]:
                    mesh_iter += 1
            diffs.append(row)
        gradient = as_matrix(diffs)
        
        # christoffel part
        i, j, k = indices(3)
        gradient += as_tensor(T.ufl_tens[k]*T.CoordSys.ch[i,k,j], (i,j))
    
    elif T.order == 2:
        # partial derivatives part
        diffs =  []
        for i in range(T.dim):
            row = []
            for j in range(T.dim):
                column = []
                mesh_iter = 0
                for k in range(T.dim):
                    if T.sym[i][j][k] == 1:
                        column.append(T.ufl_tens[i,j].dx(mesh_iter))
                    else:
                        if mLocal==0:
                            column.append(0.0)
                        else:
                            column.append(1j*mLocal*T.ufl_tens[i,j]) # Added wavenumber on homogeneous direction
                    if T.CoordSys.mesh_dims[k] == 1:
                        mesh_iter += 1
                row.append(column)
            diffs.append(row)
            
        gradient = as_matrix(diffs)
        #return Tensor(gradient, T.CoordSys, basis = T.basis + [True])
        i, j, k, l = indices(4)
        gradient += \
            as_tensor(T.ufl_tens[l,j] * T.CoordSys.ch[i,l,k], (i,j,k)) + \
            as_tensor(T.ufl_tens[i,l] * T.CoordSys.ch[j,l,k], (i,j,k))
    
    else:
        raise ValueError("Gradients of tensors of order > 3 are not " \
                         "implemented.")
        
    return Tensor(
                  gradient, 
                  T.CoordSys, 
                  basis = T.basis + [True],
                  containsTestFunction = T.containsTestFunction,  
                  )
    

def iDiv(
         tensor: Tensor,
         m=0,
         ):
    # If the expression contains a test function, the wave numbers in mean flow homogeneous directions must be multiplied by minus one for the following reason: if a gradient operator contains a test function, integration by parts has been applied. With the wave numbers, however an analytical expression for the gradient is found instead of a numerical one. The easily readable application of the tensorial framework, however will apply integration by parts 'falsely' also to the terms where derivations in mean flow homogeneous directions are applied analytically. The inversion of the sign in front of the wave number cancels this effect.

    if tensor.order == 0:
        raise ValueError("Divergence of scalars not defined.")
        tensor = convertBasis(tensor, tensor.order*[False])
        
    elif tensor.order == 1:
        gradient = iGrad(tensor, m) # Added wavenumber on homogeneous direction
        i = indices(1)
        Div = tr(gradient.ufl_tens)
    
    elif tensor.order == 2:
        gradient = iGrad(tensor, m) # Added wavenumber on homogeneous direction
        i,j = indices(2)
        Div = as_tensor(gradient.ufl_tens[i,j,j], (i))
        
    return Tensor(
                  Div, 
                  tensor.CoordSys, 
                  basis = tensor.basis[:-1],
                  containsTestFunction = tensor.containsTestFunction
                  )


def iT(tensor: Tensor):
    if tensor.order != 2:
        raise ValueError("Transpose only unambiously defined for tensors of " \
                         "order 2.")
    else:
        # transposing by permuting the basis must not be performed (e.g. when
        # order=2: tensor.basis[::-1]), because iDot cant handle the change in
        # the order of contracted indices:
        # A_ij g^j dyade g^i dot v^k g_k = A_kj v^k g^j. But iDot can only do:
        # A_jk v^k --> neighboring contracted indices.
        return Tensor(
                      tensor.ufl_tens.T, 
                      tensor.CoordSys, 
                      basis =[tensor.basis[1], tensor.basis[0]],
                      containsTestFunction = tensor.containsTestFunction,
                      )


def iTr(tensor: Tensor):
    ''' Function returning trace '''
    # tr(tensor.ufl_tens) is invariant. This function is only for when you
    # need a Tensor object returned
    if tensor.order != 2:
        raise ValueError("Trace only working for tensors of order 2.")
    else:
        return Tensor(
                      tr(tensor.ufl_tens),
                      tensor.CoordSys, 
                      basis = [],
                      containsTestFunction = tensor.containsTestFunction,
                      )
    
    
def iDev(tensor: Tensor):
    spherical = Tensor(1/3 * tr(tensor.ufl_tens) * Identity(tensor.dim), \
                       tensor.CoordSys, basis = tensor.basis)
    return tensor - spherical


def convertBasis(
                 tensor: Tensor, 
                 goal_basis: list,
                 ):
    # You have to contract the tensor with the dual basis to the goal_basis:
    # A^i_j g_i dyade g^j dot dot g^k dyade g^l = A^k,l = covariant components
    if tensor.order == 1:
        metric = getMetric(tensor.basis[0], not(goal_basis[0]), tensor.CoordSys)
        
        i,j = indices(2)
        components = as_tensor(metric[i,j] * tensor.ufl_tens[j], (i))
    
    if tensor.order == 2:
        metric1 = getMetric(tensor.basis[-2], not(goal_basis[0]), tensor.CoordSys)
        metric2 = getMetric(tensor.basis[-1], not(goal_basis[1]), tensor.CoordSys)
        
        i,j,k,l = indices(4)
        components = as_tensor(tensor.ufl_tens[k,l] * metric1[k,i] * \
                               metric2[l,j], (i,j))
    
    return Tensor(
                  components, 
                  tensor.CoordSys, 
                  basis = goal_basis,
                  containsTestFunction = tensor.containsTestFunction,
                  )


def iIdentity(tensor: Tensor):
    return Tensor(
                  Identity(tensor.dim), 
                  tensor.CoordSys,
                  basis = tensor.basis,
                  containsTestFunction = False,
                  )

def iConj(tensor: Tensor):
    return Tensor(
            conj(tensor.ufl_tens),
            tensor.CoordSys,
            basis = tensor.basis,
            containsTestFunction = tensor.containsTestFunction,
            )


def iOuter(tensorA: Tensor, tensorB: Tensor):
    if tensorA.containsTestFunction and tensorB.containsTestFunction:
        raise ValueError("iOuter product at least second order in test functions.")
    if tensorA.order == 1 and tensorB.order == 1:
        i,j = indices(2)
        outered = as_tensor(tensorA.ufl_tens[i]*tensorB.ufl_tens[j], (i,j))

    # the below should not be necessary
    elif tensorA.order == 2 and tensorB.order == 1:
        i,j,k = indices(3)
        outered = as_tensor(tensorA.ufl_tens[i,j]*tensorB.ufl_tens[k], (i,j,k))
    elif tensorA.order == 1 and tensorB.order == 2:
        i,j,k = indices(3)
        outered = as_tensor(tensorA.ufl_tens[i]*tensorB.ufl_tens[i,k], (i,j,k))
    elif tensorA.order == 1 and tensorB.order == 2:
        i,j,k,l = indices(4)
        outered = as_tensor(tensorA.ufl_tens[i,j]*tensorB.ufl_tens[k,l], (i,j,k,l))
    
    elif tensorA.order == 0 or tensorB.order == 0:
        raise Exception('Use standard multiplication with the asterisk symbol *.')
        
    return Tensor(
                  outered, 
                  tensorA.CoordSys, 
                  basis = tensorA.basis + tensorB.basis,
                  containsTestFunction = tensorA.containsTestFunction or tensorB.containsTestFunction,
                  )
