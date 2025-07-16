#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Wed Dec 14 17:59:50 2022
# @author: kai hildebrandt
# Modifications: 
# - Thomas L. Kaiser
# - Simon Demange

from typing import Union
from ufl import TrialFunction, TestFunction
from ufl import indices
from ufl import Identity, as_vector, as_matrix, as_tensor
from ufl import sin, cos, tan, sqrt
from ufl import det, tr
from ufl import conj
from dolfinx.fem import (
    Constant,
)

### CLASSES
class CoordinateSystem():
    """
    The `CoordinateSystem` class supports several coordinate systems,
    initializing key geometric quantities such as the metric tensors
    and Christoffel symbols.

    Supports various systems such as Cartesian, polar, cylindrical, and spherical.
    Initializes geometric quantities like the metric tensors and Christoffel symbols
    required for tensor analysis.

    Parameters
    ----------
    SpatialCoordinateObj : ufl.SpatialCoordinate
        Coordinate vector of the mesh.
    name : str
        Name of the coordinate system ("cartesian", "polar", etc.).
    m : int, optional
        Wave number for mean-flow homogeneous directions.
    **kwargs : dict
        Optional keyword arguments, e.g. 'mesh_dims' to reduce dimensionality.
        mesh_dims : list
            This must be specified if a symmetry direction is solved for that
            is not meshed. Then mesh_dims is a list of booleans containing
            True if the direction of the chosen coordinate system is meshed
            and False if it is not meshed.
            
        


    Raises
    ------
    ValueError
        When object is initialized: if the coordinate system name is not recognized.


    Attributes
    -------

    dim : int
        Spatial dimension of the coordinate system. Not necessarily geometric
        dimension of the mesh.
    mesh_dims : list of int
        List indicating which coordinate directions are part of the mesh.
    x : ufl.Vector
        Coordinate vector in the system.
    ch : ufl.Tensor
        Christoffel symbols.
    cov_metric : ufl.Matrix
        Covariant metric tensor.
    con_metric : ufl.Matrix
        Contravariant metric tensor.
    g : ufl.Expr
        Determinant of the covariant metric tensor.
    J_hat : ufl.Expr
        Square root of the metric determinant.


    Notes
    -----

    Below are the supported systems and their conventions:

    **Polar Coordinates** (2D):
    Defined by:

    (x, y) = (r cos(ϕ), r sin(ϕ))

    ordering: (r, ϕ)

    Christoffel Symbols:

        ch1 = [[0, 0], [0, -r]]

        ch2 = [[0, 1/r], [1/r, 0]]

    Covariant Metric:

        [[1,     0],

        [0,   r²]]

    Contravariant Metric:

        [[1,     0],

        [0, 1/r²]]

    **Cylindrical Coordinates** (3D):
    Defined by:

    (x, y, z) = (r cos(ϕ), r sin(ϕ), z)

    ordering: (r, ϕ, z)

    Tangent Basis:

        g1 = cosϕ eₓ + sinϕ e_y

        g2 = r(-sinϕ eₓ + cosϕ e_y)

        g3 = e_z

    Christoffel Symbols:

        ch1 = [[0, 0, 0], [0, -r, 0], [0, 0, 0]]

        ch2 = [[0, 1/r, 0], [1/r, 0, 0], [0, 0, 0]]

        ch3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    Covariant Metric:

        [[1,     0, 0],
        [0,   r², 0],
        [0,     0, 1]]

    Contravariant Metric:

        [[1,     0,   0],
        [0, 1/r²,   0],
        [0,     0, 1]]

    **Cylindrical (FELiCS Convention)** (3D):
    Same coordinate mapping as cylindrical, but uses a different ordering:
    
    ordering: (z, r, ϕ)

    Tangent Basis:

        g1 = e_z

        g2 = cosϕ eₓ + sinϕ e_y

        g3 = r(-sinϕ eₓ + cosϕ e_y)

    Christoffel Symbols:

        ch1 = 0

        ch2 = [[0, 0, 0], [0, 0, 0], [0, 0, -r]]

        ch3 = [[0, 0, 0], [0, 0, 1/r], [0, 1/r, 0]]

    Covariant Metric:

        [[1, 0,     0],

        [0, 1,     0],

        [0, 0,   r²]]

    Contravariant Metric:

        [[1, 0,     0],

        [0, 1,     0],

        [0, 0, 1/r²]]

    **Spherical Coordinates** (3D):
    Defined by:

    (x, y, z) = (r sinθ cosϕ, r sinθ sinϕ, r cosθ)

    ordering: (r, θ, ϕ)

    Christoffel Symbols:

        ch1 = [[0, 0, 0], [0, -r, 0], [0, 0, -r sin²θ]]

        ch2 = [[0, 1/r, 0], [1/r, 0, 0], [0, 0, -sinθ cosθ]]

        ch3 = [[0, 0, 1/r], [0, 0, 1/tanθ], [1/r, 1/tanθ, 0]]

    Covariant Metric:

        [[1,     0,           0],

        [0,   r²,           0],

        [0,     0, r² sin²θ]]

    Contravariant Metric:

        [[1,      0,              0],

        [0,   1/r²,              0],

        [0,      0, 1/(r² sin²θ)]]


    Example
    -------

    **Metric Determinant and Jacobian Weight**
    For any coordinate system, the metric determinant and Jacobian weight
    are defined as:

    >>> g = det(cov_metric)
    >>> J_hat = sqrt(g)

    These are used to properly scale integrals in variational forms.

    **Accessing Christoffel Symbols in UFL**
    UFL supports index notation for tensors:
    
    >>> ch[2,:,0]  # First column of the third Christoffel tensor

    **Example for index operations on Matrices in UFL.**

    >>> ch1 = [[1, 2, 3],[4,5,6], [7,8,9]]
    >>> ch2 = [[10, 11, 12],[13,14,15], [16,17,18]]
    >>> ch3 = [[19, 20, 21],[22,23,24], [25,26,27]]
    >>> ch = as_matrix([ch1, ch2, ch3])
    >>> print(ch3)
    >>> print(ch[2,:,0]) # first column of the last christoffel symbol

    """


    def __init__(
                self, 
                SpatialCoordinateObj, 
                name: str,
                m = 0, 
                **kwargs,
                ):
        """
        Initialize the CoordinateSystem object.

        Parameters
        ----------
        SpatialCoordinateObj : ufl.SpatialCoordinate
            Coordinate vector of the mesh.
        name : str
            Name of the coordinate system.
        m : int, optional
            Wave number for mean-flow homogeneous directions.
        **kwargs : dict
            Optional keyword arguments, e.g. 'mesh_dims' to reduce dimensionality.

        Raises
        ------
        ValueError
            If the coordinate system name is not recognized.
        """

        x = SpatialCoordinateObj
        
        # mesh_dims is a variable to handle that the ufl vector x[i] might re-
        # present for example (r,z) in cylindrical coordinates but
        # tensor_utils can only operate with (r, theta, z).
        self._m = m
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
            

    @property
    def m(self):
        """
        Returns the wave number `m`, used in homogeneous directions.

        Returns
        -------
        int
            The wave number.
        """
        return self._m

class Tensor():
    """
    Tensor object that extends the UFL tensors by bases vectors, such
    that they are valid / definable in coordinate systems other than
    the cartesian one.

    Supports scalar, vector, and matrix-valued tensors in both physical
    and tangent bases.

    Attributes
    ----------
    ufl_tens : ufl.Expr
        The underlying UFL tensor expression.
    CoordSys : CoordinateSystem
        Coordinate system in which this tensor is defined.
    order : int
        Order of the tensor (0: scalar, 1: vector, 2: matrix).
    basis : list of booleans
        A list of booleans in which False represents a covariante basis
        vector and True a contravariant basis vector.
    sym : list
        A list of order self.order + 1, containing booleans. When computing
        the gradient of the tensor this information is used to determine
        if the partial derivative in this direction exists.
    hasSpectralDimension : bool
        Whether this tensor has spectral dimension sdim, i.e. tensor = tensor_coefficients * exp(i*m*sdim) with wave number m
    m : int 
        wave number in spectral dimension
        

    Example
    -------
    >>> T = Tensor(u, CoordSys)
    >>> grad_T = iGrad(T)
    """ 

    def __init__(
        self, 
        ufl_tens, 
        CoordSys: CoordinateSystem, 
        hasSpectralDimension = False, 
        m = None,
        **kwargs,
        ):
        """
        Initialize a Tensor object.

        Parameters
        ----------
        ufl_tens : ufl.Expr
            The UFL expression representing the tensor.
        CoordSys : CoordinateSystem
            The coordinate system in which the tensor is defined.
        hasSpectralDimension : bool, optional
            Whether this tensor has spectral dimension sdim, i.e. tensor = tensor_coefficients * exp(i*m*sdim) with wave number m.
        m : int, optional
            Optional wave number.
        **kwargs : dict
            Optional arguments: 'basis' (covariant or contravariant), 'sym' (custom symmetry).

        Notes
        -----
        - Automatically converts physical basis to tangent basis for vectors.
        - Derives symbolic symmetry information if not explicitly provided.
        - Only tensors of order 0, 1, and 2 are currently supported.

        Raises
        ------
        ValueError
            If the tensor order is >= 2 and no transformation rule is provided.
        """

        self.CoordSys = CoordSys
        self.x = CoordSys.x
        self.dim = CoordSys.dim
        self.order = len(ufl_tens.ufl_shape) # scalar --> order = 0
        self.hasSpectralDimension = hasSpectralDimension
        if self.hasSpectralDimension:
            if m == None:
                self.m = self.CoordSys.m
            else:
                self.m = m
        else:
            self.m = 0
         
        # Sophie: We could check here if ufl_tens contains a test- or trialfunction and make the error messages clearer. 
        #         E.g. we could do something in the direction of: 
        #         for arg in ufl_tens.arguments():
        #             if isinstance(arg, TestFunction):
        #                 self.containsTestFunction = True
        #             if isinstance(arg, TrialFunction):
        #                 self.containsTestFunction = True
 

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
        """
        Add two tensors or a tensor and a scalar.

        Parameters
        ----------
        other : Tensor or scalar
            The tensor or scalar to add.

        Returns
        -------
        Tensor
            Sum of self and other.

        Raises
        ------
        ValueError
            If operand types are not compatible.

        Example
        -------
        >>> A + B  # where A is a UFL or Tensor object and B is a Tensor object.
        """
        if type(other) == Tensor:
            if (self.hasSpectralDimension and other.hasSpectralDimension) and (self.m != other.m):
                ValueError("In a tensor sum, both tensors must have the same wave number m, if they contain a spectral dimension.")
            if self.basis == other.basis:
                added = self.ufl_tens + other.ufl_tens
            else:
                #new_self = convertBasis(self, other.basis)
                #added = new_self.ufl_tens + other.ufl_tens
                ### changed by Sophie to make the iT (transpose) functionality work properly
                new_other = convertBasis(other, self.basis)
                added     = self.ufl_tens + new_other.ufl_tens
        if type(other) in [float,complex,int,Constant]:
            added = self.ufl_tens + other
        else:
            ValueError("Tensor addition only defined for Tensors,Constant, float, complex and integers")
        return Tensor(
                    added, 
                    self.CoordSys, 
                    basis = self.basis, 
                    hasSpectralDimension = self.hasSpectralDimension,
                    m = self.m
                    )

    # division, 
    def __truediv__(self, other): # Tensor object to the left
        """
        Divide this tensor by another tensor or scalar.

        Parameters
        ----------
        other : Tensor or scalar
            The denominator.

        Returns
        -------
        Tensor
            Result of the division.

        Raises
        ------
        ValueError
            If division is not defined for the operand types.

        Example
        -------
        >>> A / 2.0
        """
        if type(other) == Tensor:
            if other.order == 0:
                return Tensor(
                            self.ufl_tens / other.ufl_tens,
                            self.CoordSys, 
                            basis = self.basis, 
                            hasSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                            m = self.m - other.m
                            )
            elif self.order == 0:
                return Tensor(
                            self.ufl_tens / other.ufl_tens,
                            self.CoordSys,
                            basis = other.basis,
                            hasSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                            m = self.m - other.m
                            )
            else:
                raise ValueError("Division operation between Tensors only " \
                                 "defined, if one is a scalar.")
        if type(other) in [float,complex,int,Constant]:
            return Tensor(
                        self.ufl_tens / other,
                        self.CoordSys,
                        basis = self.basis,
                        hasSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                        m = self.m - other.m
                        )
        else:
            ValueError("Tensor division only defined for divisors of type Tensor, Constant, float, complex and integer")
    
    # subtraction: A - B is the same as A.__sub__(B)
    def __sub__(self, other):
        """
        Subtract one tensor from another or from a scalar.

        Parameters
        ----------
        other : Tensor or scalar
            The tensor or scalar to subtract.

        Returns
        -------
        Tensor
            Result of subtraction.

        Raises
        ------
        ValueError
            If operand types are not compatible.

        Example
        -------
        >>> A - B # where A is a UFL or Tensor object and B is a Tensor object.
        """
        if type(other) == Tensor:
            if (self.hasSpectralDimension and other.hasSpectralDimension) and (self.m != other.m):
                ValueError("In a tensor difference, both tensors must have the same wave number m, if they contain a spectral dimension.")

            if self.basis == other.basis:
                subtracted = self.ufl_tens - other.ufl_tens
            else:
                #new_self = convertBasis(self, other.basis)
                #subtracted = new_self.ufl_tens - other.ufl_tens
                ### changed by Sophie to make the iT (transpose) functionality work properly
                new_other = convertBasis(other, self.basis)
                subtracted = self.ufl_tens - new_other.ufl_tens
        if type(other) in [float,complex,int,Constant]:
            subtracted = self.ufl_tens - other
        else:
            ValueError("Tensor subtraction only defined for Tensors, float, complex and integer")
        
        return Tensor(
                    subtracted, 
                    self.CoordSys, 
                    basis = self.basis, 
                    hasSpectralDimension = self.hasSpectralDimension,
                    m = self.m
                    )
    # muliplication, both ways, because matrix mul not commutative
    def __mul__(self, other): # Tensor object to the left
        """
        Multiply this tensor with another tensor or scalar.

        Parameters
        ----------
        other : Tensor or scalar
            The operand on the right-hand side.

        Returns
        -------
        Tensor
            The product tensor.

        Raises
        ------
        ValueError
            If multiplication is undefined for operand types.

        Example
        -------
        >>> A * 2.0
        >>> A * B # where A is a UFL or Tensor object and B is a Tensor object.
        """
        if type(other) == Tensor:
            if other.order == 0:
                return Tensor(
                            self.ufl_tens * other.ufl_tens,
                            self.CoordSys, 
                            basis = self.basis, 
                            hasSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                            m = self.m + other.m 
                            )
            elif self.order == 0:
                return Tensor(
                            self.ufl_tens * other.ufl_tens,
                            self.CoordSys,
                            basis = other.basis,
                            hasSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                            m = self.m + other.m 
                            )
            else:
                raise ValueError("Multiplication between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) in [float,complex,int,Constant]:
            return Tensor(
                        self.ufl_tens * other,
                        self.CoordSys,
                        basis = self.basis,
                        hasSpectralDimension = self.hasSpectralDimension,
                        m = self.m
                        )
        else:
            ValueError("Tensor multiplication only defined for Tensors, Constant, float, complex, and integer")

    def __pow__(self, exponent): # Tensor object to the left
        """
        Returns this tensor to the power of the exponent.

        Parameters
        ----------
        exponent : scalar
            The exponent. 

        Returns
        -------
        Tensor
            The tensor, has to be of order 0.

        Raises
        ------
        ValueError
            If the tensor is not order 0 or if the exponent is not a scalar number.

        Example
        -------
        >>> A ** 2
        """
        import numpy as np
        if self.order ==0 and np.isscalar(exponent):
             return Tensor(
                       self.ufl_tens ** exponent,
                       self.CoordSys, 
                       basis = self.basis, 
                       hasSpectralDimension = self.hasSpectralDimension,
                       m = self.m * exponent
                       )
        else:
            raise ValueError("Taking exponents is only defined if the exponent is a scalar number "\
                                 "and the tensor has order 0.")
 



    def __rmul__(self, other): # Tensor object to the right
        """
        Multiply scalar or tensor from the left.

        Parameters
        ----------
        other : Tensor or scalar
            The operand on the left-hand side.

        Returns
        -------
        Tensor
            The product tensor.

        Raises
        ------
        ValueError
            If multiplication is undefined for operand types.

        Example
        -------
        >>> 3.0 * A
        """
        if type(other) == Tensor:
            if other.order == 0:
                return Tensor(
                            other.ufl_tens * self.ufl_tens, 
                            self.CoordSys,
                            basis = self.basis,
                            hasSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                            m = self.m + other.m
                            )
            elif self.order == 0:
                return Tensor(
                            other.ufl_tens * self.ufl_tens, 
                            self.CoordSys,
                            basis = other.basis,
                            hasSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                            m = self.m + other.m
                            )
            else:
                raise ValueError("Multiplication between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) in [float,complex,int,Constant]:
            return Tensor(
                        other * self.ufl_tens,
                        self.CoordSys,
                        basis = self.basis,
                        hasSpectralDimension = self.hasSpectralDimension,
                        m = self.m
                        )
            ValueError("Tensor multiplication only defined for Tensors, Constant, float, complex and integerr")
    
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
    """
    Performs a single contraction between two tensors using the appropriate metric.

    The function supports dot products between tensors of order 1 or 2 and uses 
    the metric determined by the coordinate system and the tensor bases.

    Parameters
    ----------
    tensorA : Tensor
        First tensor operand.
    tensorB : Tensor
        Second tensor operand.

    Returns
    -------
    Tensor
        Result of the dot product, with updated basis and metadata.

    Raises
    ------
    ValueError
        If both tensors contain test functions or fluctuations.
    """
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
            hasSpectralDimension = tensorA.hasSpectralDimension or tensorB.hasSpectralDimension, 
            m = tensorA.m + tensorB.m
            )


# Kai TODO: Is this validated? Which of the two tensors is transposed? What is the benefit over iDot(iT(TensorA), TensorB)?
# It is not part of my version of tensor_utils, just FYI for the others.
def iDotT(tensorA: Tensor, tensorB: Tensor):
    """
    Performs an intrinsic dot product between the transpose of a second-order tensor and a first-order tensor.

    Parameters
    ----------
    tensorA : Tensor
        A second-order tensor (not explicitly transposed).
    tensorB : Tensor
        A first-order tensor.

    Returns
    -------
    Tensor
        Resulting tensor from the intrinsic dot product.

    Raises
    ------
    ValueError
        If both tensors contain test functions or fluctuations.
    Exception
        If the combination of tensor orders is not implemented.
    """
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
            hasSpectralDimension = tensorA.hasSpectralDimension or tensorB.hasSpectralDimension, 
            m = tensorA.m + tensorB.m
            )


def iInner(tensorA: Tensor, tensorB: Tensor):
    """
    Computes the inner product between two second-order tensors.

    Parameters
    ----------
    tensorA : Tensor
        First tensor operand of order 2.
    tensorB : Tensor
        Second tensor operand of order 2.

    Returns
    -------
    Tensor
        Scalar-valued (zero order) tensor representing the inner product.

    Raises
    ------
    ValueError
        If the tensors are not of order 2 or contain invalid combinations of test functions or fluctuations.
    """
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
            hasSpectralDimension = tensorA.hasSpectralDimension or tensorB.hasSpectralDimension, 
            m = tensorA.m + tensorB.m
            )


def getMetric(basisA: Union[list, bool, int], basisB: Union[list, bool, int], \
              CoordSys: CoordinateSystem):
    """
    Determines the appropriate metric based on basis vectors.

    Parameters
    ----------
    basisA : list, bool, or int
        Basis type of the first tensor.
    basisB : list, bool, or int
        Basis type of the second tensor.
    CoordSys : CoordinateSystem
        Coordinate system providing the metric tensors.

    Returns
    -------
    ufl.Expression
        Appropriate metric tensor for contracting the basis vectors.
    """
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


def iGrad(T: Tensor):
    """
    Computes the gradient of a tensor.

    Converts bases to tangent basis and then applies the gradient including terms
    arising from the Christoffel symbols. Adjusts signs for spectral derivatives
    in presence of test functions or fluctuations.

    Parameters
    ----------
    T : Tensor
        The tensor to differentiate.

    Returns
    -------
    Tensor
        Gradient of the input tensor, with one additional contravariant basis.

    Raises
    ------
    ValueError
        If the tensor is of order > 2 (not implemented).
    """
    # implementing all kind of derivatives is tideous. --> convert every tensor
    # to tangent basis and only implement grads in tangent basis.
    if sum(T.basis) != 0:
        T = convertBasis(T, T.order*[False])

    # If the expression contains a test function, the wave numbers in mean flow homogeneous directions must be multiplied by minus one for the following reason: if a gradient operator contains a test function, integration by parts has been applied. With the wave numbers, however an analytical expression for the gradient is found instead of a numerical one. The easily readable application of the tensorial framework, however will apply integration by parts 'falsely' also to the terms where derivations in mean flow homogeneous directions are applied analytically. The inversion of the sign in front of the wave number cancels this effect.

    # compute gradient: first check symmetry, then add christoffel parts
    if T.order == 0:
        diffs =  []
        mesh_iter = 0
        for i in range(T.dim):
            if T.sym[i] == 1:
                diffs.append(T.ufl_tens.dx(mesh_iter))
            else:
                if T.m==0:
                    diffs.append(0.0)
                else:
                    diffs.append(1j * T.m * T.ufl_tens) # Added wavenumber on homogeneous direction
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
                    if T.m==0:
                        row.append(0.0)
                    else:
                        row.append(1j*T.m*T.ufl_tens[i]) # Added wavenumber on homogeneous direction
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
                        if T.m==0:
                            column.append(0.0)
                        else:
                            column.append(1j*T.m*T.ufl_tens[i,j]) # Added wavenumber on homogeneous direction
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
                  hasSpectralDimension = T.hasSpectralDimension,
                  m = T.m,
                  )
    

def iDiv(tensor: Tensor):
    """
    Computes the divergence of a tensor.

    Computes the gradient and sums of the last two indices of a tensor.
    Spectral treatment is applied if the tensor includes test functions.
    # TODO Kai: "spectral treatment is applied" is ambiguous to me.

    Parameters
    ----------
    tensor : Tensor
        Tensor of order 1 or 2.

    Returns
    -------
    Tensor
        The divergence result as a tensor of order reduced by one reduced.

    Raises
    ------
    ValueError
        If called on a scalar tensor (order 0).
    """
    # If the expression contains a test function, the wave numbers in mean flow homogeneous directions must be multiplied by minus one for the following reason: if a gradient operator contains a test function, integration by parts has been applied. With the wave numbers, however an analytical expression for the gradient is found instead of a numerical one. The easily readable application of the tensorial framework, however will apply integration by parts 'falsely' also to the terms where derivations in mean flow homogeneous directions are applied analytically. The inversion of the sign in front of the wave number cancels this effect.
    if tensor.order == 0:
        raise ValueError("Divergence of scalars not defined.")
        tensor = convertBasis(tensor, tensor.order*[False])
        
    elif tensor.order == 1:
        gradient = iGrad(tensor) # Added wavenumber on homogeneous direction
        i = indices(1)
        Div = tr(gradient.ufl_tens)
    
    elif tensor.order == 2:
        gradient = iGrad(tensor) # Added wavenumber on homogeneous direction
        i,j = indices(2)
        Div = as_tensor(gradient.ufl_tens[i,j,j], (i))
        
    return Tensor(
                  Div, 
                  tensor.CoordSys, 
                  basis = tensor.basis[:-1],
                  hasSpectralDimension = tensor.hasSpectralDimension,
                  m = tensor.m
                  )


def iT(tensor: Tensor):
    """
    Returns the transpose of a second-order tensor.

    Parameters
    ----------
    tensor : Tensor
        A tensor of order 2.

    Returns
    -------
    Tensor
        Transposed tensor with permuted bases.

    Raises
    ------
    ValueError
        If tensor order is not 2.
    """
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
                      hasSpectralDimension = tensor.hasSpectralDimension,
                      m = tensor.m
                      )


def iTr(tensor: Tensor):
    """
    Computes the trace of a second-order tensor.

    Parameters
    ----------
    tensor : Tensor
        Tensor of order 2.

    Returns
    -------
    Tensor
        Scalar-valued tensor (trace result).

    Raises
    ------
    ValueError
        If the input tensor is not of order 2.
    """
    # tr(tensor.ufl_tens) is invariant. This function is only for when you
    # need a Tensor object urned
    if tensor.order != 2:
        raise ValueError("Trace only working for tensors of order 2.")
    else:
        return Tensor(
                      tr(tensor.ufl_tens),
                      tensor.CoordSys, 
                      basis = [],
                      hasSpectralDimension = tensor.hasSpectralDimension,
                      m = tensor.m
                      )
    
    
def iDev(tensor: Tensor):
    """
    Computes the deviatoric (spherical-free) part of a second-order tensor.

    Parameters
    ----------
    tensor : Tensor
        A second-order tensor.

    Returns
    -------
    Tensor
        The deviatoric component of the input tensor.
    """
    spherical = Tensor(1/3 * tr(tensor.ufl_tens) * Identity(tensor.dim), \
                       tensor.CoordSys, basis = tensor.basis)
    return tensor - spherical


def convertBasis(
                 tensor: Tensor, 
                 goal_basis: list,
                 ):
    """
    Converts a tensor to a desired basis by contracting with the appropriate metric.
    
    A^i_j g_i dyade g^j : g^k dyade g^l = A^k,l = covariant components

    Parameters
    ----------
    tensor : Tensor
        Tensor to convert.
    goal_basis : list
        Desired basis for the tensor components.

    Returns
    -------
    Tensor
        The tensor expressed in the new basis.
    """
    # You have to contract the tensor with the dual basis to the goal_basis:
    # A^i_j g_i dyade g^j : g^k dyade g^l = A^k,l = covariant components
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
                  hasSpectralDimension = tensor.hasSpectralDimension,
                  m = tensor.m
                  )


def iIdentity(tensor: Tensor):
    """
    Returns the identity tensor corresponding to the tensor's dimension and coordinate system.

    Parameters
    ----------
    tensor : Tensor
        Input tensor, used to infer dimensionality and coordinate system.

    Returns
    -------
    Tensor
        Identity tensor with the same coordinate system and basis.
    """
    return Tensor(
                  Identity(tensor.dim), 
                  tensor.CoordSys,
                  basis = tensor.basis,
                  hasSpectralDimension = False,
                  m = 0
                  )

def iConj(tensor: Tensor):
    """
    Computes the complex conjugate of a tensor.

    Parameters
    ----------
    tensor : Tensor
        Tensor to conjugate.

    Returns
    -------
    Tensor
        Conjugated tensor with identical metadata.
    """
    return Tensor(
            conj(tensor.ufl_tens),
            tensor.CoordSys,
            basis = tensor.basis,
            hasSpectralDimension = tensor.hasSpectralDimension,
            m = -tensor.m
            )

# TODO Kai: is this now validated? Otherwise it should raise a warning.
def iOuter(tensorA: Tensor, tensorB: Tensor):
    """
    Computes the outer product of two tensors.

    Parameters
    ----------
    tensorA : Tensor
        Left operand tensor.
    tensorB : Tensor
        Right operand tensor.

    Returns
    -------
    Tensor
        Outer product tensor with combined bases.

    Raises
    ------
    ValueError
        If both tensors contain test functions or fluctuations.
    Exception
        If scalar multiplication is attempted using this function.
    """
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
                  hasSpectralDimension = tensorA.hasSpectralDimension or tensorB.hasSpectralDimension,
                  m = tensorA.m + tensorB.m
                  ) 
    
