#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

# This is a new tensor utils file, which is using the method and class names from Kai Hildebrandt,
# but has a very simple ("hands-on") handling of the two FELiCS coordinate systems.
# This is to eliminate the numerous bugs which the other version had in the cylindrical 
# coordinate system.
# Written by Sophie Knechtel

# Standard libraries
from enum   import Enum
from typing import Union

# Third party libraries
from    dolfinx.fem import (
    Constant,
)
import  numpy       as np
from    ufl         import (
    as_matrix, 
    as_tensor,
    as_vector,
    det,
    conj,
    dot,
    Identity,
    indices,
    inner, 
    sin, cos, tan, sqrt,
    tr,
    TrialFunction,
    TestFunction,
) 

class SpectralIndicator(Enum):
    """
    Enumeration describing the role of each dimension in a tensor or coordinate system.

    This enum is used to mark whether a given axis corresponds to a geometric,
    spectral, or inactive (zero) dimension. It is primarily used internally by
    the :class:`Tensor` class to determine how derivatives and spectral terms
    should be handled.

    Members
    -------
    SPECTRAL : int
        Marks a dimension as a spectral dimension. Derivatives along this axis
        are replaced by ``i * m / r`` terms.
    NOTSPECTRAL : int
        Marks a dimension as a regular geometric dimension.
    ZERO : int
        Marks a dimension that is neither geometric nor spectral (e.g. padding
        dimensions for tensors whose true dimension exceeds the geometric one).
    """
    # for handling spectral dimensions
    SPECTRAL     = 1 # this is a spectral dimension
    NOTSPECTRAL  = 2 # this is not a spectral dimension but still a geometric dimension
    ZERO         = 3 # this is not a spectral dimension and not a geometric dimension

    

class CoordinateSystem():
    """
    The `CoordinateSystem` class supports two coordinate systems,
    'cartesian' and 'cylindricalfelics'.

    Parameters
    ----------
    SpatialCoordinateObj : ufl.SpatialCoordinate
        Coordinate vector of the mesh.
    name : str
        Name of the coordinate system ("cartesian", "cylindricalfelics").
    gdim : int
        Geometrical dimension of the mesh.
    trueDim : int, optional
        True dimension of the system, including a possible spectral dimension.
    m : int, optional
        Wave number for spectral direction.
        

    Raises
    ------
    ValueError
        When object is initialized: if the coordinate system name is not recognized.


    """


    def __init__(
        self,
        SpatialCoordinateObj,
        name: str,
        gdim: int,
        trueDim = None,
        m = 0
    ):
        """
        Initialize the CoordinateSystem object.

        Parameters
        ----------
        SpatialCoordinateObj : ufl.SpatialCoordinate
            Coordinate vector of the mesh.
        name : str
            Name of the coordinate system
            (currently ``"cartesian"`` or ``"cylindricalfelics"``).
        gdim : int
            Geometric dimension of the mesh (excluding any spectral
            dimensions).
        trueDim : int, optional
            Total dimension of the system, including a possible spectral
            dimension. If ``None`` and ``m == 0``, the true dimension is
            initially set to ``gdim`` and can be corrected later using
            :meth:`setTrueDimension`. If ``None`` and ``m != 0``, the
            true dimension is set to ``gdim + 1``.
        m : int, optional
            Wave number for the spectral direction. Default is 0.

        Raises
        ------
        ValueError
            If the coordinate system name is not recognized.
        """

        x = SpatialCoordinateObj
        
        self._m   = m
        self.gdim = gdim
        if trueDim == None and m==0:
            self.dim  = gdim #should be updated later by using "setTrueDim" if there are any spectral dimensions
        elif trueDim == None and m!=0:
            self.dim  = gdim + 1
        else:
            self.dim  = trueDim
        self.x    = x

        self.name = name
        if name == "cartesian":
            self.J_hat = 1.
        
        elif name == "cylindricalfelics":

            #Given for {x,y,z} = {r\cos\phi, r\sin\phi, z}, where ordering is
            #{z, r, \phi}.
            self.J_hat = self.x[1]

        else:
            raise ValueError("The specified coordinate system isn't " \
                             "implemented.")

    def set_true_dimension(
        self,
        dim
    ):
        """
        Corrects the true dimension of the system. This has to be called if a spectral dimension is included,
        in which case the true dimension is higher than the geometrical dimension of the mesh.
         
        Parameters
        ----------
        dim: int 
            True dimension of the system, including a possible spectral dimension. 

        """

        self.dim = dim

    @property
    def m(self):
        """
        Returns the wave number `m`, used in the spectral dimension.

        Returns
        -------
        int
            The wave number.
        """
        return self._m


class Tensor():
    """
    Tensor object that extends UFL tensors.

    Supports scalar, vector, and matrix-valued tensors in both physical
    and spectral settings, including an optional spectral dimension.

    Attributes
    ----------
    ufl_tens : ufl.Expr
        The underlying UFL tensor expression.
    CoordSys : CoordinateSystem
        Coordinate system in which this tensor is defined.
    hasSpectralDimension : bool
        Whether this tensor actually carries a spectral dependence
        (i.e. a factor of ``exp(i * m * sdim)``), as determined from
        ``mayHaveSpectralDimension`` and the coordinate system.
    isSpectralDimension : list of SpectralIndicator
        Per-dimension indicator specifying whether each axis is spectral,
        geometric, or zero (non-geometric, non-spectral).
    m : int
        Wave number associated with the spectral dimension (if present).
    dim : int
        Total dimension of the system (geometric + spectral).
    order : int
        Tensor order, inferred from ``ufl_tens.ufl_shape`` (0, 1, or 2).

    Example
    -------
    >>> T = Tensor(u, CoordSys)
    >>> grad_T = iGrad(T)
    """

    def __init__(
        self,
        ufl_tens,
        CoordSys: CoordinateSystem,
        mayHaveSpectralDimension = False,
        m = None
    ):
        """
        Initialize a Tensor object.

        Parameters
        ----------
        ufl_tens : ufl.Expr
            The UFL expression representing the tensor.
        CoordSys : CoordinateSystem
            The coordinate system in which the tensor is defined.
        mayHaveSpectralDimension : bool, optional
            Whether this tensor could have a spectral dimension sdim, i.e. tensor = tensor_coefficients * exp(i*m*sdim) with wave number m (should be set to "True" for all trial and test functions)
        m : int, optional
            Optional wave number.

        Notes
        -----
        - Only tensors of order 0, 1, and 2 are currently supported.

        """

        self.CoordSys = CoordSys
        self.x        = CoordSys.x
        self.dim      = CoordSys.dim 
        self.order    = len(ufl_tens.ufl_shape) # scalar --> order = 0

        # handling a possible spectral dimension
        self.hasSpectralDimension = (mayHaveSpectralDimension and CoordSys.gdim < CoordSys.dim) #this is also true for m=0 
        self.isSpectralDimension  = [SpectralIndicator.NOTSPECTRAL]*self.dim
        if self.hasSpectralDimension:
            self.isSpectralDimension[-1] = SpectralIndicator.SPECTRAL
            if m == None:
                self.m = self.CoordSys.m
            else:
                self.m = m
        else:
            if CoordSys.gdim < CoordSys.dim:
                self.isSpectralDimension[-1] = SpectralIndicator.ZERO
            self.m = 0

        # CAUTION: the following is valid for cartesian and cylindrical felics; for any new coordinatesystem this should be adjusted
        self.r                    = CoordSys.J_hat          

        self.ufl_tens = ufl_tens

        # this is for the case that the true dimension of the system is bigger than the geometric dimension: 
        # if a 2-dimensional vector is given but the true dimension is 3, then a third dimension is added with zero values
        if self.order ==1:
            length = ufl_tens.ufl_shape[0]
            if length < self.dim:
                self.ufl_tens = as_vector((ufl_tens[0],ufl_tens[1],0.0))
           

    # addition
    def __add__(
        self,
        other
    ):
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
            added = self.ufl_tens + other.ufl_tens
        if type(other) in [float,complex,int,Constant]:
            added = self.ufl_tens + other
        else:
            ValueError("Tensor addition only defined for Tensors,Constant, float, complex and integers")
        return Tensor(
            added,
            self.CoordSys,
            mayHaveSpectralDimension = self.hasSpectralDimension,
            m = self.m,
        )

    # division, 
    def __truediv__(
        self,
        other,
    ): # Tensor object to the left
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
            if other.order == 0 or self.order == 0:
                return Tensor(
                    self.ufl_tens / other.ufl_tens,
                    self.CoordSys,
                    mayHaveSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                    m = self.m - other.m,
                )
            else:
                raise ValueError("Division operation between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) in [float,complex,int,Constant]:
            return Tensor(
                self.ufl_tens / other,
                self.CoordSys,
                mayHaveSpectralDimension = self.hasSpectralDimension,
                m = self.m,
            )
        else:
            ValueError("Tensor division only defined for divisors of type Tensor, Constant, float, complex and integer")


    # division, 
    def __rtruediv__(
        self,
        other
    ): # Tensor object to the left
        """
        Divide a scalar or tensor by this tensor.

        This implements the reflected division operator ``other / self``,
        i.e. it is called when the left-hand operand does not know how to
        divide by a :class:`Tensor`.

        Parameters
        ----------
        other : Tensor or scalar
            The numerator in the division ``other / self``.

        Returns
        -------
        Tensor
            Result of the division.

        Raises
        ------
        ValueError
            If division is not defined for the operand types or if both
            operands are non-scalar tensors.

        Example
        -------
        >>> 2.0 / A
        >>> B / A  # where A and B are Tensor objects and one is scalar.
        """
        if type(other) == Tensor:
            if other.order == 0 or self.order == 0:
                return Tensor(
                    other.ufl_tens / self.ufl_tens,
                    self.CoordSys,
                    mayHaveSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                    m = self.m - other.m,
                )
            else:
                raise ValueError("Division operation between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) in [float,complex,int,Constant]:
            return Tensor(
                other / self.ufl_tens,
                self.CoordSys,
                mayHaveSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                m = -self.m,
            )
        else:
            ValueError("Tensor division only defined for divisors of type Tensor, Constant, float, complex and integer")
    




    # subtraction: A - B is the same as A.__sub__(B)
    def __sub__(self,
    other
    ):
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

            subtracted = self.ufl_tens - other.ufl_tens
        if type(other) in [float,complex,int,Constant]:
            subtracted = self.ufl_tens - other
        else:
            ValueError("Tensor subtraction only defined for Tensors, float, complex and integer")
        
        return Tensor(
        subtracted,
        self.CoordSys,
        mayHaveSpectralDimension = self.hasSpectralDimension,
        m = self.m,
        )


    # muliplication, both ways, because matrix mul not commutative
    def __mul__(
        self,
        other
    ): # Tensor object to the left
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
                mayHaveSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                m = self.m + other.m,
                )
            elif self.order == 0:
                return Tensor(
                self.ufl_tens * other.ufl_tens,
                self.CoordSys,
                mayHaveSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                m = self.m + other.m,
                )
            else:
                raise ValueError("Multiplication between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) in [float,complex,int,Constant]:
            return Tensor(
            self.ufl_tens * other,
            self.CoordSys,
            mayHaveSpectralDimension = self.hasSpectralDimension,
            m = self.m,
            )
        else:
            ValueError("Tensor multiplication only defined for Tensors, Constant, float, complex, and integer")

    def __pow__(
        self,
        exponent
    ): # Tensor object to the left
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
        if self.order ==0 and np.isscalar(exponent):
             return Tensor(
             self.ufl_tens ** exponent,
             self.CoordSys,
             mayHaveSpectralDimension = self.hasSpectralDimension,
             m = self.m * exponent,
             )
        else:
            raise ValueError("Taking exponents is only defined if the exponent is a scalar number "\
                                 "and the tensor has order 0.")
 



    def __rmul__(
        self,
        other
    ): # Tensor object to the right
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
                mayHaveSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                m = self.m + other.m,
                )
            elif self.order == 0:
                return Tensor(
                other.ufl_tens * self.ufl_tens,
                self.CoordSys,
                mayHaveSpectralDimension = self.hasSpectralDimension or other.hasSpectralDimension,
                m = self.m + other.m,
                )
            else:
                raise ValueError("Multiplication between Tensors only " \
                                 "defined, if one is a scalar.")
        elif type(other) in [float,complex,int,Constant]:
            return Tensor(
            other * self.ufl_tens,
            self.CoordSys,
            mayHaveSpectralDimension = self.hasSpectralDimension,
            m = self.m,
            )
            ValueError("Tensor multiplication only defined for Tensors, Constant, float, complex and integerr")
    


### TENSOR OBJECT FUNCTIONS
def i_dot(
         tensorA: Tensor,
         tensorB: Tensor,
):
    """
    Performs a single contraction between two tensors.

    The function supports dot products between tensors of order 1 or 2.

    Parameters
    ----------
    tensorA : Tensor
        First tensor operand.
    tensorB : Tensor
        Second tensor operand.

    Returns
    -------
    Tensor
        Result of the dot product, with updated metadata.

    """ 
    dotted = dot(
    tensorA.ufl_tens,
    tensorB.ufl_tens,
    ) 
    
    return Tensor(
    dotted,
    tensorA.CoordSys,
    mayHaveSpectralDimension = tensorA.hasSpectralDimension or tensorB.hasSpectralDimension,
    m = tensorA.m + tensorB.m,
    )



def i_inner(tensorA: Tensor,
            tensorB: Tensor,
):
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
        Scalar-valued (order-0) tensor representing the inner product.

    Raises
    ------
    ValueError
        If the tensors are not both of order 2.
    """
    if tensorA.order != tensorB.order or tensorA.order != 2:
        raise ValueError("The order of both tensors must be two.")
    
    
    innered = inner(
    tensorA.ufl_tens,
    tensorB.ufl_tens,
    ) 
    
    return Tensor(
    innered,
    tensorA.CoordSys,
    mayHaveSpectralDimension = tensorA.hasSpectralDimension or tensorB.hasSpectralDimension,
    m = tensorA.m + tensorB.m,
    )


def i_grad(
    T: Tensor
 ):
    """
    Computes the gradient of a tensor.

    Parameters
    ----------
    T : Tensor
        The tensor to differentiate.

    Returns
    -------
    Tensor
        Gradient of the input tensor.

    Raises
    ------
    ValueError
        If the tensor is of order > 2 (not implemented).
    """

    # compute gradient: first check symmetry, then add christoffel parts
    if T.order == 0:
        diffs = []
        for i in range(T.dim):
            if T.isSpectralDimension[i] == SpectralIndicator.NOTSPECTRAL:
                diffs.append(T.ufl_tens.dx(i))
            else:
                diffs.append(1j*T.m*T.ufl_tens/T.r)
        gradient = as_vector(diffs)

    elif T.order == 1:

        # partial derivatives part
        diffs =  []
        for i in range(T.dim):
            row = []
            for j in range(T.dim):
                if T.isSpectralDimension[j] == SpectralIndicator.NOTSPECTRAL:
                    row.append(T.ufl_tens[i].dx(j))
                else:
                    term = 1j*T.m*T.ufl_tens[i]/T.r
                    if T.CoordSys.name == "cylindricalfelics" and i==1:
                        term -= T.ufl_tens[2]/T.r
                    elif T.CoordSys.name == "cylindricalfelics" and i==2:
                        term += T.ufl_tens[1]/T.r
                    row.append(term)
            diffs.append(row)
        gradient = as_matrix(diffs)
                       

    #TODO test order 2 in unittests!!!
    elif T.order == 2:

        # partial derivatives part
        diffs =  []
        for i in range(T.dim):
            row = []
            for j in range(T.dim):
                column = []
                for k in range(T.dim):
                    if T.isSpectralDimension[k] == SpectralIndicator.NOTSPECTRAL:
                        column.append(T.ufl_tens[i,j].dx(k))
                    else:
                        term = 1j*T.m*T.ufl_tens[i,j]/T.r
                        if T.CoordSys.name == "cylindricalfelics" and j==1:
                            term -= T.ufl_tens[i,2]/T.r
                        elif T.CoordSys.name == "cylindricalfelics" and j==2:
                            term += T.ufl_tens[i,1]/T.r
                        column.append(term)
                row.append(column)
            diffs.append(row)
            
        gradient = as_matrix(diffs)
    
    else:
        raise ValueError("Gradients of tensors of order > 2 are not " \
                         "implemented.")
        
    return Tensor(
    gradient,
    T.CoordSys,
    mayHaveSpectralDimension = T.hasSpectralDimension,
    m = T.m,
    )
    

def i_div(
    tensor: Tensor
 ):
    """
    Computes the divergence of a tensor.

    Computes the gradient and sums of the last two indices of a tensor.

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
    if tensor.order == 0:
        raise ValueError("Divergence of scalars not defined.")
        
    elif tensor.order == 1:
        gradient = i_grad(tensor) # Added wavenumber on homogeneous direction
        i = indices(1)
        Div = tr(gradient.ufl_tens)

    elif tensor.order == 2:
        gradient = i_grad(tensor) # Added wavenumber on homogeneous direction
        i,j = indices(2)
        Div = as_tensor(
        gradient.ufl_tens[i,j,j],
        (i),
        )
        
    return Tensor(
    Div,
    tensor.CoordSys,
    mayHaveSpectralDimension = tensor.hasSpectralDimension,
    m = tensor.m,
    )


def i_t(
    tensor: Tensor
):
    """
    Returns the transpose of a second-order tensor.

    Parameters
    ----------
    tensor : Tensor
        A tensor of order 2.

    Returns
    -------
    Tensor
        Transposed tensor.

    Raises
    ------
    ValueError
        If tensor order is not 2.
    """
    if tensor.order != 2:
        raise ValueError("Transpose only unambiously defined for tensors of " \
                         "order 2.")
    else:
        return Tensor(
        tensor.ufl_tens.T,
        tensor.CoordSys,
        mayHaveSpectralDimension = tensor.hasSpectralDimension,
        m = tensor.m,
        )


def i_tr(
    tensor: Tensor
):
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
    if tensor.order != 2:
        raise ValueError("Trace only working for tensors of order 2.")
    else:
        return Tensor(
        tr(tensor.ufl_tens),
        tensor.CoordSys,
        mayHaveSpectralDimension = tensor.hasSpectralDimension,
        m = tensor.m,
        )
    

# not implemented anymore
#def iDev(tensor: Tensor):
#    """
#    Computes the deviatoric (spherical-free) part of a second-order tensor.
#
#    Parameters
#    ----------
#    tensor : Tensor
#        A second-order tensor.
#
#    Returns
#    -------
#    Tensor
#        The deviatoric component of the input tensor.
#    """
#    spherical = Tensor(1/3 * tr(tensor.ufl_tens) * Identity(tensor.dim), \
#                       tensor.CoordSys, basis = tensor.basis)
#    return tensor - spherical



def i_identity(
    tensor: Tensor
):
    """
    Returns the identity tensor corresponding to the tensor's dimension and coordinate system.

    Parameters
    ----------
    tensor : Tensor
        Input tensor, used to infer dimensionality and coordinate system.

    Returns
    -------
    Tensor
        Identity tensor with the same coordinate system and dimension. 
    """
    return Tensor(
    Identity(tensor.dim),
    tensor.CoordSys,
    mayHaveSpectralDimension = False,
    m = 0,
    )

def i_conj(
    tensor: Tensor
):
    """
    Computes the complex conjugate of a tensor.

    The underlying UFL expression is conjugated and the sign of the
    spectral wave number ``m`` is flipped.

    Parameters
    ----------
    tensor : Tensor
        Tensor to conjugate.

    Returns
    -------
    Tensor
        Conjugated tensor, with the same coordinate system and spectral
        flags but with ``m`` replaced by ``-m``.
    """
    return Tensor(
    conj(tensor.ufl_tens),
    tensor.CoordSys,
    mayHaveSpectralDimension = tensor.hasSpectralDimension,
    m = -tensor.m,
    )

# TODO: implement unit test in TESTS folder 
def i_outer(
    tensorA: Tensor,
    tensorB: Tensor
):
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
        Outer product tensor.

    Raises
    ------
    Exception
        If scalar multiplication is attempted using this function.
    """
    if tensorA.order == 1 and tensorB.order == 1:
        i,j = indices(2)
        outered = as_tensor(
        tensorA.ufl_tens[i]*tensorB.ufl_tens[j],
        (i,j),
        )

    # the below should not be necessary
    elif tensorA.order == 2 and tensorB.order == 1:
        i,j,k = indices(3)
        outered = as_tensor(
        tensorA.ufl_tens[i,j]*tensorB.ufl_tens[k],
        (i,j,k),
        )
    elif tensorA.order == 1 and tensorB.order == 2:
        i,j,k = indices(3)
        outered = as_tensor(
        tensorA.ufl_tens[i]*tensorB.ufl_tens[i,k],
        (i,j,k),
        )
    elif tensorA.order == 1 and tensorB.order == 2:
        i,j,k,l = indices(4)
        outered = as_tensor(
        tensorA.ufl_tens[i,j]*tensorB.ufl_tens[k,l],
        (i,j,k,l),
        )
    
    elif tensorA.order == 0 or tensorB.order == 0:
        raise Exception('Use standard multiplication with the asterisk symbol *.')
       

    return Tensor(
    outered,
    tensorA.CoordSys,
    mayHaveSpectralDimension = tensorA.hasSpectralDimension or tensorB.hasSpectralDimension,
    m = tensorA.m + tensorB.m,
    ) 
    
