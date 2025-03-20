# Explaining the tensor formalism

## Introduction
This documentation introduces the tensor formalism used in the script collection `tensor_utils` which is a subcomponent inside FELiCS. For a separate version of `tensor_utils` that relies solely on FEniCS and UFL including a more complete introduction, see [here](https://git.tu-berlin.de/kai.hildebrandt1/tensor_utils).

The implementation of the wrapper is given in UPDATE ME TO THE CORRECT THING IN FELICS.

Symmetric problems are simulated much more efficiently if the symmetry direction is not meshed and solved for. Some symmetries are easier to describe in coordinate systems (CSs) other than the Cartesian one. The most common are the curvilinear CS of cylindrical and spherical coordinates. The tensor formalism can be used to specify PDEs and other expressions in arbitrary curvilinear coordinates automatically without the need to alter the underlying equations. That is a PDE or variational form can be stated with tensor analyitcal operators like $\mathrm{grad()}$ or $\mathrm{div}()$ and remains validity in different coordinates.

In the following, the tensor notation is introduced, starting with the definition of the space, through the definition of tensors and ending with tensor-analytic operators. Instead of starting with a notation limited to Cartesian space, tensors are introduced directly in the more general case of curvilinear CSs.

## Coordinate Systems
The coordinates used in this text and in the developed wrapper follow the standard notation of ISO 80000-2. This notation is also used by [Schade et al., 2018](https://doi.org/10.1515/9783110404265), whose definitions of tensors and operators are used in the following. The FELiCS version of `tensor_utils` also includes others.

All curvilinear spaces are introduced using the Cartesian coordinates $x^i=(x,y,z)$ and basis vectors $\mathbf{e}_i$. A set of curvilinear coordinates $z^i$ is given by transformation equations $\hat{z}^i(x^i)$. The transformation back to Cartesian space is done with the inverse transformation equations $\hat{x}^i(z^i)$. To allow for the inverse mapping, the transformation equations must be bijective everywhere except at singular points, e.g. $r=0$ for polar coordinates. To illustrate this notation the cylindrical coordinates $z^i=(r, \varphi, z)$ are given below:


<script type="text/javascript" async>
  MathJax = {
    tex: {
      tags: "ams" // Enables automatic equation numbering
    }
  };
</script>
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-AMS_HTML">
</script>

$$
\begin{align}
z^1 = r &= \sqrt{x^2+y^2} & z^2 = \varphi &= \mathrm{atan2} \left( \frac{y}{x} \right) & z^3 = z &= x^3 \\
x^1 = x &= r\cos(\varphi) & x^2 = y &= r\sin(\varphi) &  x^3 = z &= z^3 \ .
\end{align}
$$

When one keeps every coordinate fixed but one, coordinate curves (2D space) or coordinate surfaces (3D space) are formed. These are shown below as dashed lines for polar coordinates. 

![Coordinate curves and contravariant basis for polar coordinates. The basis is orthogonal but not orthonormal.](/img_tensor_formalism/coordinate_curves.png)

To create a basis for the CS, the coordinate curves are linearized. This is done with the partial derivatives of the transformation equations, which give tangent vectors to those curves called the contravariant basis $\mathbf{g}^{i}$. Analogous to the Cartesian case, the basis vectors are used to assign a quantity (component) to each direction of the CS.

By starting at the differentials of the position vector, the following formula for the contravariant basis can be derived:
$$
\begin{align}
  \mathbf{g}^{i} &= \frac{\partial \hat{z}^i}{\partial x^k} \mathbf{e}_k
\end{align}
$$
where the Einstein summation convention is used. Analogous, another basis called the covariant basis is computed with the inverse transformation equations:
$$
\begin{align}
  \mathbf{g}_{i} &= \frac{\partial \hat{x}^k}{\partial z^i} \mathbf{e}_k \ .
\end{align}
$$
Both bases are called the natural bases. A main property of them is that they are dual, i.g. $\mathbf{g}_i \cdot \mathbf{g}^j = \delta_i^j$ and $\mathbf{g}^i \cdot \mathbf{g}_j = \delta^i_j$. This is very convenient and will make tensor calculations much easier.

The above picture illustrates that the contravariant basis for polar coordinates is orthogonal, i.g. $\mathbf{g}^i \cdot \mathbf{g}^j = 0$ for $i \neq j$. It is however not orthonormal. This and other properties can determined by calculating the metric of a basis:

$$
\begin{align}
  g_{ij} &= \mathbf{g}_i \cdot \mathbf{g}_j & g^{ij} &= \mathbf{g}^i \cdot \mathbf{g}^j \ .
\end{align}
$$
If the metric is an identity matrix, the basis is orthonormal or, in the case of a diagonal matrix, orthogonal. If the matrix is constant, one may refer to the basis as global, if not local.

## Tensors
A tensor is an object that is invariant w.r.t. basis transformations. For example, the velocity (a vector, i.e. a rank-1 tensor) of a particle remains the same object when an observer "looks" at it from different directions (translation and rotation of the observer's position). It remains the same object even if the observer looks at it through a stretching or shearing lens, which distorts the image of physical space.

To capture this invariance, tensors $\mathbf{A}$ of arbitrary rank $n$ defined in a $k$-dimensional space are split into two parts. First, the $n$ number of bases $\mathbf{g}_i$ and/or $\mathbf{g}^i$, where $i \in \{1,...,k\}$. The bases describe the space on which the tensor is defined. Second, the components $A^{i,j,...}$, where there are $n$ number of indices. By summing over twice occurring running indices the components then carry the amount of some quantity the tensor assigns to its bases. Below are examples up to rank-2 tensors, which is the highest rank implemented in UFL:
$$
\begin{align}
  \alpha &= \alpha \\
  \mathbf{a} &= a^i \mathbf{g}_i \\
  \mathbf{A} &= A^{ij} \mathbf{g}_i \otimes \mathbf{g}_j \ .
\end{align}
$$
In the above the dyadic or tensor product is defined through below application:
$$
\begin{align}
	\bm{a} \otimes \bm{b} \cdot \bm{c} &= \bm{a} (\bm{b} \cdot \bm{c}) &
	\bm{a} \cdot \bm{b} \otimes  \bm{c} &= (\bm{a} \cdot \bm{b}) \bm{c} \ .
\end{align}
$$
When transforming a tensor between different coordinate systems, both the basis vectors and the components of the tensor can change. This is illustrated below for a vector, where $\hat{a}_i$ are the components in a Cartesian basis, $a^{i}$ in a covariant basis of some curvilinear CS and $a_{i}$ in a contravariant basis of the same curvilinear CS. The equality between the different representations of the vector denotes the invariance:
$$
\begin{equation}
  \mathbf{a} = \hat{a}_i \mathbf{e}_i = a^{i} \mathbf{g}_i = a_{i} \mathbf{g}^i \ .
\end{equation}
$$

From the invariance of tensors (10), the metric (5) and using the duality of the bases one can find the transformation of components in the natural bases by the metric:
$$
\begin{align} % \label{cov_con_components}
  a_{i} &= g_{ij} a^j & a^{i} &= g^{ij} a_j \ .
\end{align}
$$

When specifying a PDE in a curvilinear CS, the initial and/or boundary values are often given in a physical basis. A physical basis is a normalized natural basis:
$$
\begin{align} %\label{physical_bases}
  \tilde{\mathbf{g}}_i &= \frac{\mathbf{g}_i}{\sqrt{g_{ii}}} & \tilde{\mathbf{g}}^i &= \frac{\mathbf{g}^i}{\sqrt{g^{ii}}} \ 
\end{align}
$$
where the summation convention is no longer applicable. Normalizing the basis has the advantage that the basis vectors are unitless and only the components represent the physical quantity. This makes formulating the initial or boundary conditions so easy. In fact, when curvilinear CSs are mentioned in physics books, often only physical bases are meant. Furthermore, for the important case of orthogonal natural bases the physical bases coincide and the distinction between upper and lower indices can be dropped:
$$
\begin{equation}
  \tilde{\mathbf{g}}_i = \tilde{\mathbf{g}}^i \rightarrow \mathbf{a}= \tilde{a}_i \tilde{\mathbf{g}}_i = \tilde{a}^i \tilde{\mathbf{g}}^i \ .
\end{equation}
$$
In this case the physical basis forms a locally orthonormal basis, analogous to Cartesian coordinates, only that the direction of the physical basis vectors change from point to point. For this reason, orthogonal physical coordinates are also called locally Cartesian coordinates. However, having a locally orthonormal basis comes at a price: the summation convention can no longer be applied (see 12) which breaks the powerful elegance of tensor calculus!

In practice it is far more convenient to stay in natural bases as long as possible and only convert the final result into physical coordinates. **By default, `tensor_utils` assumes that the initial and boundary conditions are given in a physical basis and converts them into a natural basis, namely the covariant basis.** All tensor operations are then performed in natural bases and only the final result, i.e. the form, is converted into the physical basis at the end. This should be the most user-friendly approach, since input and output can be specified in physical coordinates, which are most commonly used. However, it is also possible to specify initial and boundary conditions in a co- or contravariant basis.


## Tensor Operations
In the following fundamental operations as well as operators of tensor algebra and analysis in curvilinear CSs are briefly introduced. It is assumed that the reader is familiar with their basics in Cartesian coordinates. The mentioned operations are limited and lack some operations, like contractions above order 2. However, it is at least suitable for specifying forms regarding linear elasticity of solid bodies, diffusion and incompressible Newtonian fluids.

### Algebra
Tensor algebra in curvilinear CSs is basically the same as in Cartesian CSs as long as the basis vectors of the operands are the same. If the basis vectors differ, one must first transform the tensor with (11) and the invariance $a^i \mathbf{g}_i=a_i \mathbf{g}^i$. Below is an example for vector addition:
$$
\begin{equation}
  \mathbf{v} + \mathbf{u} = v^i \mathbf{g}_i + u_i \mathbf{g}^i = v^i \mathbf{g}_i + g^{ij} u_j \mathbf{g}_i = (v^i + g^{ij} u_j) \mathbf{g}_i \ .
\end{equation}
$$

The identity tensor is given by:
$$
\begin{equation}
  \mathbf{1}= \mathbf{g}_i \otimes \mathbf{g}^i = \mathbf{g}^i \otimes \mathbf{g}_i
\end{equation}
$$
where other combinations of co- and contravariant bases can be used if it is compensated by the corresponding metric.

### Analysis
In the following, tensor-analytical operators are introduced, most importantly the gradient and consecutively the divergence. Both require the Nabla operator, which is introduced here as a vector differential operator:
$$
\begin{equation}
 \nabla = \frac{\partial}{\partial x^i} \mathbf{e}_i \ .
\end{equation}
$$

When using curvilinear coordinates $z^i$, the Cartesian coordinates $x^i$ are given as a function $\hat{z}$ of these curvilinear coordinates. By applying the chain rule and \autoref{contravariant_basis}, one can arrive at:
$$
\begin{equation}%\label{nabla}
 \nabla = \frac{\partial}{\partial x^i} \mathbf{e}_i = \frac{\partial}{\partial z^k} \frac{\partial \hat{z}^k}{\partial x^i} \mathbf{e}_i = \frac{\partial}{\partial z^i} \mathbf{g}^i \ .
\end{equation}
$$
This leads to the definition of the gradient and divergence. Here, as in the UFL package, the right gradient is used, i.e. $\nabla$ is placed at the right side of the tensor. The differentiation is applied to the term to the left of the Nabla operator. When calculating the gradient of a scalar, only one change has to be accounted for - the change of the scalar itself:
$$
\begin{equation}
  \alpha \nabla = \nabla \alpha = \frac{\partial \alpha}{\partial z^i} \mathbf{g}^i \ .
\end{equation}
$$
Tensors of rank greater than zero are defined with basis vectors. Since the basis is generally a function of the coordinates, the gradient of these tensors must include two changes: one in the components and one in the change of the basis. The latter is represented by the Christoffel symbols $\Gamma^{i}_{jk}$:
$$
\begin{equation}
\begin{split}
  \mathrm{grad}\left(\mathbf{a}\right) &= \mathbf{a} \otimes \nabla = (a^i \mathbf{g}_i) \otimes \frac{\partial}{\partial z^j} \mathbf{g}^j = \frac{\partial a^i}{\partial z^j} \mathbf{g}_i \otimes \mathbf{g}^j + a^i \frac{\partial \mathbf{g}_i}{\partial z^j} \otimes \mathbf{g}^j \\ &= \left( \frac{\partial a^i}{\partial z^j} + a^k \Gamma^{i}_{jk}\right) \mathbf{g}_i \otimes \mathbf{g}^j \ .
\end{split}
\end{equation}
$$

Using the duality of the natural bases and some index renaming, one can find the formula for the Christoffel symbols from their constitutive relation (left side):
$$
\begin{equation}%\label{Christoffel}
  \frac{\partial \mathbf{g}_i}{\partial z^j} = \Gamma^{k}_{ij} \mathbf{g}_k \ \rightarrow \ \Gamma^{i}_{jk} = \frac{\partial \mathbf{g}_i}{\partial z^j} \cdot \mathbf{g}^k \ . 
\end{equation}
$$
Note that Christoffel symbols are symmetric in the lower indices due to the symmetry of second derivatives (Schwarz's theorem). Plugging (17) into the gradient of a rank-2 tensor and using (20) gives:
$$
\begin{align}
  \mathbf{A} \otimes \nabla = \left( \frac{\partial A^{ij}}{\partial z^k} + A^{rj} \Gamma^{i}_{rk} + A^{ir} \Gamma^{j}_{rk} \right) \mathbf{g}_i \otimes \mathbf{g}_j \otimes \mathbf{g}^k \ .
\end{align}
$$

Up to now, only gradients of tensors defined with covariant basis vectors were considered. If one of the bases of a tensor is contravariant, then the following partial derivatives must be calculated:
$$
\begin{equation}
  \frac{\partial \mathbf{g}^i}{\partial z^j}  = - \Gamma^{i}_{jk} \mathbf{g}^k
\end{equation}
$$
where the equality can be shown by differentiating a rank-2 identity tensor.
The divergence of a tensor is given as a contraction with the Nabla operator, i.e. for vectors it is the trace of the gradient. The formula for the divergence of a rank-2 tensor is a little harder to derive. This is done below step-by-step to illustrate the introduced concepts:
$$
\begin{equation}
\begin{split}
	\mathrm{div} (\mathbf{A}) &= \mathbf{A} \cdot \nabla = A^{ij} \mathbf{g}_i \otimes \mathbf{g}_j \cdot \frac{\partial}{\partial z^k} \mathbf{g}^k = \frac{\partial A^{ij}}{\partial z^k} \mathbf{g}_i \delta_j^k + A^{ij} \underbrace{\frac{\partial \mathbf{g}_i}{\partial z^k}}_{\Gamma^r_{ik}\mathbf{g}_r} \delta_j^k + A^{ij} \mathbf{g}_i \underbrace{\frac{\partial \mathbf{g}_j}{\partial z^k}}_{\Gamma^r_{jk}\mathbf{g}_r} \cdot \mathbf{g}^k \\
	&= \frac{\partial A^{ij}}{\partial z^j} \mathbf{g}_i + A^{ij} \Gamma^r_{ij} \mathbf{g}_r + A^{ij} \mathbf{g}_i \Gamma^r_{jr} = \left(\frac{\partial A^{ij}}{\partial z^j} +  A^{rj} \Gamma^i_{rj} + A^{ij} \Gamma^r_{jr} \right) \mathbf{g}_i \ .
\end{split}
\end{equation}
$$

Just as for scalars and vectors the divergence of a rank-2 tensor can be restated as a linear combination of the coefficients of the gradient:
$$
\begin{equation}
	\mathrm{grad}(\mathbf{A}) \cdot \cdot \mathbf{1} = \left( \frac{\partial A^{ij}}{\partial z^k} + A^{rj} \Gamma^{i}_{rk} + A^{ir} \Gamma^{j}_{rk} \right) \mathbf{g}_i \otimes \mathbf{g}_j \otimes \mathbf{g}^k \cdot \cdot \, \mathbf{g}^s \otimes \mathbf{g}_s = \mathrm{div}(\mathbf{A}) \ .
\end{equation}
$$
This is convenient for implementation in the coded wrapper, since only the gradient function has to handle different bases and ranks of tensors. The divergence function can be kept simple and use calculus functions like double contractions. In the wrapper tensor operations are implemented up to rank 2. This was taken as a cut-off rank because of the limitation of UFL to rank-2 tensors.

As mentioned above, changing the CS changes not only the tensors in the domain, but also the domain itself. For example a unit sphere in Cartesian coordinates transforms to a cuboid in spherical coordinates with $\Omega = (0,1) \times (0, \pi) \times (0, 2 \pi)$. This change must be compensated in integral computations. According to \citet[p. 201]{SchadeNeemann.2018} the volume element reads:
$$
\begin{equation}
	\mathrm{d}V = \sqrt{\mathrm{det}(g_{ij})} \mathrm{d}z^1 \mathrm{d}z^2 \mathrm{d}z^3 \ .
\end{equation}
$$