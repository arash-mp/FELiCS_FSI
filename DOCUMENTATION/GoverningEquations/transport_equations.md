# Transport equations - WORK IN PROGRESS (Thomas) 

(how they are implemented in FELiCS)

## Momentum equation

### Variable density

$$
int \omega \overline{\rho} \widehat{\bf{u}} \bf{X}^* dx =
+\int j \nabla \cdot \left(\bf{X}^* \otimes \overline{\rho} \widehat{\bf{u}} \right) \cdot \overline{\bf{u}} \, dx - \int j \overline{\rho} \left( \left(\overline{\bf{u}} \otimes \bf{X}^*\right) \cdot \widehat{\bf{u}} \right)\cdot \bf{n} \, ds
+\int j \nabla \cdot \left(\bf{X}^* \otimes \overline{\rho} \overline{\bf{u}} \right) \cdot \widehat{\bf{u}} \, dx - \int j \overline{\rho} \left( \left(\widehat{\bf{u}} \otimes \bf{X}^*\right) \cdot \overline{\bf{u}} \right)\cdot \bf{n} \, ds
+\int j \nabla \cdot\left(\bf{X}^* \otimes \widehat{\rho} \overline{\bf{u}} \right) \cdot \overline{\bf{u}} \, dx - \int j \widehat{\rho} \left( \left(\overline{\bf{u}} \otimes \bf{X}^*\right) \cdot \overline{\bf{u}} \right)\cdot \bf{n} \, ds
+ j \widehat{p} \nabla \cdot \bf{X}^* dx - j \widehat{p} \bf{X}^* \cdot \bf{n} \, ds 
-j \widehat{\bf{\tau}} : \nabla \bf{X}^*
$$

### Constant density

$$
\int \omega \widehat{\bf{u}} \bf{X}^* dx =
+\int j \nabla \cdot \left(\bf{X}^* \otimes \widehat{\bf{u}} \right) \cdot \overline{\bf{u}} \, dx - \int j \left( \left(\overline{\bf{u}} \otimes \bf{X}^*\right) \cdot \widehat{\bf{u}} \right)\cdot \bf{n} \, ds
+\int j \nabla \cdot \left(\bf{X}^* \otimes \overline{\bf{u}} \right) \cdot \widehat{\bf{u}} \, dx - \int j \left( \left(\widehat{\bf{u}} \otimes \bf{X}^*\right) \cdot \overline{\bf{u}} \right)\cdot \bf{n} \, ds
+ j \widehat{p} \nabla \cdot \bf{X}^* dx - j \widehat{p} \bf{X}^* \cdot \bf{n} \, ds 
-j \widehat{\bf{\tau}} : \nabla \bf{X}^*
$$

## Mass equation

### Variable density

$$
\int \omega \widehat{\rho} \bf{X}^* dx =
\int j \nabla \bf{X}^* \cdot \widehat{\rho \bf{u}} \, dx
- \int j \widehat{\rho \bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds,
$$
where

$$
\widehat{\rho u} = \overline{\rho} \widehat{\bf{u}} + \widehat{\rho} \overline{\bf{u}}.
$$

### Constant density 

$$ 0 =
\int j \nabla \bf{X}^* \cdot \widehat{ \bf{u}} \, dx
- \int j \widehat{ \bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds,
$$

## Species mass ratio
Transport equation for a chmeical species "s"
$$
\int \omega \overline{\rho} \widehat{Y}_s \bf{X}^* dx =
\int j \, \widehat{Y}_s \, \nabla \cdot \left( \overline{\rho} \, \overline{\bf{u}} \, \bf{X}^* \right) \, dx
+ \int j \, \overline{Y}_s \, \nabla \cdot \left( \overline{\rho} \, \widehat{\bf{u}} \, \bf{X}^* \right) \, dx
+ \int j \, \overline{Y}_s \, \nabla \cdot \left( \widehat{\rho} \, \overline{\bf{u}} \, \bf{X}^* \right) \, dx 
- \int j \, \widehat{Y}_s \, \overline{\rho} \, \overline{\bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds
- \int j \, \overline{Y}_s \, \overline{\rho} \, \widehat{\bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds
- \int j \, \overline{Y}_s \, \widehat{\rho} \, \overline{\bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds 
- \int j \, \overline{D}_s \, \nabla \widehat{Y}_s \cdot \nabla \bf{X}^* \, dx
- \int j \, \widehat{D}_s \, \nabla \overline{Y}_s \cdot \nabla \bf{X}^* \, dx 
+ \int \overline{f}_s \, \bf{X}^* \, dx
$$

## Energy equations 

### Energy Pressure Equation

$$
\int \widehat{p} \, \bf{X}^* \, dx = 
\int j \, \nabla \cdot \left( \bf{X}^* \, \overline{\bf{u}} \right) \, \widehat{p} \, dx
+ \int j \, \nabla \cdot \left( \bf{X}^* \, \widehat{\bf{u}} \right) \, \overline{p} \, dx 
+ \int j \, \overline{\gamma} \, \nabla \left( \overline{p} \, \bf{X}^* \right) \cdot \widehat{\bf{u}} \, dx
+ \int j \, \overline{\gamma} \, \nabla \left( \widehat{p} \, \bf{X}^* \right) \cdot \overline{\bf{u}} \, dx 
- \int j (\overline{\gamma} + 1) \, \widehat{p} \, \overline{\bf{u}}  \bf{X}^* \cdot \bf{n} \, ds
- \int j (\overline{\gamma} + 1) \, \overline{p} \, \widehat{\bf{u}}  \bf{X}^* \cdot \bf{n} \, ds 
- \int j (\overline{\gamma} - 1) \, \left(\nabla \cdot \overline{\boldsymbol{\tau}}\right) \cdot \left( \widehat{\bf{u}}  \bf{X}^* \right) \, dx
- \int j (\overline{\gamma} - 1) \, \left( \nabla \cdot \widehat{\boldsymbol{\tau}} \right) \cdot \left( \overline{\bf{u}}  \bf{X}^* \right) \, dx 
- \int j (\overline{\gamma} - 1) \, \nabla \bf{X}^* \cdot \left( \widehat{\kappa} \nabla  \overline{T}  \right) \, dx
- \int j (\overline{\gamma} - 1) \, \nabla \bf{X}^* \cdot \left( \overline{\kappa} \nabla  \widehat{T}  \right) \, dx
$$
