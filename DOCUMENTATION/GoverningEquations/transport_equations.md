# Transport equations how they are implemented in FELiCS
The goal of this section is to reflect the equations how they are implemented in FELiCS. The goal is to reflect what is implemented currently for the sake of having a discussion basis for corrections and new implementations. The goal is not to derive the equations or physics.

## Momentum equation
$$\int \omega \rho \widehat{\bf{u}} \bf{X}^* dx =\\
+\int j \nabla \cdot \left(\bf{X}^* \otimes \overline{\rho} \widehat{\bf{u}} \right) \cdot \overline{\bf{u}} \, dx - \int j \overline{\rho} \left( \left(\overline{\bf{u}} \otimes \bf{X}^*\right) \cdot \widehat{\bf{u}} \right)\cdot \bf{n} \, ds\\
+\int j \nabla \cdot \left(\bf{X}^* \otimes \overline{\rho} \overline{\bf{u}} \right) \cdot \widehat{\bf{u}} \, dx - \int j \overline{\rho} \left( \left(\widehat{\bf{u}} \otimes \bf{X}^*\right) \cdot \overline{\bf{u}} \right)\cdot \bf{n} \, ds\\
+\int j \nabla \cdot\left(\bf{X}^* \otimes \widehat{\rho} \overline{\bf{u}} \right) \cdot \overline{\bf{u}} \, dx - \int j \widehat{\rho} \left( \left(\overline{\bf{u}} \otimes \bf{X}^*\right) \cdot \overline{\bf{u}} \right)\cdot \bf{n} \, ds\\
+ j \widehat{p} \nabla \cdot \bf{X}^* dx - j \widehat{p} \bf{X}^* \cdot \bf{n} \, ds \\
-j \widehat{\bf{\tau}} : \nabla \bf{X}^*$$


