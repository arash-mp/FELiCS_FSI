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
# Third party libraries
import copy
from dolfinx.fem import (
    Constant,
    Function
)
from petsc4py import PETSc

# Local Libraries and methods
from FELiCS.Fields.Field        import Field
from FELiCS.Misc.logging        import Logger
from FELiCS.Misc.tensorUtils    import (
    i_div,
    i_dot,
    i_grad,
    i_identity,
    i_t,
    Tensor,
)

# Get the logger
logger = Logger.get_logger("felics")

class FieldProperties:
    """
    Provides property decorators for accessing and transforming field quantities.

    This class defines a set of property methods and accessors for field variables,
    enabling consistent access to physical quantities (e.g., velocity, pressure, temperature)
    from the internal field dictionary. It is designed to be used as a mixin for
    classes such as `fluctuationClass`, `fluctuationSolution`, `meanFlowClass`, and
    `meanFlowVertexValues`.

    **Initialize the fieldProperties object**

    Attributes
    ----------
    _fieldDict : dict
        Dictionary containing all field variables for the current object.
    _zeroField : object
        Default fallback field for missing scalar quantities.
    _oneField : object
        Default fallback field with all coefficients set to one.
    _zeroVectorField : object
        Default fallback field for missing velocity quantities.
    _mean : object
        Reference to the mean flow object, if available.
    _FEMSpaces : object
        Finite element spaces, used for projections.
    _coordinateSystem : object
        Coordinate system information for tensor operations.
    __hSpec : dict
        Additional species enthalpy information, if present.
    __reaction : object
        Reaction information, if present.
    _param : object
        Simulation parameter object, if present.
    _mesh : object
        Mesh object, if present.
    _fluc : object
        Fluctuation field object, if present.
    _ScalarFunctionSpace : object
        Scalar function space for field projections.
    _meanflowFilename : str
        Filename for the mean flow data, if present.
    _transportedQuantities : list of str
        List of transported quantity names, if present.

    alpha : object
        Field variable for alpha (e.g., thermal diffusivity or similar).
    cp : object
        Field variable for specific heat at constant pressure.
    dQ : object
        Field variable for heat release rate.
    fieldDict : dict
        Dictionary of all field variables, with special handling for mean flow classes.
    FieldNames : list of str
        List of field variable names.
    fluc : object
        Fluctuation field object.
    forcingDomain : object
        Field variable for the domain where forcing is applied.
    gamma : object
        Field variable for heat capacity ratio.
    Pr : object
        Field variable for Prandtl number.
    UnitT : object
        Real unit number in tensor form.
    h : object
        Field variable for enthalpy.
    he : object
        Field variable for total enthalpy.
    hSpec : dict
        Dictionary of species enthalpy fields.
    meanflowFilename : str
        Filename for the mean flow data.
    molarMass : object
        Field variable for molar mass.
    nulam : object
        Field variable for laminar viscosity.
    nuTot : object
        Field variable for total viscosity (laminar + turbulent + SGS).
    p : object
        Field variable for pressure.
    RR_prefactor : object
        Field variable for reaction rate prefactor.
    Q : object
        Field variable for total heat release.
    R_spe : object
        Field variable for specific gas constant.
    reaction : object
        Reaction object.
    spg : object
        Field variable for sponge region term.
    rho : object
        Field variable for density.
    rhou : object
        Field variable for momentum (density * velocity).
    T : object
        Field variable for temperature.
    tau : object
        Field variable for stress tensor.
    Tb : object
        Field variable for burnt temperature.
    Tu : object
        Field variable for unburnt temperature.
    responseDomain : object
        Field variable for the response domain mask.
    u : object
        Field variable for velocity.
    u_forcing_i : object
        Field variable for imaginary part of velocity forcing.
    u_forcing_r : object
        Field variable for real part of velocity forcing.
    u_forcing : object
        Field variable for complex velocity forcing.
    ut : object
        Field variable for transverse velocity component.

    Notes
    -----
    This class is not intended to be instantiated directly, but to be inherited by
    classes that manage field data for mean or fluctuating flow quantities.
    """

    def is_mean_flow_class(
        self
    ):
        """
        Check if the current object is an instance of meanFlowClass.

        Returns
        -------
        bool
            True if the object is a meanFlowClass instance, False otherwise.
        """
        from FELiCS.Fields.MeanFlowClass import MeanFlowClass

        return isinstance(
            self,
            MeanFlowClass,
        )

    def is_mean_flow_vertex_values_class(
        self
    ):
        """
        Check if the current object is an instance of meanFlowVertexValues.

        Returns
        -------
        bool
            True if the object is a meanFlowVertexValues instance, False otherwise.
        """
        from FELiCS.Fields.MeanFlowClass import MeanFlowVertexValues
        return isinstance(
            self,
            MeanFlowVertexValues,
        )

    @property
    def alpha(
        self
    ):
        """
        Get the alpha field variable, with fallback for missing data.

        Returns
        -------
        object
            Alpha field as tensor or fallback field.
        """
        if self.is_mean_flow_class():
            if 'alpha' in list(self._fieldDict.keys()):
                return  self._fieldDict['alpha'].get_tensor()
            else:
                return self._zeroField.get_tensor()
        if self.is_mean_flow_vertex_values_class():
            if 'alpha' in list(self._fieldDict.keys()):
                return self._fieldDict['alpha']
            else:
                return self._zeroField
        else:
            return self._fieldDict['alpha']

    @property
    def c(
        self
    ):
        # TODO: Sophie, could you delete this property? I think nobody uses it.
        return self._fieldDict['c']

    @property
    def cp(
        self
    ):
        """
        Get the specific heat at constant pressure.

        Returns
        -------
        object
            Field variable for cp.
        """
        if self.is_mean_flow_class() or self.is_mean_flow_vertex_values_class():
            return self._fieldDict['cp'].get_tensor()
        else:       
            return self._fieldDict['cp']

    def d(
        self,
        specie
    ):
        """
        Get the diffusion coefficient for a given species.

        Parameters
        ----------
        specie : str
            Name of the species.

        Returns
        -------
        object
            Diffusion coefficient field for the species.
        """
        if self.is_mean_flow_class() or self.is_mean_flow_vertex_values_class():
            if 'D_' + specie in list(self._fieldDict.keys()):
                return self._fieldDict['D_' + specie].get_tensor()
            else:
                return self._zeroField.get_tensor()
        else:
            if 'D_' + specie in list(self._fieldDict.keys()):
                return self._fieldDict['D_' + specie]
            else:
                return self._zeroField


    @property
    def d_q(
        self
    ):
        """
        Get the heat release rate field variable.

        Returns
        -------
        object
            Field variable for heat release rate.
        """
        return self._fieldDict['dQ']

    @property
    def field_dict(
        self
    ):
        """
        Get the dictionary of all field variables, with special handling for mean flow classes.

        Returns
        -------
        dict
            Dictionary of field variables.
        """
        if self.is_mean_flow_class() or self.is_mean_flow_vertex_values_class():
            Output = copy.copy(self._fieldDict)
            if '__hSpec' in dir(self):
                for key in list(self.__hSpec.keys()):
                    Output['hSpec_' + key] = self.__hSpec[key]
            return Output
        else:
            OutputDict = {}
            for key in list(self._fieldDict.keys()):
                if not isinstance(
                    self._fieldDict[key],
                    Function,
                ):
                    if key == 'u':
                        FEMSpace = self._FEMSpaces.FunctionSpaceVectorVelocity
                    else:
                        FEMSpace = self._FEMSpaces.P2
                        from fenics import project
                    OutputDict[key] = project(
                        self._fieldDict[key],
                        FEMSpace,
                    ) # NOTE: (Simon) not sure what this whould be
                else:
                    OutputDict[key] = self._fieldDict[key]
            return OutputDict

    @field_dict.setter
    def field_dict(
        self,
        value
    ):
        """
        Prevent setting the fieldDict property after initialization.

        Raises
        ------
        Exception
            Always raised to prevent modification.
        """
        raise Exception('Properties of MeanFlow are not to be set after initialization!')

    @property
    def field_names(
        self
    ):
        """
        Get the list of field variable names.

        Returns
        -------
        list of str
            List of field variable names.
        """
        return self.__FieldNames

    @field_names.setter
    def field_names(
        self,
        value
    ):
        """
        Prevent setting the FieldNames property after initialization.

        Raises
        ------
        Exception
            Always raised to prevent modification.
        """
        raise Exception('Properties of MeanFlow are not to be set after initialization!')

    @property
    def fluc(
        self
    ):
        """
        Get the fluctuation field object.

        Returns
        -------
        object
            Fluctuation field object.
        """
        return self._fluc

    def forcing_i(
        self,
        solution
    ):
        """
        Get the imaginary part of the forcing field for a given solution variable.

        Parameters
        ----------
        solution : str
            Name of the solution variable.

        Returns
        -------
        object
            Imaginary part of the forcing field.
        """
        name = solution + '_forcing_i'
        if self.is_mean_flow_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name].get_tensor()
            else:
                return self._zeroField.get_tensor()
            
        elif self.is_mean_flow_vertex_values_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField
        
        else:
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField
        
    def forcing_r(
        self,
        solution
    ):
        """
        Get the real part of the forcing field for a given solution variable.

        Parameters
        ----------
        solution : str
            Name of the solution variable.

        Returns
        -------
        object
            Real part of the forcing field.
        """
        name = solution + '_forcing_r'
        if self.is_mean_flow_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name].get_tensor()
            else:
                return self._zeroField.get_tensor()
            
        elif self.is_mean_flow_vertex_values_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField
        
        else:
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField

    def forcing(
        self,
        solution
    ):
        """
        Get the complex forcing field for a given solution variable.

        Parameters
        ----------
        solution : str
            Name of the solution variable.

        Returns
        -------
        object
            Complex forcing field as a tensor.
        """
        return self.forcing_r(solution)+self.forcing_i(solution)*1j

    @property
    def forcing_domain(
        self
    ):
        """
        Get the field variable for the domain where forcing is applied.

        Returns
        -------
        object
            Field variable for the forcing domain.
        """
        return self._fieldDict['forcingDomain']

    @property
    def gamma(
        self
    ):
        """
        Get the heat capacity ratio (gamma).

        Returns
        -------
        object
            Field variable for heat capacity ratio.
        """
        return self._fieldDict['gamma'].get_tensor()

    @property
    def pr(
        self
    ):
        """
        Get the Prandtl number field variable.

        Returns
        -------
        object
            Field variable for Prandtl number.
        """
        return self._fieldDict['Pr'].get_tensor()

    @property
    def unit_t(
        self
    ):
        """
        Get the real unit number in tensor form.

        Returns
        -------
        object
            Tensor representing the real unit number.
        """
        mesh = self._fieldDict[list(self._fieldDict.keys())[0]].space.mesh
        return Tensor(
            Constant(
                mesh,
                PETSc.ScalarType(1.0 + 0j),
            ),
            self._coordinateSystem,
        )

    @property
    def h(
        self
    ):
        """
        Get the enthalpy field variable.

        Returns
        -------
        object
            Field variable for enthalpy.
        """
        return self._fieldDict['h']

    @property
    def he(
        self
    ):
        """
        Get the total enthalpy field variable.

        Returns
        -------
        object
            Field variable for total enthalpy.
        """
        return (
            self._fieldDict['he'].get_tensor()
            + 0.5 * i_dot(
                self.u,
                self.u,
            )   
        )

    @property
    def h_spec(
        self
    ):
        """
        Get the dictionary of species enthalpy fields.

        Returns
        -------
        dict
            Dictionary of species enthalpy fields.
        """
        return self.__hSpec

    @property
    def meanflow_filename(
        self
    ):
        """
        Get the filename for the mean flow data.

        Returns
        -------
        str
            Filename for the mean flow data.
        """
        return self._meanflowFilename

    @property
    def molar_mass(
        self
    ):
        """
        Get the molar mass field variable.

        Returns
        -------
        object
            Field variable for molar mass.
        """
        return self._fieldDict['molarMass']

    @property
    def nulam(
        self
    ):
        """
        Get the laminar viscosity field variable.

        Returns
        -------
        object
            Field variable for laminar viscosity.
        """
        if self.is_mean_flow_class():
            if 'nulam' in list(self._fieldDict.keys()):
                return self._fieldDict['nulam'].get_tensor()
            else:
                return self._zeroField.get_tensor()
        elif self.is_mean_flow_vertex_values_class():
            if 'nulam' in list(self._fieldDict.keys()):
                return self._fieldDict['nulam']
            else:
                return self._zeroField
        
        else:
            if 'nulam' in list(self._fieldDict.keys()):
                return self._fieldDict['nulam']
            else:
                return self._zeroField

    @property
    def nu_tot(
        self
    ):
        """
        Get the total viscosity field variable (laminar + turbulent + SGS).

        Returns
        -------
        object
            Field variable for total viscosity.
        """
        nu_tot = Field(
            self._ScalarFunctionSpace,
            self._mesh,
        )

        if 'nulam' in list(self._fieldDict.keys()):
            nu_tot += self._fieldDict['nulam']
        if 'nuturb' in list(self._fieldDict.keys()):
            nu_tot += self._fieldDict['nuturb']
        if 'nuSGS' in list(self._fieldDict.keys()):
            nu_tot += self._fieldDict['nuSGS']
        return nu_tot.get_tensor()

    @property
    def p(
        self
    ):
        """
        Get the pressure field variable, with fallback for missing data.

        Returns
        -------
        object
            Field variable for pressure.
        """
        if self.is_mean_flow_class():
            if 'p' in list(self._fieldDict.keys()):
                return self._fieldDict['p'].get_tensor()
            else:
                return self._zeroField.get_tensor()
            
        elif self.is_mean_flow_vertex_values_class():
            if 'p' in list(self._fieldDict.keys()):
                return self._fieldDict['p']
            else:
                return self._oneField

        else:
            if 'p' in list(self._fieldDict.keys()):
                return self._fieldDict['p']
            else:
                return self._zeroField
            
    @property
    def phi(
        self
    ):
        # TODO: Sophie, could you check if we can delete this property? I think nobody uses it.
        if 'phi' in list(self._fieldDict.keys()):
            return self._fieldDict['phi']
        else:
            return self._zeroField

    @property
    def rr_prefactor(
        self
    ):
        """
        Get the reaction rate prefactor field variable.

        Returns
        -------
        object
            Field variable for reaction rate prefactor.
        """
        if self.is_mean_flow_class():
            if 'RR_prefactor' in list(self._fieldDict.keys()):
                return self._fieldDict['RR_prefactor'].get_tensor()
        else:
            return self._fieldDict['RR_prefactor']

    @property
    def q(
        self
    ):
        """
        Get the total heat release field variable.

        Returns
        -------
        object
            Field variable for total heat release.
        """
        return self._fieldDict['Q']
    
    @property
    def r_spe(
        self
    ):
        """
        Get the specific gas constant field variable.

        Returns
        -------
        object
            Field variable for specific gas constant.
        """
        return self._fieldDict['R_spe'].get_tensor()

    @property
    def reaction(
        self
    ):
        """
        Get the reaction object.

        Returns
        -------
        object
            Reaction object.
        """
        return self.__reaction
    
    @property
    def spg(
        self
    ):
        """
        Get the sponge region term field variable.

        Returns
        -------
        object
            Field variable for sponge region term.
        """
        return self._fieldDict['spg'].get_tensor()

    @property
    def rho(
        self
    ):
        """
        Get the density field variable, with fallback for missing data.

        Returns
        -------
        object
            Field variable for density.
        """
        if self.is_mean_flow_class():
            # Comment from Sophie: I added rho to the quantities to read in as default, s.t. a variable density
            # without rho as fluctuation variable is possible ("cold flow"). If rho is not given as a mean field,
            # it will be automatically initialized as a function with all coefficients equal to zero. In that
            # case, a field with all coefficients equal to one is returned.
            # TODO: redo when restructuring the initialization process.
            if 'rho' in list(self._fieldDict.keys()) and sum(self._fieldDict['rho'].get_coefficient_array() ) != 0.:
                return self._fieldDict['rho'].get_tensor()
            else:
                return self._oneField.get_tensor()
            
        elif self.is_mean_flow_vertex_values_class():
            if 'rho' in list(self._fieldDict.keys()):
                return self._fieldDict['rho']
            else:
                return self._oneField

        else:
            if 'rho' in list(self._fieldDict.keys()):
                return self._fieldDict['rho']
            else:
                return self._zeroField

    @property
    def rhou(
        self
    ):
        """
        Get the momentum field variable (density * velocity).

        Returns
        -------
        object
            Field variable for momentum.
        """
        if 'rhou' in list(self._fieldDict.keys()):
            return self._fieldDict['rhou']
        else:
            return self._mean.rho * self.u + self.rho * self._mean.u

    def rho_y(
        self,
        species
    ):
        """
        Get the product of density and species mass fraction.

        Parameters
        ----------
        species : str
            Name of the species.

        Returns
        -------
        object
            Field variable for rho * Y(species).
        """
        return self._mean.rho * self.y(species) + self.rho * self._mean.y(species)

    @property
    def T(
        self
    ):
        """
        Get the temperature field variable, with fallback for missing data.

        Returns
        -------
        object
            Field variable for temperature.
        """
        if self.is_mean_flow_class():
            if 'T' in list(self._fieldDict.keys()):
                return self._fieldDict['T'].get_tensor()
            else:
                return self._oneField
            
        elif self.is_mean_flow_vertex_values_class():
            if 'T' in list(self._fieldDict.keys()):
                return self._fieldDict['T']
            else:
                return self._oneField

        else:
            if 'T' in list(self._fieldDict.keys()):
                return self._fieldDict['T']
            else:
                return self._zeroField

    @property
    def tau(
        self
    ):
        """
        Get the stress tensor field variable.

        Returns
        -------
        object
            Field variable for stress tensor.
        """
        if self.is_mean_flow_class() or self.is_mean_flow_vertex_values_class():
            mean_nu = self.nu_tot
            mean_u = self.u
            tau_out = mean_nu * i_grad(mean_u)
            tau_out += i_t(tau_out)
            # Sophie: this if-clause if not really necessary, in the incompressible case the term is just zero
            if not self._param.Case.SetOfEquations['Energy']['Equation'] == 'None':
                tau_out += ( - 2.0 / 3.0 * mean_nu * 
                            i_div(mean_u) * i_identity(i_grad(mean_u)))
                
        else:
            mean_nu = self._mean.nu_tot
            mean_u = self._mean.u
            fluc_nu = self.nulam
            tau_out = mean_nu * i_grad(self.u) + \
                        fluc_nu * i_grad(mean_u)
            tau_out += i_t(tau_out)
            # Sophie: this if-clause if not really necessary, in the incompressible case the term is just zero
            if not self._param.Case.SetOfEquations['Energy']['Equation'] == 'None':
                tau_out += -2.0/3.0 * mean_nu * i_div(self.u) * i_identity(i_grad(self.u))
                tau_out += -2.0/3.0 * fluc_nu * i_div(mean_u) * i_identity(i_grad(self.u))
        return tau_out

    @property
    def tb(
        self
    ):
        """
        Get the burnt temperature field variable.

        Returns
        -------
        object
            Field variable for burnt temperature.
        """
        if self.is_mean_flow_class():
            if 'Tb' in list(self._fieldDict.keys()):
                return self._fieldDict['Tb'].get_tensor()
        elif self.is_mean_flow_vertex_values_class():
            if 'Tb' in list(self._fieldDict.keys()):
                return self._fieldDict['Tb']
                
        else:
            return self._fieldDict['Tb']

    @property
    def tu(
        self
    ):
        """
        Get the unburnt temperature field variable.

        Returns
        -------
        object
            Field variable for unburnt temperature.
        """
        if self.is_mean_flow_class():
            if 'Tu' in list(self._fieldDict.keys()):
                return self._fieldDict['Tu'].get_tensor()
        elif self.is_mean_flow_vertex_values_class():
            if 'Tu' in list(self._fieldDict.keys()):
                return self._fieldDict['Tu']
                
        else:
            return self._fieldDict['Tu']

    @property
    def tm(
        self
    ):
        # TODO: Sophie, could you check if we can delete this property? I think nobody uses it.
        return self._fieldDict['Tm']


    @property
    def response_domain(
        self
    ):
        """
        Get the response domain mask field variable.

        Returns
        -------
        object
            Field variable for response domain mask.
        """
        return self._fieldDict['responseDomain']

    @property
    def u(
        self
    ):
        """
        Get the velocity field variable, with fallback for missing data.

        Returns
        -------
        object
            Field variable for velocity.
        """
        if self.is_mean_flow_class() or self.is_mean_flow_vertex_values_class():
            if 'u' in list(self._fieldDict.keys()):
                return self._fieldDict['u'].get_tensor()
            else:
                return self._zeroVectorField.get_tensor()
        else:
            if 'u' in self._transportedQuantities:
                return  self._fieldDict['u']
            else:
                return self._zeroVectorField

    @property
    def u_forcing_i(
        self
    ):
        """
        Get the imaginary part of the velocity forcing field variable.

        Returns
        -------
        object
            Field variable for imaginary part of velocity forcing.
        """
        name = 'u_forcing_i'
        if self.is_mean_flow_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name].get_tensor()
            else:
                return self._zeroVectorField.get_tensor()
            
        elif self.is_mean_flow_vertex_values_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField
        
        else:
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField

    @property
    def u_forcing_r(
        self
    ):
        """
        Get the real part of the velocity forcing field variable.

        Returns
        -------
        object
            Field variable for real part of velocity forcing.
        """
        name = 'u_forcing_r'
        if self.is_mean_flow_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name].get_tensor()
            else:
                return self._zeroVectorField.get_tensor()
            
        elif self.is_mean_flow_vertex_values_class():
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField
        
        else:
            if name in list(self._fieldDict.keys()):
                return self._fieldDict[name]
            else:
                return self._zeroField

    @property
    def u_forcing(
        self
    ):
        """
        Get the complex velocity forcing field variable.

        Returns
        -------
        object
            Field variable for complex velocity forcing.
        """
        forcing = self.u_forcing_r + self.u_forcing_i * 1j
        return forcing

    @property
    def ut(
        self
    ):
        """
        TODO: Deprecated variable must be deleted
        Get the transverse velocity component field variable.

        Returns
        -------
        object
            Field variable for transverse velocity component.
        """
        if self.is_mean_flow_class() or self.is_mean_flow_vertex_values_class():
            if 'ut' in list(self._fieldDict.keys()):
                return self._fieldDict['u'][2]
            else:
                mesh = self._fieldDict[
                    list(self._fieldDict.keys())[0]].space.mesh
                return Constant(
                    mesh,
                    0.0,
                )
        else:
            if  self._param.Case.TransVelFluc:
                return self._fieldDict['u'][2]
            else:
                return Constant(0)

    def y(self,
    specie
    ):
        """
        Get the mass fraction field variable for a given species.

        Parameters
        ----------
        specie : str
            Name of the species.

        Returns
        -------
        object
            Field variable for species mass fraction.
        """
        if self.is_mean_flow_class():
            if specie in list(self._fieldDict.keys()):
                return self._fieldDict[specie].get_tensor()
            else:
                return self._zeroField.get_tensor()
            
        elif self.is_mean_flow_vertex_values_class():
            if specie in list(self._fieldDict.keys()):
                return self._fieldDict[specie]
            else:
                return self._zeroField
        
        else:
            if specie in list(self._fieldDict.keys()):
                return self._fieldDict[specie]
            else:
                return self._zeroField
                        
