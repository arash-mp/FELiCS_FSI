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
# Standard libraries
from    abc     import ABC
from    inspect import isclass
import  json
import  os
import  pdb

# Third party libraries
from    h5py    import File

# Local Libraries and methods
from    FELiCS.Equation.MixtureClass    import MixtureClass

from    FELiCS.Misc.functions           import get_last_git_commit
from 	FELiCS.Misc.logging			    import Logger, log_and_raise
from    FELiCS.SpaceDisc.FELiCSMesh     import FELiCSMesh


# Get the logger
logger = Logger.get_logger("felics")

class Dotdict(dict):
    """
    Adds possibility to use dot notation to access entries of dictionary dict
    also works for nested dictionaries

    Parameters
    ----------
    dict : dictionary
        the dictionary which is going to become a dot-dictionary
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

class Config(ABC):
    """
    Abstract base class for FELiCS configuration management.

    Handles loading, parsing, and exporting simulation parameters from JSON and HDF5 files,
    manages default settings, and provides utility methods for parameter access and validation.

    **Initialize the config object**

    Parameters
    ----------
    None
    """

    def __init__(
        self
    ):
        """
        Initialize the config object with default settings.
        """
        logger.debug("Initializing config class with defaults.")
        
        # Get an instance of the defaults settings
        self.default_config = self.get_all_settings_dict()
        
        # NOTE: deprecated?
        self.__BCIDs__      = []   # move to BC

    def get_all_settings_dict(
        self
    ):
        """
        Return the settings dictionary for the config object.

        Returns
        -------
        dict
            Dictionary of input parameters structured in subcategories.
        """
        SettingsDict={
            'BoundaryCondition':{
                'BCsFilePath':              {'datatype':str,    'default':''}
            },
            'Case':{
                'AnalysisMode':             {'datatype':str,    'default':'Modal'},
                'CalculateAdjoint':         {'datatype':bool,   'default':True},
                'CoordinateSystem':         {'datatype':str,    'default':'Cartesian'},
                'needInterpolation':        {'datatype':bool,   'default':True},
                'm':                        {'datatype':int,    'default':0},
                'MeshFilePath':             {'datatype':str,    'default':''},
                'MixtureFilePath':          {'datatype':str,    'default':''},
                'MolVisc':                  {'datatype':int,    'default':0.0},
                'MolViscModel':             {'datatype':str,    'default':'Constant'},
                'MolViscPerturbModel':      {'datatype':dict,
                    'default':{
                        'type':             'Constant',
                        'Constants':        {'Viscosity':1.0}
                    }
                },
                'nDim':                     {'datatype':int,    'default':2},
                'PrandtlNumber':            {'datatype':int,    'default':0.72},
                'Reaction':                 {'datatype':bool,   'default':False},
                'SetOfEquations':           {'datatype':dict,
                    'default':{
                        'Momentum':         {'Equation':'NSPrimitive',  'Variable':'u'},
                        'Mass':             {'Equation':'Continuity',   'Variable':'p'},
                        'Energy':           {'Equation':'None',         'Variable':'None'},
                        'Species':          {'Equation':'None',         'Variable':'None'},
                        'EquationOfState':  {'Equation':'None',         'Variable':'None'}
                    }
                },
                'SpeciesFilePath':          {'datatype':str,    'default':''},
                'TurbulenceModel':          {'datatype':str,    'default':'None'}
            },
            'Export':{
                'ExportFolder':             {'datatype':str,    'default':''},
            },
            'FlowInput':{
                'MeanFlowFilePath':         {'datatype':str,    'default':''},
            },
            'IOResolvent':{
                'ForcingBoundaryIndices':   {'datatype':list,   'default':[]},
                'ForcingMode':              {'datatype':str,    'default':'Body'},
                'ForcingNorm':              {'datatype':str,    'default':'TKE'},
                'ResponseNorm':             {'datatype':str,    'default':'TKE'},
                'Omegas':                   {'datatype':list,   'default':[]}
            },
            'Numerics':{
                'EigenValueGuess':          {'datatype':list,   'default':[1.0]},
                'nSolut':                   {'datatype':int,    'default':3},
                'PolynomialOrder':          {'datatype':dict,   'default':{'u':'2'},    'options':[1,2]}
            }
        }
        return SettingsDict

    def parse_complex_list(
        self,
        data,
    ):
        """
        Convert a list of mixed strings and floats into complex numbers.

        Parameters
        ----------
        data : list
            List containing strings, ints, or floats representing complex numbers.

        Returns
        -------
        list of complex
            List of complex numbers.
        """
        result = []
        for item in data:
            if isinstance(
                item,
                str,
            ):
                try:
                    result.append(complex(item))  # Convert string to complex
                except ValueError:

                    log_and_raise(logger, f"Invalid complex number string: {item}", ValueError)
            elif isinstance(item, (int, float)):

                result.append(item)  # Convert float/int to complex
            else:
                log_and_raise(logger, f"Unsupported type {type(item)} in list. Must be str or float.", TypeError)
        return result
    
    def calculate_parameters(
        self
    ):
        """
        Calculate and set derived parameters based on current configuration.
        """
        self.BoundaryCondition.nVelocityComponents  = len(self.get_velocity_components())
        self.Case.SolutionList                      = self.get_transported_quantity_list()
        self.BoundaryCondition.VelocityComponents   = self.get_velocity_components()
        self.__mesh__.set_true_dimension(self.BoundaryCondition.nVelocityComponents)

    def check_for_mandatory_files(
        self,
        config_dict,
        mandatory_files
    ):
        """
        Check that all mandatory files exist in the configuration.

        Parameters
        ----------
        config_dict : dict
            Configuration dictionary.
        mandatory_files : list of str
            List of required file keys in the format 'Category_Field'.

        Raises
        ------
        Exception
            If any mandatory file is missing.
        """
        for field in mandatory_files:
            filename = config_dict[field.split("_")[0]][field.split("_")[1]]
            if not os.path.isfile(filename):
                log_and_raise(logger, f"File '{field.split('_')[1]}' from '{field.split('_')[0]}' not found.", Exception)

    def import_from_file(
        self,
        configFilePath
    ):

        """
        Import parameters from a .json file and update configuration.

        If a parameter is not found in the .json file, the default is used. This method also calls
        check_for_mandatory_files() and calculate_parameters().

        Parameters
        ----------
        configFilePath : str
            Path to the .json file containing parameters.
        """

        if not configFilePath.endswith(".json"):
            logger.info(" If you are using an old file with the ending '.set', \
            run the script 'set_to_json.py', wich you can find in the folder 'PREPROC_POSTPROC'. \
            The script does need the path to the directory containing the old settings file and will create \
            recursively json-files that contain the same parameters as the old '.set' and '.bc' files.")
            logger.error("The given settings file is not a json file.")

        logger.info(f"Loading configuration from {configFilePath}")
        input_file      = open(configFilePath)
        input_data      = json.load(input_file)
        
        # Get an instance of the defaults settings
        default_config  = self.default_config

        # Loop over fields and overwrite defaults by file values
        for category in default_config:
            dict = Dotdict()
            for parameter in default_config[category]:
                if parameter in input_data[category]:
                    input_value         = input_data[category][parameter]
                    if parameter in ["EigenValueGuess","Omegas"]:
                        dict[parameter] = self.parse_complex_list(input_value)
                    else:
                        dict[parameter] = input_value
                else:
                    dict[parameter]     = default_config[category][parameter]["default"]
                    logger.warning(f'Field "{parameter}" missing from file, setting default: {dict[parameter]}')
            setattr(
                self,
                category,
                dict,
            )
        input_file.close()

        # Check if mandatory files are there [category_name]
        logger.debug("Checking mandatory files")
        mandatory_files     = ['BoundaryCondition_BCsFilePath','Case_MeshFilePath']
        tmp, extension      = os.path.splitext(input_data['FlowInput']['MeanFlowFilePath'])
        if extension == ".fel":
            mandatory_files.append('FlowInput_MeanFlowFilePath')
        self.check_for_mandatory_files(
            input_data,
            mandatory_files,
        )

        # Check if export folder exists
        if not os.path.isdir(input_data['Export']['ExportFolder']):
            os.makedirs(input_data['Export']['ExportFolder'])
        
        # Move the log file to the export folder
        current_dir = os.path.dirname(configFilePath)
        log_dir     = os.path.join(*[current_dir,input_data['Export']['ExportFolder'], "log"])
        Logger.change_log_location(log_dir)

        logger.debug("Creating Mixture class")
        self.Mixture = MixtureClass(
            self.Case["MixtureFilePath"],
            self.Case["SpeciesFilePath"],
        )
        
        # Get domain data and set BCs
        self.read_domain_data(
            self.Case["MeshFilePath"],
            self.Case["nDim"],
            self.get_extended_transported_quantity_list(),
            self.Case["CoordinateSystem"],
            self.Case["m"],
        )

        # Hardcoded parameters
        self.debug                      = True # specify here if printDebug messages should be shown
        self.Numerics.nCPU              = 1 # hardcoded for now, move to defaults later
        self.Numerics.NumericalScheme   = "Continuous Galerkin" # hardcoded for now, move to defaults later
        
        # Calculate parameters
        self.calculate_parameters()
        logger.info("Configuration loaded successfully")

    def import_from_h5_file(
        self,
        h5FileName
    ):
        """
        Import parameters from an HDF5 file.

        This function is not used in the current version of FELiCS. It used to be called when reading
        'meanflow.h5' instead of a .fel file, in which case it would overwrite the parameters.

        Parameters
        ----------
        h5FileName : str
            Path to the HDF5 file containing parameters.

        Returns
        -------
        None or updates internal configuration attributes.
        """
        # NOTE: This function is not used in the current version of FELiCS
        # NOTE: It used to be called when reading 'meanflow.h5' instead of a .fel
        # NOTE: In that case it was overwritting the parameters.
        # TODO: load the parameters into a "data" dictionnary
        # similar to what we get from loading a .json

        hf = File(
            h5FileName,
            'r',
        )
        if 'param' not in hf.keys():
            return None
        settingsDict = self.get_all_settings_dict()
        if self._settingsKind in hf['param'].keys():
            for settingsParameter in list(settingsDict.keys()):
                if settingsParameter in hf[f'param/{self._settingsKind}'].attrs.keys():
                    groupName   = f'param/{self._settingsKind}'
                    datatype    = settingsDict[settingsParameter]['datatype']
                    value       = hf[groupName].attrs[settingsParameter]
                    #pdb.set_trace()
                    if 'int' in str(datatype):
                        try:
                            value = int(float(value))
                            exec(f'self.{settingsParameter} = {value}')
                        except:
                            pdb.set_trace()
                    elif 'str' in str(datatype):
                        exec(f'self.{settingsParameter} = "{value}"')
        hf.close()

        # NOTE: This will not be needed anymore
        self.Mixture = MixtureClass(
            self.Case.mixtureFilePath,
            self.SpeciesFilePath,
        )


    def export(
        self,
        filestring
    ):
        """
        Export the current configuration parameters to a file.

        Parameters
        ----------
        filestring : str
            Path of the parameter file to export to. Supports .h5 and text files.

        Notes
        -----
        This function is currently not fully working for all export types.
        """

        # Writing the parameters to a file
        if '.h5' in filestring:
            file = File(
                filestring,
                'a',
            )

            if 'parameters' not in file.keys():
                paramGroup = file.create_group('parameters')
            else:
                paramGroup = file['parameters']
            #dt = string_dtype()
            # add the fenics-version:
            if 'FELiCSVersion' not in paramGroup.keys():
                felicsVersion = paramGroup.create_group('FELiCSVersion')
            else:
                felicsVersion = paramGroup['FELiCSVersion']
            felicsVersion.attrs.create(
                'FELiCSVersion',
                data=get_last_git_commit(),
            )
        else:
            file = open(
                filestring,
                'w',
            )

        # Iterate over all attributes of the object
        for group in dir(self):
            # Only not internal (marked by _) and not callable attributes should be exported
            if not group.startswith('_') and not callable(eval('self.'+group)):
                # Iterate over all attributes of the object
                for parameter in dir(eval('self.'+group)):
                    # Only not internal (marked by _) and not callable or class attributes should be exported
                    if not parameter.startswith('_') and not callable(eval('self.'+group+'.'+parameter)) and not isclass(parameter):

                        # Check if the attributes are string, float, int, bool or list and export them
                        if type(eval('self.'+group+'.'+parameter)) == str:
                            if isinstance(
                                file,
                                File,
                            ):
                                if group not in list(paramGroup.keys()):
                                    currentGroup = paramGroup.create_group(group)
                                    # stringInArray = np.array([eval(f'self.{group}.{parameter}')])
                                    #

                                    currentGroup.attrs.create(
                                        parameter,
                                        data=eval(f'self.{group}.{parameter}'),
                                    )
                                else:

                                    #tringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
                                    file[f'parameters/{group}'].attrs.create(
                                        parameter,
                                        data=eval(f'self.{group}.{parameter}'),
                                    )
                            else:
                                file.writelines(parameter+'='+'\''+str(eval('self.'+group+'.'+parameter))+'\''+'\n')
                        elif type(eval('self.'+group+'.'+parameter)) in [float,int,bool,list,dict]:
                            if isinstance(
                                file,
                                File,
                            ):
                                if group not in list(paramGroup.keys()):
                                    currentGroup = paramGroup.create_group(group)
                                    #try:
                                    currentGroup.attrs.create(
                                        parameter,
                                        data=str(eval(f'self.{group}.{parameter}')),
                                    )
                                    # except:
                                    #   # stringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
                                    #
                                    #   currentGroup.attrs.create(parameter, data=str(eval(f'self.{group}.{parameter}')))
                                else:
                                    #try:
                                    file[f'parameters/{group}'].attrs.create(
                                        parameter,
                                        data=str(eval(f'self.{group}.{parameter}')),
                                    )
                                    # except:
                                    #   stringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
                                    #   file[f'param/{group}'].create_dataset(parameter, data=stringInArray)
                            else:
                                file.writelines(parameter+'='+str(eval('self.'+group+'.'+parameter))+'\n')
        file.close()


    def read_domain_data(
        self,
        MeshFile,
        gDim,
        ExtendedTransportedQuantityList,
        coordinate_system,
        m
    ):
        """
        Read all domain data from the mesh file and update mesh-related attributes.

        Parameters
        ----------
        MeshFile : str
            Path to the mesh file.
        gDim : int
            Geometric dimension of the mesh.
        ExtendedTransportedQuantityList : list
            List of transported quantities including velocity components.
        coordinateSystem : str
            Coordinate system type (e.g., 'Cartesian', 'Cylindrical').
        m : int
            Azimuthal wavenumber or mode.
        """
        logger.debug(f"Reading domain data from '{MeshFile}'")
        # Read mesh and store it in self.__mesh__
        if not MeshFile == '' and os.path.isfile(MeshFile):
            self.__mesh__ = FELiCSMesh(
                coordinate_system,
                MeshFile,
                gDim,
                m,
            )
            self.dim      = self.__mesh__.gdim

        self.Case.Equations            = self.get_equation_list()
        self.Case.StateVectorVariables = self.get_state_vector_variables(self.Case.Equations)
        

    def get_mesh(
        self
    ):
        """
        Return the mesh object associated with the configuration.

        Returns
        -------
        FELiCSMesh
            The mesh object.
        """
        return self.__mesh__

    def get_internal_velocity_components(
        self
    ):
        """
        Get a list of velocity components directed within the mesh dimensions.

        Returns
        -------
        list of str
            List of internal velocity component labels (e.g., ['x', 'y', 'z']).
        """
        if self.Case["CoordinateSystem"] == 'Cartesian':
            VelCompList = ['x','y']
            if self.Case.nDim > 2:
                VelCompList.append('z')
        elif self.Case["CoordinateSystem"] == 'Cylindrical':
            VelCompList = ['x','r']
        return VelCompList

    def get_external_velocity_components(
        self
    ):
        """
        Get a list of velocity components directed outside the mesh dimensions.

        Returns
        -------
        list of str
            List of external velocity component labels (e.g., ['z'], ['t'], or []).
        """
        if self.Case["CoordinateSystem"] == 'Cartesian':
            if self.Case.m != 0 and self.Case.nDim == 2:
                VelCompList = ['z']
            else:
                VelCompList = []
        elif self.Case["CoordinateSystem"] == 'Cylindrical':
            VelCompList = ['t']

        return VelCompList

    def get_velocity_components(
        self
    ):
        """
        Get a list of all velocity components, both internal and external.

        Returns
        -------
        list of str
            List of all velocity component labels.
        """
        templist = self.get_internal_velocity_components()
        templist.extend(self.get_external_velocity_components())

        return templist

    def get_transported_quantity_list(
        self
    ):
        """
        Get a list of all transported quantities for the current case settings.

        Returns
        -------
        list of str
            List of transported quantity variable names.
        """
        SolutionList = []
        if not self.Case["SetOfEquations"]['Momentum']['Variable'] == 'None':
            SolutionList.append(self.Case["SetOfEquations"]['Momentum']['Variable'])
        if not self.Case["SetOfEquations"]['Mass']['Variable'] == 'None':
            SolutionList.append(self.Case["SetOfEquations"]['Mass']['Variable'])
        if not self.Case["SetOfEquations"]['Energy']['Variable'] == 'None':
            SolutionList.append(self.Case["SetOfEquations"]['Energy']['Variable'])
        if "Custom1" in self.Case["SetOfEquations"]:
            if not self.Case["SetOfEquations"]['Custom1']['Variable'] == 'None':
                SolutionList.append(self.Case["SetOfEquations"]['Custom1']['Variable'])
        if "Custom2" in self.Case["SetOfEquations"]:
            if not self.Case["SetOfEquations"]['Custom2']['Variable'] == 'None':
                SolutionList.append(self.Case["SetOfEquations"]['Custom2']['Variable'])
        if 'Species' in list(self.Case["SetOfEquations"].keys()):
            if self.Case["SetOfEquations"]['Species']['Variable'] != 'None':
                for species in list(self.Mixture.get_species_list('transported')):
                    SolutionList.append(species)
        return SolutionList


    def get_equation_list(
        self
    ):
        """
        Get the list of equations with both defined equation and variable.

        Returns
        -------
        list of tuple
            List of (equation_name, equation_dict) tuples for active equations.
        """
        ### TODO Sophie: make this the central equation list and clean up
        # TODO Sophie: throw warning if Equation != None and Variable == None (before that: move equation of state out of those equations)
        EquationsList = []
        for equation in self.Case["SetOfEquations"].items():
            if not equation[1]["Equation"] == "None" and not equation[1]["Variable"] == "None":
                EquationsList.append(equation)
        return EquationsList


    def get_state_vector_variables(
        self,
        EquationList
    ):
        """
        Get the list of state vector variables and their components.

        Parameters
        ----------
        EquationList : list
            List of (equation_name, equation_dict) tuples.

        Returns
        -------
        list of tuple
            List of (variable, components) tuples.
        """
        ### TODO Sophie: make this the central list besides the equation list and clean up
        VariablesList = []
        for equation in EquationList:
            variable = equation[1]["Variable"]
            if variable in ["u", "rhou"]:
                components = self.get_velocity_components()
            else:
                components = []
            VariablesList.append((variable, components))
        return VariablesList 


    def get_extended_transported_quantity_list(
        self
    ):
        """
        Get a list of all transported quantities, including all velocity components.

        Returns
        -------
        list of str
            List of transported quantities with velocity components expanded.
        """
        SolutionList = []

        # First add all velocity components
        transportedQuantities = self.get_transported_quantity_list()
        if 'u' in transportedQuantities:
            for component in self.get_velocity_components():
                SolutionList.append('u'+component)

        # Then extend the list by the transported quantity list
        SolutionList.extend(self.get_transported_quantity_list())

        # Finally, remove the component 'u' if present
        if 'u' in SolutionList:
            SolutionList.remove('u')
        return SolutionList

    def get_mean_flow_field_names(
        self
    ):
        """
        Get the list of mean flow field names required for input.

        Returns
        -------
        list of str
            List of mean field variable names to be read from file.
        """
        MeanList=[]

        # Add velocity components
        MeanList.append('u')

        # Add pressure component
        MeanList.append('p')

        # Add density component (for cold flow)
        MeanList.append('rho')

        # If necessary, add density and enthalpy diffusion
        if 'rho' in self.get_transported_quantity_list():
            MeanList.append('rho')
            if self.Case.MolViscModel == 'File' or self.Case.MolViscPerturbModel == 'Sutherland mean':
                MeanList.append('alpha')

        # Add species which are transported
        for specie in self.Mixture.get_species_list('transported'):
            MeanList.append(specie)
            if self.Case.MolViscModel == 'File' or self.Case.MolViscPerturbModel == 'Sutherland mean':
                MeanList.append('D_'+specie)

        # If Input-Output analysis is used, the forcing must be read in (at least curently) for
        # every conservative variable ()...
        if self.Case.AnalysisMode in ['Input-Output']:
            CurrentList=MeanList.copy()
            for entry in CurrentList:

                # If the field is the velocity vector, all components must be considered
                if entry == 'u':

                    #for component in self.VelocityComponents:
                    MeanList.append(entry+'_forcing_r')#jvs
                    MeanList.append(entry+'_forcing_i')
                else:
                    MeanList.append(entry+'_forcing_r')
                    MeanList.append(entry+'_forcing_i')

        # If Resolvent analysis is used, we always read forcing and response domains
        if self.Case.AnalysisMode in ['Resolvent']:
            MeanList.append('responseDomain')
            MeanList.append('forcingDomain')

        # Always look for a sponge variable in the mean flow file
        # if not present it will be zero
        MeanList.append('spg')

        # Add species, which are not transported
        for specie in self.Mixture.get_species_list('constraint'):
            MeanList.append(specie)
        if self.Case.TurbulenceModel in ['File']:
            MeanList.append('nuturb')
            MeanList.append('muturb')
        if self.Case.MolViscModel in ['File'] or self.Case.MolViscPerturbModel in ['Sutherland mean']:
            MeanList.append('nulam')
            MeanList.append('mulam')
        if self.Case.Reaction:
            MeanList.append('dQ')
        return MeanList

    def get_n_velocity_components(
        self
    ):
        """
        Get the number of velocity components for the current configuration.

        Returns
        -------
        int
            Number of velocity components.
        """
        return len(self.get_velocity_components())
