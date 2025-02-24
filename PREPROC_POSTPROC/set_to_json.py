import json
import os
import sys

'''
This is a script to transform FELiCS input files to the .json format.

Run the script via

"python3 set_to_json.py <path_to_directory>".

The script will recursively (i.e. in the directory & in all the subdirectories) create json-files fron '.set' and '.bc' files with the same name. 
No old files will be deleted.

CAUTION:
    - Does not work if your boundary file and your settings file have the same name (e.g. case.set and case.bc).

    - Does not convert any mixture files. Unfortnuately, those have to be converted by hand. For examples see the felics-test repository.

    - Does use the command "eval" for the conversion of the boundary file, which is a "security risk", meaning if you have executable code in your 
      ".bc" file it could execute that. Since the chances are basically zero just a warning is written here:D

'''


dirname = sys.argv[1]
ext_set = ".set"
ext_bc  = ".bc"
#ext_mix = ".mix"

# remove the slash behind the directory name if necessary
if dirname[len(dirname)-1] == "/":
    dirname = dirname[0:len(dirname)-1]


print("CONVERTING OLD FELICS '.SET' AND '.BC' FILES RECURSIVELY TO JSON, NO OLD FILES WILL BE DELETED....")


# get a list of all files
list_set = []
list_bc  = []
list_mix  = []
for path, dirc, files in os.walk(dirname):
    for name in files:
        if name.endswith(ext_set):
            list_set.append(path+"/"+name)
        elif name.endswith(ext_bc):
            list_bc.append(path+"/"+name)
        elif name.endswith(ext_mix):
            list_mix.append(path+"/"+name)


#def adapt_set_to_json(set_files, set_files_directory, json_files_directory):
def adapt_set_to_json(set_files):
    # transforms an unstructured .set config file to a structured .json config file
    for set_file in set_files:
        # initialization
        input_data  = {}
        result      = {'BoundaryCondition':{}, 'Case':{}, 'Export':{}, 'FlowInput':{}, 'IOResolvent':{}, 'Numerics':{}}
        missing     = {}

        # open old .set file
        #with cd(set_files_directory):
        with open(set_file) as file:
            lines = file.readlines()
            for line in lines:
                key = line.split('=')[0].strip()
                value = line.split('=')[1].strip()
                exec('input_data[key] = '+value)

        # loop over default parameters and import parameters from .set config file
        SettingsDict = getAllSettingsDict()
        for category in SettingsDict:
            for parameter in SettingsDict[category]:
                if parameter in input_data:
                    if parameter in ["BCsFilePath"]:
                        bcFileName = input_data[parameter].removesuffix(".bc")
                        bcFileName += ".json"
                        input_data[parameter] = bcFileName
                    if parameter in ["EigenValueGuess", "Omegas"]:
                        if type(input_data[parameter]) is list:
                            for i,number in enumerate(input_data[parameter]):
                                if type(number) is complex:
                                    input_data[parameter][i] = str(number)[1:-1] # Remove parentheses from complex numbers
                                elif type(number) is float:
                                    input_data[parameter][i] = number
                        else:
                            if type(input_data[parameter]) is complex:
                                input_data[parameter] = [str(input_data[parameter])[1:-1]] # Remove parentheses from complex numbers
                            elif type(input_data[parameter]) is float:
                                input_data[parameter] = [str(input_data[parameter])]
                            
                        result[category][parameter] = str(input_data[parameter])
                    result[category][parameter] = input_data[parameter]
                else:
                    missing[category] = [parameter]
                    print("---- No input found for parameter "+parameter+", of category "+category+"! added to missing.json")

        # save new .json config file
        filename, file_extension = os.path.splitext(set_file)
        json_file = filename + ".json"
        directory, filename_pure = os.path.split(filename)
        json_file_missing= directory + "/missing.json"
        with open(json_file, 'w', encoding='utf-8') as output_file:
            json.dump(result, output_file, ensure_ascii=False, indent=4)
        with open(json_file_missing, 'w', encoding='utf-8') as missing_file:
            json.dump(missing, missing_file, ensure_ascii=False, indent=4)
        print("CONVERTED: "+set_file+"  =>  "+json_file)
        

def adapt_bc_to_json(bc_files):
    # transforms a .bc boundaries file to a .json boundaries file
    for bc_file in bc_files:
        result = {}
        # open .bc boundaries file
        with open(bc_file) as file:
            line = file.readlines()[0]
            main_dict = eval(line)
            # loop over file and invert
            for variable, ID_list in main_dict.items():
                for ID_dict in ID_list:
                    ID = ID_dict["ID"]
                    ID_dict.pop("ID", None)
                    ID_dict["variable"] = variable
                    if ID not in result:
                        result[ID] = []
                    result[ID].append(ID_dict)

        # save .json boundaries file
        filename, file_extension = os.path.splitext(bc_file)
        json_file = filename + ".json"
        with open(json_file, 'w', encoding='utf-8') as output_file:
            json.dump(result, output_file, ensure_ascii=False, indent=4)
        print("CONVERTED: "+bc_file+"  =>  "+json_file)


def getAllSettingsDict():
    """
    return SettingsDict.

    Returns
    -------
    SettingsDict : Dictionary 
        Dictionary of input parameters structured in subcategories
    """
    SettingsDict={
        'BoundaryCondition':{
            'BCsFilePath':          {'datatype':str,    'default':''}
        },
        'Case':{
            'AnalysisMode':         {'datatype':str,    'default':'Modal'},
            'CalculateAdjoint':     {'datatype':bool,   'default':True},
            'CoordinateSystem':     {'datatype':str,    'default':'Cartesian'},
            'm':                    {'datatype':int,    'default':0},
            'MeshFilePath':         {'datatype':str,    'default':''},
            'MixtureFilePath':      {'datatype':str,    'default':''},
            'MolVisc':              {'datatype':int,    'default':0.0},
            'MolViscModel':         {'datatype':str,    'default':'Constant'},
            'nDim':                 {'datatype':int,    'default':2},
            'Reaction':             {'datatype':bool,   'default':False},
            'SetOfEquations':       {'datatype':dict,
                'default':{
                    'Momentum':         {'Equation':'NSPrimitive',  'Variable':'u'},
                    'Mass':             {'Equation':'Continuity',   'Variable':'p'},
                    'Energy':           {'Equation':'None',         'Variable':'None'},
                    'Species':          {'Equation':'None',         'Variable':'None'},
                    'EquationOfState':  {'Equation':'None',         'Variable':'None'}
                }
            },
            'SpeciesFilePath':      {'datatype':str,    'default':''},
            'TransVelFluc':         {'datatype':bool,   'default':False},
            'TurbulenceModel':      {'datatype':str,    'default':'None'}
        },
        'Export':{
            'ExportFolder':         {'datatype':str,    'default':''},
            'Video':                {'datatype':bool,   'default':False},
        },
        'FlowInput':{
            'AveragingDirection':   {'datatype':str,    'default':'None'},
            'MeanFlowFilePath':     {'datatype':str,    'default':''},
        },
        'IOResolvent':{
            'ForcingBoundaryIndices':   {'datatype':list,   'default':[]},
            'ForcingMode':              {'datatype':str,    'default':'Body'},
            'ForcingNorm':              {'datatype':str,    'default':'TKE'},
            'ResponseNorm':             {'datatype':str,    'default':'TKE'},
            'Omegas':                   {'datatype':list,   'default':[]},
            'ResponseCoeff':            {'datatype':list,   'default':[]}
        },
        'Numerics':{
            'EigenValueGuess':      {'datatype':list,   'default':[1.0]},
            'nSolut':               {'datatype':int,    'default':3},
            'PolynomialOrder':      {'datatype':dict,   'default':{'u':'2'},    'options':[1,2]}
        }
    }
    return SettingsDict


class cd:
    ''' Context manager for changing the current working directory '''
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        if not os.path.exists(self.newPath):
            os.makedirs(self.newPath)
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


# function calls
adapt_set_to_json(list_set)
adapt_bc_to_json(list_bc)
