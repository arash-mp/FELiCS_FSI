from FELiCS.GUI.BCsSettingsClass import FELiCSMesh

class Mesh(FELiCSMesh):

    def __init__(self, parameters):
        #this is the mesh, initialized the old way (in the GUI):
        #mesh = parameters.BCs.getMesh()
        #

        super().__init__(parameters.Case.CoordinateSystem,
                         parameters.Case.MeshFilePath,
                         parameters.Case.nDim,
                         parameters.Case.m,
                         )
