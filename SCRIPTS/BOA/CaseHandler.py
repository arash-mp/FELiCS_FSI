import GeometryDeformer as geo

class CaseHandler:

    def __init__(self,  settingsFile, param, mesh, boundaryFacets):

        # read non-FELiCS-specific parameters
        BOACaseDict={'BOACase':{'datatype':str},\
                     'MolViscRampFactors':{'datatype':float}}
        param.Case.importSettings(settingsFile, BOACaseDict) 

        # set default values, if they were not present in the settings file
        if not hasattr(param.Case, 'BOACase'):
            param.Case.BOACase = 'CylinderBSpline2Pts'
        if not hasattr(param.Case, 'MolViscRampFactors'):
            param.Case.MolViscRampFactors = '[1.]'

        self.case = param.Case.BOACase

        if self.case == 'CylinderBSpline2Pts':
            self.targetValuesForSponge    = [1.,0.,0.] #ux, uy, p
            self.initialValuesForBaseFlow = [1.,0.,0.] #ux, uy, p
            self.geometryDeformer = geo.CylinderBSpline2Pts(mesh,boundaryFacets)

        elif self.case == 'ProfileBSpline8Pts':
            self.targetValuesForSponge    = [1.5,0.,0.] #ux, uy, p
            self.initialValuesForBaseFlow = [1.5,0.,0.] #ux, uy, p
            self.geometryDeformer = geo.ProfileBSpline8Pts(mesh,boundaryFacets)

        elif self.case == 'ProfileBSpline2Pts':
            self.targetValuesForSponge    = [1.5,0.,0.] #ux, uy, p
            self.initialValuesForBaseFlow = [1.5,0.,0.] #ux, uy, p
            self.geometryDeformer = geo.ProfileBSpline2Pts(mesh,boundaryFacets)


    def getTargetValuesForSponge(self):
        return self.targetValuesForSponge

    def getInitialValuesForBaseFlow(self):
        return self.initialValuesForBaseFlow

    def getGeometryDeformer(self):
        return self.geometryDeformer




