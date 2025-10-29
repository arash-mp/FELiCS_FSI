
import  h5py
import  numpy as np
from    pathlib             import Path
from    scipy.interpolate   import griddata
from    FELiCS.Fields.Field import Field
from    FELiCS.Misc.logging import Logger

from FELiCS.SpaceDisc.FEMSpaces import getFELiCSSpace

# Get the logger
logger = Logger.get_logger("felics")

class Writer:
    
    def __init__(self, mesh, exportFolder = "Output"):
        """
        Function arguments:
        - mesh: FELiCSMesh object 
        - exportFolder: str, optional 
               folder in which the files will be saved

        Function returns:

        """

        # Notes:
        # - the writer is specified for a mesh and an export folder
        # - when initialized, the writer creats the export folder and exports the mesh

        import os

        self.mesh         = mesh
        self.exportFolder = exportFolder
        self.exportMesh   = mesh.exportMesh
        self.exportSpace  = getFELiCSSpace(self.exportMesh, order = 1, dim = 1)  

        # create the export folder, if it does not already exist
        os.makedirs(exportFolder, exist_ok=True)

        # export the mesh into the export folder
        meshFileName = self.exportFolder + "/mesh.h5"
        meshPath     = Path(meshFileName)
        self.mesh.saveInFELiCSFormat(meshFileName)    


    def exportListOfFieldsToH5(self, listOfFields, fileName):

        #-----------------------------------------------------------------------
        ## conversions to export fields 
        #-----------------------------------------------------------------------
        # get a list of purely scalar fields from the field list
        listOfSubFields = []
        for field in listOfFields:
            listOfSubFields.extend(self._getListOfScalarFields(field))

        # convert the list of (sub-)fields to export fields defined on the export mesh
        listOfExportFields = []
        for subField in listOfSubFields:
            listOfExportFields.append(self._getExportField(subField))

        return self._export(listOfExportFields, fileName) 



    def exportFieldToH5(self, field, fileName):

        #-----------------------------------------------------------------------
        ## conversions to export fields
        #-----------------------------------------------------------------------
        # get a list of purely scalar fields from the field 
        # (which could be scalar/vector/mixed)
        listOfSubFields = self._getListOfScalarFields(field)

        # convert the list of (sub-)fields to export fields defined on the export mesh
        listOfExportFields = []
        for subField in listOfSubFields:
            listOfExportFields.append(self._getExportField(subField))

        return self._export(listOfExportFields, fileName)


    def _export(self, listOfExportFields, fileName):



        #-----------------------------------------------------------------------
        ## write in  h5 file and simultaneously create xmf text 
        #-----------------------------------------------------------------------
        xmfText  = self._createXMFHeader()

        if len(self.exportFolder)>0:
            fileName = self.exportFolder + "/" + fileName 

        # h5 file path
        fileName.removesuffix(".h5") # this makes sure that also fileNames with the correct ending can be given
        h5FilePath = Path(fileName + ".h5")

        with h5py.File(h5FilePath, 'w') as h5:
            for exportField in listOfExportFields:
                if exportField.isReal():
                    h5.create_dataset(exportField.getName(), data = exportField.getRealCoefficientArray(),  dtype = np.float64)
                    xmfText += self._createXMFForScalarField(h5FilePath.name, exportField.getName(), exportField.getSize())
                    print(exportField.getName(), h5FilePath.name)
                else:
                    realName = exportField.getName()+"_real"
                    imagName = exportField.getName()+"_imag"
                    h5.create_dataset(realName, data = exportField.getRealCoefficientArray(),  dtype = np.float64)
                    h5.create_dataset(imagName, data = exportField.getImagCoefficientArray(),  dtype = np.float64)
                    xmfText += self._createXMFForScalarField(h5FilePath.name, realName, exportField.getSize())
                    xmfText += self._createXMFForScalarField(h5FilePath.name, imagName, exportField.getSize())


        #-----------------------------------------------------------------------
        ## create xmf file and return xmf text without footer (s.t. other fields can be extended)
        #-----------------------------------------------------------------------
        xmfTextFinal  = xmfText + self._createXMFFooter()
        Path(fileName+".xmf").write_text(xmfTextFinal, encoding="utf-8")

        return xmfText


    def _getListOfScalarFields(self, field):

        # convert vector fields and mixed fields into a list of scalar fields
        listOfSubFields = []

        if field.info['type'] == 'scalar':
            listOfSubFields.append(field)

        elif field.info['type'] == 'vector':
            names        = field.getNamesOfSubFields()
            degree       = field.space.ufl_element().degree
            space_scalar = getFELiCSSpace(field.mesh, order = degree, dim=1)
            for i in range(field.info['num_subspaces']):
                indices_mapping = field.space.sub(i).collapse()[1]
                field_scalar    = Field(space_scalar, field.mesh, name=names[i])
                field_scalar.setCoefficientArray(field.getCoefficientArray()[indices_mapping])
                listOfSubFields.append(field_scalar)

        elif field.info['type'] == 'mixed':
            listOfFields    = field.getListOfSubFields()
            for fieldMixed in listOfFields:
                listOfSubFields.extend(self._getListOfScalarFields(fieldMixed))

        else:
            raise Exception(f"Unknown field type: {field.info['type']}")


        return listOfSubFields


    def _getExportField(self, field):
        #only works for scalar fields
        # TODO: check real quick if the field mesh is the same as the export mesh and throw error?
        exportField = Field(self.exportSpace, field.mesh.exportMesh, name=field.getName())

        degree      = field.space.ufl_element().degree

        if degree == 2:
            # check if there is already a mapping
            # TODO: change this after export mesh updates from Simon
            from FELiCS.IO.Mapping import Mapping
            if not hasattr(self, "mappingP2ToExport"): 
                self.mappingP2ToExport = Mapping.calculateMappingFromSpaces(
                                             field.space, self.exportSpace)
            exportField.setCoefficientArray( field.getCoefficientArray()[self.mappingP2ToExport])

        else:
            # TODO: write warning that this may be slow and that data may be lost, if used on degree > 2
            self._interpolateWithGridData(field, exportField)

        return exportField


    def _interpolateWithGridData(self, field, exportField):
        ## Interpolating a P1 FELiCS field to a P1 export field using griddata:
        # field: input field 
        # exportField: export field for writing in file

        # Dimension of the source and export meshes
        source_dim          = field.mesh.dolfinxMesh.geometry.dim
        # Get the coordinates of the source and export meshes
        source_coords       =       field.space.tabulate_dof_coordinates()[:, :source_dim]
        export_coords       = exportField.space.tabulate_dof_coordinates()[:, :source_dim]
        # Interpolate the source field values to the export mesh coordinates
        interpolated_values = griddata(
                source_coords, 
                field.getCoefficientArray(), 
                export_coords, 
                method='linear'
            )
        # Set the interpolated values to the export field
        exportField.setCoefficientArray(interpolated_values)


    def _createXMFHeader(self):
        from textwrap import dedent
        import h5py

        meshPath = self.exportFolder + "/mesh.h5"
        meshPath_rel = "mesh.h5"

        if self.mesh.gdim == 2:
            cellStyle = "Triangle"
        else:
            cellStyle = "Tetrahedron"

        meshh5 = h5py.File(meshPath, 'r')
        tri = meshh5['cells']['triangles']

        xmf = dedent(f"""\
              <?xml version="1.0" ?>
              <Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
                <Domain>
                  <Grid Collection="grid"  Name="grid">
                    <Topology Type="{cellStyle}" NumberOfElements="{len(tri[:,0])}">
                       <DataItem ItemType="Function" Dimensions="{len(tri[0,:])*len(tri[:,0])}" Function="$0 - 0">
                          <DataItem Format="HDF" DataType="Int" Dimensions="{len(tri[0,:])*len(tri[:,0])}">
                              {meshPath_rel}:/cells/triangles
                          </DataItem>
                       </DataItem>
                    </Topology>
                    <Geometry Type="X_Y"> """)
         
        for key in list(meshh5['coordinates'].keys()):
            xmf = xmf + dedent(f"""  
                    <DataItem Format="HDF" ItemType="Uniform" Precision="8" NumberType="Float" Dimensions="{len(meshh5['coordinates'][key])}">
                    {meshPath_rel}:/coordinates/{key}
                    </DataItem> """)
        xmf = xmf + dedent(f"""
                </Geometry>""")
        return xmf

    def _createXMFForScalarField(self, h5FileName, fieldName, size, fieldNameInFile=None):
        from textwrap import dedent

        if fieldNameInFile==None:
            fieldNameInFile = fieldName

        xmf = dedent(f""" 
            <Attribute Name="{fieldName}" Center="Node" AttributeType="Scalar">                  
               <DataItem ItemType="Function" Dimensions="{size}" Function="$0">                      
                   <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="{size}"> 
                       {h5FileName}:{fieldNameInFile}                                           
                   </DataItem>                                                                      
               </DataItem>                                                                          
            </Attribute>  
            """)
        return xmf


    def _createXMFFooter(self):
        from textwrap import dedent
        return dedent(f"""\
                    </Grid>
                  </Domain>
                </Xdmf>
                """)

