
import  h5py
import  numpy as np
from    pathlib             import Path
from    scipy.interpolate   import griddata
from    FELiCS.Fields.Field import Field
from    FELiCS.Misc.logging import Logger
from    FELiCS.IO.Mapping   import Mapping
from    FELiCS.SpaceDisc.FEMSpaces import createFunctionSpace


# Get the logger
logger = Logger.get_logger("felics")

class Writer:
    
    def __init__(self, mesh, exportDir = "Output"):
        """
        Function arguments:
        - mesh: FELiCSMesh object 
        - exportDir: str, optional 
               directory in which the files will be saved

        Function returns:

        """

        # Notes:
        # - the writer is specified for a mesh and an export directory
        # - when initialized, the writer creates the export directory and exports the mesh

        import os

        self.mesh         = mesh
        self.exportFolder = exportDir
        self.exportMesh   = mesh.exportMesh
        self.exportSpace  = createFunctionSpace(self.exportMesh, degree = 1, dim = 1)  

        # create the export folder, if it does not already exist
        os.makedirs(exportFolder, exist_ok=True)

        # export the mesh into the export folder
        meshFileName = self.exportFolder + "/mesh.h5"
        meshPath     = Path(meshFileName)
        self.mesh.saveInFELiCSFormat(meshFileName)    


    def exportListOfFieldsToH5(self, listOfFields, fileName, attributes = None):
        """
        Function arguments:
        - listOfFields: list of Field (or child class) objects
               fields that are to be written into one H5 file
        - fileName: str
               filename of H5 file (and of corresponding xmf file)

        Function returns:

        """

        # get a list of purely scalar fields from the field list
        # (which could be scalar/vector/mixed fields)
        listOfSubFields = []
        for field in listOfFields:
            listOfSubFields.extend(self._getListOfScalarFields(field))

        return self._export(listOfSubFields, fileName, attributes) 


    def exportFieldToH5(self, field, fileName, attributes = None):
        """
        Function arguments:
        - fields: Field (or child class) objects
               field that is to be written as H5 file
        - fileName: str
               filename of H5 file (and of corresponding xmf file)

        Function returns:

        """
        # get a list of purely scalar fields from the field 
        # (which could be a scalar/vector/mixed field)
        listOfSubFields = self._getListOfScalarFields(field)

        return self._export(listOfSubFields, fileName, attributes)


    def _export(self, listOfScalarFields, fileName, attributes):

        #-----------------------------------------------------------------------
        ## conversions to export fields
        #-----------------------------------------------------------------------
        # convert the list of (sub-)fields to export fields defined on the export mesh
        listOfExportFields = []
        for field in listOfScalarFields:
            listOfExportFields.append(self._getExportField(field))

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
                    h5.create_dataset(exportField.name, data = exportField.getRealCoefficientArray(),  dtype = np.float64)
                    xmfText += self._createXMFForScalarField(h5FilePath.name, exportField.name, exportField.getSize())
                else:
                    realName = exportField.name+"_real"
                    imagName = exportField.name+"_imag"
                    h5.create_dataset(realName, data = exportField.getRealCoefficientArray(),  dtype = np.float64)
                    h5.create_dataset(imagName, data = exportField.getImagCoefficientArray(),  dtype = np.float64)
                    xmfText += self._createXMFForScalarField(h5FilePath.name, realName, exportField.getSize())
                    xmfText += self._createXMFForScalarField(h5FilePath.name, imagName, exportField.getSize())
            if attributes !=None:
                for attr in attributes:
                    h5.create_dataset(attr.name, data = attr.value)



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
            space_scalar = createFunctionSpace(field.mesh, degree = degree, dim=1)
            for i in range(field.info['num_subspaces']):
                indices_mapping = field.space.sub(i).collapse()[1]
                field_scalar    = Field(space_scalar, field.mesh, name=names[i])
                field_scalar.setCoefficientArray(field.getCoefficientArray()[indices_mapping])
                listOfSubFields.append(field_scalar)

        elif field.info['type'] == 'mixed':
            names        = field.getNamesOfSubFields()
            for i in range(field.info['num_subspaces']):
                degree          = field.space.sub(i).ufl_element().degree
                num_subspaces   = field.space.sub(i).num_sub_spaces
                space_scalar    = createFunctionSpace(field.mesh, degree = degree, dim=1)
                if num_subspaces == 0: # mapping of scalar field
                    indices_mapping = field.space.sub(i).collapse()[1]
                    field_scalar    = Field(space_scalar, field.mesh, name=names[i])
                    field_scalar.setCoefficientArray(field.getCoefficientArray()[indices_mapping])
                    listOfSubFields.append(field_scalar)
                else:
                    axisNames = field.mesh.axisNames
                    for j in range(num_subspaces): # mapping of vector field
                        indices_mapping = field.space.sub(i).sub(j).collapse()[1]
                        field_scalar    = Field(space_scalar, field.mesh, name=names[i]+axisNames[j])
                        field_scalar.setCoefficientArray(field.getCoefficientArray()[indices_mapping])
                        listOfSubFields.append(field_scalar)

        else:
            raise Exception(f"Unknown field type: {field.info['type']}")

        return listOfSubFields


    def _getExportField(self, field):
        #only works for scalar fields
        # TODO: check real quick if the field mesh is the same as the export mesh and throw error?
        exportField = Field(self.exportSpace, field.mesh.exportMesh, name=field.name)

        degree      = field.space.ufl_element().degree

        if degree == 2:
            # check if there is already a mapping
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

