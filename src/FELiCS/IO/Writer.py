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
    """
    Utility class for exporting FELiCS fields and meshes to HDF5/XDMF.

    A writer instance is associated with a mesh and an export directory.
    When initialized with a mesh, the mesh is immediately exported to
    ``mesh.h5`` in the target directory. If no mesh is given at
    initialization, the mesh from the first exported field is used and
    written when data is exported.

    Multiple HDF5/XDMF files can be created within the same export
    directory. All fields written by a single writer are assumed to be
    defined on the same mesh.

    **Initialize the Writer object**

    Parameters
    ----------
    mesh : FELiCSMesh or None, optional
        Mesh associated with this writer. If provided, the mesh is exported
        immediately to ``mesh.h5`` in the export directory.
    exportDir : str, optional
        Directory where all HDF5/XDMF output files and the mesh file are
        written. The directory is created if it does not exist.

    Attributes
    ----------
    mesh : FELiCSMesh or None
        Mesh associated with the writer. If initially ``None``, it is set
        from the first field passed to an export method.
    exportFolder : str
        Path to the directory used for all exported files.
    meshFileName : str
        Path to the exported mesh HDF5 file, typically
        ``<exportFolder>/mesh.h5``.
    exportMesh : object
        Mesh object used specifically for export, usually
        ``mesh.exportMesh`` of the associated FELiCS mesh.
    exportSpace : dolfinx.fem.FunctionSpace
        Finite element function space associated with the export mesh, used
        to represent exported scalar fields.
    mappingP2ToExport : ndarray, optional
        Mapping from P2 degrees of freedom of a source space to the export
        space. Created on demand during the first P2 export.

    Notes
    -----
    All fields exported by a given writer are assumed to be compatible with
    the same mesh. If fields are defined on different meshes, a separate
    :class:`Writer` instance with a different export directory must be used
    for each mesh; otherwise the generated XMF files will not be readable
    correctly.
    """
    def __init__(self, mesh = None, exportDir = "Output"):
        """
        Initialize a Writer instance.

        Parameters
        ----------
        mesh : FELiCSMesh or None, optional
            Mesh to be exported and associated with this writer. If ``None``,
            the mesh is taken from the first field provided to an export
            method.
        exportDir : str, optional
            Directory where all exported files (mesh, HDF5, XMF) are stored.
            The directory is created if it does not already exist.
        """
        # Notes for docstring:
        # - the writer is specified for a mesh and an export directory
        # - when initialized, the writer creates the export directory and exports the mesh, if given
        # - if no mesh is given, the writer exports the mesh from the first field that is exported
        # - there can be as many exports in one export folder as is needed
        # - CAUTION: if any of the exported fields are defined on a different mesh, a separate Writer object
        #            has to be created, which has a different export directory, else reading the xmf files will 
        #            not work

        import os

        self.mesh         = mesh
        self.exportFolder = exportDir

        # create the export folder, if it does not already exist
        os.makedirs(exportDir, exist_ok=True)

        # export the mesh into the export folder
        self.meshFileName = self.exportFolder + "/mesh.h5"
        if self.mesh != None:
            self.exportMesh   = mesh.exportMesh
            self.exportSpace  = createFunctionSpace(self.exportMesh, degree = 1, dim = 1)  
            self.mesh.saveInFELiCSFormat(self.meshFileName)    



    def exportListOfFieldsToH5(self, listOfFields, fileName, attributes = None):
        """
        Export a list of fields to a single HDF5 file and generate XMF content.

        Each field in the list may be scalar, vector, or mixed. Internally,
        the fields are decomposed into scalar subfields, mapped or
        interpolated to the export space, written to `<fileName>.h5`, and
        referenced in an associated XMF file `<fileName>.xmf`.

        Parameters
        ----------
        listOfFields : list of Field
            List of `Field` (or subclasses) to be written to a single HDF5
            file.
        fileName : str
            Base name of the HDF5/XMF files. The export directory is
            automatically prefixed. If the name ends with ``.h5``, the
            suffix is stripped before adding the correct extensions.
        attributes : list of objects, optional
            Additional attributes to store as datasets in the HDF5 file. Each
            element is expected to provide ``name`` and ``value`` attributes.

        Returns
        -------
        str
            XMF header/body text corresponding to the exported fields but
            without the final XMF footer. This can be reused to append
            additional attributes before closing the XMF document.
        """
        # get a list of purely scalar fields from the field list
        # (which could be scalar/vector/mixed fields)
        listOfSubFields = []
        for field in listOfFields:
            listOfSubFields.extend(self._getListOfScalarFields(field))

        return self._export(listOfSubFields, fileName, attributes) 


    def exportFieldToH5(self, field, fileName, attributes = None):
        """
        Export a single field to an HDF5 file and generate XMF content.

        The field may be scalar, vector, or mixed. It is internally
        decomposed into scalar subfields, mapped or interpolated to the
        export space, written to `<fileName>.h5`, and referenced in an
        associated XMF file `<fileName>.xmf`.

        Parameters
        ----------
        field : Field
            Field (or subclass) to be written to an HDF5 file.
        fileName : str
            Base name of the HDF5/XMF files. The export directory is
            automatically prefixed. If the name ends with ``.h5``, the
            suffix is stripped before adding the correct extensions.
        attributes : list of objects, optional
            Additional attributes to store as datasets in the HDF5 file. Each
            element is expected to provide ``name`` and ``value`` attributes.

        Returns
        -------
        str
            XMF header/body text corresponding to the exported field but
            without the final XMF footer.
        """
        # get a list of purely scalar fields from the field 
        # (which could be a scalar/vector/mixed field)
        listOfSubFields = self._getListOfScalarFields(field)

        return self._export(listOfSubFields, fileName, attributes)


    def _export(self, listOfScalarFields, fileName, attributes):
        """
        Core export routine for scalar fields.

        This internal method performs the following steps:

        1. If no mesh is set on the writer, it is taken from the first
           scalar field and exported to ``mesh.h5``.
        2. Each scalar field is mapped or interpolated to the export space,
           creating an export field.
        3. All export fields are written as datasets into `<fileName>.h5`.
           Complex-valued fields are split into real and imaginary datasets.
        4. XMF content is constructed for each dataset and written to
           `<fileName>.xmf`.

        Parameters
        ----------
        listOfScalarFields : list of Field
            Scalar fields to be exported. All fields are assumed to be
            defined on a mesh compatible with the writer's export mesh.
        fileName : str
            Base name of the HDF5/XMF files. The export directory is
            automatically prefixed. If the name ends with ``.h5``, the
            suffix is stripped before adding the correct extensions.
        attributes : list of objects or None
            Optional additional attributes to store in the HDF5 file as
            datasets. Each attribute object should provide ``name`` and
            ``value`` attributes.

        Returns
        -------
        str
            XMF header/body text corresponding to the exported fields but
            without the final XMF footer, allowing further extension.
        """
        #-----------------------------------------------------------------------
        ## save the mesh of the first field if no mesh has been given at initialization 
        #-----------------------------------------------------------------------
        if self.mesh == None:
            self.mesh         = listOfScalarFields[0].mesh
            self.exportMesh   = self.mesh.exportMesh
            self.exportSpace  = createFunctionSpace(self.exportMesh, degree = 1, dim = 1)  
            self.mesh.saveInFELiCSFormat(self.meshFileName)    

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
        """
        Decompose a field into a list of scalar subfields.

        Depending on the field type, this method behaves as follows:

        * **scalar**: returns a list containing the field itself.
        * **vector**: creates one scalar field per component.
        * **mixed**: iterates over all subspaces; subspaces with
          ``num_sub_spaces == 0`` are treated as scalar, while subspaces
          with multiple components (e.g. vector-valued) are split into
          scalar component fields.

        The scalar subfields reuse the original mesh but are defined on
        newly created scalar function spaces.

        Parameters
        ----------
        field : Field
            Field to decompose. Its type is inferred from
            ``field.info['type']``.

        Returns
        -------
        list of Field
            List of scalar `Field` objects representing all components of
            the original field.

        Raises
        ------
        Exception
            If ``field.info['type']`` is unknown or not supported.
        """
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
        """
        Map or interpolate a scalar field to the export mesh and space.

        For P2 fields, a precomputed index mapping
        (``mappingP2ToExport``) is used for an efficient transfer of
        coefficients. For other polynomial degrees, the data is
        interpolated using grid-based interpolation.

        Parameters
        ----------
        field : Field
            Scalar field defined on the original mesh and space.

        Returns
        -------
        Field
            New scalar field defined on the export mesh and export space,
            containing the mapped or interpolated coefficients.

        Notes
        -----
        This method assumes that the provided field is scalar. If the field
        lives on a mesh that is not compatible with the export mesh,
        interpolation is used and some loss of accuracy may occur.
        """
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
        """
        Interpolate a scalar field from its original mesh to the export mesh.

        This method uses :func:`scipy.interpolate.griddata` to interpolate
        the degrees of freedom of a scalar FELiCS field onto the degrees of
        freedom of the export field.

        Parameters
        ----------
        field : Field
            Input scalar field defined on the source mesh and function space.
        exportField : Field
            Target scalar field defined on the export mesh and export space.
            Its coefficient array is overwritten with the interpolated values.

        Notes
        -----
        The interpolation is performed in physical coordinates using linear
        interpolation. It may be slower and less accurate for higher-order
        fields or strongly distorted meshes. A warning is advisable when
        using this on fields with polynomial degree greater than 2.
        """
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
        """
        Create the XMF header and mesh description.

        The header includes:

        * XML prolog and XDMF root tags.
        * A grid collection with topology referencing the cell connectivity
          stored in ``mesh.h5:/cells/triangles``.
        * Geometry data items referencing the coordinate datasets in
          ``mesh.h5:/coordinates/*``.

        The cell type is automatically set to ``Triangle`` for 2D grids and
        ``Tetrahedron`` otherwise.

        Returns
        -------
        str
            XMF header string including topology and geometry sections but
            not the closing footer.
        """
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
        """
        Create an XMF attribute block for a scalar field.

        Parameters
        ----------
        h5FileName : str
            Name of the HDF5 file containing the scalar dataset (without
            path).
        fieldName : str
            Logical name of the field as it should appear in the XMF
            attribute.
        size : int
            Number of nodal values in the scalar dataset.
        fieldNameInFile : str or None, optional
            Name of the dataset inside the HDF5 file. If ``None``, the
            value of ``fieldName`` is used.

        Returns
        -------
        str
            XMF snippet describing the scalar attribute and its link to the
            HDF5 dataset.
        """
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
        """
        Create the XMF footer.

        The footer closes the grid, domain, and XDMF tags started in the
        header.

        Returns
        -------
        str
            XMF footer string that finalizes the XDMF document.
        """
        from textwrap import dedent
        return dedent(f"""\
                    </Grid>
                  </Domain>
                </Xdmf>
                """)

