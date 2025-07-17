
import  dolfinx
from    basix.ufl           import mixed_element
from    scipy.interpolate   import griddata
from    mpi4py              import MPI
from    FELiCS.Fields.Field import Field
from    FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class Writer:
    
    def __init__(self, param, FEMSpaces):
        """
        Function arguments:
        - FEMSpaces: Object containing the Cacluation FEM-Spaces, the export
            FEM-Spaces and the corresponding meshes.
        - param: Object, containing the parameters of the calculation.

        Function returns:

        """
        self._FEMSpaces     = FEMSpaces
        self._param         = param
        self._exportMesh    = FEMSpaces.exportMesh


    # Method to return fields on the export mesh
    def getFieldsOnExportMesh(self, sourceField):
        """
        Returns a dictionary of fields on the export mesh.
        TODO: Include field names

        Parameters
        ----------
        fieldDict : dict
            Dictionary containing fields to be exported.

        Returns
        -------
        dict
            Dictionary of fields on the export mesh.
        """
        
        # Type of FunctionSpace in the Field
        infoSpaceField          = sourceField.describeFunctionSpace()
        
        # Mappings from original to export mesh
        dofsMappingP2P1         = self._FEMSpaces.mappingObj.P2CalcToP1ExportIndecies # NOTE: This mappings seems to work for scalars on single fields and velocity components from mean flow
        VectorCalcToP1ExportIndecies    = self._FEMSpaces.mappingObj.VectorCalcToP1ExportIndecies # NOTE: this mapping seems to work only for vectors from the mixed space
        
        # Prepare the export FunctionSpace
        exportSpaceScalar       = self._FEMSpaces.P1Export
        exportSpaceVector       = self._FEMSpaces.FunctionSpaceVectorVelocityExport
        exportVMixed            = []
        
        def add_export_space(space_info, exportVMixed):
            """Helper function to add export spaces based on space type."""
            if space_info['type'] == 'scalar':
                exportVMixed.append(exportSpaceScalar)
            elif space_info['type'] == 'vector':
                exportVMixed.append(exportSpaceVector)
            else:
                raise Exception(f"Unknown space type: {space_info['type']}")
            return exportVMixed
        
        # Create an export field with the same functionSpace structure as the original field
        # At this point the functions are empty, but the structure is defined
        if infoSpaceField['type'] == 'mixed':
            # For mixed type, we must create a VMixed space with P1 versions of subspaces
            for subspace in infoSpaceField['subspaces']:
                exportVMixed    = add_export_space(subspace, exportVMixed)
            # Create the mixed function space from the export spaces
            exportVMixed_ufl    = mixed_element([space.ufl_element() for space in exportVMixed])
            exportMixedSpace    = dolfinx.fem.functionspace(self._exportMesh.dolfinxMesh, exportVMixed_ufl)
            # Create the export field
            exportField         = Field(exportMixedSpace, self._exportMesh, name=sourceField.name)
        else:
            exportVMixed        = add_export_space(infoSpaceField, exportVMixed)
            exportField         = Field(exportVMixed[0], self._exportMesh, name=sourceField.name)
            
        def interpolate_ScalarP1felics_to_P1export(sourceField, exportField):
            """ Interpolates a P1 FELiCS field to a P1 export field using griddata."""
            # Dimension of the source and export meshes
            source_dim          = sourceField.function.function_space.mesh.geometry.dim
            # Get the coordinates of the source and export meshes
            source_coords       = sourceField.function.function_space.tabulate_dof_coordinates()[:, :source_dim]
            export_coords       = exportField.function.function_space.tabulate_dof_coordinates()[:, :source_dim]
            # Interpolate the source field values to the export mesh coordinates
            interpolated_values = griddata(
                source_coords, 
                sourceField.getCoefficientArray(), 
                export_coords, 
                method='linear'
            )
            # Set the interpolated values to the export field
            exportField.setCoefficientArray(interpolated_values)
            return exportField
            
        # Function to map the sub-fields of a vector field
        def map_vector_field(originField, destinationField, dofsmapping):
            """Map a vector field from the source mesh to the export mesh."""
            # Split the fields into sub-fields
            sub_destinationFields   = destinationField.getListOfSingleFields()
            sub_originField         = originField.getListOfSingleFields()
            # Loop over sub-fields and set arrays with mapping
            for iSub in range(len(sub_destinationFields)):
                sub_destinationFields[iSub].setCoefficientArray(sub_originField[iSub].getCoefficientArray()[dofsmapping])
            # Assemble the export field from sub-fields
            destinationField.setListOfSingleFields(sub_destinationFields)
            return destinationField
                
        # Populate the export field based on the type of source field
        if infoSpaceField['type'] == 'mixed':
            list_ofExportFields     = exportField.getListOfSingleFields()
            list_ofSourceFields     = sourceField.getListOfSingleFields()
            
            # Loop over sub-fields and set arrays
            for iField in range(len(list_ofExportFields)):
                subFieldSource      = list_ofSourceFields[iField]
                degreeSource        = subFieldSource.function.function_space._ufl_element.degree
                
                if infoSpaceField['subspaces'][iField]['type'] == 'vector':
                    if degreeSource == 2:
                        list_ofExportFields[iField] = map_vector_field(subFieldSource, list_ofExportFields[iField], VectorCalcToP1ExportIndecies)
                    else:
                        raise Exception(f"Element degree {degreeSource} not supported for export vector.")
                
                elif infoSpaceField['subspaces'][iField]['type'] == 'scalar':
                    if degreeSource == 1:
                        # Interpolate the function values on the export mesh TODO: check if this works
                        list_ofExportFields[iField] = interpolate_ScalarP1felics_to_P1export(subFieldSource, list_ofExportFields[iField])
                    elif degreeSource == 2:
                        list_ofExportFields[iField].setCoefficientArray(subFieldSource.getCoefficientArray()[dofsMappingP2P1])
                    else:
                        raise Exception(f"Element degree {degreeSource} not supported for export.")
                    
            # Assemble all the interpolated/mapped subfields into the export field
            exportField.setListOfSingleFields(list_ofExportFields)    
        
        elif infoSpaceField['type'] == 'scalar':
            degree = sourceField.function.function_space._ufl_element.degree
            if degree == 1:
                # Interpolate the function values on the export mesh TODO: check if this works
                exportField = interpolate_ScalarP1felics_to_P1export(sourceField, exportField)
            elif degree == 2:
                exportField.setCoefficientArray(sourceField.getCoefficientArray()[dofsMappingP2P1])
            else:
                raise Exception(f"Element degree {degree} not supported for export.")

        elif infoSpaceField['type'] == 'vector':
            degree = sourceField.function.function_space._ufl_element.degree
            if degree == 2:
                # NOTE: For some reason I need to use "dofsMappingP2P1" when using a vector from meanflow
                # NOTE: but then I need to use "VectorCalcToP1ExportIndecies" when using a vector from the fluctuations
                exportField = map_vector_field(sourceField, exportField, VectorCalcToP1ExportIndecies)
            else:
                raise Exception(f"Element degree {degree} not supported for export vector.")
        
        return exportField
    
    # Method to export the fields to xdmf
    def writeFieldToXDMF(self, field, filename):
        
        # DolfinX export mesh
        exportMesh          = self._exportMesh.dolfinxMesh
        
        # Check the info from the field
        infoSpaceField      = field.describeFunctionSpace()
        
        # Set the names in the fem functions and save
        if infoSpaceField['type'] == 'vector' or infoSpaceField['type'] == 'scalar':
            field.function.name = field.name[0][0]
            with dolfinx.io.XDMFFile(MPI.COMM_WORLD, filename+".xdmf", "w") as xdmf:
                xdmf.write_mesh(exportMesh)
                xdmf.write_function(field.function)
            
        # For mixed spaces we need to iterate over the sub-fields
        elif infoSpaceField['type'] == 'mixed':
            list_ofFields   = field.getListOfSingleFields()
            with dolfinx.io.XDMFFile(MPI.COMM_WORLD, filename+".xdmf", "w") as xdmf:
                xdmf.write_mesh(exportMesh)
                for i in range(len(list_ofFields)):
                    list_ofFields[i].function.name = list_ofFields[i].name[0][0]
                    xdmf.write_function(list_ofFields[i].function)
        else:
            raise Exception(f"Unknown field type: {infoSpaceField['type']}")
        
        
        
        