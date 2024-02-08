# Third party libraries
import numpy as np


class laminarDiffusionHandler:
    """
    This class is used to initialize the laminar diffusion.

    Parent classes:

    Child classes:
    - fluctuationClass
    - fluctuationSolution

    Private attributes:

    Protected attributes:

    Public attributes:

    """

    def __init__(
            self,
            param,
            mean
    ):

        from dependentVariables.getAlpha import getAlpha
        '''
        Adding laminar diffusions to the fieldDict

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        '''
        if (param.Case.MolViscPerturbModel == 'Sutherland') \
                and ('T' in list(self._fieldDict.keys())):

            self._fieldDict['nulam'] = self.getSutherlandViscFluctuationClass(
                mean,
                self._fieldDict,
                self._meanfieldDict,
                self.rho
            )
            self._fieldDict['alpha'] = getAlpha(
                param,
                mean,
                self._fieldDict,
                self._meanfieldDict,
                self.rho,
                isMeanFlowClass=False
            )
            for specie in param.Case.Mixture.getSpeciesList('transported'):
                Sc = param.Case.Mixture.species[specie]['Sc']
                self._fieldDict['D_' + specie] = self._fieldDict['nulam'] / Sc

        elif param.Case.MolViscPerturbModel == 'Sutherland mean':
            #self._fieldDict['nulam'], fluct = self.getSutherlandMeanVisc(
            foobar, fluct = self.getSutherlandMeanVisc(
                mean,
                self._meanfieldDict if hasattr(self, "_meanfieldDict") else None,
                self.rho
            )
            #self._fieldDict['alpha'] = getAlpha(
            #    param,
            #    mean,
            #    self._fieldDict,
            #    self._meanfieldDict if hasattr(self, "_meanfieldDict") else None,
            #    self.rho,
            #    isMeanFlowClass=False
            #)
            for specie in param.Case.Mixture.getSpeciesList('transported'):
                if isinstance(mean.T, np.ndarray):
                    if mean.T.shape[0] != self.rho.shape[0]:
                        local_mean_D = self._meanfieldDict[f'D_{specie}'].compute_vertex_values()
                    else:
                        local_mean_D = mean.D(specie)
                else:
                    local_mean_D = mean.D(specie)
                self._fieldDict['D_' + specie] = local_mean_D * fluct
