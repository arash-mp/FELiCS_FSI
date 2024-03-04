# Third party libraries
import numpy as np

class reactionClass:
    def __init__(
                 self,
                 reactionData,
                 reactionType,
                 mean,
                 fluc,
                ):
        self._educts = reactionData['educts']
        self._stochiometricCoefficientsEducts = reactionData['stochiometricCoefficientsEducts']
        self._products = reactionData['products']
        self._stochiometricCoefficientsProducts = reactionData['stochiometricCoefficientsProducts']
        self._reactionType = reactionType
        if reactionType == 'KaiserCnF2023':
            specie = self._products[0]
            if fluc._isSolution:
                mean_RR_prefactor = mean.fieldDict['RR_prefactor']
                # The following if clause is a quick fix. Should be removed once the mean flow class has been changed to Tensor notation
                mean_Yspecie = mean.fieldDict[specie]
                if 'rho' in list(mean.fieldDict.keys()):
                    mean_rho = mean.fieldDict['rho']
                else:
                    #mean_rho = mean.oneField
                    mean_rho = mean.rho
            # else it is a tensor object and will be used to build a ufl formulation
            else:
                mean_RR_prefactor = mean.RR_prefactor
                mean_rho = mean.rho
                mean_Yspecie = mean.Y(specie)
        
            self._RR = mean_RR_prefactor * mean_rho * (fluc.Y(specie) - 2 * fluc.Y(specie) * mean_Yspecie)\
                     + mean_RR_prefactor * fluc.rho * (mean_Yspecie - mean_Yspecie * mean_Yspecie)
        else:
            raise Exception('Reaction mechanism ' + reaction['type'] + ' unknown!')
        self._determineReactionFormula

        

    def _determineReactionFormula(self):
        self._reactionFormula = ''
        for i_educt,educt in enumerate(self._educts):
            self._reactionFormula += str(self._stochiometricCoefficientsEducts) + educt + ' + '
        if not self._reactionFormula == '':
            self._reactionFormula = self._reactionFormula[:-3]
        self._reactionFormula += ' --> '
        for i_product,product in enumerate(self._products):
            self._reactionFormula += str(self._stochiometricCoefficientsEducts) + product + ' + '
        if not self._reactionFormula == '':
            self._reactionFormula = self._reactionFormula[:-3]
        
    @property
    def RR(self):
        return self._RR

    @property
    def species(self):
        return self._educts + self_products
        
    def sourceTerm(
            self,
            specie,
            ):
        if specie in self._educts:
            index = self._educts.index(specie)
            stoch = self._stochiometricCoefficientsEducts[index]
            return - self.RR * stoch
        elif specie in self._products:
            index = self._products.index(specie)
            stoch = self._stochiometricCoefficientsProducts[index]
            return self.RR * stoch
        else: 
            raise Exception(specie +' not a part of reaction!') 

class reactionHandler:
    """
    This class is used to manage the variables for the linear reactions. 
    It is appleid last when all other variables are already calculated.

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
                ):
        pass

    def _additionalFieldsToBeReadReaction(self):
        if self._param.Case.Mixture.reactionMechanism['type'] == 'KaiserCnF2023': 
            outList = ['RR_prefactor']
        elif self._param.Case.Mixture.reactionMechanism['type'] == 'None':
            outList = []
        else:
            raise Exception('Reaction mechanism ' + self._param.Case.Mixture.reactionMechanism['type'] + ' unknown!')
        return outList
    
    def _initializeReactions(
                        self,
                        mean = 'None',
                        ):
        """
        Initializing the class and link the reactions to the previsously
        determined state and secondary variables.

        Function arguments:

        Function returns:
        """
        if mean == 'None':
            mean = self._mean
        reactionMechanism = self._param.Case.Mixture.reactionMechanism
        self._reactionMechanismType = reactionMechanism['type']
        self._reactions = []
        for i_reaction, reaction in enumerate(reactionMechanism['reactions']):
            reaction = reactionClass(
                                     reaction, 
                                     self._reactionMechanismType,
                                     mean,
                                     self,
                                    )
            self._reactions.append(reaction)
            self._fieldDict['RR'+str(i_reaction)] = reaction.RR

    def omega(
            self,
            specie,
            ):
        omegaTotal = 'not initialized'
        for i_reaction, reaction in enumerate(self._reactions):
            if specie in reaction._educts:
                index = reaction._educts.index(specie)
                stoch = reaction._stochiometricCoefficientsEducts[index]
                omegaReaction = - stoch * self.RR(index)
                if omegaTotal == 'not initialized':
                    omegaTotal = omegaReaction
                else:
                    omegaTotal += omegaReaction
            elif specie in reaction._products:
                index = reaction._products.index(specie)
                stoch = reaction._stochiometricCoefficientsProducts[index]
                omegaReaction = stoch * self.RR(index)
                if omegaTotal == 'not initialized':
                    omegaTotal = omegaReaction
                else:
                    omegaTotal += omegaReaction
                input(omegaTotal.ufl_tens)
        return omegaTotal
        
    def RR(
            self,
            index,
          ):
        return self._fieldDict['RR'+str(index)]
        
