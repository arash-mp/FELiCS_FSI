import numpy as np

# Third party libraries

class reactionClass:
    """
    Represents a reaction in the system.

    Parameters:
    ----------
    reactionData : dict
        A dictionary containing reaction data.
    reactionType : str
        The type of reaction.
    mean : object
        The mean flow object.
    fluc : object
        The fluctuation object.

    Attributes:
    ----------
    _educts : list
        List of educts in the reaction.
    _stochiometricCoefficientsEducts : list
        List of stoichiometric coefficients for the educts.
    _products : list
        List of products in the reaction.
    _stochiometricCoefficientsProducts : list
        List of stoichiometric coefficients for the products.
    _reactionType : str
        The type of reaction.
    _RR : float
        The reaction rate.

    Methods:
    -------
    _determineReactionFormula()
        Determines the reaction formula.
    RR()
        Returns the reaction rate.
    species()
        Returns the species involved in the reaction.
    sourceTerm(specie)
        Returns the source term for a given species.
    """

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
        self._determineReactionFormula()

    def _determineReactionFormula(self):
        """
        Determines the reaction formula.
        """
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
        
    def RR(self):
        """
        Returns the reaction rate.

        Returns:
        -------
        float
            The reaction rate.
        """
        return self._RR

    def species(self):
        """
        Returns the species involved in the reaction.

        Returns:
        -------
        list
            List of species involved in the reaction.
        """
        return self._educts + self._products
        
    def sourceTerm(self, specie):
        """
        Returns the source term for a given species.

        Parameters:
        ----------
        specie : str
            The species for which the source term is calculated.

        Returns:
        -------
        float
            The source term for the given species.
        """
        if specie in self._educts:
            index = self._educts.index(specie)
            stoch = self._stochiometricCoefficientsEducts[index]
            return - self.RR() * stoch
        elif specie in self._products:
            index = self._products.index(specie)
            stoch = self._stochiometricCoefficientsProducts[index]
            return self.RR() * stoch
        else: 
            raise Exception(specie +' not a part of reaction!') 

class reactionHandler:
    """
    This class is used to manage the variables for the linear reactions. 
    It is applied last when all other variables are already calculated.

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
        """
        Returns a list of additional fields to be read for the reaction.

        Returns:
        -------
        list
            List of additional fields to be read for the reaction.
        """
        if self._param.Mixture.reactionMechanism['type'] == 'KaiserCnF2023': 
            outList = ['RR_prefactor']
        elif self._param.Mixture.reactionMechanism['type'] == 'None':
            outList = []
        else:
            raise Exception('Reaction mechanism ' + self._param.Mixture.reactionMechanism['type'] + ' unknown!')
        return outList
    
    def _initializeReactions(self, mean='None'):
        """
        Initializes the reactions and links them to the previously determined state and secondary variables.

        Parameters:
        ----------
        mean : object, optional
            The mean flow object. Default is 'None'.

        Returns:
        -------
        None
        """
        if mean == 'None':
            mean = self._mean
        reactionMechanism = self._param.Mixture.reactionMechanism
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
            self._fieldDict['RR'+str(i_reaction)] = reaction.RR()

    def omega(self, specie):
        """
        Calculates the omega value for a given species.

        Parameters:
        ----------
        specie : str
            The species for which the omega value is calculated.

        Returns:
        -------
        float
            The omega value for the given species.
        """
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
        
    def RR(self, index):
        """
        Returns the reaction rate for a given index.

        Parameters:
        ----------
        index : int
            The index of the reaction.

        Returns:
        -------
        float
            The reaction rate for the given index.
        """
        return self._fieldDict['RR'+str(index)]
