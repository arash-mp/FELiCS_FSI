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
import numpy as np

# Third party libraries

class ReactionClass:
    """
    Represents a chemical reaction and computes associated reaction rates and source terms.

    The class supports the KaiserCnF2023 reaction mechanism and uses mean and fluctuating fields
    to evaluate the reaction rate. It also determines the reaction formula and provides access
    to species involved in the reaction and source terms.

    **Initialize the reactionClass object**

    Parameters
    ----------
    reactionData : dict
        A dictionary containing the educts, products, and their stoichiometric coefficients.
    reactionType : str
        The type of the reaction mechanism, e.g., 'KaiserCnF2023'.
    mean : object
        The mean flow field object.
    fluc : object
        The fluctuation field object.

    Attributes
    ----------
    _educts : list
        List of educt species in the reaction.
    _stochiometricCoefficientsEducts : list
        Stoichiometric coefficients corresponding to the educts.
    _products : list
        List of product species in the reaction.
    _stochiometricCoefficientsProducts : list
        Stoichiometric coefficients corresponding to the products.
    _reactionType : str
        Type of the reaction mechanism.
    _RR : float
        Computed reaction rate based on mean and fluctuating fields.
    """

    def __init__(
                 self,
    reactionData,
    reactionType,
    mean,
    fluc,
    ):
        """
        Initializes the reactionClass instance.

        Parameters
        ----------
        reactionData : dict
            A dictionary containing educts, products, and their stoichiometric coefficients.
        reactionType : str
            The type of the reaction mechanism.
        mean : object
            The mean flow object.
        fluc : object
            The fluctuation object.
        """

        self._educts = reactionData['educts']
        self._stochiometricCoefficientsEducts = reactionData['stochiometricCoefficientsEducts']
        self._products = reactionData['products']
        self._stochiometricCoefficientsProducts = reactionData['stochiometricCoefficientsProducts']
        self._reactionType = reactionType
        if reactionType == 'KaiserCnF2023':
            specie = self._products[0]
            if fluc._isSolution:
                mean_RR_prefactor = mean.field_dict['RR_prefactor']
                # The following if clause is a quick fix. Should be removed once the mean flow class has been changed to Tensor notation
                mean_Yspecie = mean.field_dict[specie]
                if 'rho' in list(mean.field_dict.keys()):
                    mean_rho = mean.field_dict['rho']
                else:
                    #mean_rho = mean.oneField
                    mean_rho = mean.rho
            # else it is a tensor object and will be used to build a ufl formulation
            else:
                mean_RR_prefactor = mean.rr_prefactor
                mean_rho = mean.rho
                mean_Yspecie = mean.y(specie)
       
            self._RR = mean_RR_prefactor * mean_rho * (fluc.y(specie) - 2 * fluc.y(specie) * mean_Yspecie)\
                     + mean_RR_prefactor * fluc.rho * (mean_Yspecie - mean_Yspecie * mean_Yspecie)
            # self._RR = 860 * mean_rho * (fluc.Y(specie) - 2 * fluc.Y(specie) * mean_Yspecie)\
            #          + 860 * fluc.rho * (mean_Yspecie - mean_Yspecie * mean_Yspecie)
        else:
            raise Exception('Reaction mechanism ' + reaction['type'] + ' unknown!')
        self._determine_reaction_formula()

    def _determine_reaction_formula(self):
        """
        Determines the reaction formula from educts and products.
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
        
    def rr(self):
        """
        Returns the computed reaction rate.

        Returns:
        -------
        float
            The reaction rate.
        """
        return self._RR

    def species(self):
        """
        Returns the list of species involved in the reaction.

        Returns:
        -------
        list
            Combined list of educt and product species.
        """

        return self._educts + self._products
        
    def source_term(self,
    specie,
    ):
        """
        Computes the source term for a specific species based on its role in the reaction.

        Parameters
        ----------
        specie : str
            The species for which the source term is calculated.

        Returns
        -------
        float
            The source term corresponding to the input species.

        Raises
        ------
        Exception
            If the species is not part of the reaction.
        """

        if specie in self._educts:
            index = self._educts.index(specie)
            stoch = self._stochiometricCoefficientsEducts[index]
            return - self.rr() * stoch
        elif specie in self._products:
            index = self._products.index(specie)
            stoch = self._stochiometricCoefficientsProducts[index]
            return self.rr() * stoch
        else: 
            raise Exception(specie +' not a part of reaction!') 

class ReactionHandler:
    """
    Manages and initializes chemical reactions and computes source terms for species.

    This class is responsible for setting up reaction objects, managing reaction data,
    and computing omega values for species. It is intended to be used after all other
    flow-related variables have been established.

    **Initialize the reactionHandler object**

    Parameters
    ----------
    None
    """

    def __init__(
                self,
                ):
        pass

    def _additional_fields_to_be_read_reaction(self):
        """
        Determines which additional fields are needed for the reaction mechanism.

        Returns
        -------
        list
            A list of additional field names required for the reaction mechanism.

        Raises
        ------
        Exception
            If the reaction mechanism type is unknown.
        """

        if self._param.Mixture.reaction_mechanism['type'] == 'KaiserCnF2023': 
            outList = ['RR_prefactor']
        elif self._param.Mixture.reaction_mechanism['type'] == 'None':
            outList = []
        else:
            raise Exception('Reaction mechanism ' + self._param.Mixture.reaction_mechanism['type'] + ' unknown!')
        return outList
    
    def _initialize_reactions(self,
    mean='None',
    ):
        """
        Creates and initializes reactionClass objects using the defined mechanism and mean fields.

        Parameters
        ----------
        mean : object, optional
            The mean flow field object to be used for initializing reactions. Defaults to 'None'.

        Returns
        -------
        None
        """

        if mean == 'None':
            mean = self._mean
        reaction_mechanism = self._param.Mixture.reaction_mechanism
        self._reactionMechanismType = reaction_mechanism['type']
        self._reactions = []
        for i_reaction, reaction in enumerate(reaction_mechanism['reactions']):
            reaction = ReactionClass(
            reaction,
            self._reactionMechanismType,
            mean,
            self,
            )
            self._reactions.append(reaction)
            self._fieldDict['RR'+str(i_reaction)] = reaction.rr()

    def omega(self,
    specie,
    ):
        """
        Computes the net source term (omega) for a given species by aggregating contributions
        from all reactions.

        Parameters
        ----------
        specie : str
            The species for which omega is computed.

        Returns
        -------
        float
            The total omega value for the given species.
        """

        omegaTotal = 'not initialized'
        for i_reaction, reaction in enumerate(self._reactions):
            if specie in reaction._educts:
                index = reaction._educts.index(specie)
                stoch = reaction._stochiometricCoefficientsEducts[index]
                omegaReaction = - stoch * self.rr(index)
                if omegaTotal == 'not initialized':
                    omegaTotal = omegaReaction
                else:
                    omegaTotal += omegaReaction
            elif specie in reaction._products:
                index = reaction._products.index(specie)
                stoch = reaction._stochiometricCoefficientsProducts[index]
                omegaReaction = stoch * self.rr(index)
                if omegaTotal == 'not initialized':
                    omegaTotal = omegaReaction
                else:
                    omegaTotal += omegaReaction
                input(omegaTotal.ufl_tens)
        return omegaTotal
        
    def rr(self,
    index,
    ):
        """
        Retrieves the reaction rate for a specified reaction.

        Parameters
        ----------
        index : int
            Index of the reaction.

        Returns
        -------
        float
            The reaction rate for the specified reaction.
        """
        
        return self._fieldDict['RR'+str(index)]
