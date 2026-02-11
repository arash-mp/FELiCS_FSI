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
class ReactionClass:
    """
    Represents a chemical reaction with educts, products, and mechanism.

    This class manages the definition, discretization, and evaluation of a chemical reaction, including educt and product species, and provides access to reaction rates and equations.

    **Initialize the reactionClass object**

    Parameters
    ----------
    educts : list or str
        List of educt species.
    products : list or str
        List of product species.
    reactionMechanism : str
        Name of the reaction mechanism.
    fluc : object
        Fluctuation field (optional).
    mean : object
        Mean field (optional).

    Attributes
    ----------
    __educts : list or str
        Educt species.
    __products : list or str
        Product species.
    __isDiscretized : bool
        Flag indicating if the reaction is discretized.
    __name : str
        Reaction equation string.
    __rr : object
        Reaction rate (if discretized).
    """

    def __init__(self,
    educts,
    products,
    reaction_mechanism,
    fluc,
    mean,
    ):
        """
        Initializes the reactionClass instance.

        Parameters
        ----------
        educts : list or str
            List of educt species.
        products : list or str
            List of product species.
        reactionMechanism : str
            Name of the reaction mechanism.
        fluc : object
            Fluctuation field (optional).
        mean : object
            Mean field (optional).
        """
        educts_name = ''
        self.__educts = educts
        self.__products = products
        self.__isDiscretized = False
        for educt in educts:
            educts_name += educts + '+'
        educts = educts[:-1]
        products_name = ''
        for product in products:
            products_name += products + '+'
        products = products[:-1]
        self.__name = educts + '-->' + products
        if not (fluc == None or mean == None):
            self.discretize_reaction()

    def consumption(self,
    specie,
    ):
        """
        Returns the reaction rate for a consumed educt species.

        Parameters
        ----------
        specie : str
            Name of the educt species.

        Returns
        -------
        rr : object
            Reaction rate for the educt.

        Raises
        ------
        Prints error if the species is not an educt.
        """
        if specie in educts:
            return self.__rr
        else: 
            printError('Specie ' + specie + ' not a educt of reaction ' + self.__name + '.')

    def discretize_reaction(self,
    mean,
    fluc,
    ):
        """
        Discretizes the reaction using the provided mean and fluctuation fields.

        Parameters
        ----------
        mean : object
            Mean field.
        fluc : object
            Fluctuation field.

        Notes
        -----
        Sets the reaction rate according to the specified mechanism.
        """
        self.__isDisretized == True
        if reaction_mechanism in ['EBU_CnF_Kaiser2023']:
            self.__rr = mean.rr_prefactor * mean.rho * (fluc.y('progress') - 2 * fluc.y('progress') * mean.y('progress'))\
                      + mean.rr_prefactor * fluc.rho * (mean.y('progress') -     mean.y('progress') * mean.y('progress'))
            # self.__rr = 860 * mean.rho * (fluc.Y('progress') - 2 * fluc.Y('progress') * mean.Y('progress'))\
            #           + 860 * fluc.rho * (mean.Y('progress') -     mean.Y('progress') * mean.Y('progress'))
        else:
            printWarning('Reaction ' + self.__name + ' could not be discretized.') 
            self.__isDisretized == False

    @property
    def educts(self):
        """
        Returns the educt species.

        Returns
        -------
        educts : list or str
            Educt species.
        """
        return self.__educts

    def production(self,
    specie,
    ):
        """
        Returns the reaction rate for a produced product species.

        Parameters
        ----------
        specie : str
            Name of the product species.

        Returns
        -------
        rr : object
            Reaction rate for the product.

        Raises
        ------
        Prints error if the species is not a product.
        """
        # NOTE / TODO (LUKAS) the following code will not work, products is not defined
        if specie in products:
            return self.__rr
        else: 
            printError('Specie ' + specie + 'not a product of reaction ' + self.__name)

    @property
    def products(self):
        """
        Returns the product species.

        Returns
        -------
        products : list or str
            Product species.
        """
        return self.__products

    @property
    def reaction_equation(self):
        """
        Returns the reaction equation string.

        Returns
        -------
        equation : str
            Reaction equation.
        """
        return self.__name

    @property
    def rr(self):
        """
        Returns the reaction rate.

        Returns
        -------
        rr : object
            Reaction rate.
        """
        return self.__rr


class ReactionMechanismClass:
    """
    Represents a collection of chemical reactions for a given mechanism.

    This class manages a list of reactions, species, and additional mean field quantities for a specified reaction mechanism.

    **Initialize the reactionMechanismClass object**

    Parameters
    ----------
    reactionMechanism : str
        Name of the reaction mechanism.
    fluc : object, optional
        Fluctuation field.
    mean : object, optional
        Mean field.

    Attributes
    ----------
    __numberOfSpecies : int
        Number of species in the mechanism.
    __numberOfReactions : int
        Number of reactions in the mechanism.
    __reactionList : list
        List of reactionClass instances.
    __additionalMeanFieldQuantities : list
        List of additional mean field quantities.
    __reactionMechanism : str
        Name of the reaction mechanism.
    """

    def __init__(self,
    reaction_mechanism,
    fluc=None,
    mean=None,
    ):
        """
        Initializes the reactionMechanismClass instance.

        Parameters
        ----------
        reactionMechanism : str
            Name of the reaction mechanism.
        fluc : object, optional
            Fluctuation field.
        mean : object, optional
            Mean field.
        """
        self.__numberOfSpecies = 0
        self.__numberOfReactions = 0
        self.__reactionList = []
        self.__additionalMeanFieldQuantities = []
        self.__reactionMechanism = reaction_mechanism
        if self.__reactionMechanism in ['EBU_CnF_Kaiser2023']:            
            self.__numberOfSpecies = 1
            self.__numberOfReactions = 1
            self.__reactionList.append(ReactionClass(
            '',
            'progress',
            reaction_mechanism,
            fluc,
            mean,
            ))
            self.__additionalMeanFieldQuantities.append('rr_prefactor')
            self.__additionalMeanFieldQuantities.append('T_u') # unburnt temperature
            self.__additionalMeanFieldQuantities.append('T_b') # burnt temperature


    @property
    def number_of_reactions(self):
        """
        Returns the number of reactions in the mechanism.

        Returns
        -------
        numberOfReactions : int
            Number of reactions.
        """
        return self._numberOfReactions

    @property
    def number_of_species(self):
        """
        Returns the number of species in the mechanism.

        Returns
        -------
        numberOfSpecies : int
            Number of species.
        """
        return self._numberOfSpecies

    @property
    def reaction_mechanism(self):
        """
        Returns the name of the reaction mechanism.

        Returns
        -------
        reactionMechanism : str
            Name of the reaction mechanism.
        """
        return self.__reactionMechanism

    @property
    def reactions(self):
        """
        Returns the list of reactions.

        Returns
        -------
        reactions : list
            List of reactionClass instances.
        """
        return self.__reactionList

    @property
    def additional_mean_field_quantities(self):
        """
        Returns the list of additional mean field quantities.

        Returns
        -------
        additionalMeanFieldQuantities : list
            List of additional mean field quantities.
        """
        return self.__additionalMeanFieldQuantities

