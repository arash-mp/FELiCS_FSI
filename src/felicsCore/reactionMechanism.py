class reactionClass:
    def __init__(
            self,
            educts, 
            products,
            reactionMechanism, 
            fluc, 
            mean,
            ):
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
            self.discretizeReaction()

    def consumption(
            self, 
            specie
            ):
        if specie in educts:
            return self.__rr
        else: 
            printError('Specie ' + specie + ' not a educt of reaction ' + self.__name + '.')

    def discretizeReaction(
            self,
            mean,
            fluc,
            ):
        self.__isDisretized == True
        if reactionMechanism in ['EBU_CnF_Kaiser2023']:
            self.__rr = mean.RR_prefactor * mean.rho * (fluc.Y('progress') - 2 * fluc.Y('progress') * mean.Y('progress'))\
                      + mean.RR_prefactor * fluc.rho * (mean.Y('progress') -     mean.Y('progress') * mean.Y('progress'))
        else:
            printWarning('Reaction ' + self.__name + ' could not be discretized.') 
            self.__isDisretized == False

    @property
    def educts(self):
        return self.__educts

    def production(
            self, 
            specie,
            ):
        if specie in products:
            return self.__rr
        else: 
            printError('Specie ' + specie + 'not a product of reaction ' + self.__name)

    @property
    def products(self):
        return self.__products

    @property
    def reactionEquation(self):
        return self.__name

    @property
    def rr(self):
        return self.__rr


class reactionMechanismClass:
    def __init__(
            self,reactionMechanism, 
            fluc = None, 
            mean = None,
            ):
        self.__numberOfSpecies = 0
        self.__numberOfReactions = 0
        self.__reactionList = []
        self.__additionalMeanFieldQuantities = []
        self.__reactionMechanism = reactionMechanism
        if self.__reactionMechanism in ['EBU_CnF_Kaiser2023']:            
            self.__numberOfSpecies = 1
            self.__numberOfReactions = 1
            self.__reactionList.append(reactionClass('','progress',reactionMechanism, fluc, mean))
            self.__additionalMeanFieldQuantities.append('rr_prefactor')
            self.__additionalMeanFieldQuantities.append('T_u') # unburnt temperature
            self.__additionalMeanFieldQuantities.append('T_b') # burnt temperature


    @property
    def numberOfReactions(self):
        return self._numberOfReactions

    @property
    def numberOfSpecies(self):
        return self._numberOfSpecies

    @property
    def reactionMechanism(self):
        return self.__reactionMechanism

    @property
    def reactions(self):
        return self.__reactionList

    @property
    def additionalMeanFieldQuantities(self):
        return self.__additionalMeanFieldQuantities

