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
# Third party libraries
from colorama import Fore,Style

# Local Libraries and methods
from FELiCS.Misc.functions import get_last_git_commit

print(Fore.RED+"(     "+Style.RESET_ALL+"    "+Fore.RED+"(     "+Style.RESET_ALL+"          "+Fore.RED+"(     "+Style.RESET_ALL)
print(Fore.RED+")"+Fore.YELLOW+"\\ "+Fore.RED+")      "+Fore.RED+")"+Fore.YELLOW+" "+Fore.RED+")       "+Fore.RED+"(    "+Fore.RED+")"+Fore.YELLOW+"\\ "+Fore.RED+")  "+Style.RESET_ALL)
print(Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+")"+Fore.YELLOW+"/"+Fore.RED+"(  (  "+Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+")"+Fore.YELLOW+"/"+Fore.RED+"( (    "+Fore.RED+")\\   "+Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+")"+Fore.YELLOW+"/"+Fore.RED+"(  "+Style.RESET_ALL)
print(Fore.RED+"/"+Fore.YELLOW+"("+Style.RESET_ALL+"_"+Fore.YELLOW+")"+Fore.RED+") )\\  /"+Fore.YELLOW+"("+Style.RESET_ALL+"_"+Fore.YELLOW+")"+Fore.RED+"))\\  "+Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+"("+Fore.YELLOW+"_"+Fore.RED+")  "+Fore.RED+"/"+Fore.YELLOW+"("+Style.RESET_ALL+"_"+Fore.YELLOW+")"+Fore.RED+") "+Style.RESET_ALL)
print(Fore.CYAN+"("+Style.RESET_ALL+"_"+Fore.CYAN+")"+Style.RESET_ALL+"_"+Fore.CYAN+")"+Style.RESET_ALL+""+Fore.CYAN+"(("+Style.RESET_ALL+"_"+Fore.CYAN+") ("+Style.RESET_ALL+"_"+Fore.CYAN+"))"+Style.RESET_ALL+" "+Fore.CYAN+"(("+Style.RESET_ALL+"_"+Fore.CYAN+") "+Fore.RED+")"+Fore.CYAN+"\\"+Style.RESET_ALL+"___ "+Fore.CYAN+"("+Style.RESET_ALL+"_"+Fore.CYAN+"))   "+Style.RESET_ALL)
print("| __|| __|| |   (_)"+Fore.RED+"("+Fore.CYAN+"("+Style.RESET_ALL+"/ __|/ __|  "+Style.RESET_ALL)
print("| _| | _| | |__ | | | (__ \\__ \\  "+Style.RESET_ALL)
print("|_|  |___||____||_|  \\___||___/  "+Style.RESET_ALL)
#
#   from colorama import Fore,Style
#   print(Fore.RED+'ERROR!!! '+string+' Aborting Program...'+Style.RESET_ALL)   
#                 "+Fore.YELLOW+"
#         "+Fore.RED+"                
#         "+Fore.CYAN+"
#         "+Style.RESET_ALL+"

label = get_last_git_commit()
# NOTE: The logger is not initialized at this point
print('Git commit: '+str(label))

