from colorama import Fore,Style
print(Fore.RED+"(     "+Style.RESET_ALL+"    "+Fore.RED+"(     "+Style.RESET_ALL+"          "+Fore.RED+"(     "+Style.RESET_ALL)
print(Fore.RED+")"+Fore.YELLOW+"\ "+Fore.RED+")      "+Fore.RED+")"+Fore.YELLOW+"\ "+Fore.RED+")       "+Fore.RED+"(    "+Fore.RED+")"+Fore.YELLOW+"\ "+Fore.RED+")  "+Style.RESET_ALL)
print(Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+")"+Fore.YELLOW+"/"+Fore.RED+"(  (  "+Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+")"+Fore.YELLOW+"/"+Fore.RED+"( (    "+Fore.RED+")\   "+Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+")"+Fore.YELLOW+"/"+Fore.RED+"(  "+Style.RESET_ALL)
print(Fore.RED+"/"+Fore.YELLOW+"("+Style.RESET_ALL+"_"+Fore.YELLOW+")"+Fore.RED+") )\  /"+Fore.YELLOW+"("+Style.RESET_ALL+"_"+Fore.YELLOW+")"+Fore.RED+"))\  "+Fore.RED+"("+Fore.YELLOW+"("+Style.RESET_ALL+"("+Fore.YELLOW+"_"+Fore.RED+")  "+Fore.RED+"/"+Fore.YELLOW+"("+Style.RESET_ALL+"_"+Fore.YELLOW+")"+Fore.RED+") "+Style.RESET_ALL)
print(Fore.CYAN+"("+Style.RESET_ALL+"_"+Fore.CYAN+")"+Style.RESET_ALL+"_"+Fore.CYAN+")"+Style.RESET_ALL+""+Fore.CYAN+"(("+Style.RESET_ALL+"_"+Fore.CYAN+") ("+Style.RESET_ALL+"_"+Fore.CYAN+"))"+Style.RESET_ALL+" "+Fore.CYAN+"(("+Style.RESET_ALL+"_"+Fore.CYAN+") "+Fore.RED+")"+Fore.CYAN+"\\"+Style.RESET_ALL+"___ "+Fore.CYAN+"("+Style.RESET_ALL+"_"+Fore.CYAN+"))   "+Style.RESET_ALL)
print("| __|| __|| |   (_)"+Fore.RED+"("+Fore.CYAN+"("+Style.RESET_ALL+"/ __|/ __|  "+Style.RESET_ALL)
print("| _| | _| | |__ | | | (__ \__ \  "+Style.RESET_ALL)
print("|_|  |___||____||_|  \___||___/  "+Style.RESET_ALL)
#
#	from colorama import Fore,Style
#	print(Fore.RED+'ERROR!!! '+string+' Aborting Program...'+Style.RESET_ALL)	
#                 "+Fore.YELLOW+"
#		  "+Fore.RED+"                
#		  "+Fore.CYAN+"
#		  "+Style.RESET_ALL+"
from FELiCS.functions import getLastGitCommit
label=getLastGitCommit()
print('Git commit: '+str(label))

