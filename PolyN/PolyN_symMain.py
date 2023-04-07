

import PolyN_V3F as polyn
import PolyN_symb as spolyn




####-------------INPUT----------------------------------
###Place input data here 

if(0):## example files for Yld89
    fName='symb_Yld89.txt'
    subDir='AA6022T4' 
    
if(True):## example files for the Yld2000_2D model 
    #fName,subDir='sym_Yld2002_Tian.txt','AA6022T4'
    fName='symb_AA6016T4_Yld2000_Postech.txt'
    subDir='AA6016T4' 
    
if(0):### files for Yld2004-18p
    fName,subDir='symb_AA6016T4_Yld2004_Postech.txt','AA6016T4' 


if(0):##file for BBC2005 (parameters and input experimental data shared by Prof. D. Banabic)
    fName,subDir='symb_AA6016T4_BBC2005.txt','AA6016T4'

if(0):##file for  FACET (model reported in Habraken et al)
    fName,subDir='symb_AA6016T4_FACET.txt','AA6016T4'

if(0):##file for  Caz2018-Ort (model reported in Habraken et al)
    fName,subDir='symb_AA6016T4_Caz.txt','AA6016T4'

####----------END of INPUT----------------------------------



##Execute task 
ddata=spolyn.readInputFile(fName,subDir,errReport=True,echo=True)
spolyn.yldFunc_To_PolyN(ddata)



