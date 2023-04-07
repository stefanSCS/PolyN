
import PolyN_V3F as polyn

fName,subDir='',''


####-------------INPUT----------------------------------
###Place input data here 

if(0):##example of fitting Poly6 to the TUAT data 
    fName='polyn_AA6016T4_TUAT_P6.txt'
    subDir='AA6016T4' 

if(0):## example of testing a Poly6 model  
    fName='polyn_AA6016T4_TUAT_P6_test.txt'
    subDir='AA6016T4' 

if(1):##example of testing convexity on a denser grid 
    fName='polyn_AA6016T4_TUAT_P6_cvxCheck.txt'
    subDir='AA6016T4' 

####----------END of INPUT----------------------------------

##Execute task 
ddata=polyn.readInputFile(fName,subDir,echo=True)
polyn.mainPolyN(ddata)

