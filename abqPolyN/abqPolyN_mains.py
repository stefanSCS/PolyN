

import abqPolyNshell as abqp

##--- input data block 
fName='abq_uniaxial_AA6022T4.txt'
subDir='AA6022T4'

fName='abq_cpDraw_AA6016T4.txt'
subDir='AA6016T4'

###---execution block
inpD=abqp.readInpData(fName,subDir)
abqp.zexec(inpD)


