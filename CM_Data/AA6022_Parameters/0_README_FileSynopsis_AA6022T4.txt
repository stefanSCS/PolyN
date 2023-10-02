

This file maps Figures to models (except for the PolyN clones of Hill'48 (see below)).
The parameters of a model are in the specified file. 



-----------------------------------Main ARTICLE--------------------------------------------------------------

Fig.10(dashed):    AA6022T4_Yld2000_Tian_deg8_Err_and_Coeff.txt  
Fig.10(solid):   AA6022T4_R_deg8_Err_and_Coeff.txt 



-----------------------------------APPENDIX--------------------------------------------------------------


Fig.4:   AA6022T4_Yld2000_Tian_deg8_Err_and_Coeff.txt  

Fig.5:     AA6022T4_H_deg4_Err_and_Coeff.txt  

Fig.6:   AA6022T4_H_deg6_Err_and_Coeff.txt 

Fig.7:     AA6022T4_H_ExptP6_deg8_Err_and_Coeff.txt  

Fig.8:  AA6022T4_H_ExptP6_ExptP8_deg10_Err_and_Coeff.txt 

Fig.9:     AA6022T4_H_ExptP6_ExptP8_ExptP10_deg12_Err_and_Coeff.txt 


   
-------------------------------PolyN clones of Hill'48--------------------------------------------------------

Hill'48 function: 
GH*sx*sx+FH*sy*sy-2.0*H*sx*sy+N2*sxy*sxy

Hill'48 model of AA6022-T4 constructed using: 
s0=1.0, r0=0.8,  r45=0.37,  r90=0.54

Hill'48 parameters for AA6022-T4:
GH,FH,H,N2 =  1.0, 1.2674897119341564, 0.4444444444444445, 2.3987654320987657  


Note: Below are FE input files (as implemented by the fast evaluation algorithm).
The actual parameters can be extracted as follow: 
--The 7-th row stores the number M of parameters (for example, 9, in case of Poly4)
--Then the PolyN parameters are stored in rows 8 to 8+M


Poly4:  symb_AA6022T4_HillPoly_deg4_FEdata.txt 
Poly6:  symb_AA6022T4_HillPoly_deg6_FEdata.txt
Poly8:  symb_AA6022T4_HillPoly_deg8_FEdata.txt
Poly10:  symb_AA6022T4_HillPoly_deg10_FEdata.txt  
Poly12:  symb_AA6022T4_HillPoly_deg12_FEdata.txt









