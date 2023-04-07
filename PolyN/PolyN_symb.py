

##Developed and released by SCS 

import numpy as np
##from matplotlib import pyplot as plt
import sympy as sp
import PolyN_V3F as polyn
sp.init_printing(use_unicode=True,wrap_line=False)



def readInputFile(fName,subDir,errReport=False,echo=False):
    if(polyn.osn=='nt'):drsep='\\'
    else:drsep='/'
    try:
        ff=open(polyn.figDirData+drsep+subDir+drsep+fName)
    except IOError as err:
        print(err);exit()
    ddata={'func':'','degree':0,'alpha_1':[],'alpha_2':[],'alpha':0,
    'fileInpData':'','subDir':subDir,'fileFACET':'','fileFACETsep':' ','errReport':errReport}
    msg='Calculations aborted'   
    for line in ff:
        line=line.strip()
        if(line=='' or line[0]=='#'):continue
        if(':' in line):line=line.split(':');lnz=line[0].strip()
        else:print(fName+': Unknown format on line');print(line);print(msg);exit()
        if(lnz=='func'):
            tval=line[1].strip()
            if(tval in ['Yld89', 'Yld2000_2D', 'BBC2005', 'Yld2004_18p', 'FACET', 'Caz2018_Ort']):
                ddata['func']=tval
            else:print('func: Unknown option');print(msg);exit()
        if(lnz=='degree'):
            try:tval=int(float(line[1]))
            except ValueError as err:print(err);print(msg);exit()
            if(tval%2):print('degree: must be even');print(msg);exit()
            ddata['degree']=tval
        if(lnz=='fileInpData'):ddata['fileInpData']=line[1].strip()
        if(lnz=='fileFACET'):ddata['fileFACET']=line[1].strip()
        if(lnz=='fileFACETsep'):
            sep=line[1].strip()
            if(sep in ['|',',']):ddata['fileFACETsep']=sep
            elif(sep=='b'):ddata['fileFACETsep']=' '            
            else:print('Unknown file separator');print(msg);exit()
        if(lnz=='alpha_1'):
            tval=line[1].strip().split(',')
            for kval in tval:
                try:ddata['alpha_1'].append(float(kval))
                except ValueError as err:print(err);print(msg);exit()
        if(lnz=='alpha_2'):
            tval=line[1].strip().split(',')
            for kval in tval:
                try:ddata['alpha_2'].append(float(kval))
                except ValueError as err:print(err);print(msg);exit() 
        if(lnz=='alpha'):
            try:ddata['alpha']=float(line[1].strip())
            except ValueError as err:print(err);print(msg);exit()            
    if(echo):
        for key in ddata:print(key,': ',ddata[key])
    return ddata    
    
        
def paramHill(rVals):
    r0,r45,r90=rVals
    G=1/(1+r0)
    H=r0/(1+r0)
    F=H/r90
    N2=(2*r45+1)*(F+G)
    GH,FH=G+H,F+H
    return GH,FH,H,N2

def fHill(vsx,vsy,vsxy,GH,FH,H,N2):
    vf=np.sqrt(GH*vsx*vsx+FH*vsy*vsy-2.0*H*vsx*vsy+N2*vsxy*vsxy)
    vg1=(GH*vsx-H*vsy)/vf
    vg2=(FH*vsy-H*vsx)/vf
    vg3=(N2*vsxy)/vf
    vh11=(GH-vg1*vg1)/vf
    vh22=(FH-vg2*vg2)/vf
    vh12=(-H-vg1*vg2)/vf
    vh13=(-vg1*vg3)/vf
    vh23=(-vg2*vg3)/vf
    vh33=(N2-vg3*vg3)/vf
    return vf,[vg1,vg2,vg3],[vh11,vh12,vh13,vh22,vh23,vh33]


def sHill48(degree,GH,FH,H,N2):
    if(degree%2):print('degree must be even');exit()
    deg2=int(degree/2)
    ###GH,FH,H,N2=paramHill(rVals)
    sx,sy,txy=sp.symbols("sx,sy,txy")
    vHill=GH*sx*sx+FH*sy*sy-2.0*H*sx*sy+N2*txy*txy
    vQ=[]
    jj=0
    while(jj<=degree):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degree+1-jj)]
        jj+=2  
    #print(vQ);exit()
    PP=vHill**deg2
    PP=sp.poly(PP,sx,sy,txy)
    kk,vCoeff=0,[]
    #for mon in PP.monoms():
    #    pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
    #    cff=sp.simplify(PP.coeff_monomial(pMon))
    #    print(kk,mon,cff);kk+=1
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(cff)
        ##print(kk,mon,cff);kk+=1
    #print(vCoeff)
    #print('DONE')
    return vCoeff


def plotHill48(GH,FH,H,N2,degree,zvcf,astheta,vsexp,artheta,vrexp,name,figSave=False):
    rad,strDeg=180/np.pi,str(degree)
    vtheta=np.linspace(0,np.pi/2,200)
    ct,st=np.cos(vtheta),np.sin(vtheta)
    vsxy,vsx,vsy=ct*st,ct*ct,st*st
    aVtheta=rad*vtheta
    vstheta,[vdx,vdy,vdxy],[vh11,vh12,vh13,vh22,vh23,vh33]=fHill(vsx,vsy,vsxy,GH,FH,H,N2)
    vstheta=1.0/vstheta
    vrtheta=(vdxy*vsxy-(vdx*vsy+vdy*vsx))/(vdx+vdy)
    nbp=vtheta.shape[0]
    zstheta,zrtheta,zdx,zdy,zdxy=np.zeros(nbp),np.zeros(nbp),np.zeros(nbp),np.zeros(nbp),np.zeros(nbp)
    for k in range(nbp):
        zstheta[k],[zdx[k],zdy[k],zdxy[k]]=polyn.fGYF(vsx[k],vsy[k],vsxy[k],zvcf)
    zstheta=1.0/zstheta
    zrtheta[:]=(zdxy*vsxy-(zdx*vsy+zdy*vsx))/(zdx+zdy)  
    fg1=polyn.plt.figure()
    ax1=fg1.add_subplot()
    ax1.plot(aVtheta,zstheta,linewidth=4,linestyle='-',color='r',label='Poly{s}'.format(s=strDeg))
    ax1.plot(aVtheta,vstheta,linewidth=2,linestyle='dashed',color='k',label='Hill48 fitted')
    ##vsexp=np.array(mMat['sExp'])
    xtk=[0,15,30,45,60,75,90]
    #lxtk=[str(x) for x in xtk]
    ax1.plot(astheta,vsexp,linestyle='',marker='o',markersize=9,markeredgewidth=2,
    color='g',markerfacecolor=(1,1,1,0),markeredgecolor='g',label='Exp data')
    ax1.set_xticks(xtk) ##;ax1.set_xticklabels(lxtk,fontsize=12)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ##ax1.tick_params(axis='both', which='minor', labelsize=8)
    ylim=ax1.get_ylim()
    ax1.text(0.1,0.5*(ylim[0]+ylim[1]),r"$\overline{\sigma}_{\theta}$",fontsize=16)
    ax1.text(88,ylim[0]+0.01*(ylim[1]-ylim[0]),r"$\theta$",fontsize=15)
    ax1.set_xlim(-1.25,91.25)
    ax1.grid()
    ax1.legend(fontsize=14)     
    fg2=polyn.plt.figure()
    ax2=fg2.add_subplot()
    ##ax2.plot(aVtheta,vrtheta,linewidth=1,linestyle='dashed',color='k',label='Hill fitted r-value')
    ax2.plot(aVtheta,zrtheta,linewidth=4,linestyle='-',color='r',label='Poly{s}'.format(s=strDeg))
    ax2.plot(aVtheta,vrtheta,linewidth=2,linestyle='dashed',color='k',label='Hill48 fitted')
    ax2.plot(artheta,vrexp,linestyle='',marker='o',markersize=9,markeredgewidth=2,
    color='g',markerfacecolor=(1,1,1,0),markeredgecolor='g',label='Exp data')
    ##vrexp=np.array(mMat['rValsExp2'])
    ax2.set_xticks(xtk)
    ylim=ax2.get_ylim()
    ##ax2.text(0.1,vrexp[0]-0.1*(ylim[1]-ylim[0]),r"$r_{\theta}$",fontsize=15)
    ax2.text(0.1,0.52*(ylim[1]+ylim[0]),r"$r_{\theta}$",fontsize=16)
    ax2.text(88,ylim[0]+0.01*(ylim[1]-ylim[0]),r"$\theta$",fontsize=15)
    ax2.set_xlim(-1.25,91.25)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.grid()
    ax2.legend(fontsize=14)
    if(0):
        fg3=plt.figure()
        ax3=fg3.add_subplot()
        ax3.plot(aVtheta,zHH[:,0],linewidth=3,linestyle='-',color='r',label='$H_{ij}$'+'(Poly_{s}-Hill'.format(s=strDeg))
        ax3.plot(aVtheta,vh11,linewidth=2,linestyle='dashed',color='k',label='$H_{ij}$ (Hill)')
        ax3.plot(aVtheta,zHH[:,1],linewidth=3,linestyle='-',color='r')
        ax3.plot(aVtheta,vh12,linewidth=2,linestyle='dashed',color='k')
        ax3.plot(aVtheta,zHH[:,2],linewidth=3,linestyle='-',color='r')
        ax3.plot(aVtheta,vh13,linewidth=2,linestyle='dashed',color='k')
        ax3.plot(aVtheta,zHH[:,3],linewidth=3,linestyle='-',color='r')
        ax3.plot(aVtheta,vh22,linewidth=2,linestyle='dashed',color='k')
        ax3.plot(aVtheta,zHH[:,4],linewidth=3,linestyle='-',color='r')
        ax3.plot(aVtheta,vh23,linewidth=2,linestyle='dashed',color='k')
        ax3.plot(aVtheta,zHH[:,5],linewidth=3,linestyle='-',color='r')
        ax3.plot(aVtheta,vh33,linewidth=2,linestyle='dashed',color='k')
        ylim=ax3.get_ylim()
        ax3.text(88,ylim[0]+0.01*(ylim[1]-ylim[0]),r"$\theta$",fontsize=14)
        ax3.set_xlim(-1.25,91.25)
        ax3.set_xticks([k*15 for k in range(0,7)])
        ax3.legend()
    #figSave=False
    if(figSave):
        ffz='./';ffz=polyn.figDirPlot
        fg1.savefig(ffz+name+'_Poly_{s}_Hill48_Stheta.png'.format(s=strDeg),dpi=300,bbox_inches='tight')
        fg2.savefig(ffz+name+'_Poly_{s}_Hill48_Rtheta.png'.format(s=strDeg),dpi=300,bbox_inches='tight')
        ##fg3.savefig(ffz+'Poly_{s}_Hill48_HTheta.png'.format(s=strDeg),dpi=300,bbox_inches='tight')
    polyn.plt.show() 


def binomialCoeff(degree):
    if(type(degree)!=int):print('binomialCoeff: degree must be an integer');exit()
    if(degree%2):print('binomialCoeff: degree must be an even integer');exit()
    deg2=int(degree/2)
    vBin=[0 for k in range(deg2+1)]
    vBin[0],vBin[-1]=1.0,1.0
    for k in range(1,deg2+1):
        if(k%2):continue
        nom,denom=1,1
        for j in range(1,k+1):
            denom*=j
            nom*=degree-k+j
        kk=int(k/2)    ##;print(kk)
        vBin[kk]=nom/denom
        vBin[deg2-kk]=vBin[kk]
    return vBin
    
def binomialCoeff2(degree):
    if(type(degree)!=int):print('binomialCoeff: degree must be an integer');exit()
    if(degree%2):print('binomialCoeff: degree must be an even integer');exit()
    deg2=int(degree/2)
    vBin=[0 for k in range(deg2+1)]
    vBin[0]=1.0
    for k in range(1,deg2+1):
        nom,denom=1,1
        for j in range(1,k+1):
            denom*=j
            nom*=degree-k+j
        vBin[k]=nom/denom
    return vBin    

    
def Yld2004_To_Poly8(zAlpha_1,zAlpha_2):
    sx,sy,txy=sp.symbols("sx,sy,txy")
    deg=8;degOne=deg+1;deg2=int(deg/2)
    vW=[1.0, -8.0, 28.0, -56.0, 70.0]
    vQ,jj=[],0
    while(jj<=deg):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(deg-jj-k,k,jj) for k in range(degOne-jj)]
        jj+=2 
    vAlpha_1=[0,zAlpha_1[0]+zAlpha_1[1],zAlpha_1[1]-2*zAlpha_1[0],
    zAlpha_1[3]-2*zAlpha_1[2],zAlpha_1[2]+zAlpha_1[3],
    zAlpha_1[5]-2*zAlpha_1[4],zAlpha_1[4]-2*zAlpha_1[5],zAlpha_1[6]]  
    vAlpha_1=[x/3.0 for x in vAlpha_1[0:-1]]+[vAlpha_1[-1]]    
    SX=vAlpha_1[1]*sx+vAlpha_1[2]*sy
    SY=vAlpha_1[3]*sx+vAlpha_1[4]*sy
    SZ=vAlpha_1[5]*sx+vAlpha_1[6]*sy
    TXY=vAlpha_1[7]*txy
    L11,L12,SZ_1 = SX+SY, SX*SY-TXY*TXY, SZ
    vQ_1=[0.0 for k in range(degOne)]
    vP_1=[0.0 for k in range(degOne)]
    vP_1[0],vP_1[1]=2,L11
    vQ_1[0],vQ_1[1]=vP_1[0]+1,vP_1[1]+SZ_1
    vAlpha_2=[0,zAlpha_2[0]+zAlpha_2[1],zAlpha_2[1]-2*zAlpha_2[0],
    zAlpha_2[3]-2*zAlpha_2[2],zAlpha_2[2]+zAlpha_2[3],
    zAlpha_2[5]-2*zAlpha_2[4],zAlpha_2[4]-2*zAlpha_2[5],zAlpha_2[6]]
    vAlpha_2=[x/3.0 for x in vAlpha_2[0:-1]]+[vAlpha_2[-1]]     
    SX=vAlpha_2[1]*sx+vAlpha_2[2]*sy
    SY=vAlpha_2[3]*sx+vAlpha_2[4]*sy
    SZ=vAlpha_2[5]*sx+vAlpha_2[6]*sy
    TXY=vAlpha_2[7]*txy
    L21,L22,SZ_2 = SX+SY, SX*SY-TXY*TXY, SZ
    vQ_2=[0.0 for k in range(degOne)]
    vP_2=[0.0 for k in range(degOne)]
    vP_2[0],vP_2[1]=2,L21
    vQ_2[0],vQ_2[1]=vP_2[0]+1,vP_2[1]+SZ_2
    for k in range(2,degOne):
        vP_1[k]=L11*vP_1[k-1]-L12*vP_1[k-2]
        #SZ_1*=SZ_1
        vQ_1[k]=vP_1[k]+SZ_1**k
        vP_2[k]=L21*vP_2[k-1]-L22*vP_2[k-2]
        #SZ_2*=SZ_2
        vQ_2[k]=vP_2[k]+SZ_2**k
    ff=vW[deg2]*vQ_1[deg2]*vQ_2[deg2]
    for k in range(deg2):
        #print(k)
        ff+=vW[k]*(vQ_1[k]*vQ_2[deg-k]+vQ_1[deg-k]*vQ_2[k])
    ff=sp.expand(ff)    
    PP=sp.poly(ff,sx,sy,txy)
    #PP=sp.poly(ff,sx)
    kk,vCoeff=0,[]
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(cff/4)
        #print(kk,mon,cff);kk+=1
    return vCoeff 


def Yld2004_To_PolyN(zAlpha_1,zAlpha_2,degree):
    sx,sy,txy=sp.symbols("sx,sy,txy")
    deg=degree;degOne=deg+1;deg2=int(deg/2)
    vW=binomialCoeff2(deg)
    vW=[((-1)**k)*vW[k] for k in range(0,len(vW))]##;print(vW);exit()
    vQ,jj=[],0
    while(jj<=deg):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(deg-jj-k,k,jj) for k in range(degOne-jj)]
        jj+=2 
    vAlpha_1=[0,zAlpha_1[0]+zAlpha_1[1],zAlpha_1[1]-2*zAlpha_1[0],
    zAlpha_1[3]-2*zAlpha_1[2],zAlpha_1[2]+zAlpha_1[3],
    zAlpha_1[5]-2*zAlpha_1[4],zAlpha_1[4]-2*zAlpha_1[5],zAlpha_1[6]]  
    vAlpha_1=[x/3.0 for x in vAlpha_1[0:-1]]+[vAlpha_1[-1]]    
    SX=vAlpha_1[1]*sx+vAlpha_1[2]*sy
    SY=vAlpha_1[3]*sx+vAlpha_1[4]*sy
    SZ=vAlpha_1[5]*sx+vAlpha_1[6]*sy
    TXY=vAlpha_1[7]*txy
    L11,L12,SZ_1 = SX+SY, SX*SY-TXY*TXY, SZ
    vQ_1=[0.0 for k in range(degOne)]
    vP_1=[0.0 for k in range(degOne)]
    vP_1[0],vP_1[1]=2,L11
    vQ_1[0],vQ_1[1]=vP_1[0]+1,vP_1[1]+SZ_1
    vAlpha_2=[0,zAlpha_2[0]+zAlpha_2[1],zAlpha_2[1]-2*zAlpha_2[0],
    zAlpha_2[3]-2*zAlpha_2[2],zAlpha_2[2]+zAlpha_2[3],
    zAlpha_2[5]-2*zAlpha_2[4],zAlpha_2[4]-2*zAlpha_2[5],zAlpha_2[6]]
    vAlpha_2=[x/3.0 for x in vAlpha_2[0:-1]]+[vAlpha_2[-1]]     
    SX=vAlpha_2[1]*sx+vAlpha_2[2]*sy
    SY=vAlpha_2[3]*sx+vAlpha_2[4]*sy
    SZ=vAlpha_2[5]*sx+vAlpha_2[6]*sy
    TXY=vAlpha_2[7]*txy
    L21,L22,SZ_2 = SX+SY, SX*SY-TXY*TXY, SZ
    vQ_2=[0.0 for k in range(degOne)]
    vP_2=[0.0 for k in range(degOne)]
    vP_2[0],vP_2[1]=2,L21
    vQ_2[0],vQ_2[1]=vP_2[0]+1,vP_2[1]+SZ_2
    S31,S32=SZ_1,SZ_2
    for k in range(2,degOne):
        vP_1[k]=L11*vP_1[k-1]-L12*vP_1[k-2]
        #S31=sp.expand(SZ_1*S31)
        vQ_1[k]=sp.expand(vP_1[k]+SZ_1**k)
        vP_2[k]=L21*vP_2[k-1]-L22*vP_2[k-2]
        #S32=sp.expand(SZ_2*S32)
        vQ_2[k]=sp.expand(vP_2[k]+SZ_2**k)
    ff=vW[deg2]*vQ_1[deg2]*vQ_2[deg2]
    for k in range(deg2):
        #print(k)
        ff+=vW[k]*(vQ_1[k]*vQ_2[deg-k]+vQ_1[deg-k]*vQ_2[k])
    ff=sp.expand(ff)    
    PP=sp.poly(ff,sx,sy,txy)
    kk,vCoeff=0,[]
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(cff/4)
        #print(kk,mon,cff);kk+=1
    return vCoeff 


    
    
def Yld2000_2D_To_PolyN(vAlpha,degree=8):
    if(type(degree)!=int):print('Yld2000_2D_To_PolyN: degree must be an integer');exit()
    if(degree%2):print('Yld2000_2D_To_PolyN: degree must be an even integer');exit()
    if(len(vAlpha)!=8):print('Yld2000_2D_To_PolyN: vAlpha must have length = 8');exit()
    print('Calculating Poly{} coeffs of Yld2000_2D...'.format(degree))
    deg2=int(degree/2)
    sx,sy,txy=sp.symbols("sx,sy,txy")
    vBin=binomialCoeff(degree)
    alpha=[1]+list(vAlpha)
    L11,L12,L21,L22,L33=2.0*alpha[1]/3.0,-alpha[1]/3.0,-alpha[2]/3.0,2.0*alpha[2]/3.0,alpha[7]
    S11=L11*sx+L12*sy
    S22=L21*sx+L22*sy
    S12=L33*txy
    Ione=S11+S22
    Itwo=S11*S22-S12*S12
    pp=((Ione*Ione-4.0*Itwo)**deg2)/2.0
    L11=2.0*(-alpha[3]+alpha[4]+4.0*alpha[5]-alpha[6])/9.0
    L12=(alpha[3]-4.0*alpha[4]-4.0*alpha[5]+4.0*alpha[6])/9.0
    L21=(4.0*alpha[3]-4.0*alpha[4]-4.0*alpha[5]+alpha[6])/9.0
    L22=2.0*(-alpha[3]+4.0*alpha[4]+alpha[5]-alpha[6])/9.0
    L33=alpha[8]
    S11=L11*sx+L12*sy
    S22=L21*sx+L22*sy
    S12=L33*txy
    Ione=(S11+S22)/2;Ione=Ione*Ione
    Itwo=S11*S22-S12*S12
    for j in range(deg2+1):
        pp+=vBin[j]*((9.0*Ione)**(deg2-j))*((Ione-Itwo)**j)
    degOne=degree+1
    vQ=[]
    jj=0
    while(jj<=degree):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degOne-jj)]
        jj+=2
    pp=sp.expand(pp)
    PP=sp.poly(pp,sx,sy,txy)
    #PP=sp.poly(ff,sx)
    kk,vCoeff=0,[]
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(cff)
        #print(kk,mon,cff);kk+=1
    return vCoeff    


def Yld89_To_PolyN(vAlpha,degree=8):
    if(type(degree)!=int):print('Yld89_To_PolyN: degree must be an integer');exit()
    if(degree%2):print('Yld89_To_PolyN: degree must be an even integer');exit()
    if(len(vAlpha)!=4):print('Yld89_To_PolyN: vAlpha must have length = 4');exit()
    if(abs(vAlpha[0]+vAlpha[1]-2.0)>1.0e-9):print('Yld89_To_PolyN:(Warning) a+c != 2')
    print('Calculating Poly{} coeffs of Yld89...'.format(degree))
    deg2=int(degree/2)
    sx,sy,txy=sp.symbols("sx,sy,txy")
    vBin=binomialCoeff(degree)
    aa,cc,hh,pp=vAlpha[0],vAlpha[1],vAlpha[2],vAlpha[3]
    pp=4.0*pp*pp
    KKone=(sx+hh*sy)*(sx+hh*sy)
    KKtwo=(sx-hh*sy)*(sx-hh*sy)+pp*txy*txy
    vpp1=(cc/2.0)*(KKtwo**deg2)
    vpp2=0
    for j in range(deg2+1):
        vpp2+=vBin[j]*(KKone**(deg2-j))*(KKtwo**j)
    vpp2=vpp1+(aa/2.0**degree)*vpp2    
    degOne=degree+1
    vQ=[]
    jj=0
    while(jj<=degree):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degOne-jj)]
        jj+=2
    vpp2=sp.expand(vpp2)    
    PP=sp.poly(vpp2,sx,sy,txy)
    #PP=sp.poly(ff,sx)
    kk,vCoeff=0,[]
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(cff)
        #print(kk,mon,cff);kk+=1
    return vCoeff     
    
    
def BBC2005_To_PolyN(vAlpha,degree=8):
    if(type(degree)!=int):print('BBC2005_To_PolyN: degree must be an integer');exit()
    if(degree%2):print('BBC2005_To_PolyN: degree must be an even integer');exit()
    if(len(vAlpha)!=8):print('BBC2005_To_PolyN: vAlpha must have length = 8');exit()
    print('Calculating Poly{} coeffs of BBC2005...'.format(degree))
    deg2=int(degree/2)
    sx,sy,txy=sp.symbols("sx,sy,txy")
    vBin=binomialCoeff(degree)
    aa,bb,L,M,N,P,Q,R=vAlpha[:]
    zLAMBDA=N*sx-P*sy;zLAMBDA=zLAMBDA*zLAMBDA+txy*txy
    zGAMMA=L*sx+M*sy;zGAMMA*=zGAMMA
    zPSI=Q*sx-R*sy;zPSI=zPSI*zPSI+txy*txy
    vpp1=0
    for j in range(deg2+1):
        vpp1+=vBin[j]*(zLAMBDA**(deg2-j))*(zGAMMA**j)
    vpp2=0
    for j in range(deg2+1):
        vpp2+=vBin[j]*(zLAMBDA**(deg2-j))*(zPSI**j)    
    vpp2=sp.expand(2*(aa*vpp1+bb*vpp2))  
    degOne=degree+1
    vQ=[]
    jj=0
    while(jj<=degree):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degOne-jj)]
        jj+=2
    PP=sp.poly(vpp2,sx,sy,txy)
    #PP=sp.poly(ff,sx)
    kk,vCoeff=0,[]
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(cff)
        #print(kk,mon,cff);kk+=1
    return vCoeff     


def Caz2018_TO_P8(vAlpha_1,vAlpha_2,alpha):
    #print(vAlpha_1);print(vAlpha_2);print(alpha)
    a1,a2,a3,a4=vAlpha_1
    b1,b2,b3,b4,b5,b10=vAlpha_2
    #print(a1,a2,a3,a4);print(b1,b2,b3,b4,b5,b10)
    sx,sy,txy=sp.symbols("sx,sy,txy")
    sx2,sy2,txy2=sx*sx,sy*sy,txy*txy
    P2=(a1/6.)*(sx-sy)**2+(a2/6.)*sy2+(a3/6.)*sx2+a4*txy2
    P3=(b1+b2)*sx2*sx/27+(b3+b4)*sy2*sy/27-b1*sx2*sy/9-b4*sx*sy2/9+((2*b10-b5)*sx+b5*sy)*txy2/3
    ff=sp.expand(P2**4-alpha*P2*P3*P3)
    #ff=P2**4-alpha*P2*P3*P3
    PP=sp.poly(ff,sx,sy,txy)
    B=(3*np.sqrt(2.0))**8/(3*(a1+a3)*(27*(a1+a3)**3-8*alpha*(b1+b2)**2))
    degree,vQ=8,[]
    jj=0;degOne=degree+1
    while(jj<=degree):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degOne-jj)]
        jj+=2
    kk,vCoeff=0,[]    
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(B*cff)
        #print(kk,mon,cff);kk+=1
    return vCoeff

def FACET_To_PolyN(inpFile,sep,degree=8,norm=True):
    vdd=[]
    try:
        ff=open(inpFile,'r')
        for line in ff:
            line=line.strip()
            if(line==''):break
            line=line.split(sep)
            vdd.append((float(line[1]),float(line[2]),float(line[3]),float(line[4])))
        ff.close()
    except (IOError,ValueError) as err:
        print(err);exit()
    #sq2=np.sqrt(2.0)
    sx,sy,txy=sp.symbols("sx,sy,txy")
    vs=((sx-sy),(sx+sy)/np.sqrt(3.0),2*txy)
    #vs=((sx-sy)/np.sqrt(2.0),(sx+sy)/np.sqrt(6.0),txy*np.sqrt(2.0))
    vpp=0
    for dd in vdd:
        vpp+=dd[0]*((dd[1]*vs[0]+dd[2]*vs[1]+dd[3]*vs[2])**degree)
    vpp=sp.expand(vpp/(2**(degree/2)))
    vQ=[]
    jj=0;degOne=degree+1
    while(jj<=degree):
        lv=len(vQ)
        ##vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degOne-jj)]
        jj+=2
    PP=sp.poly(vpp,sx,sy,txy)
    #PP=sp.poly(ff,sx)
    kk,vCoeff=0,[]
    for mon in vQ:
        pMon=sx**mon[0]*sy**mon[1]*txy**mon[2]
        cff=sp.simplify(PP.coeff_monomial(pMon))
        vCoeff.append(cff)
    if(norm):    
        vCoeff=[cc/vCoeff[0] for cc in vCoeff]    
    return vCoeff    


    

'''
Block 'Yld_To_Poly'
Calculates PolyN coefficients of the polynomial form of various yield functions 
Current options:'Yld89', 'Yld2000_2D', 'BBC2005', 'Yld2004_18p', 'FACET', 'Caz2018_Ort'
'''
def yldFunc_To_PolyN(ddata):
    yldFunc,degree=ddata['func'],ddata['degree']
    inpData,subDir=ddata['fileInpData'],ddata['subDir']
    vAlpha_1,vAlpha_2=ddata['alpha_1'],ddata['alpha_2']
    if(111):
        if(yldFunc=='Yld2000_2D'):##----- Yld2000_2D
            vCoeff=Yld2000_2D_To_PolyN(vAlpha_1,degree=degree)
            data,pdata=polyn.readData(inpData,subDir)
            vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            if(ddata['errReport']):minKG,negKG=polyn.PolyN_predictions(vcf,pdata,export=True)
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
            ##exit()
        if(yldFunc=='Yld2004_18p'):##------   Yld2004-18p
            ## test parameters from Habraken et al(2022) paper
            #if(0):##Postech (experiments)
                #vAlpha_1=[1.0, 1.0, -0.0870, 0.4496, 1.0192, 1.1724, 0.0]
                #vAlpha_2=[1.3329, 1.2510, 1.0948, 1.0855, -0.0065, 0.4812, 1.5461]
                #inpData,subDir='matCustom_AA6016T4_Yld2004_Postech.txt','AA6016T4'
            #if(0):##NTNU (experiments)
                #vAlpha_1=[1.1998, 1.2289, 0.1315, 0.8081, 1.1386, 0.2870,  1.3739]
                #vAlpha_2=[1.1510, -0.2576, 0.4873, 0.9934, 0.9110, 0.7965, 0.2679]
                #inpData,subDir='matCustom_AA6016T4_Yld2004_NTNU.txt','AA6016T4'
            #if(0):##NTNU (cp)
                #vAlpha_1=[1.1718, 1.1284, -0.0403, 0.6912, 1.1030, 0.4439, 0.9922]
                #vAlpha_2=[1.2910, -0.1740, 0.7416, 1.0787, 1.0278, 0.9367, 0.8146]   
                #inpData,subDir='matCustom_AA6016T4_Yld2004_NTNU_cp.txt','AA6016T4'            
            data,pdata=polyn.readData(inpData,subDir)
            #vCoeff=Yld2004_To_Poly8(vAlpha_1,vAlpha_2)
            vCoeff=Yld2004_To_PolyN(vAlpha_1,vAlpha_2,degree)
            vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            if(ddata['errReport']):minKG,negKG=polyn.PolyN_predictions(vcf,pdata,export=True)
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
            return        
        if(yldFunc=='Yld89'):##----- Yld89
            ##parameters from Tian_et_al(2016) paper
            #vAlpha=[1.211,0.789,1.125,0.958]
            #degree=8
            vCoeff=Yld89_To_PolyN(vAlpha,degree=degree)
            #inpData='customPlot_Yld89_AA6022T4.txt'
            data,pdata=polyn.readData(inpData,subDir)
            vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            if(ddata['errReport']):minKG,negKG=polyn.PolyN_predictions(vcf,pdata,export=True)
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
            return    
        if(yldFunc=='BBC2005'):##----- BBC2005
            ##Parameters are such that: a,b,L,M,N,P,Q,R=vAlpha[0,1,2,3,4,5,6,7] 
            vCoeff=BBC2005_To_PolyN(vAlpha_1,degree=degree)
            data,pdata=polyn.readData(inpData,subDir)
            vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            if(ddata['errReport']):minKG,negKG=polyn.PolyN_predictions(vcf,pdata,export=True)
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
            return
        if(yldFunc=='Caz2018_Ort'):##----Caz2018_Ort
            ## AA6016-T4 parameters of plane stress Cazacu2018-Orth in Habraken et al(2022)
            vCoeff=Caz2018_TO_P8(vAlpha_1,vAlpha_2,ddata['alpha'])
            data,pdata=polyn.readData(inpData,subDir)
            vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            if(ddata['errReport']):minKG,negKG=polyn.PolyN_predictions(vcf,pdata,export=True)
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
            return
        if(yldFunc=='FACET'):##----- FACET
            ##Parameters from Table-14 in Habraken et al(2022) paper
            #inpFile,subDir,sep=fData+'AA6016T4/FACETparams.txt','AA6016T4',' '
            #degree=8
            subDir,sep=ddata['subDir'],ddata['fileFACETsep']
            inpFile=polyn.figDirData+subDir+'/'+ddata['fileFACET']
            vCoeff=FACET_To_PolyN(inpFile,sep,degree=degree,norm=True)
            #inpData='matAA6016T4_TUAT.txt'
            data,pdata=polyn.readData(inpData,subDir)
            vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            if(ddata['errReport']):minKG,negKG=polyn.PolyN_predictions(vcf,pdata,export=True)
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
            return
        #if(yldFunc=='fYld2000_2D'): ### raw Yld2000_2D: used to illustrate the non-integer exponent case (export data)
        #    data,pdata=polyn.readData(inpData,subDir)
        #    fYld2000_2D(vAlpha_1,degree,pdata,export=True,saveFig=True)            
        if(0):##------  Poly6
            ##parameters from Habraken et al(2022) paper
            vCoeff=[1.00000, -2.06815, 3.53177, -3.93091, 3.91107, -2.49639, 1.10835, 16.41230,
                    -7.21657, 11.92569, -5.02031, 18.22988, 22.46969, 17.98091, 24.43169, 13.93065]  
            inpData,subDir='matCustom_AA6016T4_Poly6.txt','AA6016T4'                
            data,pdata=polyn.readData(inpData,subDir)
            degree=6
            vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            minKG,negKG=polyn.PolyN_predictions(vcf,pdata,export=True)
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
            return        
    '''
    Block 'paramHill'
    Calculates PolyN coefficients corresponding to 
    plane stress Hill48 reformulated as PolyN 
    '''
    if(False):
        eps=1.0e-12
        matName='AA6022T4' ##set the data export name of the material  
        inpData,subDir='matAA6022T4_P4.txt','AA6022T4_dig_bx'###input data file and its subDir 
        data,pdata=polyn.readData(inpData,subDir)
        ##Hill48 parameters are calculated based on r0,r45 and r90
        ##Extract r0,r45 and r90 r-values from input file
        kk,rVals=0,np.ones(3)
        for theta in pdata['thetaR']:
            #ang=np.pi*(theta/180.0)
            if(np.abs(theta-0.0)<eps):rVals[0]=pdata['rValue'][kk]
            if(np.abs(theta-45.0)<eps):rVals[1]=pdata['rValue'][kk]
            if(np.abs(theta-90.0)<eps):rVals[2]=pdata['rValue'][kk]
            kk+=1
        GH,FH,H,N2=spolyn.paramHill(rVals)        
        degree=12 ##set the degree N of the equivalent PolyN representation    
        vCoeff=spolyn.sHill48(degree,GH,FH,H,N2)
        exportName='symb_'+matName+'_HillPoly'
        vcf=polyn.PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=exportName,fPrint=True)
        spolyn.plotHill48(GH,FH,H,N2,degree,vcf,pdata['thetaS'],pdata['sigma'],pdata['thetaR'],pdata['rValue'],name=matName,figSave=True)
        exit()    
        