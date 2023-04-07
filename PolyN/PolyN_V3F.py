

###Developed and released by SCS 

import numpy as np
#import quadprog as qpg
import cvxopt
from matplotlib import pyplot as plt
from os import path as osp
from os import name as osn
from os import makedirs as osmaked
from time import time

figDir="./figsPolyN/"
figDirData="./figsPolyN/DATA/"
figDirPlot="./figsPolyN/PLOTS/"
if(osn=='nt'):
    figDir=".\\figsPolyN\\"
    figDirData=".\\figsPolyN\\DATA\\"
    figDirPlot=".\\figsPolyN\\PLOTS\\"
zro=0.0
zroTol=1.0e-10
gTol=1.0e-9

def readData(fName,subDir=''):
    global figDirData,figDirPlot
    if(not (osp.exists(figDir) and osp.isdir(figDir))):
        print("The local folder for saving reports and figures was not found")
        print("Make sure the folder \'figsPolyN\' exists at the same location as this script")
        print("Calculations aborted");exit()
    if(not (osp.exists(figDirData) and osp.isdir(figDirData))):
        print("The local folder for saving reports was not found")
        print("Make sure the folder \'DATA\' exists in \'figsPolyN\'")
        print("Calculations aborted");exit()
    if(not (osp.exists(figDirPlot) and osp.isdir(figDirPlot))):
        print("The local folder for saving plots was not found")
        print("Make sure a folder \'PLOTS\' exists in \'figsPolyN\'")
        print("Calculations aborted");exit()
    if(subDir):
        if(osn=='nt'):figDirData+=subDir+'\\'
        else:figDirData+=subDir+'/'
        if(not (osp.exists(figDirData) and osp.isdir(figDirData))):
            print("The local folder for saving reports was not found")
            print("Make sure the folder {} exists in \'figsPolyN\'".format(figDirData))
            print("Calculations aborted");exit()
        if(osn=='nt'):figDirPlot+=subDir+'\\'
        else:figDirPlot+=subDir+'/'
        if(not(osp.exists(figDirPlot) and osp.isdir(figDirPlot))):
            try:osmaked(figDirPlot)
            except OSError as err:printMsg('readData:');printMsg(err);exit()    
    data={'name':'','degree':6,'data':[],'wa':0.9,'war':0.6,'wvr':0.9,'nas':0,'nar':0,'nvs':0,'nvr':0,
           'dataAct':[],'wAct':0.9999,'wActR':0.25,'dataHill':{},'qHill':1}
    pdata={'fullName':'','degree':6,'thetaS':[],'sigma':[],'thetaR':[],'rValue':[],
           'sbiax':0.0, 'rbiax':0.0,'pstr':[],'FEdata':['*','*','*','*','*']}       
    tempData=[]
    cc,eq='#','='
    mexit,mline='\nCalculations aborted','Line in input file with incorrect format:\n'
    nas,nar,nvs,nvr,rad=0,0,0,0,np.pi/180
    try:ff=open(figDirData+fName,'r')
    except IOError as err:print(err);exit()
    for line in ff:
        line=line.strip()
        if((not line) or (line[0]==cc)):continue
        if(line[0:5]=='name='):data['name']=line[5:].strip()
        if(line[0:7]=='degree='):
            try:data['degree']=int(line[7:])
            except ValueError as err:print(mline+line+'\n'+str(err)+mexit);exit()
            if(data['degree']%2):print('Degree must be even'+mexit);exit()
        if(line[0:3]=='EE='):
            try:
                pdata['FEdata'][0]=float(line[3:].strip())
            except ValueError: pass
        if(line[0:3]=='NU='):
            try:
                pdata['FEdata'][1]=float(line[3:].strip())
            except ValueError: pass
        if(line[0:3]=='AA='):
            try:
                pdata['FEdata'][2]=float(line[3:].strip())
            except ValueError: pass
        if(line[0:3]=='BB='):
            try:
                pdata['FEdata'][3]=float(line[3:].strip())
            except ValueError: pass
        if(line[0:3]=='CC='):
            try:
                pdata['FEdata'][4]=float(line[3:].strip())
            except ValueError: pass            
        if(line[0:2]=='d='):
            vline=line[2:].split(',')
            if(len(vline) != 5):
                print(mline+line+mexit);exit()
            vline[4]=vline[4].strip()
            if(vline[4] not in ['a','v']):
                print(mline+line+mexit);exit()
            try:
                v1,v2=float(vline[0]),float(vline[1])
            except ValueError as err:print(mline+line+'\n'+str(err)+mexit);exit()
            if(np.abs(v1)>=1.000000001):print('Stress ratio q must be in [-1,+1]'+mexit);exit()
            try:
                v3=float(vline[2])
                if(vline[4]=='a'):nas+=1
                else:nvs+=1
            except ValueError as err:
                vline[2]=vline[2].strip()
                if(vline[2] != '*'):print(mline+line+'\n'+str(err)+mexit);exit()
                v3='*'
            try: 
                v4=float(vline[3])
                if(vline[4]=='a'):nar+=1
                else:nvr+=1
            except ValueError as err:
                vline[3]=vline[3].strip()
                if(vline[3] != '*'):print(mline+line+'\n'+str(err)+mexit);exit()
                v4='*'
            tempData.append((v1,v2*rad,v3,v4,vline[4]))
            if((v4!='*') and (-1.0e-6<v4<1.0e-6)):##plane strain 
                pdata['pstr'].append((v1,v2*rad,v3,vline[4]))            
        if(line[0:3]=='wa='):
            try: data['wa']=float(line[3:])
            except ValueError: pass
            if(data['wa']>1.0):print('wa > 1 (must be < 1)'+mexit);exit()
        if(line[0:4]=='war='):
            try: data['war']=float(line[4:])
            except ValueError: pass
            if(data['war']>1.0):print('war > 1 (must be < 1)'+mexit);exit()
        if(line[0:4]=='wvr='):
            try: data['wvr']=float(line[4:])
            except ValueError: pass
            if(data['wvr']>1.0):print('wvr > 1 (must be < 1)'+mexit);exit()
        if(line[0:6]=='qHill='):
            try: data['qHill']=int(line[6:])
            except ValueError: pass
            if(data['qHill']<0):print('qHill must be >=0'+mexit);exit()
            if(data['qHill']>1):data['qHill']=1
    ff.close()
    if(nas==0):print('Input file: Zero number of actual stresses'+mexit);exit()
    if(nar==0):print('Input file: Zero number of actual r-values'+mexit);exit()
    data['nas'],data['nar'],data['nvs'],data['nvr']=nas,nar,nvs,nvr
    wa,was,war,wv,wvs,wvr=data['wa'],0.0,0.0,0.0,0.0,0.0
    pwas,pwar,pwvs,pwvr=1.0,1.0,1.0,1.0
    if(nvs+nvr==0):##no data from virtual tests
        wa=1.0
    elif(nvr==0):wvs=1-wa;pwvs=wvs/nvs
    elif(nvs==0):wvr=1-wa;pwvr=wvr/nvr
    else:wv=1.0-wa;wvr=data['wvr']*wv;wvs=wv-wvr;pwvs,pwvr=wvs/nvs,wvr/nvr
    war=data['war']*wa;was=wa-war;pwas,pwar=was/nas,war/nar
    wActR=data['wActR']*data['wAct'];wActS=(1.0-data['wActR'])*data['wAct'];pwActS,pwActR=wActS/nas,wActR/nar
    nTotal,nParam=nas+nar+nvs+nvr,(int(data['degree']/2+1))**2
    if(nTotal<nParam-1):
        print("\nWarning:\n Total number of data points = {} < {} = Total number of parameters\n".format(nTotal,nParam))
    weight=0.0
    for item in tempData:
        if(item[2] !='*'):##yield stress data point
            if(item[4]=='a'):
                weight=pwas
                data['dataAct'].append((0,item[0],item[1],item[2],pwActS))
            else:weight=pwvs
            data['data'].append((0,item[0],item[1],item[2],weight))
        if(item[3] !='*'):##r-value data point
            if(item[4]=='a'):
                weight=pwar
                data['dataAct'].append((1,item[0],item[1],item[3],pwActR))
            else:weight=pwvr
            data['data'].append((1,item[0],item[1],item[3],weight))
        if(-1.0e-9<item[0]<1.0e-9):
            theta=item[1]
            if((item[4]=='a') and (item[2] != '*')):
                pdata['thetaS'].append(theta/rad);pdata['sigma'].append(item[2])
                #print('theta:',theta)
                if(45-1.0e-5<theta/rad<45+1.0e-5):
                    data['dataHill']['s45']=item[2]
                if(90-1.0e-5<theta/rad<90+1.0e-5):
                    data['dataHill']['s90']=item[2]
            if((item[4]=='a') and (item[3] != '*')):
                pdata['thetaR'].append(theta/rad);pdata['rValue'].append(item[3])
                if(45-1.0e-5<theta/rad<45+1.0e-5):
                    data['dataHill']['r45']=item[3]
        if((1.0-1.0e-9<item[0]<1.0+1.0e-9) and (-1.0e-9<item[1]<1.0e-9)):##balanced-biaxial
            pdata['sbiax'],pdata['rbiax']=item[2],item[3]  
    pdata['fullName'],pdata['degree']=data['name'],data['degree']            
    if(0):###verify weights
        ss=0.0
        for item in data['data']:
            ss+=item[-1]
        print('Weights verification: Sum of all weights = ',ss)
    ##testData(data,pdata);exit()   
    return data,pdata    


def testData(data,pdata):
    for kk in data:
        print('data['+kk+']\n',data[kk])
    print('--------------PlotData')    
    for kk in pdata:
        print('pdata['+kk+']\n',pdata[kk])
    return        

def nMonoms(degree):
    if(int(degree)%2):
        print("PolyN degree must be an even integer\nCalculations aborted")
        exit()    
    return (int(degree/2)+1)**2

def vPoly(degree):
    dd={'nQ':degree,'nDim':nMonoms(degree)}
    vQ=[];vQidx=[];jj=0
    while(jj<=degree):
        lv=len(vQ)
        vQidx.append([jj,[k for k in range(lv,lv+degree+1-jj)]])
        vQ+=[(degree-jj-k,k,jj) for k in range(degree+1-jj)]
        jj+=2 
    dd['vQ']=vQ;dd['vQidx']=vQidx
    dd['vD1Q']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in vQ]
    dd['vC1Q']=np.array([mon[0] if(mon[0]>0) else 0 for mon in vQ])
    dd['vD2Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in vQ]
    dd['vC2Q']=np.array([mon[1] if(mon[1]>0) else 0 for mon in vQ])
    dd['vD3Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in vQ]
    dd['vC3Q']=np.array([mon[2] if(mon[2]>0) else 0 for mon in vQ])
    dd['vH11Q']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH11Q']=np.array([cf*mon[0] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH12Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH12Q']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH13Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH13Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH22Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD2Q']]
    dd['vCH22Q']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC2Q'],dd['vD2Q'])])
    dd['vH23Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD2Q']]
    dd['vCH23Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC2Q'],dd['vD2Q'])])
    dd['vH33Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD3Q']]
    dd['vCH33Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC3Q'],dd['vD3Q'])])     
    return dd


'''    
'''
def Bezier5YSSeg(Bs,Ts,Be,Te,LshapeS,LshapeE,plot=False,nTT=11):
    if(plot):
        nTT=101;vTT=np.linspace(0.0,1.0,nTT)
    else:    
        vTT=np.linspace(0.2,0.8,nTT)
    TTs=LshapeS*Ts;TTe=LshapeE*Te
    B1=Bs+TTs;B2=B1+TTs;B4=Be-TTe;B3=B4-TTe
    A1=5.0*(B1-Bs)
    A2=10.0*(Bs+B2-2.0*B1)
    A3=10.0*(3.0*B1-Bs+B3-3.0*B2)
    A4=5.0*(Bs-4.0*B1+6.0*B2-4.0*B3+B4)
    A5=5.0*(B1-2.0*B2+2.0*B3-B4)-Bs+Be
    return Bs+vTT*(A1+vTT*(A2+vTT*(A3+vTT*(A4+vTT*A5))))


'''
'''
def lambdaMax(bOne,tOne,bTwo,tTwo):
    bb=bTwo-bOne
    mdot=np.sum(tOne*tTwo)
    ddet=1.0-mdot*mdot
    t1=0.5*np.sum(bb*(tOne-mdot*tTwo))/ddet
    t2=0.5*np.sum(bb*(tTwo-mdot*tOne))/ddet
    #return min(t1,t2)
    if(t1<0): print('(lambdaMax)Warning: negative t1 = ',t1)
    if(t1<0): print('(lambdaMax)Warning: negative t2 = ',t2)
    return t1,t2


'''
wS=0.7,wE=1.0
'''
def biaxData(data,wS=1.0,wE=1.0):
    vPoints={'s0':np.zeros((2,1)), 's90':np.zeros((2,1)), 'sb':np.zeros((2,1)), 
             't0':np.zeros((2,1)), 't90':np.zeros((2,1)), 'tb':np.zeros((2,1)),
             'spsX':np.zeros((2,1)),'spsY':np.zeros((2,1)),
             'tpsX':np.array([0.0,1.0]).reshape((2,1)),'tpsY':np.array([-1.0,0.0]).reshape((2,1))}
    eps,fGrad,degree,vFidx,vFcoeff,r90,s90=1.0e-8,1.0,data['degree'],np.zeros(4,dtype=int),np.ones(4),1.0,1.0
    psX,psY,qqX,qqY,spsX,spsY,psXtheta,psYtheta=False,False,0.0,0.0,np.zeros((2,1)),np.zeros((2,1)),0.0,0.0
    ang=180/np.pi
    #print(data['dataAct'])
    for item in data['dataAct']:
        #print(item)
        qq,theta=item[1],item[2]*ang
        if(-eps<qq<eps): ###uniaxial 
            if(-eps<theta<eps):##theta=0
                if(item[0]):##r-value
                    fGrad=-item[3]/(1.0+item[3]);fN=np.sqrt(1.0+fGrad*fGrad)
                    vPoints['t0'][:,0]=-fGrad/fN,1.0/fN
                    vFidx[1],vFcoeff[1]=1,(-degree*item[3])/(1+item[3])
                else:vPoints['s0'][:,0]=1.0,0.0
            elif(90.0-eps<theta<90.0+eps):##theta=90
                if(item[0]):##r-value
                    r90=item[3]
                    fGrad=-r90/(1.0+r90);fN=np.sqrt(1.0+fGrad*fGrad)
                    vPoints['t90'][:,0]=-1.0/fN,fGrad/fN
                else:s90=item[3];vPoints['s90'][:,0]=0.0,s90
        if(eps<qq<1-eps): 
            if(item[0] and (-eps<theta<eps) and (-eps<item[3]<eps)): ## plane strain along RD
                psX,qqX,psXtheta=True,qq,theta
            if(item[0] and (90-eps<theta<90+eps) and (-eps<item[3]<eps)): ## plane strain along TD
                psY,qqY,psYtheta=True,qq,theta
            #print(psX,psY)                
        if((1.0-eps<qq<1.0+eps) and (-eps<theta<eps)):##balanced-biaxial
            if(item[0]):
               fGrad=-item[3]/(1.0+item[3]);fN=np.sqrt(1.0+fGrad*fGrad)
               vPoints['tb'][:,0]=-fGrad/fN,1.0/fN
            else:vPoints['sb'][:,0]=item[3],item[3]            
    vPoints['s0'][:,0]=1.0,0.0
    s90=s90**degree               
    vFidx[2],vFcoeff[2]=degree-1,(-degree*r90)/((1+r90)*s90)
    vFidx[3],vFcoeff[3]=degree,1.0/s90
    vSegments=[('s0','t0')]
    if(psX):
        for item in data['dataAct']:
            if(item[0]):continue
            if(abs(item[1]-qqX)>eps):continue
            if(abs(item[2]-psXtheta)>eps): continue
            vPoints['spsX'][:,0]=item[3],qqX*item[3]
            vSegments.append(('spsX','tpsX'))
            wS,wE=1.0,1.0
    vSegments.append(('sb','tb'))        
    if(psY):
        #print(data['dataAct'])
        for item in data['dataAct']:
            #print(item)
            if(item[0]):continue
            if(abs(item[1]-qqY)>eps):continue
            if(abs(item[2]*ang-psYtheta)>eps): continue
            vPoints['spsY'][:,0]=qqY*item[3],item[3]
            vSegments.append(('spsY','tpsY'))
            wS,wE=1.0,1.0
    #print(vSegments);exit()        
    vSegments.append(('s90','t90'))
    nP,nSeg=degree+1,len(vSegments)-1
    #for kk in vPoints:
    #    print(kk,': ', vPoints[kk])
    #print(psX,psY);exit()    
    vv2=np.zeros((2,nSeg*nP))
    LBDS,LBDE=100.0,100.0
    for kk in range(nSeg):
        pt1,pt2=vSegments[kk],vSegments[kk+1]
        lbd1,lbd2=lambdaMax(vPoints[pt1[0]],vPoints[pt1[1]],vPoints[pt2[0]],vPoints[pt2[1]])
        LBDS,LBDE=min(LBDS,lbd1),min(LBDE,lbd2)
    LBDS,LBDE=wS*LBDS,wE*LBDE
    kVV=0
    for kk in range(nSeg):
        pt1,pt2=vSegments[kk],vSegments[kk+1]
        vv=Bezier5YSSeg(vPoints[pt1[0]],vPoints[pt1[1]],vPoints[pt2[0]],vPoints[pt2[1]],LBDS,LBDE,plot=False,nTT=nP)
        vv2[0,kVV:kVV+nP],vv2[1,kVV:kVV+nP]=vv[0,:],vv[1,:]
        kVV+=nP        
    #lbd11,lbd12=lambdaMax(vPoints['s0'],vPoints['t0'],vPoints['sb'],vPoints['tb'])
    #lbd21,lbd22=lambdaMax(vPoints['sb'],vPoints['tb'],vPoints['s90'],vPoints['t90'])
    #lbdS,lbdE=wS*min(lbd11,lbd21),wE*min(lbd12,lbd22)
    #vv=Bezier5YSSeg(vPoints['s0'],vPoints['t0'],vPoints['sb'],vPoints['tb'],lbdS,lbdE,plot=False)
    ####print(vv[0,:]);exit()
    #vv2=np.zeros((2,2*vv.shape[1]))
    #vv2[0,:vv.shape[1]],vv2[1,:vv.shape[1]]=vv[0,:],vv[1,:]
    #vv=Bezier5YSSeg(vPoints['sb'],vPoints['tb'],vPoints['s90'],vPoints['t90'],lbdS,lbdE,plot=False)
    #vv2[0,vv.shape[1]:],vv2[1,vv.shape[1]:]=vv[0,:],vv[1,:]
    #print(vv2[0,:]);exit()
    return vv2,vFidx,vFcoeff



def tstPlotBx(data,pdata,vCoeff,vMonoms,saveFig=False):
    vv,vFidx,vFcoeff=biaxData(data)
    plotBiax(vv,pdata,vCoeff,vMonoms,saveFig)
    plt.show()
  
def plotBiax(vpoints,data,vCoeff,vMonoms,saveFig=False):
    fig=plt.figure()
    ax=fig.add_subplot()
    ax.plot(vpoints[0,:],vpoints[1,:],linewidth=2,linestyle='dashed',color='r',label='Data')
    ax.plot(-vpoints[0,:],-vpoints[1,:],linewidth=2,linestyle='dashed',color='r')
    nomega=401
    vomega=np.linspace(0.0,2*np.pi,nomega)
    vsx,vsy,vrho=np.cos(vomega),np.sin(vomega),np.zeros(nomega)
    deg=vCoeff.shape[0]
    vPow=vMonoms['vQ'][0:deg] ###; print(vPow)
    for j in range(deg):
        vrho[:]+=vCoeff[j]*(vsx**vPow[j][0])*(vsy**vPow[j][1])
    vrho=1.0/vrho**(1/(deg-1))
    ax.plot(vrho*vsx,vrho*vsy,linewidth=1,color='k',label='Poly{}-fit'.format(deg-1))
    ax.set_aspect('equal');ax.grid()
    ax.legend(loc='center',title='')
    if(saveFig):
        name=data['fullName']+'_P{}'.format(str(data['degree']))
        fig.savefig(figDirPlot+name+'_Biax.png',bbox_inches='tight',dpi=300)

def genConstraintsPoints2DOpt(nPoints=250,number=False):
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1[1:]);vs1=np.sin(vt1[1:])
    vN2=np.int_(nPoints*vs1)
    tt=[(np.cos(t),np.sin(t)) for t in np.linspace(0,np.pi,72)[1:-1]]        
    nVec=len(tt)
    npt=np.sum(vN2)+1
    if(number): return (npt,3*(nVec+3))
    vP=np.zeros((npt,3*(nVec+3)))
    nConstr=npt*(nVec+2)
    vN2[:]+=1
    vP[0,2]=1.0
    vP[0,3]=1.0  #g1[0]
    vP[0,7]=1.0  ##g2[1]
    i=9
    for t in tt:
        gg=t[0]*vP[0,3:6]+t[1]*vP[0,6:9]
        gg[:]/=np.sqrt(np.sum(gg*gg,axis=0))
        vP[0,i],vP[0,i+1],vP[0,i+2]=gg
        i+=3
    jj=1
    for kk in range(0,len(vt1)-1):
        ct1=vc1[kk];st1=vs1[kk]
        N2=vN2[kk];jN2=jj+N2-1
        vt2=np.linspace(0,np.pi,N2)
        vc2,vs2=np.cos(vt2[0:N2-1]),np.sin(vt2[0:N2-1])
        ##print("kk={}, N2={}, jN2={}, jj={}".format(kk,N2,jN2,jj))
        vP[jj:jN2,0]=st1*vc2
        vP[jj:jN2,1]=st1*vs2
        vP[jj:jN2,2]=ct1
        vP[jj:jN2,3]=-vs2*st1  ##g1[0]
        vP[jj:jN2,4]=vc2*st1
        vP[jj:jN2,6]=ct1*vc2  #g2[0]
        vP[jj:jN2,7]=ct1*vs2
        vP[jj:jN2,8]=-st1
        i=9
        for t in tt:
            gg=t[0]*vP[jj:jN2,3:6]+t[1]*vP[jj:jN2,6:9]
            vMod=np.sqrt(np.sum(gg*gg,axis=1))
            gg[:,0]/=vMod;gg[:,1]/=vMod;gg[:,2]/=vMod
            vP[jj:jN2,i]=gg[:,0];vP[jj:jN2,i+1]=gg[:,1];vP[jj:jN2,i+2]=gg[:,2]
            i+=3
        jj=jN2 
    return nConstr,vP

def genConstraintsPoints2DOptB(pNegKG):
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1[1:]);vs1=np.sin(vt1[1:])
    vN2=np.int_(nPoints*vs1)
    tt=[(np.cos(t),np.sin(t)) for t in np.linspace(0,np.pi,72)[1:-1]]        
    nVec=len(tt)
    pNegKG=np.array(pNegKG)
    nKG=pNegKG.shape[0]
    npt=nKG+np.sum(vN2)+1
    if(number): return (npt,3*(nVec+3))
    vP=np.zeros((npt,3*(nVec+3)))
    nConstr=npt*(nVec+2)
    vN2[:]+=1
    vP[0,2]=1.0
    vP[0,3]=1.0  #g1[0]
    vP[0,7]=1.0  ##g2[1]
    i=9
    for t in tt:
        gg=t[0]*vP[0,3:6]+t[1]*vP[0,6:9]
        gg[:]/=np.sqrt(np.sum(gg*gg,axis=0))
        vP[0,i],vP[0,i+1],vP[0,i+2]=gg
        i+=3
    jj=1
    for kk in range(0,len(vt1)-1):
        ct1=vc1[kk];st1=vs1[kk]
        N2=vN2[kk];jN2=jj+N2-1
        vt2=np.linspace(0,np.pi,N2)
        vc2,vs2=np.cos(vt2[0:N2-1]),np.sin(vt2[0:N2-1])
        ##print("kk={}, N2={}, jN2={}, jj={}".format(kk,N2,jN2,jj))
        vP[jj:jN2,0]=st1*vc2
        vP[jj:jN2,1]=st1*vs2
        vP[jj:jN2,2]=ct1
        vP[jj:jN2,3]=-vs2*st1  ##g1[0]
        vP[jj:jN2,4]=vc2*st1
        vP[jj:jN2,6]=ct1*vc2  #g2[0]
        vP[jj:jN2,7]=ct1*vs2
        vP[jj:jN2,8]=-st1
        i=9
        for t in tt:
            gg=t[0]*vP[jj:jN2,3:6]+t[1]*vP[jj:jN2,6:9]
            vMod=np.sqrt(np.sum(gg*gg,axis=1))
            gg[:,0]/=vMod;gg[:,1]/=vMod;gg[:,2]/=vMod
            vP[jj:jN2,i]=gg[:,0];vP[jj:jN2,i+1]=gg[:,1];vP[jj:jN2,i+2]=gg[:,2]
            i+=3
        jj=jN2
    if(nKG):    
        st1=np.sin(np.arccos(pNegKG[:,2]))    
        vt2=np.arctan2(pNegKG[:,1],pNegKG[:,0])
        vc2,vs2=np.cos(vt2),np.sin(vt2)
        vP[jj:,0]=pNegKG[:,0];vP[jj:,1]=pNegKG[:,1];vP[jj:,2]=pNegKG[:,2]
        vP[jj:,3]=-vs2*st1  ##g1[0]
        vP[jj:,4]=vc2*st1
        vP[jj:,6]=pNegKG[:,2]*vc2  #g2[0]
        vP[jj:,7]=pNegKG[:,2]*vs2
        vP[jj:,8]=-st1
        i=9
        for t in tt:
            gg=t[0]*vP[jj:,3:6]+t[1]*vP[jj:,6:9]
            vMod=np.sqrt(np.sum(gg*gg,axis=1))
            gg[:,0]/=vMod;gg[:,1]/=vMod;gg[:,2]/=vMod
            vP[jj:,i]=gg[:,0];vP[jj:,i+1]=gg[:,1];vP[jj:,i+2]=gg[:,2]
            i+=3        
    return nConstr,vP    


def genConstraintsPoints2DOptZ(pNegKG):
    tt=[(np.cos(t),np.sin(t)) for t in np.linspace(0,np.pi,72)[1:-1]]        
    nVec=len(tt)
    ppNegKG=np.array(pNegKG)
    nKG=ppNegKG.shape[0]
    vP=np.zeros((nKG,3*(nVec+3)))
    nConstr=nKG*(nVec+2)
    st1=np.sin(np.arccos(ppNegKG[:,2]))    
    vt2=np.arctan2(ppNegKG[:,1],ppNegKG[:,0])
    vc2,vs2=np.cos(vt2),np.sin(vt2)
    vP[:,0]=ppNegKG[:,0];vP[:,1]=ppNegKG[:,1];vP[:,2]=ppNegKG[:,2]
    vP[:,3]=-vs2*st1  ##g1[0]
    vP[:,4]=vc2*st1
    vP[:,6]=ppNegKG[:,2]*vc2  #g2[0]
    vP[:,7]=ppNegKG[:,2]*vs2
    vP[:,8]=-st1
    i=9
    for t in tt:
        gg=t[0]*vP[:,3:6]+t[1]*vP[:,6:9]
        vMod=np.sqrt(np.sum(gg*gg,axis=1))
        gg[:,0]/=vMod;gg[:,1]/=vMod;gg[:,2]/=vMod
        vP[:,i]=gg[:,0];vP[:,i+1]=gg[:,1];vP[:,i+2]=gg[:,2]
        i+=3        
    return nConstr,vP
    
def genTanDir2D(pPoint,nDirs):
    tt=[(np.cos(t),np.sin(t)) for t in np.linspace(0,2*np.pi,nDirs+1)[0:-1]]        
    vP=np.zeros((nDirs,3)) #; print(nDirs);print(vP.shape)
    st1=np.sin(np.arccos(pPoint[2]))    
    vt2=np.arctan2(pPoint[1],pPoint[0])
    vc2,vs2=np.cos(vt2),np.sin(vt2)
    g1,g2=np.zeros(3),np.zeros(3)
    g1[0],g1[1]=-vs2*st1,vc2*st1
    g2[0],g2[1],g2[2]=pPoint[2]*vc2,pPoint[2]*vs2,-st1
    i=0
    for t in tt:
        gg=t[0]*g1+t[1]*g2
        vMod=np.sqrt(np.sum(gg*gg))
        gg[:]/=vMod
        vP[i,:]=gg[:]
        i+=1   #;print(i)  
    ##print(vP)        
    return vP    
    
def genConstraintsTG(vMonoms,pNegKG=[]):
    print('Generating constraints.......')
    degree,ndim=vMonoms['nQ'],vMonoms['nDim']
    oDegree=1.0-1.0/degree
    if(pNegKG):
        nTotal,vvp=genConstraintsPoints2DOptZ(pNegKG)
    else:
        if(degree<11):
            nTotal,vvp=genConstraintsPoints2DOpt()
        else:
            nTotal,vvp=genConstraintsPoints2DOpt(275)    
    nvOne=vvp.shape[0]    
    QQ=np.zeros((nTotal,ndim,ndim))
    vQ=vMonoms['vQ']
    vD1,vD2,vD3=vMonoms['vD1Q'],vMonoms['vD2Q'],vMonoms['vD3Q']
    vC1,vC2,vC3=vMonoms['vC1Q'],vMonoms['vC2Q'],vMonoms['vC3Q']
    vH11,vH12,vH13=vMonoms['vH11Q'],vMonoms['vH12Q'],vMonoms['vH13Q']
    vH22,vH23,vH33=vMonoms['vH22Q'],vMonoms['vH23Q'],vMonoms['vH33Q']
    vCH11,vCH12,vCH13=vMonoms['vCH11Q'],vMonoms['vCH12Q'],vMonoms['vCH13Q']
    vCH22,vCH23,vCH33=vMonoms['vCH22Q'],vMonoms['vCH23Q'],vMonoms['vCH33Q']
    ######
    D1,D2,D3=np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim))
    H11,H12,H13=np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim))
    H22,H23,H33=np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim))
    vAlpha,vBeta,vGamma=np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim)),np.zeros((nvOne,ndim))
    vsx,vsy,vsxy=vvp[:,0],vvp[:,1],vvp[:,2]
    vdsx,vdsy,vdsxy=np.zeros(nvOne),np.zeros(nvOne),np.zeros(nvOne)
    idx,idk=0,0
    #print((vvp.shape[1]-3)/3);print(vvp.shape);print(nvOne);print(vvp[:,3].shape)
    for kt in range(int((vvp.shape[1]-3)/3)):
        idx+=3
        vdsx[:],vdsy[:],vdsxy[:]=vvp[:,idx],vvp[:,idx+1],vvp[:,idx+2]
        #v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
        for j in range(ndim):
            vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
            D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
            D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
            D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
            vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
            H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
            H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
            H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
            H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
            H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
            H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
            vGamma[:,j]=H11[:,j]*vdsx*vdsx+H22[:,j]*vdsy*vdsy+H33[:,j]*vdsxy*vdsxy+\
            2.0*(H12[:,j]*vdsx*vdsy+H13[:,j]*vdsx*vdsxy+H23[:,j]*vdsy*vdsxy)-(D1[:,j]*vsx+D2[:,j]*vsy+D3[:,j]*vsxy)
        for ii in range(nvOne):
            QQ[idk+ii,:,:]=  degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
            oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
            0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
            ####QQ[idk+ii,:,:]/=np.max(np.abs(QQ[idk+ii,:,:]))
        idk+=nvOne           
    return nTotal,QQ

def genConstraints7B():
    N1,nPoints=71,143
    #dpi=np.pi/144
    #vt1=np.linspace(8*dpi,np.pi/2-dpi,N1)
    vt1=np.linspace((0.25/180.0)*np.pi,(89.75/180.0)*np.pi,N1)
    vc1,vs1=np.cos(vt1),np.sin(vt1)
    vN2=np.int_(nPoints*vs1)     
    nTotal=np.sum(vN2)
    vE1,vE2=np.zeros((nTotal,2)),np.zeros((nTotal,3))
    jj=0
    for kk in range(0,len(vt1)):
        ct1=vc1[kk];st1=vs1[kk]
        jN2=jj+vN2[kk]
        vt2=np.linspace(0,np.pi,vN2[kk])##half
        vc2,vs2=np.cos(vt2),np.sin(vt2)
        vE1[jj:jN2,0]=np.cos(vt2-np.pi/2)
        vE1[jj:jN2,1]=np.sin(vt2-np.pi/2)
        vE2[jj:jN2,0]=-ct1*vE1[jj:jN2,1]
        vE2[jj:jN2,1]=ct1*vE1[jj:jN2,0]
        vE2[jj:jN2,2]=st1*vc2*vE1[jj:jN2,1]-st1*vs2*vE1[jj:jN2,0]
        jj=jN2
    return vE1,vE2 

'''
nPsiZ=61,nOmegaZ=111
nPsiZ=71,nOmegaZ=121
'''
def genConstraints(vMonoms,nPsiZ=501,nOmegaZ=601):
    print('Generating constraints.......')
    degree,ndim=vMonoms['nQ'],vMonoms['nDim']
    oDegree=1.0-1.0/degree
    nSet1,nSet2,nSet3=2*nPsiZ,nPsiZ,2*nPsiZ
    nOmega1,nOmega2,nOmega3=nOmegaZ,2*nOmegaZ,2*nOmegaZ
    nTotal=(nSet1-2)*nOmega1+(nSet2-2)*nOmega2+(nSet3-1)*nOmega3
    nSet4,nSet5=nPsiZ,nPsiZ
    nTotal+=(nSet4-2)*nOmega3+(nSet5-2)*nOmega3
    #nSet4,nSet5=nPsi-2,nPsi-2
    #nOmega3=3*nOmega
    #nOmega3o=nOmega3-1    
    ##nTotal=int((nSet1+nSet2+nSet4+nSet5)*nOmega+nSet3*nOmega3o)
    #nTotal=int((nSet2+nSet4+nSet5)*nOmega+(nSet1-2+nSet3)*nOmega3o)
    nUp,afa,afaB=7,False,False
    if(degree>nUp):
        afa=False
        #vafa=(np.pi/12,np.pi/6,np.pi/3,5*np.pi/12,-np.pi/12,-np.pi/6,-np.pi/3,-5*np.pi/12)
        #vafa=((50.0/180.0)*np.pi,(40.0/180.0)*np.pi, (55.0/180.0)*np.pi,(35.0/180.0)*np.pi)
        vafa=((50.0/180.0)*np.pi,(40.0/180.0)*np.pi, (55.0/180.0)*np.pi,(35.0/180.0)*np.pi)
        #vafa=(np.pi/8,)
        #vafa=(np.pi/12,np.pi/6,np.pi/3,5*np.pi/12,7*np.pi/12,4*np.pi/6,5*np.pi/6,11*np.pi/12, np.pi/5, 3*np.pi/10)
        #nafa=17
        #vafa=[k*np.pi/nafa for k in range(1,nafa)]
        kafa,nSet6=len(vafa),nSet1
        if(afa):
            nTotal+=kafa*(nSet6-2)*nOmega1
        afaB=False    
        vafaB=((50.0/180.0)*np.pi, (55.0/180.0)*np.pi, (60.0/180.0)*np.pi,
               (65.0/180.0)*np.pi, (70.0/180.0)*np.pi, (75.0/180.0)*np.pi,
               (80.0/180.0)*np.pi, (85.0/180.0)*np.pi,               
               (40.0/180.0)*np.pi, (35.0/180.0)*np.pi, (30.0/180.0)*np.pi,
               (25.0/180.0)*np.pi, (20.0/180.0)*np.pi, (15.0/180.0)*np.pi,
               (10.0/180.0)*np.pi, (5.0/180.0)*np.pi)
        kafaB,nSet2B=len(vafaB),nSet2
        if(afaB):
            nTotal+=kafaB*(nSet2B-2)*nOmega2        
        afa2=0
        if(afa2):        
            vE1,vE2=genConstraints7B()
            nTotal+=vE1.shape[0]*nOmega3 
        #print('vE1.shape = ',vE1.shape)
        #print(vE1);print(vE2);exit()        
    QQ=np.zeros((nTotal,ndim,ndim))
    vQ=vMonoms['vQ']
    vD1,vD2,vD3=vMonoms['vD1Q'],vMonoms['vD2Q'],vMonoms['vD3Q']
    vC1,vC2,vC3=vMonoms['vC1Q'],vMonoms['vC2Q'],vMonoms['vC3Q']
    vH11,vH12,vH13=vMonoms['vH11Q'],vMonoms['vH12Q'],vMonoms['vH13Q']
    vH22,vH23,vH33=vMonoms['vH22Q'],vMonoms['vH23Q'],vMonoms['vH33Q']
    vCH11,vCH12,vCH13=vMonoms['vCH11Q'],vMonoms['vCH12Q'],vMonoms['vCH13Q']
    vCH22,vCH23,vCH33=vMonoms['vCH22Q'],vMonoms['vCH23Q'],vMonoms['vCH33Q']
    R11,R21,R31,R12,R22,R32,q0=0.0,0.0,0.0,0.0,0.0,0.0,0.0
    sq2=np.sqrt(2.0)
    idx1=0
    #######first family of constraints
    nPsi,nOmega=nSet1,nOmega1
    D1,D2,D3=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    H11,H12,H13=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    H22,H23,H33=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    vAlpha,vBeta,vGamma=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    vsx,vsy,vsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    vdsx,vdsy,vdsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    v11,v12,v13=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    v22,v23,v33=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    dOmega=(np.pi/2.0)/nOmega
    vOmega=np.array([k*dOmega for k in range(nOmega)])
    comega,somega=np.cos(vOmega),np.sin(vOmega)
    vPsi=np.linspace(0.0,np.pi,nPsi)
    cpsi,spsi=np.cos(vPsi),np.sin(vPsi)
    for k in range(1,nSet1-1):
        R11,R21,R31,R12,R22,R32=1/sq2,1/sq2,0.0,-cpsi[k]/sq2,cpsi[k]/sq2,-spsi[k]
        vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
        vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
        v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
        for j in range(ndim):
            vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
            D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
            D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
            D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
            vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
            H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
            H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
            H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
            H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
            H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
            H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
            vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
        for ii in range(nOmega):
            QQ[idx1+ii,:,:]=  degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                   oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                   0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
            #       
        idx1+=nOmega
    if(afa):
    #if(0):    
        for alpha in vafa:    
            calpha,salpha=np.cos(alpha),np.sin(alpha)
            for k in range(1,nSet6-1):##6-th family of constraints
                R11,R21,R31,R12,R22,R32=calpha,salpha,0.0,-salpha*cpsi[k],calpha*cpsi[k],spsi[k]
                vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
                vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
                v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
                for j in range(ndim):
                    vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
                    D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
                    D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
                    D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
                    vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
                    H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
                    H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
                    H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
                    H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
                    H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
                    H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
                    vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
                for ii in range(nOmega):
                    QQ[idx1+ii,:,:]=  degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                           oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                           0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
                idx1+=nOmega    
    #######second family of constraints
    nPsi,nOmega=nSet2,nOmega2
    D1,D2,D3=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    H11,H12,H13=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    H22,H23,H33=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    vAlpha,vBeta,vGamma=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    vsx,vsy,vsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    vdsx,vdsy,vdsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    v11,v12,v13=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    v22,v23,v33=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    dOmega=np.pi/nOmega
    vOmega=np.array([k*dOmega for k in range(nOmega)])
    comega,somega=np.cos(vOmega),np.sin(vOmega)
    vPsi=np.linspace(0.0,np.pi/2.0,nPsi)
    cpsi,spsi=np.cos(vPsi),np.sin(vPsi)        
    for k in range(1,nSet2-1):
        R11,R21,R31,R12,R22,R32=1/sq2,-1/sq2,0.0,cpsi[k]/sq2,cpsi[k]/sq2,-spsi[k]
        vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
        vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
        v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
        for j in range(ndim):
            vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
            D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
            D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
            D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
            vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
            H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
            H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
            H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
            H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
            H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
            H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
            vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
        for ii in range(nOmega):
            QQ[idx1+ii,:,:]=degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                       oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                       0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
        idx1+=nOmega
    if(afaB):    
        for alpha in vafaB:    
            calpha,salpha=np.cos(alpha),np.sin(alpha)
            for k in range(1,nSet2B-1):## 2-nd family extended 
                R11,R21,R31,R12,R22,R32=-salpha,calpha,0.0,-calpha*cpsi[k],-salpha*cpsi[k],-spsi[k]
                vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
                vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
                v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
                for j in range(ndim):
                    vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
                    D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
                    D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
                    D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
                    vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
                    H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
                    H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
                    H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
                    H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
                    H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
                    H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
                    vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
                for ii in range(nOmega):
                    QQ[idx1+ii,:,:]=  degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                           oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                           0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
                idx1+=nOmega        
    #######third family of constraints
    nPsi,nOmega=nSet3,nOmega3
    D1,D2,D3=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    H11,H12,H13=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    H22,H23,H33=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    vAlpha,vBeta,vGamma=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    vsx,vsy,vsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    vdsx,vdsy,vdsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    v11,v12,v13=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    v22,v23,v33=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    dOmega=np.pi/nOmega
    vOmega=np.array([k*dOmega for k in range(nOmega)])
    comega,somega=np.cos(vOmega),np.sin(vOmega)
    vPsi=np.linspace(0.0,np.pi,nPsi)
    cpsi,spsi=np.cos(vPsi),np.sin(vPsi)        
    for k in range(0,nSet3-1):
        R11,R21,R31,R12,R22,R32=cpsi[k],spsi[k],0.0,0.0,0.0,1.0
        vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
        vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
        v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
        for j in range(ndim):
            vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
            D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
            D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
            D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
            vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
            H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
            H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
            H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
            H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
            H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
            H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
            vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
        for ii in range(nOmega):
            QQ[idx1+ii,:,:]=degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                       oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                       0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
        idx1+=nOmega
    if(1):        
        for k in range(1,nSet4-1):##fourth family of constraints
            R11,R21,R31,R12,R22,R32=1.0,0.0,0.0,0.0,spsi[k],cpsi[k]
            vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
            vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
            v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
            for j in range(ndim):
                vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
                D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
                D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
                D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
                vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
                H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
                H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
                H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
                H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
                H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
                H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
                vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
            for ii in range(nOmega):
                QQ[idx1+ii,:,:]=  degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                       oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                       0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
                #       
            idx1+=nOmega        
        for k in range(1,nSet5-1):##fifth family of constraints
            R11,R21,R31,R12,R22,R32=0.0,1.0,0.0,-spsi[k],0.0,cpsi[k]
            vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
            vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
            v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
            for j in range(ndim):
                vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
                D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
                D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
                D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
                vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
                H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
                H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
                H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
                H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
                H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
                H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
                vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
            for ii in range(nOmega):
                QQ[idx1+ii,:,:]=  degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                       oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                       0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
                #       
            idx1+=nOmega
    if(degree<nUp):
        return nTotal,QQ
    ####
    if(afa2):                
        for k in range(vE1.shape[0]):##7-th family of constraints
            R11,R21,R31,R12,R22,R32=vE1[k,0],vE1[k,1],0.0,vE2[k,0],vE2[k,1],vE2[k,2]
            vsx[:],vsy[:],vsxy[:]=R11*comega+R12*somega,R21*comega+R22*somega,R31*comega+R32*somega
            vdsx[:],vdsy[:],vdsxy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega,-R31*somega+R32*comega
            v11[:],v12[:],v13[:],v22[:],v23[:],v33[:]=vdsx*vdsx,vdsx*vdsy,vdsx*vdsxy,vdsy*vdsy,vdsy*vdsxy,vdsxy*vdsxy
            for j in range(ndim):
                vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])*(vsxy**vQ[j][2])
                D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])*(vsxy**vD1[j][2])
                D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])*(vsxy**vD2[j][2])
                D3[:,j]=vC3[j]*(vsx**vD3[j][0])*(vsy**vD3[j][1])*(vsxy**vD3[j][2])
                vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy+D3[:,j]*vdsxy
                H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])*(vsxy**vH11[j][2])
                H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])*(vsxy**vH12[j][2])
                H13[:,j]=vCH13[j]*(vsx**vH13[j][0])*(vsy**vH13[j][1])*(vsxy**vH13[j][2])
                H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])*(vsxy**vH22[j][2])
                H23[:,j]=vCH23[j]*(vsx**vH23[j][0])*(vsy**vH23[j][1])*(vsxy**vH23[j][2])
                H33[:,j]=vCH33[j]*(vsx**vH33[j][0])*(vsy**vH33[j][1])*(vsxy**vH33[j][2])
                vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+H33[:,j]*v33+2.0*(H12[:,j]*v12+H13[:,j]*v13+H23[:,j]*v23)-D1[:,j]*vsx-D2[:,j]*vsy-D3[:,j]*vsxy
            for ii in range(nOmega):
                QQ[idx1+ii,:,:]=  degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                       oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                       0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
            idx1+=nOmega                    
    return nTotal,QQ

 
def genLinConstraints(vQQ,va,vFidx,vFcoeff,vhh,vGG,eLambda=0.999):
    nTotal,nDIM=vQQ.shape[0:2]
    ndim=nDIM-vFidx.shape[0]
    last=vhh.shape[0]-1
    GG,hh=np.zeros((nTotal,ndim)),np.zeros((nTotal,1))
    vhh=np.concatenate((vhh,hh),axis=0)
    vGG=np.concatenate((vGG,GG),axis=0)
    vv,pa,vu=np.zeros((nDIM,1)),np.zeros((nDIM,1)),np.zeros((nDIM,1))
    QQ,qM,qu,qMu=np.zeros((nDIM,nDIM)),0.0,0.0,0.0
    qGrad,Delta,lbd,zidx=np.zeros(nDIM),1.0,0.0,last
    if(nDIM==16):##Mises Poly6
        vm=np.array([1, -3, 6, -7, 6, -3, 1, 9, -18, 27, -18, 9, 27, -27, 27, 27]).reshape((nDIM,1))
    elif(nDIM==25):##Mises Poly8
        vm=np.array([1, -4, 10, -16, 19, -16, 10, -4, 1, 12, -36, 72, -84, 72, -36, 12, 54, -108, 162, -108, 54, 108, -108, 108, 81]).reshape((nDIM,1))
    elif(nDIM==36):##Mises Poly10
        vm=np.array([1, -5, 15, -30, 45, -51, 45, -30, 15, -5, 1, 15, -60, 150, -240, 285, -240, 150, -60, 15, 90, -270, 540, -630, 540, -270, 90, 270,
                     -540, 810, -540, 270, 405, -405, 405, 243]).reshape((nDIM,1))
    elif(nDIM==49): ###Mises Poly12
        vm=np.array([1, -6, 21, -50, 90, -126, 141, -126, 90, -50, 21, -6, 1, 18, -90, 270, -540, 810, -918, 810, -540, 270, -90, 18, 135, -540, 1350, 
                     -2160, 2565, -2160, 1350, -540, 135, 540, -1620, 3240, -3780, 3240, -1620, 540, 1215, -2430, 3645, -2430, 1215, 1458, -1458, 1458, 
                     729]).reshape((nDIM,1))
    elif(nDIM==9):##Mises Poly4
        vm=np.array([1, -2, 3, -2, 1, 6, -6, 6, 9]).reshape((nDIM,1))                     
    else:print('func:genLinConstraints: PolyN degree must be one of {4,6,8,10,12}.\nCalculations aborted');exit()
    vu[:,0]=va[:,0]-vm[:,0]
    for kk in range(nTotal):
        QQ[:,:]=vQQ[kk,:,:]
        vv[:,0]=np.matmul(QQ,vm)[:,0]
        #qM=np.dot(vm[:,0],vv[:,0])
        #if(qM<0):print('(exit from)func:genLinConstraints:negative Mises quad (epsCVX too large): qM = ', qM);exit()
        qM=np.dot(vm[:,0],vv[:,0])
        qMu=np.dot(vu[:,0],vv[:,0])
        qu=np.dot(vu[:,0],np.matmul(QQ,vu)[:,0])
        #Delta=qMu*qMu-qM*qu
        #if(Delta<0):print('func:genLinConstraints:Negative Delta: Delta,qM,qu,qMu',Delta,qM,qu,qMu);exit()
        Delta=np.sqrt(qMu*qMu-qM*qu)
        #if(qu>=0):print('qu:',qu);exit()
        #lbd=(Delta+qMu)/(-qu)
        lbd=eLambda*((Delta+qMu)/(-qu))
        ###lbd=min((Delta+qMu)/(-qu),(Delta-qMu)/qu)
        if(lbd<0 or lbd>1):print('func:genLinConstraints:lambda NOK: lbd,qM,qu,qMu\nCalculations aborted',lbd,qM,qu,qMu);exit()
        pa[:,0]=vm[:,0]+lbd*vu[:,0]
        qGrad[:]=2.0*np.matmul(QQ,pa)[:,0]
        qu=np.dot(qGrad,pa[:,0])
        for j in range(vFidx.shape[0]):
            qu-=vFcoeff[j]*qGrad[vFidx[j]]
        vhh[zidx,0]=-qu
        #vhh[zidx,0]=-qu #
        ii=0
        for j in range(nDIM):
            if(j in vFidx):continue
            vGG[zidx,ii]=-qGrad[j];ii+=1
        zidx+=1
    #vGG[zidx,:]=-vu[:,0]
    #vhh[zidx,0]=np.dot(vm[:,0],vGG[zidx,:])
    #print(QQ.shape,qq.shape,cc)
    return vhh,vGG

def gen2BiaxConstraints(vMonoms,nOmega=1001):
    degree=vMonoms['nQ'];oDegree=1.0-1.0/degree;ndim=degree+1
    QQ,vv,q0=np.zeros((nOmega,ndim,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,1))
    vQ=vMonoms['vQ']
    vD1,vD2=vMonoms['vD1Q'],vMonoms['vD2Q']
    vC1,vC2=vMonoms['vC1Q'],vMonoms['vC2Q']
    vH11,vH12=vMonoms['vH11Q'],vMonoms['vH12Q']
    vH22,vH23=vMonoms['vH22Q'],vMonoms['vH23Q']
    vCH11,vCH12,vCH22=vMonoms['vCH11Q'],vMonoms['vCH12Q'],vMonoms['vCH22Q']
    D1,D2=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    H11,H12,H22=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    vAlpha,vBeta,vGamma=np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim)),np.zeros((nOmega,ndim))
    R11,R21,R12,R22=0.0,0.0,0.0,0.0
    vsx,vsy,vsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    vdsx,vdsy,vdsxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    #vd2sx,vd2sy,vd2sxy=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    v11,v12,v22=np.zeros(nOmega),np.zeros(nOmega),np.zeros(nOmega)
    dOmega=np.pi/nOmega
    vOmega=np.array([k*dOmega for k in range(nOmega)])
    comega,somega=np.cos(vOmega),np.sin(vOmega)
    sq2=1.0/np.sqrt(2.0)
    R11,R21,R12,R22=sq2,sq2,-sq2,sq2
    vsx[:],vsy[:]=R11*comega+R12*somega,R21*comega+R22*somega
    vdsx[:],vdsy[:]=-R11*somega+R12*comega,-R21*somega+R22*comega
    v11[:],v12[:],v22[:]=vdsx*vdsx,vdsx*vdsy,vdsy*vdsy
    for j in range(ndim):
        vAlpha[:,j]=(vsx**vQ[j][0])*(vsy**vQ[j][1])
        D1[:,j]=vC1[j]*(vsx**vD1[j][0])*(vsy**vD1[j][1])
        D2[:,j]=vC2[j]*(vsx**vD2[j][0])*(vsy**vD2[j][1])
        vBeta[:,j]=D1[:,j]*vdsx+D2[:,j]*vdsy
        H11[:,j]=vCH11[j]*(vsx**vH11[j][0])*(vsy**vH11[j][1])
        H12[:,j]=vCH12[j]*(vsx**vH12[j][0])*(vsy**vH12[j][1])
        H22[:,j]=vCH22[j]*(vsx**vH22[j][0])*(vsy**vH22[j][1])
        vGamma[:,j]=H11[:,j]*v11+H22[:,j]*v22+2.0*H12[:,j]*v12-D1[:,j]*vsx-D2[:,j]*vsy
    for ii in range(nOmega):
        QQ[ii,:,:]=degree*(vAlpha[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))-\
                   oDegree*(vBeta[ii,:].reshape((ndim,1))*vBeta[ii,:].reshape((1,ndim)))+\
                   0.5*(vAlpha[ii,:].reshape((ndim,1))*vGamma[ii,:].reshape((1,ndim))+vGamma[ii,:].reshape((ndim,1))*vAlpha[ii,:].reshape((1,ndim)))
    return nOmega,QQ



def gen2BiaxLinConstraints(vQQ,va,vFidx,vFcoeff,vhh,vGG,eLambda=0.9999):
    nTotal,degree=vQQ.shape[0:2] #;print(vQQ.shape)
    ndim=degree-vFidx.shape[0]
    last=vhh.shape[0]-1
    GG,hh=np.zeros((nTotal,ndim)),np.zeros((nTotal,1))
    vhh=np.concatenate((vhh,hh),axis=0)
    vGG=np.concatenate((vGG,GG),axis=0)
    vv,pa,vu=np.zeros((degree,1)),np.zeros((degree,1)),np.zeros((degree,1))
    QQ,qM,qu,qMu=np.zeros((degree,degree)),0.0,0.0,0.0
    qGrad,Delta,lbd,zidx=np.zeros(degree),1.0,0.0,last
    if(degree==7):
        vm=np.array([1,-3, 6, -7, 6, -3, 1]).reshape((degree,1))
    elif(degree==9):
        vm=np.array([1,-4, 10, -16, 19, -16, 10, -4, 1]).reshape((degree,1))
    elif(degree==11):
        vm=np.array([1, -5, 15, -30, 45, -51, 45, -30, 15, -5, 1]).reshape((degree,1))
    elif(degree==13):
        vm=np.array([1, -6, 21, -50, 90, -126, 141, -126, 90, -50, 21, -6, 1]).reshape((degree,1))
    else:print('func:genLinConstraints: PolyN degree must be one of {4,6,8,10,12}.\nCalculations aborted');exit()        
    vu[:,0]=va[:,0]-vm[:,0]
    for kk in range(nTotal):
        QQ[:,:]=vQQ[kk,:,:]
        vv[:,0]=np.matmul(QQ,vm)[:,0]
        qM=np.dot(vm[:,0],vv[:,0])
        qMu=np.dot(vu[:,0],vv[:,0])
        qu=np.dot(vu[:,0],np.matmul(QQ,vu)[:,0])
        #Delta=qMu*qMu-qM*qu
        #if(Delta<0):print('NEGATIVE DELTA');exit()
        Delta=np.sqrt(qMu*qMu-qM*qu)
        lbd=eLambda*((Delta+qMu)/(-qu))
        if(lbd<0 or lbd>1):print('lbd NOK: lbd,qM,qu,qMu',lbd,qM,qu,qMu);exit()
        pa[:,0]=vm[:,0]+lbd*vu[:,0]
        qGrad[:]=2.0*np.matmul(QQ,pa)[:,0]
        qu=np.dot(qGrad,pa[:,0])
        for j in range(vFidx.shape[0]):
            qu-=vFcoeff[j]*qGrad[vFidx[j]]
        vhh[zidx,0]=-qu
        ii=0
        for j in range(degree):
            if(j in vFidx):continue
            vGG[zidx,ii]=-qGrad[j];ii+=1
        zidx+=1
    #vGG[zidx,:]=-vu[:,0]
    #vhh[zidx,0]=np.dot(vm[:,0],vGG[zidx,:])
    #print(QQ.shape,qq.shape,cc)
    return vhh,vGG


def fCostBiax(data,vMonoms,vpoints,vFidx,vFcoeff,bxScale=1.0):
    degree,nFcoeff=vMonoms['nQ']+1,vFidx.shape[0];ndim=degree-nFcoeff
    MM=np.zeros((ndim,ndim))
    VV=np.zeros((1,ndim))
    vIdx=np.zeros(ndim,dtype=int)
    ii=0
    for k in range(degree):
        if(k in vFidx):continue
        vIdx[ii]=k;ii+=1
    #print(vFidx);print(vIdx);exit()
    di,vvi=0.0,np.zeros((1,degree))
    vQ=vMonoms['vQ']
    vD1,vD2=vMonoms['vD1Q'],vMonoms['vD2Q']#,vMonoms['vD3Q']
    vC1,vC2=vMonoms['vC1Q'],vMonoms['vC2Q']#,vMonoms['vC3Q']
    D1,D2=np.zeros(degree),np.zeros(degree)
    eps=1.0e-6
    for item in data['data']:####theta=item[2] must me 0 or 90 (for biax)
        if((eps<item[2]<90-eps) or (item[2]<-eps) or (item[2]>90+eps)):continue
        if(-eps<item[2]<eps):##theta=0
            sx,sy,sxy=1.0,item[1],0.0
            if(item[0]):##r-value equation
                for k in range(0,degree):
                    D1[k]=vC1[k]*(sx**vD1[k][0])*(sy**vD1[k][1])#*(sxy**vD1[k][2])
                    D2[k]=vC2[k]*(sx**vD2[k][0])*(sy**vD2[k][1])#*(sxy**vD2[k][2])
                    #D3[k-1]=vC3[k]*(sx**vD3[k][0])*(sy**vD3[k][1])*(sxy**vD3[k][2])
                vvi[:]=item[3]*D1[:]+(1.0+item[3])*D2[:]
                #di=-item[3]*vC1[0]*(sx**vD1[0][0])
                for k in range(nFcoeff):
                    #print(k,vFidx,vFcoeff)
                    di-=vFcoeff[k]*vvi[0,vFidx[k]]
                MM[:,:]+=(item[4]*vvi[0,vIdx].reshape((ndim,1)))*vvi[0,vIdx].reshape((1,ndim))
                VV[:]+=item[4]*di*vvi[0,vIdx].reshape((1,ndim))
                di=0.0
            else:
                for k in range(0,degree):
                    vvi[0,k]=(sx**vQ[k][0])*(sy**vQ[k][1])#*(sxy**vQ[k][2])
                di=1.0/(item[3]**(degree-1))
                for k in range(nFcoeff):
                    di-=vFcoeff[k]*vvi[0,vFidx[k]]
                MM[:,:]+=(item[4]*vvi[0,vIdx].reshape((ndim,1)))*vvi[0,vIdx].reshape((1,ndim))
                VV[:]+=item[4]*di*vvi[0,vIdx].reshape((1,ndim))
                di=0.0
    #Bezier5YS points
    MM2=np.zeros((ndim,ndim))
    VV2=np.zeros((1,ndim))
    npoints=vpoints.shape[1]-1
    for ii in range(1,npoints):
        sx,sy,rho=vpoints[0,ii],vpoints[1,ii],1.0
        if(sx>=sy):rho=sx;sy=sy/sx;sx=1.0
        else:rho=sy;sx=sx/sy;sy=1.0
        for k in range(0,degree):
            vvi[0,k]=(sx**vQ[k][0])*(sy**vQ[k][1])#*(vsxy[ii]**vQ[k][2])
        di=1/rho**(degree-1)
        for k in range(nFcoeff):
            di-=vFcoeff[k]*vvi[0,vFidx[k]]
        MM2[:,:]+=vvi[0,vIdx].reshape((ndim,1))*vvi[0,vIdx].reshape((1,ndim))
        VV2[:]+=di*vvi[0,vIdx].reshape((1,ndim))
    wv=(1.0-bxScale*data['wa'])/npoints
    return VV+wv*VV2,MM+wv*MM2


def fCost2(data,vMonoms,vFidx,vFcoeff):
    degree,nDIM,nFcoeff=vMonoms['nQ'],vMonoms['nDim'],vFidx.shape[0];ndim=nDIM-nFcoeff
    MM,VV=np.zeros((ndim,ndim)),np.zeros((1,ndim))
    vIdx=np.zeros(ndim,dtype=int)
    ii=0
    for k in range(nDIM):
        if(k in vFidx):continue
        vIdx[ii]=k;ii+=1
    #print('vIdx:',vIdx)
    di,vvi=0.0,np.zeros((1,nDIM))
    vQ=vMonoms['vQ']
    vD1,vD2,vD3=vMonoms['vD1Q'],vMonoms['vD2Q'],vMonoms['vD3Q']
    vC1,vC2,vC3=vMonoms['vC1Q'],vMonoms['vC2Q'],vMonoms['vC3Q']
    D1,D2,D3=np.zeros(nDIM),np.zeros(nDIM),np.zeros(nDIM)
    eps=1.0e-6
    vData,qHill=data['data'],False
    if(data['qHill']):vData,qHill=data['dataAct'],True
    #print(vData)
    for item in vData:
        if(-eps<item[2]<eps or 90-eps<item[2]<90+eps):continue ##skip biaxial data
        cphi,sphi=np.cos(item[2]),np.sin(item[2])
        cphi2,sphi2,csphi=cphi*cphi,sphi*sphi,cphi*sphi
        sx,sy,sxy=cphi2+item[1]*sphi2,sphi2+item[1]*cphi2,(1-item[1])*csphi
        if(item[0]):##r-value equation
            for k in range(0,nDIM):
                D1[k]=vC1[k]*(sx**vD1[k][0])*(sy**vD1[k][1])*(sxy**vD1[k][2])
                D2[k]=vC2[k]*(sx**vD2[k][0])*(sy**vD2[k][1])*(sxy**vD2[k][2])
                D3[k]=vC3[k]*(sx**vD3[k][0])*(sy**vD3[k][1])*(sxy**vD3[k][2])
            vvi[0,:]=(item[3]+sphi2)*D1[:]+(item[3]+cphi2)*D2[:]-csphi*D3[:]
            for k in range(nFcoeff):
                di-=vFcoeff[k]*vvi[0,vFidx[k]]
            MM[:,:]+=(item[4]*vvi[0,vIdx].reshape((ndim,1)))*vvi[0,vIdx].reshape((1,ndim))
            VV[:]+=item[4]*di*vvi[0,vIdx].reshape((1,ndim))
            di=0.0
        else:
            for k in range(0,nDIM):
                vvi[0,k]=(sx**vQ[k][0])*(sy**vQ[k][1])*(sxy**vQ[k][2])
            di=1.0/(item[3]**degree)
            for k in range(nFcoeff):
                di-=vFcoeff[k]*vvi[0,vFidx[k]]
            MM[:,:]+=(item[4]*vvi[0,vIdx].reshape((ndim,1)))*vvi[0,vIdx].reshape((1,ndim))
            VV[:]+=item[4]*di*vvi[0,vIdx].reshape((1,ndim))
            di=0.0
    #print(VV);print(MM)
    if(not qHill): return VV,MM
    else:print("(Note: Calibration using Hill'48 samples)\n")
    #Hill quadratic points
    a0,a2=1.0,1.0/(data['dataHill']['s90'])**2
    s45,r45=data['dataHill']['s45'],data['dataHill']['r45']
    a3=((r45+0.5)/(r45+1.0))*(2/s45)**2
    a1=(2/s45)**2-(a0+a2+a3)##;print('a0,a1,a2,a3:',a0,a1,a2,a3)
    nT=5+int(degree/12)
    #vt=np.linspace(np.pi/72,np.pi/24,nT)
    vt=np.linspace(np.pi/288,np.pi/32,nT)
    vct,vst=np.cos(vt),np.sin(vt)
    nOmega=np.array([2**(k+2) for k in range(nT)])
    nHill=np.sum(nOmega)
    vsx,vsy,vsxy,vrho=np.zeros(nHill),np.zeros(nHill),np.zeros(nHill),np.zeros(nHill)
    ii,deg2=0,degree/2
    for k in range(nT):
        vOmega=np.linspace(0,np.pi,nOmega[k])
        for omega in vOmega:
            vsx[ii],vsy[ii],vsxy[ii]=vst[k]*np.cos(omega),vst[k]*np.sin(omega),vct[k]
            vrho[ii]=(a0*vsx[ii]**2+a1*vsx[ii]*vsy[ii]+a2*vsy[ii]**2+a3*vsxy[ii]**2)**deg2
            ii+=1
    MM2=np.zeros((ndim,ndim))
    VV2=np.zeros((1,ndim))
    for ii in range(nHill):
        for k in range(0,nDIM):
            vvi[0,k]=(vsx[ii]**vQ[k][0])*(vsy[ii]**vQ[k][1])*(vsxy[ii]**vQ[k][2])
        di=vrho[ii]
        for k in range(nFcoeff):
            di-=vFcoeff[k]*vvi[0,vFidx[k]]
        MM2[:,:]+=vvi[0,vIdx].reshape((ndim,1))*vvi[0,vIdx].reshape((1,ndim))
        VV2[:]+=di*vvi[0,vIdx].reshape((1,ndim))
    wH=(1.0-data['wAct'])/nHill
    return VV+wH*VV2,MM+wH*MM2
    


def fCoeff(ndim,vCoeff,vFidx,vFcoeff):
    vv=np.zeros(ndim)
    ii,flag=0,False
    for k in range(ndim):
        flag=False
        for j in range(vFidx.shape[0]):
            if(k==vFidx[j]):
                vv[k],flag=vFcoeff[j],True;break
        if(flag):continue
        vv[k]=vCoeff[ii,0];ii+=1
    return vv

def polyBiaxOptim(data,vMonoms,vpoints,vFidx,vFcoeff,bxscale=1.0,epsTol=1.0e-7):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 200
    degree,nFcoeff=vMonoms['nQ']+1,vFidx.shape[0];ndim=degree-nFcoeff
    VV,MM=fCostBiax(data,vMonoms,vpoints,vFidx,vFcoeff,bxScale=bxscale)
    MM,VV=cvxopt.matrix(MM),cvxopt.matrix(-VV.T)
    vSol=cvxopt.solvers.qp(MM,VV)
    vCoeff=np.array([k for k in vSol['x']]).reshape((ndim,1))
    vvv=fCoeff(degree,vCoeff,vFidx,vFcoeff)
    nConst,zQQ=gen2BiaxConstraints(vMonoms)
    ##print('zQQ.shape = ',zQQ.shape)
    hh,GG=np.zeros((1,1)),np.zeros((1,ndim))
    iz=0
    cvx=np.zeros(nConst)
    vZero=np.zeros((ndim,1));vZero[:,0]=vCoeff[:,0]
    while(1):
        for kx in range(nConst):
            cvx[kx]=np.dot(vvv,np.matmul(zQQ[kx],vvv.reshape((degree,1)))[:,0])
        mcvx=np.min(cvx)
        if(mcvx>=0):break
        idx=np.argwhere(cvx<0.0)[:,0]
        hh,GG=gen2BiaxLinConstraints(zQQ[idx],vvv.reshape((degree,1)),vFidx,vFcoeff,hh,GG)
        print('Biax: {}-----'.format(iz),hh.shape,GG.shape,mcvx)
        vSol=cvxopt.solvers.qp(MM,VV,cvxopt.matrix(GG),cvxopt.matrix(hh))
        vCoeff=np.array([k for k in vSol['x']]).reshape((ndim,1))
        vvv=fCoeff(degree,vCoeff,vFidx,vFcoeff)
        #print('vCoeff:',vCoeff[:,0])
        if(np.sqrt(np.sum((vZero-vCoeff)**2,axis=0))<epsTol):
            print('Total negative constraints = ',np.sum(cvx<0))
            print('epsToll exit at iz = ',iz);break
        vZero[:,0]=vCoeff[:,0]
        iz+=1
        if(iz>200):
            print('(MaxIter)Total negative constraints = ',np.sum(cvx<0))
            break
    ##print('Total negative constraints = ',np.sum(cvx<0))
    ##print('Min cvx = ',mcvx)
    return fCoeff(degree,vCoeff,vFidx,vFcoeff)


def readVcoeff(fName,fPath=figDirData):
    print('file: '+fName)
    try:ff=open(fPath+fName,'r')
    except IOError as err:print(err);exit()
    if('FEdata' in fName):
        print('reading FEdata file')
        vx=ff.readlines()
        vx=[float(xx.strip()) for xx in vx]
        degree,nCoeff=int(vx[5]),int(vx[6])
        #print('degree: ',degree)
        #for kk in range(7,nCoeff+7):
        #    print(kk-7,': ', vx[kk])
        vcf=vx[5:]            
    elif('Err_and_Coeff' in fName):
        print('reading Err_and_Coeff file')
        for line in ff:
            if ('Poly_N degree:' in line):
                degree=int(line.strip().split(':')[1])
            if('--Polynomial coefficients' in line):break
        vCoeff=[]
        for line in ff:
            vCoeff.append(float(line.strip()))    
        vCoeff=np.array(vCoeff)
        #print('degree: ',degree)    
        #print(vCoeff.shape)
        #for k in range(vCoeff.shape[0]):
        #    print(k,': ',vCoeff[k])    
        vcf=PolyNparam(degree,vCoeff,matProp=[],matName='',fPrint=False)
    else:print('Unknown data file');exit()    
    ff.close()
    return degree,vcf
    
def testValYF(dd,vcf):
    qq,th=dd[0],(dd[1]/180.0)*np.pi
    ct,st=np.cos(th),np.sin(th)
    ct2,st2,ctst=ct*ct,st*st,ct*st
    yf,[gx,gy,gxy]=fGYF(ct2+qq*st2,qq*ct2+st2,(1.0-qq)*ctst,vcf)
    rv=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
    yf=1.0/yf
    print('degree: ',int(vcf[0]))
    print('Yield stress: data, prediction')
    print('{:.4f}, {:.4f}'.format(dd[2],yf))
    print('R-value: data, prediction')
    print('{:.4f}, {:.4f}\n'.format(dd[3],rv))
    return yf,rv    



def polyOptim(data,vMonoms,nConst,zQQ,bxScale=1.0,epsTol=1.0e-9):
    print('Running polyOptim.......') #;  print(data['dataAct']);exit()
    vv,vFidx,vFcoeff=biaxData(data)
    vFcoeff=polyBiaxOptim(data,vMonoms,vv,vFidx,vFcoeff,bxscale=bxScale)#;print('vFcoeff:',vFcoeff);exit()
    degree=vMonoms['nQ']
    vFidx=np.array([k for k in range(degree+1)],dtype=int)
    if(degree==6):eLambda=1.0
    if(degree==8):eLambda=0.9992
    if(degree==10):eLambda=0.998
    else:eLambda=0.995
    nDIM,ndim=vMonoms['nDim'],vMonoms['nDim']-vFidx.shape[0]
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 200
    VV,MM=fCost2(data,vMonoms,vFidx,vFcoeff)
    MM+=1.0e-12*np.eye(MM.shape[0])
    MM,VV=cvxopt.matrix(MM),cvxopt.matrix(-VV.T)
    vSol=cvxopt.solvers.qp(MM,VV)
    vCoeff=np.array([k for k in vSol['x']]).reshape((ndim,1))
    vvv=fCoeff(nDIM,vCoeff,vFidx,vFcoeff)
    ####nConst,zQQ=genConstraints(vMonoms)
    ##nConst,zQQ=genConstraintsTG(vMonoms,pNegKG)
    ####print('zQQ.shape = ',zQQ.shape)
    hh,GG=np.zeros((1,1)),np.zeros((1,ndim))
    iz=0
    cvx=np.zeros(nConst)
    vZero=np.zeros((ndim,1));vZero[:,0]=vCoeff[:,0]
    while(1):
        for kx in range(nConst):
            #cvx[kx]=np.dot(vvv,np.matmul(zQQ[kx],vvv.reshape((nDIM,1)))[:,0])-epsCVX
            cvx[kx]=np.dot(vvv,np.matmul(zQQ[kx],vvv.reshape((nDIM,1)))[:,0])
        mcvx=np.min(cvx)
        if(mcvx>=0.0):break
        idx=np.argwhere(cvx<0.0)[:,0]
        if(idx.shape[0]>4*(10**5)):
            print('The number of violated constraints is too large')
            print('Convergence is unlikely to be achieved')
            print('Iterations terminated');exit()
        hh,GG=genLinConstraints(zQQ[idx],vvv.reshape((nDIM,1)),vFidx,vFcoeff,hh,GG,eLambda)
        print('{}-----'.format(iz),nConst,hh.shape,GG.shape,mcvx)
        vSol=cvxopt.solvers.qp(MM,VV,cvxopt.matrix(GG),cvxopt.matrix(hh))
        vCoeff=np.array([k for k in vSol['x']]).reshape((ndim,1))
        vvv=fCoeff(nDIM,vCoeff,vFidx,vFcoeff)
        #print('vCoeff:',vCoeff[:,0])
        if(np.sqrt(np.sum((vZero-vCoeff)**2,axis=0))<epsTol):
            print('Total negative constraints = ',np.sum(cvx<0))
            print('epsToll exit at iz = ',iz);break
        vZero[:,0]=vCoeff[:,0]
        iz+=1
        if(iz>200):
            print('(MaxIter)Total negative constraints = ',np.sum(cvx<0))
            break
    ##print('Total negative constraints = ',np.sum(cvx<0))
    print('Min cvx = ',mcvx)
    return fCoeff(nDIM,vCoeff,vFidx,vFcoeff)
    


def polyNOptim(degree,data,pdata,bxscale=1.0):
    if(type(degree)!=int):print('Degree must be an integer\nCalculations aborted');exit()
    if(degree%2):print('Degree must be an even integer\nCalculations aborted');exit()
    if(bxscale>1.0 or bxscale<0.0):print('bxscale must be in [0.0,1.0]\nCalculations aborted');exit()
    zdd=vPoly(degree)
    ###Optimization section: calculate PolyN parameters
    print('-----------Calculating parameters....')   
    #print(data['dataAct']);exit()
    nConst,zQQ=genConstraintsTG(zdd)
    vCoeff=polyOptim(data,zdd,nConst,zQQ,bxScale=bxscale)
    ###Generate predictions report and input file for FE (use 'fPrint=True')
    ###Use 'export = True' to generate data from PolyN model
    print('-----------Postprocessing....')  
    vcf=PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
    minKG,negKG=PolyN_predictions(vcf,pdata,export=True)
    if(negKG):##
        nOpt=12 ##Max number of additional overall optimization iterations 
        vepsNeg,vnDir=0.002,64
        dnegKG_old=[item for item in negKG] 
        kOpt=0       
        while(kOpt<nOpt):
            print('-----------Calculating parameters({})....'.format(kOpt+1))
            dnConst,dzQQ=genConstraintsTG(zdd,negKG)
            nConst+=dnConst
            zQQ=np.concatenate((zQQ,dzQQ),axis=0)
            vCoeff=polyOptim(data,zdd,nConst,zQQ,bxScale=bxscale)   
            print('-----------Postprocessing({})....'.format(kOpt+1))  
            vcf=PolyNparam(degree,vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
            minKG,dnegKG=PolyN_predictions(vcf,pdata,export=True)
            if(not dnegKG):break
            eqFlag,neq=0,len(dnegKG)
            if(neq==len(dnegKG_old)):
                for kk in range(neq):
                    if(not(dnegKG[kk][0]==dnegKG_old[kk][0])):break
                    if(not(dnegKG[kk][1]==dnegKG_old[kk][1])):break
                    if(not(dnegKG[kk][2]==dnegKG_old[kk][2])):break
                    eqFlag+=1
            if(eqFlag==neq): 
                print('-------------------------------------------------------')            
                print('Negative KG points are the same: iterations terminated')
                break
            dnegKG_old=[item for item in dnegKG]     
            for kk in range(len(dnegKG)):
                point=dnegKG[kk]
                print(kk,': ',point)
                vtdir=genTanDir2D(point,vnDir)
                for vr in vtdir:
                    x,y,z=point[0]+vepsNeg*vr[0],point[1]+vepsNeg*vr[1],point[2]+vepsNeg*vr[2]
                    md=np.sqrt(x*x+y*y+z*z)
                    dnegKG.append((x/md,y/md,z/md,0))                     
            negKG+=dnegKG
            kOpt+=1
    return vcf
    
def genConstraintsPoints2DOptPoints(nPoints):
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1[1:]);vs1=np.sin(vt1[1:])
    vN2=np.int_(nPoints*vs1)     
    vP=np.zeros((np.sum(vN2)+1,3))
    vN2[:]+=1
    vP[0,2]=1.0
    jj=1
    for kk in range(0,len(vt1)-1):
        ct1=vc1[kk];st1=vs1[kk]
        N2=vN2[kk];jN2=jj+N2-1
        vt2=np.linspace(0,np.pi,N2)
        vc2=np.cos(vt2[0:N2-1]);vs2=np.sin(vt2[0:N2-1])
        vP[jj:jN2,0]=st1*vc2
        vP[jj:jN2,1]=st1*vs2
        vP[jj:jN2,2]=ct1
        jj=jN2
    return vP       
    
def PolyN_GaussCheck(vcf,nRandom=10**5):
    np.random.seed(99)     
    vuRandom=np.random.normal(0,1,(nRandom,3)) ##nPoints on 2D unit sphere of 3D
    vNorm=np.sqrt(np.sum(vuRandom**2,axis=1))
    vuRandom[:,0]/=vNorm;vuRandom[:,1]/=vNorm;vuRandom[:,2]/=vNorm
    vu=genConstraintsPoints2DOptPoints(nPoints=300)
    vu=np.concatenate((vu[:,0:3],vuRandom),axis=0)
    nPoints=vu.shape[0]
    zyf,G1,G2,G3,H11,H12,H13,H22,H23,H33=0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    md1,md2,md3,mdPoints=1.0e+9,1.0e+9,1.0e+9,[]    
    mm,mmPoints,KG=1.0e+009,[],0.0
    for k in range(nPoints):
        zyf,[G1,G2,G3],[H11,H12,H13,H22,H23,H33]=HGYF(vu[k,0],vu[k,1],vu[k,2],vcf)
        md1=min(md1,H11)
        md2=min(md2,H11*H22-H12*H12)
        md3=min(md3,H33*(H11*H22-H12*H12)+H13*(H12*H23-H22*H13)+H23*(H12*H13-H11*H23))
        if(md1<0 or md2<0 or md3<0):mdPoints.append(k)
        #H11s=H22*H33-H23*H23
        #H12s=H23*H13-H12*H33
        #H13s=H12*H23-H22*H13
        #H22s=H11*H33-H13*H13
        #H23s=H12*H13-H11*H23
        #H33s=H11*H22-H12*H12
        zyf*=zyf
        H11s=zyf*(H22*H33-H23*H23)
        H12s=zyf*(H23*H13-H12*H33)
        H13s=zyf*(H12*H23-H22*H13)
        H22s=zyf*(H11*H33-H13*H13)
        H23s=zyf*(H12*H13-H11*H23)
        H33s=zyf*(H11*H22-H12*H12)
        vNorm2=G1**2+G2**2+G3**2
        vNorm=np.sqrt(vNorm2)
        G1/=vNorm;G2/=vNorm;G3/=vNorm
        KG=(H11s*G1*G1+H22s*G2*G2+H33s*G3*G3+2.0*(H12s*G1*G2+H13s*G1*G3+H23s*G2*G3))/vNorm2
        mm=min(mm,KG)
        if(KG<0):mmPoints.append(np.array((vu[k,0],vu[k,1],vu[k,2],KG)))
    return mm,mmPoints,[md1,md2,md3],nPoints  

def PolyN_GaussCheck2(vcf):
    np.random.seed(99)     
    vuRandom=np.random.normal(0,1,(3500,3)) ##nPoints on 2D unit sphere of 3D
    vNorm=np.sqrt(np.sum(vuRandom**2,axis=1))
    vuRandom[:,0]/=vNorm;vuRandom[:,1]/=vNorm;vuRandom[:,2]/=vNorm
    vu=genConstraintsPoints2DOptPoints(nPoints=200)
    vu=np.concatenate((vu[:,0:3],vuRandom),axis=0)
    nPoints=vu.shape[0]
    zyf,G1,G2,G3,H11,H12,H13,H22,H23,H33=0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    md1,md2,md3,mdPoints=1.0e+9,1.0e+9,1.0e+9,[]    
    mm,mmPoints,KG=1.0e+009,[],0.0
    for k in range(nPoints):
        zyf=fYF(vu[k,0],vu[k,1],vu[k,2],vcf)
        zyf,[G1,G2,G3],[H11,H12,H13,H22,H23,H33]=HGYF(vu[k,0]/zyf,vu[k,1]/zyf,vu[k,2]/zyf,vcf)
        md1=min(md1,H11)
        md2=min(md2,H11*H22-H12*H12)
        md3=min(md3,H33*(H11*H22-H12*H12)+H13*(H12*H23-H22*H13)+H23*(H12*H13-H11*H23))
        if(md1<0 or md2<0 or md3<0):mdPoints.append(k)
        H11s=H22*H33-H23*H23
        H12s=H23*H13-H12*H33
        H13s=H12*H23-H22*H13
        H22s=H11*H33-H13*H13
        H23s=H12*H13-H11*H23
        H33s=H11*H22-H12*H12
        vNorm2=G1**2+G2**2+G3**2
        vNorm=np.sqrt(vNorm2)
        G1/=vNorm;G2/=vNorm;G3/=vNorm
        KG=(H11s*G1*G1+H22s*G2*G2+H33s*G3*G3+2.0*(H12s*G1*G2+H13s*G1*G3+H23s*G2*G3))/vNorm2
        mm=min(mm,KG)
        if(mm<0):mmPoints.append(k)
    return mm,mmPoints,[md1,md2,md3],nPoints     

###---------------------------------------------------------------------------------------------------------

'''

'''
def PolyNparam(deg,vvvCoeff,matProp=[],matName='',fPrint=False):
    nMon=int(deg/2)
    NN=nMon+1
    nCoeff=NN*NN
    vcf=np.zeros(2+12*nCoeff)
    vcf[0]=deg
    mxV=max([abs(xx) for xx in vvvCoeff])
    if(mxV<zroTol):print('All parameters are = 0\n Computations aborted');exit()
    if(False):
        xvcf[1]=mxV
        xvCoeff=[xx/mxV for xx in vvvCoeff]
    vcf[1]=nCoeff
    vCoeff=[xx for xx in vvvCoeff]    
    d1,d2=np.zeros(nCoeff),np.zeros(nCoeff)
    d11,d22,d12=np.zeros(nCoeff),np.zeros(nCoeff),np.zeros(nCoeff)
    kcf=0
    for powerXY in range(nMon):
        kdeg=deg-2*powerXY
        for jj in range(kdeg+1):
            ##calculate derivatives 
            d1[kcf]=(kdeg-jj)*vCoeff[kcf]
            d2[kcf]=jj*vCoeff[kcf]
            d11[kcf]=(kdeg-jj-1)*d1[kcf]
            d22[kcf]=(jj-1)*d2[kcf]
            d12[kcf]=jj*d1[kcf]
            kcf+=1
        ##shift backwards d/dSigmaY    
        d2[kcf-kdeg-1:kcf-1]=d2[kcf-kdeg:kcf]
        d2[kcf-1]=zro         
        ##shift backwards d2/(dSigmaY,dSigmaY)    
        d22[kcf-kdeg-1:kcf-2]=d22[kcf-kdeg+1:kcf]
        d22[kcf-1],d22[kcf-2]=zro,zro
        ##shift backwards d2/(dSigmaX,dSigmaY)    
        d12[kcf-kdeg-1:kcf-2]=d12[kcf-kdeg:kcf-1]
        d12[kcf-1],d12[kcf-2]=zro,zro
    if(1): ##make pure zero some awkward -0.0        
        for kcf in range(nCoeff):
            if(abs(d1[kcf])<zroTol): d1[kcf]=zro
            if(abs(d2[kcf])<zroTol): d2[kcf]=zro
            if(abs(d11[kcf])<zroTol): d11[kcf]=zro
            if(abs(d22[kcf])<zroTol): d22[kcf]=zro
            if(abs(d12[kcf])<zroTol): d12[kcf]=zro
    ###setup the reversed sets of coefficients 
    ###use indexing starting from 1 
    ###(to avoid the problem of vec[-1] which Python interprets as last position in vec)
    nOne=nCoeff+1
    zvcf=np.zeros(nOne)
    zd1,zd2=np.zeros(nOne),np.zeros(nOne)
    zd11,zd22,zd12=np.zeros(nOne),np.zeros(nOne),np.zeros(nOne)
    zvcf[1:],zd1[1:],zd2[1:],zd11[1:],zd22[1:],zd12[1:]=vCoeff,d1,d2,d11,d22,d12
    pvcf=np.zeros(nCoeff)
    pd1,pd2=np.zeros(nCoeff),np.zeros(nCoeff)
    pd11,pd22,pd12=np.zeros(nCoeff),np.zeros(nCoeff),np.zeros(nCoeff)
    kcf=1 ##index starts from 1
    for ii in range(0,nMon):
        kdeg=deg-2*ii
        kcf2=kcf+kdeg+1
        pvcf[kcf-1:kcf2-1:1]=zvcf[kcf2-1:kcf-1:-1] ##if kcf were 0 =>>Problem
        pd1[kcf-1:kcf2-2:1]=zd1[kcf2-2:kcf-1:-1]
        pd2[kcf-1:kcf2-2:1]=zd2[kcf2-2:kcf-1:-1]
        pd11[kcf-1:kcf2-3:1]=zd11[kcf2-3:kcf-1:-1]
        pd22[kcf-1:kcf2-3:1]=zd22[kcf2-3:kcf-1:-1]
        pd12[kcf-1:kcf2-3:1]=zd12[kcf2-3:kcf-1:-1]
        kcf=kcf2
    pvcf[-1]=vCoeff[-1]        
    vcf[2:]=np.concatenate((vCoeff,d1,d2,d11,d22,d12,pvcf,pd1,pd2,pd11,pd22,pd12),axis=0)    
    #return d1,d2,d11,d22,d12,pvcf,pd1,pd2,pd11,pd22,pd12
    if(fPrint):
        fName=matName+'_deg'+str(deg)+'_FEdata.txt'
        ff=open(figDirData+fName,'w')
        for item in matProp:
            ff.write('{:.12f}\n'.format(item))
        for item in vcf:
            ff.write('{:.12f}\n'.format(item))
        ff.close()        
    return vcf



'''
fYF: Not vectorized yet
'''
def fYF(sx,sy,sxy,vcf):
    deg=int(vcf[0])
    nMon=int(deg/2)
    nCoeff=int(vcf[1])
    KKo=nCoeff
    KKo-=1
    msx,msy=abs(sx),abs(sy)
    if(msx+msy<gTol):
        #KKo+=2
        #tau=abs(sxy)
        #yf=vcf[KKo]*tau       
        return (vcf[KKo+2]**(1.0/float(deg)))*abs(sxy)    
    if(msx>=msy):
        tt=sx;rho=sy/sx;gamm=sxy/sx;gamm=gamm*gamm
        #kcfA=2;kcfB=kcfA+nCoeff 
        cc=vcf[2:nCoeff+2]
    else: 
        tt=sy;rho=sx/sy;gamm=sxy/sy;gamm=gamm*gamm
        kcfA=2+6*nCoeff;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]  
    yf=cc[KKo]
    jj,jj2=nMon,nMon
    while(jj>0):
        KKo-=1
        bb=cc[KKo]
        jj-=1
        mm=deg-2*jj ##degree of multiplying polynomial in sx and sy
        while(mm>0):
            mm-=1
            KKo-=1
            bb=cc[KKo]+bb*rho
        yf=bb+yf*gamm  
    #zyf=(vcf[1]*yf)**odeg  ##in case of coefficients normalization 
    #zyf=yf**odeg
    #zyf=att*zyf    
    return abs(tt)*(yf**(1.0/float(deg)))


'''
fGYF: Not vectorized yet
'''
def fGYF(sx,sy,sxy,vcf):
    deg=int(vcf[0])
    odeg=1.0/float(deg)
    nMon=int(deg/2)
    nCoeff=int(vcf[1])
    KKo=nCoeff
    KKo-=1
    msx,msy=abs(sx),abs(sy)
    if(msx+msy<gTol):
        KKo+=2
        tau=abs(sxy)
        yf=(vcf[KKo]**odeg)*tau 
        tau=sxy/tau ##( = +/- 1)
        ##gyf=[d/dsx,d/dsy,d/dsxy]     
        return yf,[0.0,0.0,tau]
    if(msx>=msy):
        tt=sx;rho=sy/sx;gamm=sxy/sx;gamm=gamm*gamm
        kcfA=2;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
    else: 
        tt=sy;rho=sx/sy;gamm=sxy/sy;gamm=gamm*gamm
        kcfA=2+6*nCoeff;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
    d3=np.zeros(nMon)   
    yf,g1yf,g2yf=cc[KKo],d1[KKo],d2[KKo]
    jj,jj2=nMon,nMon
    while(jj>0):
        KKo-=1
        bb=cc[KKo]
        d1bb,d2bb=d1[KKo],d2[KKo]
        jj-=1
        mm=deg-2*jj ##degree of multiplying polynomial in sx and sy
        while(mm>0):
            mm-=1
            KKo-=1
            bb=cc[KKo]+bb*rho
            d1bb=d1[KKo]+d1bb*rho
            d2bb=d2[KKo]+d2bb*rho           
        d3[jj-1]=jj*bb
        yf=bb+yf*gamm
        g1yf=d1bb+g1yf*gamm
        g2yf=d2bb+g2yf*gamm  
    #zyf=(vcf[1]*yf)**odeg  ##in case of coefficients normalization 
    zyf=yf**odeg
    yval,att,y2val=zyf/(deg*yf),abs(tt),float(deg-1)/zyf
    sgntt=(tt/att)*yval
    g1yf=sgntt*g1yf
    g2yf=sgntt*g2yf  
    jj2-=1    
    g3yf=float(nMon)*cc[-1]
    while(jj2>0):
        jj2-=1
        g3yf=d3[jj2]+g3yf*gamm
    tt=(2.0*sxy)/tt    
    g3yf=sgntt*tt*g3yf
    zyf=att*zyf    
    return zyf,[g1yf,g2yf,g3yf]


'''
HGYF: Not vectorized yet
'''
def HGYF(sx,sy,sxy,vcf):
    deg=int(vcf[0])
    odeg=1.0/float(deg)
    pdeg=1.0-odeg
    nMon=int(deg/2)
    nCoeff=int(vcf[1])
    KKo=nCoeff
    KKo-=1
    msx,msy=abs(sx),abs(sy)
    if(msx+msy<gTol):
        KKo+=2
        tau=abs(sxy)
        yf=(vcf[KKo]**odeg)*tau 
        tau=sxy/tau ##( = +/- 1)
        ##gyf=[d/dsx,d/dsy,d/dsxy]
        gyf=[0.0,0.0,tau]
        ##print('deg,vcf[KKo],pdeg,tau: ', deg,vcf[KKo],pdeg,tau)
        tau=deg*(vcf[KKo]**pdeg)*tau
        ##print('tau: ', tau)
        ##hyf=[d2/dsxdsx,d2/dsxdsy,d2/dsxdsxy,d2/dsydsy,d2/dsydsxy,d2/dsxydsxy]
        hyf=[(2.0*vcf[KKo-3])/tau,vcf[KKo-2]/tau,0.0,(2.0*vcf[KKo-1])/tau,0.0,0.0]        
        return yf,gyf,hyf   
    if(msx>=msy):
        tt=sx;rho=sy/sx;gamm=sxy/sx;gamm=gamm*gamm
        kcfA=2;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d11=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d22=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d12=vcf[kcfA:kcfB]
    else: 
        tt=sy;rho=sx/sy;gamm=sxy/sy;gamm=gamm*gamm
        kcfA=2+6*nCoeff;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d11=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d22=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d12=vcf[kcfA:kcfB]
    d3=np.zeros(nMon)
    d13,d23,d33=np.zeros(nMon),np.zeros(nMon),np.zeros(nMon)     
    yf,g1yf,g2yf=cc[KKo],d1[KKo],d2[KKo]
    h11,h12,h22=d11[KKo],d12[KKo],d22[KKo]
    jj,jj2=nMon,nMon
    while(jj>0):
        KKo-=1
        bb=cc[KKo]
        d1bb,d2bb,d11bb,d22bb,d12bb=d1[KKo],d2[KKo],d11[KKo],d22[KKo],d12[KKo]
        jj-=1
        mm=deg-2*jj ##degree of multiplying polynomial in sx and sy
        while(mm>0):
            mm-=1
            KKo-=1
            bb=cc[KKo]+bb*rho
            d1bb=d1[KKo]+d1bb*rho
            d2bb=d2[KKo]+d2bb*rho
            d11bb=d11[KKo]+d11bb*rho
            d22bb=d22[KKo]+d22bb*rho
            d12bb=d12[KKo]+d12bb*rho            
        d3[jj-1]=jj*bb
        d13[jj-1]=jj*d1bb
        d23[jj-1]=jj*d2bb
        d33[jj-1]=jj*bb*(2*jj-1)
        yf=bb+yf*gamm
        g1yf=d1bb+g1yf*gamm
        g2yf=d2bb+g2yf*gamm
        h11=d11bb+h11*gamm
        h22=d22bb+h22*gamm
        h12=d12bb+h12*gamm   
    #zyf=(vcf[1]*yf)**odeg  ##in case of coefficients normalization 
    zyf=yf**odeg
    yval,att,y2val=zyf/(deg*yf),abs(tt),float(deg-1)/zyf
    sgntt=(tt/att)*yval
    g1yf=sgntt*g1yf
    g2yf=sgntt*g2yf
    h11=(yval*h11-y2val*g1yf*g1yf)/att
    h22=(yval*h22-y2val*g2yf*g2yf)/att
    h12=(yval*h12-y2val*g1yf*g2yf)/att    
    jj2-=1    
    g3yf=float(nMon)*cc[-1]
    h13,h23,h33=d13[-1],d23[-1],float(deg-1)*g3yf
    while(jj2>0):
        jj2-=1
        g3yf=d3[jj2]+g3yf*gamm
        h13=d13[jj2]+h13*gamm
        h23=d23[jj2]+h23*gamm
        h33=d33[jj2]+h33*gamm
    tt=(2.0*sxy)/tt    
    g3yf=sgntt*tt*g3yf
    tt=yval*tt
    h13=(tt*h13-y2val*g1yf*g3yf)/att
    h23=(tt*h23-y2val*g2yf*g3yf)/att
    h33=(2.0*yval*h33-y2val*g3yf*g3yf)/att
    zyf=att*zyf    
    return zyf,[g1yf,g2yf,g3yf],[h11,h12,h13,h22,h23,h33]






def plotPoly(vCoeff,data,vMonoms={},saveFig=True):
    degree=data['degree'];oDeg=1.0/degree
    if(not vMonoms):
        vMonoms=vPoly(degree)
    NN=nMonoms(vMonoms['nQ'])
    Nq=11
    vq=np.linspace(0,-0.4,Nq) ###; print(vq)
    vLW=np.ones(Nq);vLW[0]=2
    Nphi=101
    vphi=np.linspace(0,np.pi/2,Nphi);vDegs=(180/np.pi)*vphi
    vcos=np.cos(vphi);vsin=np.sin(vphi)
    vcos2,vsin2,vsc=vcos*vcos,vsin*vsin,vcos*vsin
    vsx,vsy,vsxy=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    vS,vR, mxS,mnS = np.zeros(Nphi),np.zeros(Nphi), -1.0, 100.0
    vQ=vMonoms['vQ']
    vD1=vMonoms['vD1Q'];vD2=vMonoms['vD2Q'];vD3=vMonoms['vD3Q']
    vC1=vMonoms['vC1Q'];vC2=vMonoms['vC2Q'];vC3=vMonoms['vC3Q']
    vDX,vDY,vDXY=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    fig,fig2=plt.figure(),plt.figure()
    ax,ax2=fig.add_subplot(),fig2.add_subplot()
    for kq in range(Nq):
        q=vq[kq]
        vsx[:]=vcos2+q*vsin2
        vsy[:]=vsin2+q*vcos2
        vsxy[:]=(1.0-q)*vsc
        for k in range(NN):
            vS[:]+=vCoeff[k]*(vsx**vQ[k][0])*(vsy**vQ[k][1])*(vsxy**vQ[k][2])
            vDX[:]+=vCoeff[k]*vC1[k]*(vsx**vD1[k][0])*(vsy**vD1[k][1])*(vsxy**vD1[k][2])
            vDY[:]+=vCoeff[k]*vC2[k]*(vsx**vD2[k][0])*(vsy**vD2[k][1])*(vsxy**vD2[k][2])
            vDXY[:]+=vCoeff[k]*vC3[k]*(vsx**vD3[k][0])*(vsy**vD3[k][1])*(vsxy**vD3[k][2])
        vR[:]=(vDXY*vsc-vDX*vsin2-vDY*vcos2)/(vDX+vDY)
        vS[:]=1.0/vS**oDeg;mxS=max(max(vS),mxS);mnS=min(min(vS),mnS)
        ax.plot(vDegs,vR,linewidth=vLW[kq],color='k')
        ax2.plot(vDegs,vS,linewidth=vLW[kq],color='k')
        vS[:],vDX[:],vDY[:],vDXY[:]=0.0,0.0,0.0,0.0
    ax.plot(data['thetaR'],data['rValue'],linestyle='',marker='o',markersize=6,markerfacecolor='b')
    majorXticks=[0,15,30,45,60,75,90]
    minorXticks=np.linspace(0,90,4*6+1)
    ax.set_xticks(minorXticks,minor=True)
    ax.set_xticks(majorXticks,minor=False)
    ax.set_xticklabels([r'$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'])
    #ax.grid()
    yMin=int(min(data['rValue']))
    yMax=int(max(data['rValue']))+1
    majorYticks=np.linspace(yMin,yMax,yMax-yMin+1)
    #minorYticks=np.linspace(yMin,yMax,4*(yMax-yMin)+1)
    #ax.set_yticks(majorYticks,minor=False)
    #ax.set_yticks(minorYticks,minor=True)
    #ax.set_yticklabels([str(round(k,2)) for k in minorYticks],minor=True,fontsize=10)
    #ax.set_yticklabels([str(round(k,2)) for k in ax.get_yticks(minor=False)],minor=False,fontsize=10)
    #ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.65)
    ax.set_xlim(-0.85,90.85)
    ax.set_ylim(yMin-0.0275,0.9*yMax)
    fName=data['fullName']+'_P{}'.format(str(data['degree']))
    if(saveFig):
        fig.savefig(figDirPlot+fName+'_rValue.png',bbox_inches='tight',dpi=300)
    ax2.plot(data['thetaS'],data['sigma'],linestyle='',marker='o',markersize=6,markerfacecolor='b')
    x0,x1=-0.25,90.25
    ax2.set_xlim(x0,x1)
    y0,y1=ax2.get_ylim()
    ax2.text(x1-0.05*(x1-x0),y0+0.01*(y1-y0),r'$\theta$',fontsize=16)
    ax2.text(x0+0.01*(x1-x0),y1-0.05*(y1-y0),r'$\overline{\sigma}(q,\theta)$',fontsize=16)
    ax2.set_xticks([0,15,30,45,60,75,90])
    ax2.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'])
    ax2.grid()
    if(saveFig):
        fig2.savefig(figDirPlot+fName+'_Sigma.png',bbox_inches='tight',dpi=300)    
    fig.clf()
    ax=fig.add_subplot()
    NX,NY,maxX,maxY=201,201,1.25,1.25
    vsx=np.linspace(-maxX,maxX,NX)
    vsy=np.linspace(-maxY,maxY,NY)
    maxSXY=1.0/(vCoeff[-1]**(1/degree))
    #vsxy=[0.0]
    vsxy=maxSXY*np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.975,0.99,0.998])
    nIdx=int(degree/2+1)
    vTerms=np.zeros((nIdx,NY,NX))
    for kidx in range(nIdx):
        vPowers=vMonoms['vQidx'][kidx][1]#; print(vPowers)
        #nPowers=len(vPowers)   
        for ky in range(NY):
            for i in vPowers:
                vTerms[kidx,ky,:]+=vCoeff[i]*(vsx**(vMonoms['vQ'][i][0]))*(vsy[ky]**(vMonoms['vQ'][i][1]))
    for sxy in vsxy:
        sxy2=sxy*sxy
        vz=vTerms[nIdx-1]
        for i in range(nIdx-2,-1,-1):
            vz=vTerms[i]+sxy2*vz#;print(sxy2,i)
        #vz=vz**(1.0/degree)
        ax.contour(vsx,vsy,vz,levels=[1.0],linewidths=1,colors=['k'])
        #vz2=vTerms[0]**(1.0/degree);print(np.max(vz-vz2),np.min(vz-vz2))
        #ax.contour(vsx,vsy,vz2,levels=[1.0],linewidths=1,colors=['r'])
    ax.grid()
    ax.set_aspect('equal')
    x0,x1=ax.get_xlim()
    y0,y1=ax.get_ylim()
    ax.text(x1-0.075*(x1-x0),y0+0.02*(y1-y0),r'$\overline{\sigma}_x$',fontsize=14)
    ax.text(x0+0.01*(x1-x0),y1-0.075*(y1-y0),r'$\overline{\sigma}_y$',fontsize=14)
    if(saveFig):
        fig.savefig(figDirPlot+fName+'_sxySections.png',bbox_inches='tight',dpi=300)
    else:
        plt.show()    



def plotPoly22(vcf,data,saveFig=True,actualPts=False,zaxes=False,zfg=None,zfg2=None,zfg3=None,zstyle=False,preName=''):
    print('Generating {} plots...'.format(data['fullName']))
    degree,NN=int(vcf[0]),int(vcf[1])
    oDeg=1.0/degree
    #NN=nMonoms(vMonoms['nQ'])
    Nq=6
    vq=np.linspace(0,-0.2,Nq) ###; print(vq)
    vLW=np.ones(Nq);vLW[0]=2.5
    Nphi=101
    vphi=np.linspace(0,np.pi/2,Nphi);vDegs=180.0*(vphi/np.pi)
    vcos=np.cos(vphi);vsin=np.sin(vphi)
    vcos2,vsin2,vsc=vcos*vcos,vsin*vsin,vcos*vsin
    #vsx,vsy,vsxy=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    vS,vR, mxS,mnS = np.zeros(Nphi),np.zeros(Nphi), -1.0, 100.0
    if(zaxes):
        if((not zfg) or (not zfg2) or (not zfg3)):print('missing axes args to plot func: abort');exit()
        fig,fig2,fig3=zfg,zfg2,zfg3
        ax,ax2,ax3=zfg.axes[0],zfg2.axes[0],zfg3.axes[0]
    else:
        fig,fig2,fig3=plt.figure(),plt.figure(),plt.figure()
        ax,ax2,ax3=fig.add_subplot(),fig2.add_subplot(),fig3.add_subplot()
    if(zstyle):linestyle='--'
    else:linestyle='-'        
    for kq in range(Nq):
        q=vq[kq]
        for k in range(Nphi):
            vsx=vcos2[k]+q*vsin2[k]
            vsy=vsin2[k]+q*vcos2[k]
            vsxy=(1.0-q)*vsc[k]
            yf,[vDX,vDY,vDXY]=fGYF(vsx,vsy,vsxy,vcf)            
            vR[k]=(vDXY*vsc[k]-vDX*vsin2[k]-vDY*vcos2[k])/(vDX+vDY)
            vS[k]=1.0/yf
        ax.plot(vDegs,vR,linestyle=linestyle,linewidth=vLW[kq],color='k')
        ax2.plot(vDegs,vS,linestyle=linestyle,linewidth=vLW[kq],color='k')
    thetaTicks=[0,15,30,45,60,75,90]
    if(not actualPts):
        for k in range(len(data['thetaR'])):
            if(int(round(data['thetaR'][k],1)) in thetaTicks):    
                ax.plot(data['thetaR'][k],data['rValue'][k],linestyle='',marker='o',markersize=8,
                markerfacecolor='b',markeredgecolor='b')
    else:
        ax.plot(data['thetaR'],data['rValue'],linestyle='',marker='o',markersize=8,
                markerfacecolor='b',markeredgecolor='b')    
    #thetaTicks=[0,15,30,45,60,75,90]
    thetaLabels=['0','15','30','45','60','75','90']
    ax.set_xticks(thetaTicks,minor=False)
    #ax.set_xticklabels(thetaLabels,fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    x0,x1=-1,91
    ax.set_xlim(x0,x1)
    y0,y1=ax.get_ylim()
    ax.text(x1-0.05*(x1-x0),y0+0.01*(y1-y0),r'$\theta$',fontsize=16)
    ax.text(x0+0.01*(x1-x0),y1-0.05*(y1-y0),r'$r(q,\theta)$',fontsize=16)
    ax.grid(True)
    fName=preName+data['fullName']+'_P{}'.format(str(degree))
    if(saveFig):
        fig.savefig(figDirPlot+fName+'_rValue.png',bbox_inches='tight',dpi=300)
    if(not actualPts):    
        for k in range(len(data['thetaS'])):
            if(int(round(data['thetaS'][k],1)) in thetaTicks):    
                ax2.plot(data['thetaS'][k],data['sigma'][k],linestyle='',marker='o',markersize=8,
                markerfacecolor='b',markeredgecolor='b')
    else:
        ax2.plot(data['thetaS'],data['sigma'],linestyle='',marker='o',markersize=8,
                markerfacecolor='b',markeredgecolor='b')    
    x0,x1=-1,91
    ax2.set_xlim(x0,x1)
    y0,y1=ax2.get_ylim()
    ax2.text(x1-0.05*(x1-x0),y0+0.01*(y1-y0),r'$\theta$',fontsize=16)
    ax2.text(x0+0.01*(x1-x0),y1-0.15*(y1-y0),r'$\overline{\sigma}(q,\theta)$',fontsize=15)
    ax2.set_xticks(thetaTicks)
    #ax2.set_xticklabels(thetaLabels,fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.grid(True)
    if(saveFig):
        fig2.savefig(figDirPlot+fName+'_Sigma.png',bbox_inches='tight',dpi=300)    
    #fig3=plt.figure()
    #ax3=fig3.add_subplot()
    NX,NY,maxX,maxY=201,201,1.25,1.25
    vsx=np.linspace(-maxX,maxX,NX)
    vsy=np.linspace(-maxY,maxY,NY)
    vz=np.zeros((NY,NX))
    maxSXY=1.0/(vcf[1+NN]**(1/degree))
    #vsxy=[0.0]
    vsxy=maxSXY*np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.975,0.99,0.998])
    for sxy in vsxy:
        for k in range(NY):
            sy=vsy[k]
            for j in range(NX):
                vz[k,j]=fYF(vsx[j],sy,sxy,vcf)
        ax3.contour(vsx,vsy,vz,levels=[1.0],linestyles=linestyle,linewidths=1,colors=['k'])
        #vz2=vTerms[0]**(1.0/degree);print(np.max(vz-vz2),np.min(vz-vz2))
        #ax.contour(vsx,vsy,vz2,levels=[1.0],linewidths=1,colors=['r'])
    ax3.grid(True)
    ax3.set_aspect('equal')
    x0,x1=ax3.get_xlim()
    y0,y1=ax3.get_ylim()
    ax3.text(x1-0.075*(x1-x0),y0+0.02*(y1-y0),r'$\overline{\sigma}_x$',fontsize=14)
    ax3.text(x0+0.01*(x1-x0),y1-0.075*(y1-y0),r'$\overline{\sigma}_y$',fontsize=14)
    ax3.tick_params(axis='both', which='major', labelsize=12)
    if(saveFig):
        fig3.savefig(figDirPlot+fName+'_sxySections.png',bbox_inches='tight',dpi=300)
    #else:
    #    plt.show()
    return fig,fig2,fig3


def plotCmPoly():
    pass    
        

def funcYld2000_2D(vAlpha,degree,sx,sy,txy):
    deg2=degree/2.0
    alpha=[1]+list(vAlpha)
    L11,L12,L21,L22,L33=2.0*alpha[1]/3.0,-alpha[1]/3.0,-alpha[2]/3.0,2.0*alpha[2]/3.0,alpha[7]
    S11=L11*sx+L12*sy
    S22=L21*sx+L22*sy
    S12=L33*txy
    Ione=(S11+S22)/2; Itwo=S11*S22-S12*S12
    pp=(2.0**(degree-1))*((Ione*Ione-Itwo)**deg2)
    #Ione=(S11+S22); Itwo=S11*S22-S12*S12
    #pp=0.5*((Ione*Ione-4*Itwo)**deg2)
    L11=2.0*(-alpha[3]+alpha[4]+4.0*alpha[5]-alpha[6])/9.0
    L12=(alpha[3]-4.0*alpha[4]-4.0*alpha[5]+4.0*alpha[6])/9.0
    L21=(4.0*alpha[3]-4.0*alpha[4]-4.0*alpha[5]+alpha[6])/9.0
    L22=2.0*(-alpha[3]+4.0*alpha[4]+alpha[5]-alpha[6])/9.0
    L33=alpha[8]
    S11=L11*sx+L12*sy
    S22=L21*sx+L22*sy
    S12=L33*txy
    Ione=(S11+S22)/2.0; Itwo=S11*S22-S12*S12
    tta=10*Ione*Ione-Itwo;ttb=6*Ione*np.sqrt(Ione*Ione-Itwo)
    pp2=0.5*((tta+ttb)**deg2+(tta-ttb)**deg2)
    return (pp+pp2)**(1.0/degree)
 
def gfuncYld2000_2D(vAlpha,degree,sx,sy,txy):
    deg2=degree/2.0
    alpha=[1]+list(vAlpha)
    L11,L12,L21,L22,L33=2.0*alpha[1]/3.0,-alpha[1]/3.0,-alpha[2]/3.0,2.0*alpha[2]/3.0,alpha[7]
    S11=L11*sx+L12*sy
    S22=L21*sx+L22*sy
    S12=L33*txy
    Ione=(S11+S22)/2; Itwo=S11*S22-S12*S12
    pp=(2.0**(degree-1))*((Ione*Ione-Itwo)**deg2)
    dI1dx=0.5*(L11+L21);dI1dy=0.5*(L12+L22);dI1dxy=0.0
    dI2dx=L11*S22+L21*S11;dI2dy=L12*S22+L22*S11;dI2dxy=-2.0*L33*S12
    tt=2**degree*(Ione*Ione-Itwo)**(deg2-1)
    df1dx,df1dy,df1dxy=tt*(2*Ione*dI1dx-dI2dx),tt*(2*Ione*dI1dy-dI2dy),tt*(2*Ione*dI1dxy-dI2dxy)
    L11=2.0*(-alpha[3]+alpha[4]+4.0*alpha[5]-alpha[6])/9.0
    L12=(alpha[3]-4.0*alpha[4]-4.0*alpha[5]+4.0*alpha[6])/9.0
    L21=(4.0*alpha[3]-4.0*alpha[4]-4.0*alpha[5]+alpha[6])/9.0
    L22=2.0*(-alpha[3]+4.0*alpha[4]+alpha[5]-alpha[6])/9.0
    L33=alpha[8]
    S11=L11*sx+L12*sy
    S22=L21*sx+L22*sy
    S12=L33*txy
    Ione=(S11+S22)/2.0; Itwo=S11*S22-S12*S12
    tt=np.sqrt(Ione*Ione-Itwo)
    tta=10*Ione*Ione-Itwo;ttb=6*Ione*tt
    pp2=0.5*((tta+ttb)**deg2+(tta-ttb)**deg2)
    yf=(pp+pp2)**(1.0/degree)
    dI1dx=0.5*(L11+L21);dI1dy=0.5*(L12+L22);dI1dxy=0.0
    dI2dx=L11*S22+L21*S11;dI2dy=L12*S22+L22*S11;dI2dxy=-2.0*L33*S12
    df2d1=(tta+ttb)**(deg2-1)*(20*Ione+6*tt+(6*Ione*Ione)/tt)+(tta-ttb)**(deg2-1)*(20*Ione-6*tt-(6*Ione*Ione)/tt)
    df2d2=(tta+ttb)**(deg2-1)*(-1-(3*Ione)/tt)+(tta-ttb)**(deg2-1)*(-1+(3*Ione)/tt)
    df2dx,df2dy,df2dxy=df2d1*dI1dx+df2d2*dI2dx,df2d1*dI1dy+df2d2*dI2dy,df2d1*dI1dxy+df2d2*dI2dxy
    return yf, [df1dx+df2dx,df1dy+df2dy,df1dxy+df2dxy] 
    
    
def plotYld2000(degree,vAlpha,data,saveFig=True,dataPlot=True,export=True,uniaxA=False,actualPts=False,zaxes=False,zfg=None,zfg2=None,zfg3=None,zstyle=False,preName=''):
    print('Generating {} plots...'.format(data['fullName']))
    Nq=6
    vq=np.linspace(0,-0.2,Nq) ###; print(vq)
    vLW=np.ones(Nq);vLW[0]=2.5
    Nphi=101
    vphi=np.linspace(0,np.pi/2,Nphi);vDegs=180.0*(vphi/np.pi)
    vcos=np.cos(vphi);vsin=np.sin(vphi)
    vcos2,vsin2,vsc=vcos*vcos,vsin*vsin,vcos*vsin
    #vsx,vsy,vsxy=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    vS,vR, mxS,mnS = np.zeros(Nphi),np.zeros(Nphi), -1.0, 100.0
    if(zaxes):
        if((not zfg) or (not zfg2) or (not zfg3)):print('missing axes args to plot func: abort');exit()
        fig,fig2,fig3=zfg,zfg2,zfg3
        ax,ax2,ax3=zfg.axes[0],zfg2.axes[0],zfg3.axes[0]
    else:
        fig,fig2,fig3=plt.figure(),plt.figure(),plt.figure()
        ax,ax2,ax3=fig.add_subplot(),fig2.add_subplot(),fig3.add_subplot()
    if(zstyle):linestyle='--'
    else:linestyle='-'        
    for kq in range(Nq):
        q=vq[kq]
        for k in range(Nphi):
            vsx=vcos2[k]+q*vsin2[k]
            vsy=vsin2[k]+q*vcos2[k]
            vsxy=(1.0-q)*vsc[k]
            yf,[vDX,vDY,vDXY]=gfuncYld2000_2D(vAlpha,degree,vsx,vsy,vsxy)
            #yf,[vDX,vDY,vDXY]=fGYF(vsx,vsy,vsxy,vcf)            
            vR[k]=(vDXY*vsc[k]-vDX*vsin2[k]-vDY*vcos2[k])/(vDX+vDY)
            vS[k]=1.0/yf
        ax.plot(vDegs,vR,linestyle=linestyle,linewidth=vLW[kq],color='k')
        ax2.plot(vDegs,vS,linestyle=linestyle,linewidth=vLW[kq],color='k')
    thetaTicks=[0,15,30,45,60,75,90]
    if(dataPlot):
        if(not actualPts):
            for k in range(len(data['thetaR'])):
                if(int(round(data['thetaR'][k],1)) in thetaTicks):    
                    ax.plot(data['thetaR'][k],data['rValue'][k],linestyle='',marker='o',markersize=8,
                    markerfacecolor='b',markeredgecolor='b')
        else:
            ax.plot(data['thetaR'],data['rValue'],linestyle='',marker='o',markersize=8,
                    markerfacecolor='b',markeredgecolor='b')    
        #thetaTicks=[0,15,30,45,60,75,90]
        thetaLabels=['0','15','30','45','60','75','90']
        ax.set_xticks(thetaTicks,minor=False)
        #ax.set_xticklabels(thetaLabels,fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)
        x0,x1=-1,91
        ax.set_xlim(x0,x1)
        y0,y1=ax.get_ylim()
        ax.text(x1-0.05*(x1-x0),y0+0.01*(y1-y0),r'$\theta$',fontsize=16)
        ax.text(x0+0.01*(x1-x0),y1-0.05*(y1-y0),r'$r(q,\theta)$',fontsize=16)
        ax.grid(True)
    fName=preName+data['fullName']+'_P{}'.format(str(degree))
    if(saveFig):
        fig.savefig(figDirPlot+fName+'_rValue.png',bbox_inches='tight',dpi=300)
    if(dataPlot):    
        if(not actualPts):    
            for k in range(len(data['thetaS'])):
                if(int(round(data['thetaS'][k],1)) in thetaTicks):    
                    ax2.plot(data['thetaS'][k],data['sigma'][k],linestyle='',marker='o',markersize=8,
                    markerfacecolor='b',markeredgecolor='b')
            if(0):
                ff=open(r'.\figsPolyN\DATA\AA6016T4\Yld2000_Gent_Sigma_Fig58.txt')
                vvx,vvy=[],[]
                for line in ff:
                    line=line.strip()
                    if(line):line=line.split(',')
                    vvx.append(float(line[0]));vvy.append(float(line[1]))
                ff.close()
                ax2.plot(vvx,vvy,linestyle='',marker='+',markersize=8,
                    markerfacecolor='r',markeredgecolor='k')        
        else:
            ax2.plot(data['thetaS'],data['sigma'],linestyle='',marker='o',markersize=8,
                    markerfacecolor='b',markeredgecolor='b')    
        x0,x1=-1,91
        ax2.set_xlim(x0,x1)
        y0,y1=ax2.get_ylim()
        ax2.text(x1-0.05*(x1-x0),y0+0.01*(y1-y0),r'$\theta$',fontsize=16)
        ax2.text(x0+0.01*(x1-x0),y1-0.15*(y1-y0),r'$\overline{\sigma}(q,\theta)$',fontsize=15)
        ax2.set_xticks(thetaTicks)
        #ax2.set_xticklabels(thetaLabels,fontsize=12)
        ax2.tick_params(axis='both', which='major', labelsize=12)
        ax2.grid(True)
    if(saveFig):
        fig2.savefig(figDirPlot+fName+'_Sigma.png',bbox_inches='tight',dpi=300)    
    #fig3=plt.figure()
    #ax3=fig3.add_subplot()
    NX,NY,maxX,maxY=201,201,1.25,1.25
    vsx=np.linspace(-maxX,maxX,NX)
    vsy=np.linspace(-maxY,maxY,NY)
    vz=np.zeros((NY,NX))
    maxSXY=(2.0/((2*vAlpha[6])**degree+2*vAlpha[7]**degree))**(1/degree)
    #vsxy=[0.0]
    vsxy=maxSXY*np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.975,0.99,0.998])
    for sxy in vsxy:
        for k in range(NY):
            sy=vsy[k]
            for j in range(NX):
                vz[k,j]=funcYld2000_2D(vAlpha,degree,vsx[j],sy,sxy)
        ax3.contour(vsx,vsy,vz,levels=[1.0],linestyles=linestyle,linewidths=1,colors=['k'])
        #vz2=vTerms[0]**(1.0/degree);print(np.max(vz-vz2),np.min(vz-vz2))
        #ax.contour(vsx,vsy,vz2,levels=[1.0],linewidths=1,colors=['r'])
    ax3.grid(True)
    ax3.set_aspect('equal')
    x0,x1=ax3.get_xlim()
    y0,y1=ax3.get_ylim()
    ax3.text(x1-0.075*(x1-x0),y0+0.02*(y1-y0),r'$\overline{\sigma}_x$',fontsize=14)
    ax3.text(x0+0.01*(x1-x0),y1-0.075*(y1-y0),r'$\overline{\sigma}_y$',fontsize=14)
    ax3.tick_params(axis='both', which='major', labelsize=12)
    if(saveFig):
        fig3.savefig(figDirPlot+fName+'_sxySections.png',bbox_inches='tight',dpi=300)
    #else:
    #    plt.show()
    if(not export):return fig,fig2,fig3
    print('Exporting data....')
    trdata,tsdata=np.array(data['thetaR']),np.array(data['thetaS'])
    nt,nq=25,15
    ##if(degree>9):nt,nq=49,36
    vq=np.linspace(0.0,-0.35,nq)
    dtheta=(0.5*np.pi)/(nt-1)
    vt=np.array([k*dtheta for k in range(0,nt)])
    vthet=180.0*(np.array(vt)/np.pi)
    vrdata,vsdata=np.zeros((nq,nt)),np.zeros((nq,nt))
    if(uniaxA):##type of uniaxial data: 'a' or 'v'
        dstype=['a' for jj in range(nt)]
    else:
        dstype=['v' for jj in range(nt)]    
    for jj in range(nt):###uniaxial data 
        ct,st=np.cos(vt[jj]),np.sin(vt[jj])
        ct2,st2,ctst=ct*ct,st*st,ct*st
        thet=vthet[jj]
        ath=np.abs(tsdata-thet)
        jk=np.argmin(ath) ##; ##print('---idx: jj, jk: ', jj,jk)
        if(ath[jk]<1.0e-6):
            #vsdata[0,jj]=data['sigma'][jk]
            vsdata[0,jj]=1.0/funcYld2000_2D(vAlpha,degree,ct2,st2,ctst)
            dstype[jj]='a'
        else:    
            vsdata[0,jj]=1.0/funcYld2000_2D(vAlpha,degree,ct2,st2,ctst)
            ###dstype.append('v')
        ath=np.abs(trdata-thet)
        jk=np.argmin(ath)
        if(ath[jk]<1.0e-6):
            #vrdata[0,jj]=data['rValue'][jk]
            yf,[gx,gy,gxy]=gfuncYld2000_2D(vAlpha,degree,ct2,st2,ctst)        
            vrdata[0,jj]=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
        else:
            yf,[gx,gy,gxy]=gfuncYld2000_2D(vAlpha,degree,ct2,st2,ctst)        
            vrdata[0,jj]=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
    for kk in range(1,nq):
        qq=vq[kk]
        for jj in range(nt):
            ct,st=np.cos(vt[jj]),np.sin(vt[jj])
            ct2,st2,ctst=ct*ct,st*st,ct*st
            yf,[gx,gy,gxy]=gfuncYld2000_2D(vAlpha,degree,ct2+qq*st2,qq*ct2+st2,(1-qq)*ctst)
            vsdata[kk,jj]=1.0/yf            
            vrdata[kk,jj]=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
    ff=open(figDirData+data['fullName']+'_deg'+str(degree)+'_Export.txt','w')
    ff.write('### Exported data from PolyN fit with N = '+str(degree)+'\n')
    ff.write('### Can be used as input data for higher order PolyN (default degree: N+2)\n')
    ff.write('name= '+data['fullName']+'_ExptP'+str(degree)+'\n')
    ff.write('degree= {}\n'.format(degree+2))
    ff.write('### FEA section---------\n')
    ff.write('EE= {}\n'.format(data['FEdata'][0]))
    ff.write('NU= {}\n'.format(data['FEdata'][1]))
    ff.write('AA= {}\n'.format(data['FEdata'][2]))
    ff.write('BB= {}\n'.format(data['FEdata'][3]))
    ff.write('CC= {}\n'.format(data['FEdata'][4]))
    ff.write('### Uniaxial data-----------\n')
    for jj in range(nt):##write uniaxial data
        ff.write('d= 0.0, {:.2f}, {:.4f}, {:.4f}, {}\n'.format(vthet[jj],vsdata[0,jj],vrdata[0,jj],dstype[jj]))
    ff.write('### Balanced-Biaxial data------\n')
    vvdd=[(1.0,0.0),(0.5,0.0),(0.5,np.pi/2)]
    for item in vvdd:
        qq,th=item[:]
        ct,st=np.cos(th),np.sin(th)
        ct2,st2,ctst=ct*ct,st*st,ct*st
        yf,[gx,gy,gxy]=gfuncYld2000_2D(vAlpha,degree,ct2+qq*st2,qq*ct2+st2,(1-qq)*ctst)
        vss=1.0/yf            
        vrr=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
        ff.write('d= {:.4f}, {:.2f}, {:.4f}, {:.4f}, a\n'.format(item[0],180.0*(item[1]/np.pi),vss,vrr))
    #####ff.write('d= 1.0, 0.0, {:.4f}, {:.4f}, a\n'.format(data['sbiax'],data['rbiax']))
    if(data['pstr']):
        ff.write('### Plane strain data -----------\n')
        for item in data['pstr']:
            ff.write('d= {:.4f}, {:.2f}, {:.4f}, 0.0, a\n'.format(item[0],180.0*(item[1]/np.pi),item[2]))
    ff.write('### Biaxial r-value data (drawing zone)-------\n')
    for kk in range(1,nq):
        for jj in range(nt):
            ff.write('d= {:.4f}, {:.2f}, {:.4f}, {:.4f}, v\n'.format(vq[kk],vthet[jj],vsdata[kk,jj],vrdata[kk,jj]))
    ff.write('### Weights--------------\n')
    ff.write('wa= 0.975\n')
    ff.write('war= 0.1\n')
    ff.write('wvr= 0.07\n')
    ff.write('### Hill48 samples (no)--------\n')
    ff.write('qHill= 0')    
    ff.close()    
    return fig,fig2,fig3    
    

    
def PolyN_predictions(vcf,pdata,export=False,uniaxA=False):
    print('Caclulating predictions: Err_and_Coeff file...')
    degree,name=pdata['degree'],pdata['fullName']
    ff=open(figDirData+name+'_deg'+str(degree)+'_Err_and_Coeff.txt','w')
    ff.write(name+'\n')
    ff.write('Poly_N degree: '+str(degree)+'\n')
    err=0.0
    ff.write('------Directional stresses: Theta, Data, Predicted\n')
    for kk in range(len(pdata['thetaS'])):
        th=np.pi*(pdata['thetaS'][kk]/180.0)
        ct,st=np.cos(th),np.sin(th)
        yf=1.0/fYF(ct*ct,st*st,ct*st,vcf)
        ff.write('{:4.1f}, {:.4f}, {:.4f}\n'.format(pdata['thetaS'][kk],pdata['sigma'][kk],yf))
        err+=(pdata['sigma'][kk]-yf)*(pdata['sigma'][kk]-yf)
    ff.write('------Balanced biaxial stress: Data, Predicted\n')
    yf=1.0/fYF(1.0,1.0,0.0,vcf)
    ff.write('{:.4f}, {:.4f}\n'.format(pdata['sbiax'],yf))
    err+=(pdata['sbiax']-yf)*(pdata['sbiax']-yf)
    ff.write('---All stresses square root error: {:.4f}\n'.format(np.sqrt(err)))
    err=0.0    
    ff.write('------Directional r-Values: Theta, Data, Predicted\n')
    for kk in range(len(pdata['thetaR'])):
        th=np.pi*(pdata['thetaR'][kk]/180.0)
        ct,st=np.cos(th),np.sin(th)
        ct2,st2,ctst=ct*ct,st*st,ct*st
        yf,[gx,gy,gxy]=fGYF(ct2,st2,ctst,vcf)
        rv=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
        ff.write('{:4.1f}, {:.4f}, {:.4f}\n'.format(pdata['thetaR'][kk],pdata['rValue'][kk],rv))    
        err+=(pdata['rValue'][kk]-rv)*(pdata['rValue'][kk]-rv)
    ff.write('------Balanced biaxial r-value: Data, Predicted\n')
    yf,[gx,gy,gxy]=fGYF(1.0,1.0,0.0,vcf)
    rv=-gy/(gx+gy)
    ff.write('{:.4f}, {:.4f}\n'.format(pdata['rbiax'],rv))
    err+=(pdata['rbiax']-rv)*(pdata['rbiax']-rv)   
    ff.write('---All r-values square root error: {:.4f}\n'.format(np.sqrt(err)))
    if(pdata['pstr']):
        ff.write('------Plane strain: qRatio, Theta, Data s_L, Predicted s_L, Data r-value, Predicted r-value\n')
        for item in pdata['pstr']:
            if(item[3]=='v'):continue
            qq,th=item[0],item[1]
            ct,st=np.cos(th),np.sin(th)
            ct2,st2,ctst=ct*ct,st*st,ct*st
            yf,[gx,gy,gxy]=fGYF(ct2+qq*st2,qq*ct2+st2,(1.0-qq)*ctst,vcf)
            rv=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
            ff.write('{:.4f}, {:.1f}, {:.4f}, {:.4f}, 0.0, {:.4f}\n'.format(qq,180.0*(th/np.pi),item[2],1.0/yf,rv))
    minKG,negKG,[md1,md2,md3],cvxPoints=PolyN_GaussCheck(vcf)
    print('min Gauss KG = ',minKG)
    print('min det1, det2, det3: {}, {}, {}'.format(md1,md2,md3))
    ff.write('------Convexity check with Hessian leading principal minors: Min(det_1), Min(det_2), Min(det_3)\n')
    ff.write('{}, {}, {}\n'.format(md1,md2,md3))
    ff.write('------Convexity check with Gaussian curvature: Min(KG) = {}\n'.format(minKG))
    if(minKG<0):
        msg='---Negative KG at {} out of {} points\n'.format((str(len(negKG))), str(cvxPoints))
        print(msg)
        ff.write(msg)
        ff.write('---Negative KG locations and values: x, y, z, KG\n')        
        for point in negKG:
            #print(point,'   norm = ', np.sqrt(np.sum(point[0:3]*point[0:3])))
            ff.write('{}, {}, {}, {}\n'.format(point[0],point[1],point[2],point[3]))
    ff.write('------Polynomial coefficients\n')
    for cc in vcf[2:int(vcf[1])+2]:
        ff.write('{}\n'.format(cc))         
    ff.close()
    if(not export):return minKG,negKG
    trdata,tsdata=np.array(pdata['thetaR']),np.array(pdata['thetaS'])
    nt,nq=25,15
    ##if(degree>9):nt,nq=49,36
    vq=np.linspace(0.0,-0.35,nq)
    dtheta=(0.5*np.pi)/(nt-1)
    vt=np.array([k*dtheta for k in range(0,nt)])
    vthet=180.0*(np.array(vt)/np.pi)
    vrdata,vsdata=np.zeros((nq,nt)),np.zeros((nq,nt))
    if(uniaxA):##type of uniaxial data: 'a' or 'v'
        dstype=['a' for jj in range(nt)]
    else:
        dstype=['v' for jj in range(nt)]    
    for jj in range(nt):###uniaxial data 
        ct,st=np.cos(vt[jj]),np.sin(vt[jj])
        ct2,st2,ctst=ct*ct,st*st,ct*st
        thet=vthet[jj]
        ath=np.abs(tsdata-thet)
        jk=np.argmin(ath) ##; ##print('---idx: jj, jk: ', jj,jk)
        if(ath[jk]<1.0e-6):
            vsdata[0,jj]=pdata['sigma'][jk]
            dstype[jj]='a'
        else:    
            vsdata[0,jj]=1.0/fYF(ct2,st2,ctst,vcf)
            ###dstype.append('v')
        ath=np.abs(trdata-thet)
        jk=np.argmin(ath)
        if(ath[jk]<1.0e-6):
            vrdata[0,jj]=pdata['rValue'][jk]
        else:
            yf,[gx,gy,gxy]=fGYF(ct2,st2,ctst,vcf)        
            vrdata[0,jj]=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
    for kk in range(1,nq):
        qq=vq[kk]
        for jj in range(nt):
            ct,st=np.cos(vt[jj]),np.sin(vt[jj])
            ct2,st2,ctst=ct*ct,st*st,ct*st
            yf,[gx,gy,gxy]=fGYF(ct2+qq*st2,qq*ct2+st2,(1-qq)*ctst,vcf)
            vsdata[kk,jj]=1.0/yf            
            vrdata[kk,jj]=(gxy*ctst-(gx*st2+gy*ct2))/(gx+gy)
    ff=open(figDirData+name+'_deg'+str(degree)+'_Export.txt','w')
    ff.write('### Exported data from PolyN fit with N = '+str(degree)+'\n')
    ff.write('### Can be used as input data for higher order PolyN (default degree: N+2)\n')
    ff.write('name= '+pdata['fullName']+'_ExptP'+str(degree)+'\n')
    ff.write('degree= {}\n'.format(degree+2))
    ff.write('### FEA section---------\n')
    ff.write('EE= {}\n'.format(pdata['FEdata'][0]))
    ff.write('NU= {}\n'.format(pdata['FEdata'][1]))
    ff.write('AA= {}\n'.format(pdata['FEdata'][2]))
    ff.write('BB= {}\n'.format(pdata['FEdata'][3]))
    ff.write('CC= {}\n'.format(pdata['FEdata'][4]))
    ff.write('### Uniaxial data-----------\n')
    for jj in range(nt):##write uniaxial data
        ff.write('d= 0.0, {:.2f}, {:.4f}, {:.4f}, {}\n'.format(vthet[jj],vsdata[0,jj],vrdata[0,jj],dstype[jj]))
    ff.write('### Balanced-Biaxial data------\n')
    ff.write('d= 1.0, 0.0, {:.4f}, {:.4f}, a\n'.format(pdata['sbiax'],pdata['rbiax']))
    if(pdata['pstr']):
        ff.write('### Plane strain data -----------\n')
        for item in pdata['pstr']:
            ff.write('d= {:.4f}, {:.2f}, {:.4f}, 0.0, a\n'.format(item[0],180.0*(item[1]/np.pi),item[2]))
    ff.write('### Biaxial r-value data (drawing zone)-------\n')
    for kk in range(1,nq):
        for jj in range(nt):
            ff.write('d= {:.4f}, {:.2f}, {:.4f}, {:.4f}, v\n'.format(vq[kk],vthet[jj],vsdata[kk,jj],vrdata[kk,jj]))
    ff.write('### Weights--------------\n')
    ff.write('wa= 0.975\n')
    ff.write('war= 0.1\n')
    ff.write('wvr= 0.07\n')
    ff.write('### Hill48 samples (no)--------\n')
    ff.write('qHill= 0')    
    ff.close()    
    return minKG,negKG





def readInputFile(fName,subDir,echo=False):
    if(osn=='nt'):drsep='\\'
    else:drsep='/'
    try:
        ff=open(figDirData+drsep+subDir+drsep+fName)
    except IOError as err:
        print(err);exit()
    ddata={'task':'', 'fileInpData':'','subDir':subDir,'alpha':[],'fileCVXData':'','degree':0}
    msg='Calculations aborted'   
    for line in ff:
        line=line.strip()
        if(line=='' or line[0]=='#'):continue
        if(':' in line):line=line.split(':');lnz=line[0].strip()
        else:print(fName+': Unknown format on line');print(line);print(msg);exit()
        if(lnz=='task'):
            tval=line[1].strip()
            if(tval in ['modelFit', 'convexCheck', 'testPoly']):
                ddata['task']=tval
            else:print('task: Unknown option');print(msg);exit()
        if(lnz=='fileInpData'):ddata['fileInpData']=line[1].strip()
        if(lnz=='alpha'):
            valpha=[]
            tval=line[1].strip().split(',')
            try:
                for item in tval:valpha.append(float(item))
            except ValueError as err:print(err);print(line[1]);print(msg);exit()
            ddata['alpha']=valpha
        if(lnz=='degree'):
            try:ddata['degree']=int(float(line[1].strip()))
            except ValueError as err:print(err);print(line[1]);print(msg);exit()            
        if(lnz=='fileCVXData'):
            ddata['fileCVXData']=line[1].strip()        
    if(echo):
        for key in ddata:print(key,': ',ddata[key])
    return ddata    

def mainPolyN(dData):
    task,inpData,subDir=dData['task'],dData['fileInpData'],dData['subDir']
    if(task=='modelFit'):
        tStart=time()
        ###Read mechanical data and other global parameters from text file  
        ###'data' is the whole data structure
        ###'pdata' contains info relevant only to plotting and reporting    
        data,pdata=readData(inpData,subDir)
        degree=data['degree']
        print('PolyN degree: ',degree)
        t1=time()
        #vcf=polyn.polyNOptim(degree,data,pdata)
        vcf=polyNOptim(degree,data,pdata,bxscale=0.9)        
        t2=time()
        ### Generate plots (use 'figSave = False' if only want to view (not save))
        plotPoly22(vcf,pdata,saveFig=True)
        t3=time()
        print('-----------Elapsed times (approxs to minutes:seconds)')    
        dt=int(t1-tStart)
        print('preProcessing time (data input and monoms) = {}:{}'.format(dt//60,dt%60))
        dt=int(t2-t1)
        print('optimization time (Convexity and Parameters) = {}:{} '.format(dt//60,dt%60))
        dt=int(t3-t2)
        print('graphics time (plots) = {}:{}'.format(dt//60,dt%60))
        dt=int(t3-tStart)
        print('Overall elapsed time= {}:{}'.format(dt//60,dt%60))
        plt.show()
        exit()
    #################END of modelFit block
    if(task=='testPoly'):
        vCoeff=dData['alpha']                
        data,pdata=readData(inpData,subDir)
        vcf=PolyNparam(data['degree'],vCoeff,matProp=pdata['FEdata'],matName=data['name'],fPrint=True)
        minKG,negKG=PolyN_predictions(vcf,pdata,export=True)
        plotPoly22(vcf,pdata,saveFig=True);plt.show()
        exit()
    if(task=='convexCheck'):
        if(dData['fileCVXData']):
            ##print(figDirData)
            degree,vcf=readVcoeff(dData['fileCVXData'],figDirData+'/'+subDir+'/')
        else:
            degree,vcf=dData['degree'],dData['alpha']
            if(degree not in [4,6,8]):print('degree must be in [4,6,8]\nCalculations aborted');exit()
            vcf=PolyNparam(degree,vcf,fPrint=False)
        print('Checking convexity.......')    
        minKG,negKG,[md1,md2,md3],cvxPoints=PolyN_GaussCheck(vcf,nRandom=2*10**5)
        print('Convexity checked at {} locations:'.format(cvxPoints))
        print('min Gauss KG = ',minKG)
        print('min det1, det2, det3: {}, {}, {}'.format(md1,md2,md3))
        if(minKG<0):
            msg='---Negative KG at {} out of {} points'.format((str(len(negKG))), str(cvxPoints))
            print(msg)
            #print('---Negative KG locations and values: x, y, z, KG')        
            #for point in negKG:
            #    print('{}, {}, {}, {}'.format(point[0],point[1],point[2],point[3]))
        exit()        
    if(0):##(task=='plotFigs'):
        if(0):##FACET plot 
            subDir='AA6016T4' 
            fName='AA6016T4_FACET_RV_deg8_FEdata.txt'
            fPath=fData+subDir+'/'    
            degree,vcf=polyn.readVcoeff(fName,fPath)
            inpData='matAA6016T4_TUAT_replot.txt' ##actual data for AA6016T4
            data,pdata=polyn.readData(inpData,subDir)    
            polyn.plotPoly22(vcf,pdata,saveFig=True);polyn.plt.show()
        if(0):##double plot (of two UA-models)
            subDir='AA6016T4' 
            fName='AA6016T4_UA_RV5_deg8_FEdata.txt'
            fPath=fData+subDir+'/'    
            degree,vcf=polyn.readVcoeff(fName,fPath)
            inpData='matAA6016T4_UA_RV5.txt' ##actual data (and save name) for UA-characterization of AA6016T4
            data,pdata=polyn.readData(inpData,subDir)    
            zfg,zfg2,zfg3=polyn.plotPoly22(vcf,pdata,saveFig=False)
            fName='AA6016T4_UA_RV8_deg8_FEdata.txt'
            degree,vcf=polyn.readVcoeff(fName,fPath)
            polyn.plotPoly22(vcf,pdata,saveFig=True,zaxes=True,zfg=zfg,zfg2=zfg2,zfg3=zfg3,zstyle=True,preName='dbPlot_')
            polyn.plt.show()
        if(0):##double plot (of two TUAT-models)
            subDir='AA6016T4' 
            fName='AA6016T4_TUAT_Yld2004_Postech_deg8_FEdata.txt'
            fPath=fData+subDir+'/'    
            degree,vcf=polyn.readVcoeff(fName,fPath)
            inpData='matAA6016T4_TUAT_Postech.txt' ##actual data (and save name) for TUAT-characterization of AA6016T4
            data,pdata=polyn.readData(inpData,subDir)    
            zfg,zfg2,zfg3=polyn.plotPoly22(vcf,pdata,saveFig=False)
            fName='AA6016T4_Yld2000_Siegen_deg8_FEdata.txt'
            degree,vcf=polyn.readVcoeff(fName,fPath)
            polyn.plotPoly22(vcf,pdata,saveFig=True,zaxes=True,zfg=zfg,zfg2=zfg2,zfg3=zfg3,zstyle=True,preName='dbPlot_')
            polyn.plt.show()
        if(0):##
            subDir='AA6016T4' 
            #fName='AA6016T4_FACET_RV_deg8_FEdata.txt'
            #fPath=fData+subDir+'/'    
            #degree,vcf=polyn.readVcoeff(fName,fPath)
            inpData='matAA6016T4_TUAT_UGent.txt' ##actual data for AA6016T4
            data,pdata=polyn.readData(inpData,subDir)
            degree,vAlpha=6.11,[0.6062, 1.2813, 1.2030, 1.0214, 1.0333, 0.8686, 0.7258, 1.3764]        
            polyn.plotYld2000(degree,vAlpha,pdata,saveFig=True);polyn.plt.show()
        if(1):
            subDir='AA6016T4' 
            inpData='matAA6016T4_TUAT_UGent.txt' # ##actual data (and save name) for TUAT-characterization of AA6016T4
            data,pdata=polyn.readData(inpData,subDir)
            degree,vAlpha=6.11,[0.6062, 1.2813, 1.2030, 1.0214, 1.0333, 0.8686, 0.7258, 1.3764]        
            zfg,zfg2,zfg3=polyn.plotYld2000(degree,vAlpha,pdata,saveFig=False,dataPlot=False,export=False)
            fPath=fData+subDir+'/'    
            fName='AA6016T4_Yld2000_UGent_ExptP6.11_deg8_FEdata.txt'
            degree,vcf=polyn.readVcoeff(fName,fPath)
            polyn.plotPoly22(vcf,pdata,saveFig=True,zaxes=True,zfg=zfg,zfg2=zfg2,zfg3=zfg3,zstyle=True,preName='dbPlot_')
            polyn.plt.show()    
        exit()
    '''
    Block 'convexCheck' allows for custom convexity check 
    One can read parameters from '*Err_and_Coeff.txt' or '*FEdata.txt' files  
    '''
    if(task=='convexCheck'):
        fPath=fData ##Default search path 
        fName='AA6022T4_H_ExptP6_deg10_Err_and_Coeff.txt'
        #fName='AA6022T4_H_ExptP6_deg10_FEdata.txt'
        fPath='.\\AA6022T4files\\'
        degree,vcf=polyn.readVcoeff(fName,fPath)
        ##By default, convexity is checked at a regular grid of about 14000 points
        ##plus nRandom=10**5 points spread randomly over the unit sphere.    
        ##We increase nRandom here to, say, 2*10**5
        print('Checking convexity.......')    
        minKG,negKG,[md1,md2,md3],cvxPoints=polyn.PolyN_GaussCheck(vcf,nRandom=2*10**5)
        print('Convexity checked at {} locations:'.format(cvxPoints))
        print('min Gauss KG = ',minKG)
        print('min det1, det2, det3: {}, {}, {}'.format(md1,md2,md3))
        if(minKG<0):
            msg='---Negative KG at {} out of {} points'.format((str(len(negKG))), str(cvxPoints))
            print(msg)
            print('---Negative KG locations and values: x, y, z, KG')        
            for point in negKG:
                print('{}, {}, {}, {}'.format(point[0],point[1],point[2],point[3]))
        exit()        
    ################END of convexCheck block 
    '''
    Predictions are included in the report file '*Err_and_Coeff.txt' only for data tagged as 'a'.
    One may want to check the yield surface model for other values. 
    Also, one might be interested in the predictions of other yield functions 
    that are particular cases of polynomial yield functions. 
    One example is Yld2000_2D 
    (see 'PolyN_symMain.py' for an example showing how to generate the Poly8 coefficients from Yld2000_2D parameters). 
    We want to check its prediction for the plane strain data point at 45 degs from RD 
    and for the balanced-biaxial state.
    This is illustrated below.  
    Note:reported r-values are generalized r-values 
    (they coincide with the classical r-values only for pure uniaxial stress states)
    For example, the classical balanced-biaxial r-value (rB) can be obtained from the generalized r-value (r)
    by using the formula:
    rB = -r/(1+r)
    '''
    if(task=='testVal'):
        fName='AA6022T4_H_Yld2000_Tian_deg8_FEdata.txt'
        fPath=fData+'AA6022T4_dig_bx/'    
        degree,vcf=polyn.readVcoeff(fName,fPath) 
        ##Calculate the overall error for uniaxial directional properties
        ##d=(q,theta,s_L,r-value)
        vd=[(0.0, 0.0, 1.000, 0.800),
        (0.0, 15, 1.000, 0.74),
        (0.0, 30, 1.001, 0.54),
        (0.0, 45, 0.973, 0.37),
        (0.0, 60, 0.962, 0.41),
        (0.0, 75, 0.954, 0.47),
        (0.0, 90, 0.955, 0.54)]
        derrr,derrs=0.0,0.0
        for d in vd:
            yf,rv=polyn.testValYF(d,vcf)
            derrs+=(yf-d[2])**2
            derrr+=(rv-d[3])**2
        print('overall LSQ error directional yield: ', round(np.sqrt(derrs),3))
        print('overall LSQ error directional r-value: ', round(np.sqrt(derrr),3))    
        ##The data row of the plane strain test at 45 degs from RD
        d=(0.5000, 45.0, 1.0190, 0.0)
        yf,rv=polyn.testValYF(d,vcf)
        ##The data row of the plane strain test along RD
        d=(0.476,0.0,1.096,0.0)
        yf,rv=polyn.testValYF(d,vcf)
        ##The data row of the plane strain test along TD
        d=(0.500,90.0,0.978,0.0)
        yf,rv=polyn.testValYF(d,vcf)
        ##The data row of balanced-biaxial
        d=(1.0, 0.0, 0.95, -0.519)
        yf,rv=polyn.testValYF(d,vcf)
        ### Classical balanced-biaxial r-value
        print("Classical r-value: data, prediction")
        print(-d[3]/(1+d[3]),', ',-rv/(1+rv))
        exit()
    ################END of testVal block    
    if(task=='testBiax'):
        fName="AA6022T4_H_Yld2000_Tian_deg8_FEdata.txt"
        fPath=fData+'AA6022T4_dig_bx\\'   
        degree,vcf=polyn.readVcoeff(fName,fPath)
        inpData,subDir='AA6022T4_H_deg6_Export.txt','AA6022T4_dig_bx'
        data,pdata=polyn.readData(inpData,subDir)
        zdd=polyn.vPoly(degree)
        polyn.tstPlotBx(data,pdata,np.array(vcf[2:3+degree]),zdd,saveFig=False)
        #######npt,nvec=polyn.genConstraintsPoints2DOpt(number=1);print('nconst: ',npt,nvec)
        exit()
    
    
if __name__=='__main__':
    print('\n'+3*'Use the driver script, Luke !\n')


