from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs
from sklearn.cross_validation import train_test_split
from scipy.optimize import curve_fit
from sys import argv
import _velcalc
# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'


def overlay(xt,ys,coef):
    a = coef[0]
    b = coef[1]
    c = coef[2]
    d = coef[3]
    e = coef[4]
    f = coef[5]
    g = coef[6]
    h = coef[7]
    i = coef[8]
    j = coef[9]

    return a + b*xt + c*ys + d*xt**2 + e*xt*ys + f*ys**2 + g*xt**3 + h*xt**2*ys + i*xt*ys**2 + j*ys**3



def veldist(dn,lat,men,sdv1,sdv2,sdv3,sdv4,rat,tns,spr1,spr2,spr3,spr4,scl1,scl2,scl3,sdv_gom,spr_gom):

    if sdv_gom == 0:
        sdv_v = sdv3*sdv2*sdv1*exp(sdv2*dn)*exp(-sdv1*exp(sdv2*dn))#+sdv4
    elif sdv_gom == 1:
        sdv_v = sdv1

    if spr_gom == 0:
        spr_v = spr3*spr2*spr1*exp(spr2*dn)*exp(-spr1*exp(spr2*dn))#+spr4
    elif spr_gom == 1:
        spr_v = 1.

    f1 = -1./(sdv_v*sqrt(2.*pi))*exp(-((lat/spr_v)-men)**2/(2.*sdv_v**2))*(1./(1.+exp(rat*fabs((lat/spr_v))-tns)))
    f2 = scl3*scl2*scl1*exp(scl2*dn)*exp(-scl1*exp(scl2*dn))

    # if sdv_gom == 0:
    #     sdv_v = sdv3*sdv2*sdv1*exp(sdv2*dn)*exp(sdv1)*exp(-sdv1*exp(sdv2*dn))+sdv4
    # elif sdv_gom == 1:
    #     sdv_v = sdv1
    #
    # if spr_gom == 0:
    #     spr_v = spr3*spr2*spr1*exp(spr2*dn)*exp(spr1)*exp(-spr1*exp(spr2*dn))+spr4
    # elif spr_gom == 1:
    #     spr_v = 1.
    #
    # f1 = -1./(sdv_v*sqrt(2.*pi))*exp(-((lat/spr_v)-men)**2/(2.*sdv_v**2))*(1./(1.+exp(rat*fabs((lat/spr_v))-tns)))
    # f2 = scl3*scl2*scl1*exp(scl2*dn)*exp(scl1)*exp(-scl1*exp(scl2*dn))

    return f1*f2 + 1.



def starccm_read(fdata,dia,windd,length,opt_print):
    start = length/30.
    lendat =  np.linspace(start,length,30)/dia

    for j in range(np.size(fdata)):
        for k in range(30):
            name = str(k+1)
            exec('pos'+name+' = np.array([])')
            exec('velo'+name+' = np.array([])')
        wind = windd[j]

        f = open(fdata[j])

        csv_f = csv.reader(f)

        i = 0
        for row in csv_f:
            if row[0] != 'null' and i != 0:
                pos1 = np.append(pos1,float(row[0]))
            if row[2] != 'null' and i != 0:
                pos2 = np.append(pos2,float(row[2]))
            if row[4] != 'null' and i != 0:
                pos3 = np.append(pos3,float(row[4]))
            if row[6] != 'null' and i != 0:
                pos4 = np.append(pos4,float(row[6]))
            if row[8] != 'null' and i != 0:
                pos5 = np.append(pos5,float(row[8]))
            if row[10] != 'null' and i != 0:
                pos6 = np.append(pos6,float(row[10]))
            if row[12] != 'null' and i != 0:
                pos7 = np.append(pos7,float(row[12]))
            if row[14] != 'null' and i != 0:
                pos8 = np.append(pos8,float(row[14]))
            if row[16] != 'null' and i != 0:
                pos9 = np.append(pos9,float(row[16]))
            if row[18] != 'null' and i != 0:
                pos10 = np.append(pos10,float(row[18]))
            if row[20] != 'null' and i != 0:
                pos11 = np.append(pos11,float(row[20]))
            if row[22] != 'null' and i != 0:
                pos12 = np.append(pos12,float(row[22]))
            if row[24] != 'null' and i != 0:
                pos13 = np.append(pos13,float(row[24]))
            if row[26] != 'null' and i != 0:
                pos14 = np.append(pos14,float(row[26]))
            if row[28] != 'null' and i != 0:
                pos15 = np.append(pos15,float(row[28]))
            if row[30] != 'null' and i != 0:
                pos16 = np.append(pos16,float(row[30]))
            if row[32] != 'null' and i != 0:
                pos17 = np.append(pos17,float(row[32]))
            if row[34] != 'null' and i != 0:
                pos18 = np.append(pos18,float(row[34]))
            if row[36] != 'null' and i != 0:
                pos19 = np.append(pos19,float(row[36]))
            if row[38] != 'null' and i != 0:
                pos20 = np.append(pos20,float(row[38]))
            if row[40] != 'null' and i != 0:
                pos21 = np.append(pos21,float(row[40]))
            if row[42] != 'null' and i != 0:
                pos22 = np.append(pos22,float(row[42]))
            if row[44] != 'null' and i != 0:
                pos23 = np.append(pos23,float(row[44]))
            if row[46] != 'null' and i != 0:
                pos24 = np.append(pos24,float(row[46]))
            if row[48] != 'null' and i != 0:
                pos25 = np.append(pos25,float(row[48]))
            if row[50] != 'null' and i != 0:
                pos26 = np.append(pos26,float(row[50]))
            if row[52] != 'null' and i != 0:
                pos27 = np.append(pos27,float(row[52]))
            if row[54] != 'null' and i != 0:
                pos28 = np.append(pos28,float(row[54]))
            if row[56] != 'null' and i != 0:
                pos29 = np.append(pos29,float(row[56]))
            if row[58] != 'null' and i != 0:
                pos30 = np.append(pos30,float(row[58]))

            if row[1] != 'null' and i != 0:
                velo1 = np.append(velo1,float(row[1]))
            if row[3] != 'null' and i != 0:
                velo2 = np.append(velo2,float(row[3]))
            if row[5] != 'null' and i != 0:
                velo3 = np.append(velo3,float(row[5]))
            if row[7] != 'null' and i != 0:
                velo4 = np.append(velo4,float(row[7]))
            if row[9] != 'null' and i != 0:
                velo5 = np.append(velo5,float(row[9]))
            if row[11] != 'null' and i != 0:
                velo6 = np.append(velo6,float(row[11]))
            if row[13] != 'null' and i != 0:
                velo7 = np.append(velo7,float(row[13]))
            if row[15] != 'null' and i != 0:
                velo8 = np.append(velo8,float(row[15]))
            if row[17] != 'null' and i != 0:
                velo9 = np.append(velo9,float(row[17]))
            if row[19] != 'null' and i != 0:
                velo10 = np.append(velo10,float(row[19]))
            if row[21] != 'null' and i != 0:
                velo11 = np.append(velo11,float(row[21]))
            if row[23] != 'null' and i != 0:
                velo12 = np.append(velo12,float(row[23]))
            if row[25] != 'null' and i != 0:
                velo13 = np.append(velo13,float(row[25]))
            if row[27] != 'null' and i != 0:
                velo14 = np.append(velo14,float(row[27]))
            if row[29] != 'null' and i != 0:
                velo15 = np.append(velo15,float(row[29]))
            if row[31] != 'null' and i != 0:
                velo16 = np.append(velo16,float(row[31]))
            if row[33] != 'null' and i != 0:
                velo17 = np.append(velo17,float(row[33]))
            if row[35] != 'null' and i != 0:
                velo18 = np.append(velo18,float(row[35]))
            if row[37] != 'null' and i != 0:
                velo19 = np.append(velo19,float(row[37]))
            if row[39] != 'null' and i != 0:
                velo20 = np.append(velo20,float(row[39]))
            if row[41] != 'null' and i != 0:
                velo21 = np.append(velo21,float(row[41]))
            if row[43] != 'null' and i != 0:
                velo22 = np.append(velo22,float(row[43]))
            if row[45] != 'null' and i != 0:
                velo23 = np.append(velo23,float(row[45]))
            if row[47] != 'null' and i != 0:
                velo24 = np.append(velo24,float(row[47]))
            if row[49] != 'null' and i != 0:
                velo25 = np.append(velo25,float(row[49]))
            if row[51] != 'null' and i != 0:
                velo26 = np.append(velo26,float(row[51]))
            if row[53] != 'null' and i != 0:
                velo27 = np.append(velo27,float(row[53]))
            if row[55] != 'null' and i != 0:
                velo28 = np.append(velo28,float(row[55]))
            if row[57] != 'null' and i != 0:
                velo29 = np.append(velo29,float(row[57]))
            if row[59] != 'null' and i != 0:
                velo30 = np.append(velo30,float(row[59]))
            i += 1

        f.close()

        if opt_print == True:
            print 'Imported Data Set',j+1
        # print 'DATA IMPORTED'

        posdn = np.array([])
        poslt = np.array([])
        velod = np.array([])
        for i in range(30):
            name = str(i+1)
            ind = str(i)

            #Ordering the data numerically by the position
            exec('pos'+name+', velo'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', velo'+name+'))))')
            #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
            exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)/dia\nvelo'+name+' = np.copy(velo'+name+'_0)/wind')
            #Deleting wall boundary data
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5. or pos'+name+'[j] < -5.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            #Deleting values greater than 1*wind
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif velo'+name+'[j] > 1. or fabs(pos'+name+'[j]) > 1.2:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            exec('lensize = np.size(pos'+name+')')
            exec('posdn = np.append(posdn,np.ones(lensize)*lendat['+ind+'])\nposlt = np.append(poslt,pos'+name+')\nvelod = np.append(velod,velo'+name+')')

    return posdn,poslt,velod


def sheet(pos,param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18,param19,param20,param21,param22,param23,param24,param25,param26,param27,param28,param29,param30,param31,param32,param33,param34,param35,param36,param37,param38,param39,param40,param41,param42,param43,param44,param45,param46,param47,param48,param49,param50,param51,param52,param53,param54,param55,param56,param57,param58,param59,param60,param61,param62,param63,param64,param65,param66,param67,param68,param69,param70,param71,param72,param73,param74,param75,param76,param77,param78,param79,param80,param81,param82,param83,param84,param85,param86,param87,param88,param89,param90,param91,param92,param93,param94,param95,param96,param97,param98,param99,param100,param101,param102,param103,param104,param105,param106,param107,param108,param109,param110,param111,param112,param113,param114,param115,param116,param117,param118,param119,param120,param121,param122,param123,param124,param125,param126,param127,param128,param129,param130,param131,param132,param133,param134,param135,param136,param137,param138,param139,param140):
    global xttr
    global ystr
    global sdv_gom
    global spr_gom


    coef0 = np.array([])
    coef1 = np.array([])
    coef2 = np.array([])
    coef3 = np.array([])
    coef4 = np.array([])
    coef5 = np.array([])
    coef6 = np.array([])
    coef7 = np.array([])
    coef8 = np.array([])
    coef9 = np.array([])
    coef10 = np.array([])
    coef11 = np.array([])
    coef12 = np.array([])
    coef13 = np.array([])

    k = 1
    for i in range(14):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',param'+str(k)+')')
            k += 1

    posdn = pos[0,:]
    poslt = pos[1,:]

    vel = _velcalc.sheet(xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,coef10,coef11,coef12,coef13,sdv_gom,spr_gom)
    # print vel

    ## Using Python instead of Fortran
    # vel = np.zeros_like(posdn)
    # for i in range(np.size(pos[0,:])):
    #     men = overlay(xttr[i],ystr[i],coef0)
    #     sdv1 = overlay(xttr[i],ystr[i],coef1)
    #     sdv2 = overlay(xttr[i],ystr[i],coef2)
    #     sdv3 = overlay(xttr[i],ystr[i],coef3)
    #     sdv4 = overlay(xttr[i],ystr[i],coef4)
    #     rat = overlay(xttr[i],ystr[i],coef5)
    #     tns = overlay(xttr[i],ystr[i],coef6)
    #     spr1 = overlay(xttr[i],ystr[i],coef7)
    #     spr2 = overlay(xttr[i],ystr[i],coef8)
    #     spr3 = overlay(xttr[i],ystr[i],coef9)
    #     spr4 = overlay(xttr[i],ystr[i],coef10)
    #     scl1 = overlay(xttr[i],ystr[i],coef11)
    #     scl2 = overlay(xttr[i],ystr[i],coef12)
    #     scl3 = overlay(xttr[i],ystr[i],coef13)
    #
    #
    #     vel[i] = veldist(posdn[i],poslt[i],men,sdv1,sdv2,sdv3,sdv4,rat,tns,spr1,spr2,spr3,spr4,scl1,scl2,scl3,sdv_gom,spr_gom)
    #     print i,'of (',np.size(posdn),')'

    # print vel

    return vel

def obj_func(xdict):
    global posdntr_opt
    global velodtr_opt

    param = xdict['param']
    funcs = {}

    coef0 = np.array([])
    coef1 = np.array([])
    coef2 = np.array([])
    coef3 = np.array([])
    coef4 = np.array([])
    coef5 = np.array([])
    coef6 = np.array([])
    coef7 = np.array([])
    coef8 = np.array([])
    coef9 = np.array([])
    coef10 = np.array([])
    coef11 = np.array([])
    coef12 = np.array([])
    coef13 = np.array([])

    k = 0
    for i in range(14):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',param['+str(k)+'])')
            k += 1

    posdn = posdntr_opt[0,:]
    poslt = posdntr_opt[1,:]

    vel = _velcalc.sheet(xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,coef10,coef11,coef12,coef13,sdv_gom,spr_gom)
    error = np.sum((vel-velodtr_opt)**2)

    ##Print
    # print error

    funcs['obj'] = error

    fail = False

    return funcs, fail

if __name__ == "__main__":

    comp = 'mac'
    # comp = 'fsl'

    global xttr
    global ystr
    global sdv_gom
    global spr_gom
    global posdntr_opt
    global velodtr_opt
    read_data = int(argv[3])#4
    print 'Read Data:',read_data

    cv_on = int(argv[4])
    fitoptn = int(argv[5])
    if fitoptn == 0:
        fit_opt = 'scipy'
    elif fitoptn == 1:
        fit_opt = 'snopt'

    if cv_on == 0:
        print '********Cross Validating Wake Model********'
    elif cv_on == 1:
        print '********Curve Fitting Wake Model********'
    print '\nReading in '+str(read_data)+' data sets'

    print '\nUsing',fit_opt

    if comp == 'mac':
        opt_print = True
        sdv_gom = int(argv[1])
        spr_gom = int(argv[2])
        # print 'TSR Order:',int(argv[1])-1
        # print 'Solidity Order:',int(argv[2])-1
        if int(argv[1]) == 0:
            print 'SDV Fit is Active'
        elif int(argv[1]) == 1:
            print 'SDV Fit is Disabled'
        if int(argv[2]) == 0:
            print 'SPR Fit is Activated\n'
        elif int(argv[2]) == 1:
            print 'SPR Fit is Disabled\n'
    elif comp == 'fsl':
        opt_print = False
        sdv_gom = int(argv[1])
        spr_gom = int(argv[2])
        # print 'TSR Order:',int(argv[1])-1
        # print 'Solidity Order:',int(argv[2])-1
        if int(argv[1]) == 0:
            print 'SDV Fit is Active'
        elif int(argv[1]) == 1:
            print 'SDV Fit is Disabled'
        if int(argv[2]) == 0:
            print 'SPR Fit is Activated\n'
        elif int(argv[2]) == 1:
            print 'SPR Fit is Disabled\n'




    s1length = np.array([210,210,205,196,185,178,170,165,160,145,140,123,115,112,108,101,101,90,85,80,78,75,70])*1.
    s2length = np.array([197,193,185,176,140,146,126,114,103,96,100,86,77,72,70,68,60,64,54,50,47,45,44])*1.
    s3length = np.array([185,150,100,95,83,76,72,63,60,49,50,41,39,36,34,33,30,31,28,30,29,28,27])*1.
    s4length = np.array([145,100,73,60,53,44,42,37,38,30,33,26,22,24,23,21,21,19,24,23,22,21,20])*1.
    s5length = np.array([78,70,52,43,37,32,29,27,26,23,20,20,23,21,20,19,19,18,18,16,16,15,14])*1.
    solidity = np.array(['s1','s2','s3','s4','s5'])
    solidity_cv = np.array([1,2,3,4,5])
    sol_cv = np.array([0.15,0.25,0.5,0.75,1.0])
    tsr_cv = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])
    tsr = np.linspace(150,700,23)
    tsrr_cv = tsr/100.

    s = np.array([])
    t = np.array([])
    scv = np.array([])
    tcv = np.array([])
    xtcv = np.array([])
    yscv = np.array([])
    for i in range(np.size(solidity)):
        for j in range(np.size(tsr)):
            s = np.append(s,solidity[i])
            t = np.append(t,str(int(tsr[j])))
            scv = np.append(scv,solidity_cv[i])
            tcv = np.append(tcv,tsr_cv[j])
            xtcv = np.append(xtcv,tsrr_cv[j])
            yscv = np.append(yscv,sol_cv[i])


    slength = np.array([])
    slength = np.append(slength,s1length)
    slength = np.append(slength,s2length)
    slength = np.append(slength,s3length)
    slength = np.append(slength,s4length)
    slength = np.append(slength,s5length)
    slength = slength*1.

    q = 0
    for i in range(5):
        for j in range(23):
            sname = str(i+1)
            tname = str(j+1)

            t2 = t[q]+'.0'

            wfit = s[q]+'_'+t2
            wfit2 = s[q]+'_'+t2
            wfit3 = s[q]+'_'+t2
            wfit4 = s[q]+'_'+t2
            wfit5 = s[q]+'_'+t2
            wfit6 = s[q]+'_'+t2

            length2 = slength[q]
            length3 = slength[q]
            length4 = slength[q]
            length5 = slength[q]
            length6 = slength[q]
            wind = 15.
            wind2 = 14.
            wind3 = 12.
            wind4 = 16.

            rad = 3.
            dia = rad*2.
            tsr = float(wfit[3]+'.'+wfit[4]+wfit[5])
            rot = tsr*wind/rad
            rot2 = tsr*wind2/rad
            rot3 = tsr*wind3/rad
            rot4 = tsr*wind4/rad
            rot5 = 17.
            rot6 = 18.
            wind5 = rot5*rad/tsr
            wind6 = rot6*rad/tsr

            if comp == 'mac':
                # fdata = '/Users/ning1/Documents/Flow Lab/STAR-CCM+/NACA0021/MoveForward/test.csv'
                fdata = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/Velocity Sections/'+wfit+'.csv'
                fdata2 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel14/Velocity/'+wfit2+'.csv'
                fdata3 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel12/Velocity/'+wfit3+'.csv'
                fdata4 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel16/Velocity/'+wfit4+'.csv'
                fdata5 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot17/Velocity/'+wfit5+'.csv'
                fdata6 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot18/Velocity/'+wfit6+'.csv'
            elif comp == 'fsl':
                fdata = '/fslhome/ebtingey/compute/moveForward/Velocity/'+wfit+'.csv'
                fdata2 = '/fslhome/ebtingey/compute/moveForward/vel14/Velocity/'+wfit2+'.csv'
                fdata3 = '/fslhome/ebtingey/compute/moveForward/vel12/Velocity/'+wfit3+'.csv'
                fdata4 = '/fslhome/ebtingey/compute/moveForward/vel16/Velocity/'+wfit4+'.csv'
                fdata5 = '/fslhome/ebtingey/compute/moveForward/rot17/Velocity/'+wfit5+'.csv'
                fdata6 = '/fslhome/ebtingey/compute/moveForward/rot18/Velocity/'+wfit6+'.csv'

            if read_data ==1:
                posdn,poslt,velod = starccm_read(np.array([fdata]),dia,np.array([wind]),slength[q],opt_print)
            if read_data ==2:
                posdn,poslt,velod = starccm_read(np.array([fdata,fdata2]),dia,np.array([wind,wind2]),slength[q],opt_print)
            if read_data ==3:
                posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3]),dia,np.array([wind,wind2,wind3]),slength[q],opt_print)
            if read_data ==4:
                posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4]),dia,np.array([wind,wind2,wind3,wind4]),slength[q],opt_print)
            if read_data ==5:
                posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5]),dia,np.array([wind,wind2,wind3,wind4,wind5]),slength[q],opt_print)
            if read_data ==6:
                posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5,fdata6]),dia,np.array([wind,wind2,wind3,wind4,wind5,wind6]),slength[q],opt_print)

            exec('posdns'+sname+'t'+tname+' = posdn')
            exec('poslts'+sname+'t'+tname+' = poslt')
            exec('velods'+sname+'t'+tname+' = velod')

            print t[q],s[q]
            q += 1

    if cv_on == 0:
        cvtest = 0.3
        scvtr,scvts,tcvtr,tcvts,xttr,xtts,ystr,ysts = train_test_split(scv,tcv,xtcv,yscv,test_size=cvtest)

        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdn'+name+'tr = posdns'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('poslt'+name+'tr = poslts'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('velod'+name+'tr = velods'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('xt'+name+'tr = xttr[i]')
            exec('ys'+name+'tr = ystr[i]')

        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdn'+name+'ts = posdns'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('poslt'+name+'ts = poslts'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('velod'+name+'ts = velods'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('xt'+name+'ts = xtts[i]')
            exec('ys'+name+'ts = ysts[i]')

        posdntr = np.array([])
        poslttr = np.array([])
        velodtr = np.array([])
        xttr = np.array([])
        ystr = np.array([])
        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdntr = np.append(posdntr,posdn'+name+'tr)')
            exec('poslttr = np.append(poslttr,poslt'+name+'tr)')
            exec('velodtr = np.append(velodtr,velod'+name+'tr)')
            exec('xttr = np.append(xttr,np.ones_like(posdn'+name+'tr)*xt'+name+'tr)')
            exec('ystr = np.append(ystr,np.ones_like(posdn'+name+'tr)*ys'+name+'tr)')
        posdntr = np.vstack([posdntr,poslttr])

        if fit_opt == 'snopt':
            posdntr_opt = posdntr
            velodtr_opt = velodtr

        posdnts = np.array([])
        posltts = np.array([])
        velodts = np.array([])
        xtts = np.array([])
        ysts = np.array([])
        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdnts = np.append(posdnts,posdn'+name+'ts)')
            exec('posltts = np.append(posltts,poslt'+name+'ts)')
            exec('velodts = np.append(velodts,velod'+name+'ts)')
            exec('xtts = np.append(xtts,np.ones_like(posdn'+name+'ts)*xt'+name+'ts)')
            exec('ysts = np.append(ysts,np.ones_like(posdn'+name+'ts)*ys'+name+'ts)')
        posdnts = np.vstack([posdnts,posltts])

    elif cv_on == 1:
        cvtest = 0.0
        scvtr,scvts,tcvtr,tcvts,xttr,xtts,ystr,ysts = train_test_split(scv,tcv,xtcv,yscv,test_size=cvtest)

        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdn'+name+'tr = posdns'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('poslt'+name+'tr = poslts'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('velod'+name+'tr = velods'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('xt'+name+'tr = xttr[i]')
            exec('ys'+name+'tr = ystr[i]')

        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdn'+name+'ts = posdns'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('poslt'+name+'ts = poslts'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('velod'+name+'ts = velods'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('xt'+name+'ts = xtts[i]')
            exec('ys'+name+'ts = ysts[i]')

        posdntr = np.array([])
        poslttr = np.array([])
        velodtr = np.array([])
        xttr = np.array([])
        ystr = np.array([])
        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdntr = np.append(posdntr,posdn'+name+'tr)')
            exec('poslttr = np.append(poslttr,poslt'+name+'tr)')
            exec('velodtr = np.append(velodtr,velod'+name+'tr)')
            exec('xttr = np.append(xttr,np.ones_like(posdn'+name+'tr)*xt'+name+'tr)')
            exec('ystr = np.append(ystr,np.ones_like(posdn'+name+'tr)*ys'+name+'tr)')
        posdntr = np.vstack([posdntr,poslttr])

        if fit_opt == 'snopt':
            posdntr_opt = posdntr
            velodtr_opt = velodtr

        posdnts = np.array([])
        posltts = np.array([])
        velodts = np.array([])
        xtts = np.array([])
        ysts = np.array([])
        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdnts = np.append(posdnts,posdn'+name+'ts)')
            exec('posltts = np.append(posltts,poslt'+name+'ts)')
            exec('velodts = np.append(velodts,velod'+name+'ts)')
            exec('xtts = np.append(xtts,np.ones_like(posdn'+name+'ts)*xt'+name+'ts)')
            exec('ysts = np.append(ysts,np.ones_like(posdn'+name+'ts)*ys'+name+'ts)')
        posdnts = np.vstack([posdnts,posltts])




    param0 = np.array([0.8948,-0.5217,5.864,0.04435,-0.4683,-11.85,0.0002901,0.03235,0.4074,6.154,   #men
                      3.509,-2.787,19.54,0.4326,-0.0563,-40.31,-0.01907,-0.08837,1.319,19.79,   #sdv1
                      2.816,-1.148,-4.651,0.1687,1.363,2.214,-0.008379,-0.08856,-0.24,-0.323,   #sdv2
                      32.51,-15.18,-40.38,4.555,-18.65,151.5,-0.4148,2.334,-1.695,-81.1,   #sdv3
                      5.486,-3.383,4.927,0.7519,-1.539,-6.293,-0.05414,0.1356,0.5925,2.937,   #sdv4
                      -14.47,6.027,19.79,-2.124,44.46,-87.36,0.2503,-3.983,-17.56,56.73,   #rat
                      19.06,-5.559,-0.782,0.952,17.13,-25.05,-0.03737,-1.634,-4.769,13.02,   #tns
                      1.311,-0.1151,-3.818,-0.007227,0.6576,5.142,0.001608,-0.04294,-0.2846,-2.276,   #spr1
                      -0.1603,0.06013,0.007054,-0.009813,0.2327,0.1205,0.000678,-0.01671,-0.1247,-0.1181,   #spr2
                      5.847,18.5,134.2,-1.814,-73.19,-93.16,-0.04043,5.52,28.66,19.76,   #spr3
                      1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,   #spr4
                      2.241,-0.8028,-4.366,0.1012,1.184,3.889,-0.004359,-0.03936,-0.6317,-0.5302,   #scl1
                      -0.1756,0.04034,0.4726,0.001038,0.1089,-0.6035,-0.0001736,-0.01091,0.009909,0.2134,   #scl2
                      155.9,-125.9,572.2,35.25,-191.3,-331.1,-2.914,16.57,34.04,102.8])   #scl3
    """
#cubic- Matlab
f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3
men
Coefficients (with 95% confidence bounds):
       p00 =      0.8948  (-1.734, 3.523)
       p10 =     -0.5217  (-2.272, 1.229)
       p01 =       5.864  (-2.626, 14.36)
       p20 =     0.04435  (-0.3685, 0.4572)
       p11 =     -0.4683  (-2.106, 1.169)
       p02 =      -11.85  (-25.83, 2.124)
       p30 =   0.0002901  (-0.03132, 0.0319)
       p21 =     0.03235  (-0.1131, 0.1778)
       p12 =      0.4074  (-0.518, 1.333)
       p03 =       6.154  (-1.481, 13.79)

       0.8948,-0.5217,5.864,0.04435,-0.4683,-11.85,0.0002901,0.03235,0.4074,6.154

Goodness of fit:
  SSE: 38.56
  R-square: 0.1815
  Adjusted R-square: 0.1114
  RMSE: 0.606

sdv1
Coefficients (with 95% confidence bounds):
       p00 =       3.509  (-2.385, 9.402)
       p10 =      -2.787  (-6.711, 1.138)
       p01 =       19.54  (0.4963, 38.57)
       p20 =      0.4326  (-0.4931, 1.358)
       p11 =     -0.0563  (-3.728, 3.616)
       p02 =      -40.31  (-71.65, -8.97)
       p30 =    -0.01907  (-0.08995, 0.0518)
       p21 =    -0.08837  (-0.4146, 0.2379)
       p12 =       1.319  (-0.7562, 3.394)
       p03 =       19.79  (2.672, 36.91)

       3.509,-2.787,19.54,0.4326,-0.0563,-40.31,-0.01907,-0.08837,1.319,19.79

Goodness of fit:
  SSE: 193.9
  R-square: 0.2112
  Adjusted R-square: 0.1436
  RMSE: 1.359

sdv2
Coefficients (with 95% confidence bounds):
       p00 =       2.816  (1.947, 3.686)
       p10 =      -1.148  (-1.727, -0.5693)
       p01 =      -4.651  (-7.46, -1.842)
       p20 =      0.1687  (0.03213, 0.3053)
       p11 =       1.363  (0.8213, 1.905)
       p02 =       2.214  (-2.41, 6.838)
       p30 =   -0.008379  (-0.01883, 0.002078)
       p21 =    -0.08856  (-0.1367, -0.04043)
       p12 =       -0.24  (-0.5462, 0.06611)
       p03 =      -0.323  (-2.849, 2.203)

       2.816,-1.148,-4.651,0.1687,1.363,2.214,-0.008379,-0.08856,-0.24,-0.323

Goodness of fit:
  SSE: 4.22
  R-square: 0.5797
  Adjusted R-square: 0.5437
  RMSE: 0.2005


sdv3
Coefficients (with 95% confidence bounds):
       p00 =       32.51  (-0.1639, 65.19)
       p10 =      -15.18  (-36.94, 6.581)
       p01 =      -40.38  (-145.9, 65.18)
       p20 =       4.555  (-0.5774, 9.687)
       p11 =      -18.65  (-39.01, 1.704)
       p02 =       151.5  (-22.24, 325.3)
       p30 =     -0.4148  (-0.8077, -0.0219)
       p21 =       2.334  (0.5251, 4.142)
       p12 =      -1.695  (-13.2, 9.809)
       p03 =       -81.1  (-176, 13.82)

       32.51,-15.18,-40.38,4.555,-18.65,151.5,-0.4148,2.334,-1.695,-81.1

Goodness of fit:
  SSE: 5960
  R-square: 0.1839
  Adjusted R-square: 0.1139
  RMSE: 7.534


sdv4
Coefficients (with 95% confidence bounds):
       p00 =       5.486  (1.818, 9.155)
       p10 =      -3.383  (-5.826, -0.9402)
       p01 =       4.927  (-6.923, 16.78)
       p20 =      0.7519  (0.1757, 1.328)
       p11 =      -1.539  (-3.824, 0.7464)
       p02 =      -6.293  (-25.8, 13.21)
       p30 =    -0.05414  (-0.09825, -0.01003)
       p21 =      0.1356  (-0.06749, 0.3386)
       p12 =      0.5925  (-0.699, 1.884)
       p03 =       2.937  (-7.718, 13.59)

       5.486,-3.383,4.927,0.7519,-1.539,-6.293,-0.05414,0.1356,0.5925,2.937

Goodness of fit:
  SSE: 75.11
  R-square: 0.3007
  Adjusted R-square: 0.2407
  RMSE: 0.8457


rat
Coefficients (with 95% confidence bounds):
       p00 =      -14.47  (-43, 14.06)
       p10 =       6.027  (-12.97, 25.03)
       p01 =       19.79  (-72.38, 112)
       p20 =      -2.124  (-6.606, 2.358)
       p11 =       44.46  (26.68, 62.23)
       p02 =      -87.36  (-239.1, 64.36)
       p30 =      0.2503  (-0.09277, 0.5935)
       p21 =      -3.983  (-5.563, -2.404)
       p12 =      -17.56  (-27.61, -7.514)
       p03 =       56.73  (-26.15, 139.6)

       -14.47,6.027,19.79,-2.124,44.46,-87.36,0.2503,-3.983,-17.56,56.73

Goodness of fit:
  SSE: 4544
  R-square: 0.6637
  Adjusted R-square: 0.6348
  RMSE: 6.579


tns
Coefficients (with 95% confidence bounds):
       p00 =       19.06  (-4.902, 43.03)
       p10 =      -5.559  (-21.52, 10.4)
       p01 =      -0.782  (-78.2, 76.63)
       p20 =       0.952  (-2.812, 4.716)
       p11 =       17.13  (2.204, 32.06)
       p02 =      -25.05  (-152.5, 102.4)
       p30 =    -0.03737  (-0.3255, 0.2508)
       p21 =      -1.634  (-2.96, -0.3072)
       p12 =      -4.769  (-13.21, 3.669)
       p03 =       13.02  (-56.6, 82.63)

       19.06,-5.559,-0.782,0.952,17.13,-25.05,-0.03737,-1.634,-4.769,13.02

Goodness of fit:
  SSE: 3206
  R-square: 0.2046
  Adjusted R-square: 0.1365
  RMSE: 5.525


spr1
Coefficients (with 95% confidence bounds):
       p00 =       1.311  (0.4849, 2.138)
       p10 =     -0.1151  (-0.6654, 0.4352)
       p01 =      -3.818  (-6.487, -1.148)
       p20 =   -0.007227  (-0.137, 0.1226)
       p11 =      0.6576  (0.1428, 1.172)
       p02 =       5.142  (0.7481, 9.536)
       p30 =    0.001608  (-0.008329, 0.01155)
       p21 =    -0.04294  (-0.08868, 0.002804)
       p12 =     -0.2846  (-0.5756, 0.006336)
       p03 =      -2.276  (-4.676, 0.1247)

       1.311,-0.1151,-3.818,-0.007227,0.6576,5.142,0.001608,-0.04294,-0.2846,-2.276

Goodness of fit:
  SSE: 3.812
  R-square: 0.1841
  Adjusted R-square: 0.1142
  RMSE: 0.1905


spr2
Coefficients (with 95% confidence bounds):
       p00 =     -0.1603  (-0.394, 0.07346)
       p10 =     0.06013  (-0.09551, 0.2158)
       p01 =    0.007054  (-0.748, 0.7621)
       p20 =   -0.009813  (-0.04652, 0.0269)
       p11 =      0.2327  (0.08711, 0.3783)
       p02 =      0.1205  (-1.122, 1.363)
       p30 =    0.000678  (-0.002133, 0.003489)
       p21 =    -0.01671  (-0.02965, -0.003773)
       p12 =     -0.1247  (-0.207, -0.04241)
       p03 =     -0.1181  (-0.797, 0.5608)

       -0.1603,0.06013,0.007054,-0.009813,0.2327,0.1205,0.000678,-0.01671,-0.1247,-0.1181

Goodness of fit:
  SSE: 0.3049
  R-square: 0.6457
  Adjusted R-square: 0.6153
  RMSE: 0.05389


spr3
Coefficients (with 95% confidence bounds):
       p00 =       5.847  (-34.36, 46.05)
       p10 =        18.5  (-8.276, 45.27)
       p01 =       134.2  (4.364, 264.1)
       p20 =      -1.814  (-8.129, 4.501)
       p11 =      -73.19  (-98.24, -48.14)
       p02 =      -93.16  (-307, 120.6)
       p30 =    -0.04043  (-0.5239, 0.4431)
       p21 =        5.52  (3.294, 7.745)
       p12 =       28.66  (14.51, 42.82)
       p03 =       19.76  (-97.03, 136.5)

       5.847,18.5,134.2,-1.814,-73.19,-93.16,-0.04043,5.52,28.66,19.76

Goodness of fit:
  SSE: 9023
  R-square: 0.5194
  Adjusted R-square: 0.4782
  RMSE: 9.27


spr4
Coefficients (with 95% confidence bounds):
       p00 =           1  (1, 1)
       p10 =    1.83e-16  (-3.365e-15, 3.731e-15)
       p01 =   3.875e-16  (-1.682e-14, 1.76e-14)
       p20 =   -3.54e-17  (-8.723e-16, 8.015e-16)
       p11 =  -5.037e-17  (-3.37e-15, 3.269e-15)
       p02 =  -1.128e-15  (-2.946e-14, 2.72e-14)
       p30 =   1.095e-18  (-6.298e-17, 6.517e-17)
       p21 =   5.041e-17  (-2.445e-16, 3.453e-16)
       p12 =  -3.634e-16  (-2.239e-15, 1.513e-15)
       p03 =   6.172e-16  (-1.486e-14, 1.609e-14)

       1.,0.,0.,0.,0.,0.,0.,0.,0.,0.

Goodness of fit:
  SSE: 1.585e-28
  R-square: NaN
  Adjusted R-square: NaN
  RMSE: 1.228e-15


scl1
Coefficients (with 95% confidence bounds):
       p00 =       2.241  (1.458, 3.023)
       p10 =     -0.8028  (-1.324, -0.2815)
       p01 =      -4.366  (-6.895, -1.837)
       p20 =      0.1012  (-0.02177, 0.2242)
       p11 =       1.184  (0.6966, 1.672)
       p02 =       3.889  (-0.2738, 8.052)
       p30 =   -0.004359  (-0.01377, 0.005055)
       p21 =    -0.03936  (-0.0827, 0.003974)
       p12 =     -0.6317  (-0.9073, -0.3561)
       p03 =     -0.5302  (-2.804, 1.744)

       2.241,-0.8028,-4.366,0.1012,1.184,3.889,-0.004359,-0.03936,-0.6317,-0.5302

Goodness of fit:
  SSE: 3.421
  R-square: 0.646
  Adjusted R-square: 0.6156
  RMSE: 0.1805


scl2
Coefficients (with 95% confidence bounds):
       p00 =     -0.1756  (-0.3572, 0.005986)
       p10 =     0.04034  (-0.08059, 0.1613)
       p01 =      0.4726  (-0.114, 1.059)
       p20 =    0.001038  (-0.02749, 0.02956)
       p11 =      0.1089  (-0.004248, 0.222)
       p02 =     -0.6035  (-1.569, 0.3621)
       p30 =  -0.0001736  (-0.002357, 0.00201)
       p21 =    -0.01091  (-0.02096, -0.0008583)
       p12 =    0.009909  (-0.05403, 0.07385)
       p03 =      0.2134  (-0.3141, 0.7409)

       -0.1756,0.04034,0.4726,0.001038,0.1089,-0.6035,-0.0001736,-0.01091,0.009909,0.2134

Goodness of fit:
  SSE: 0.1841
  R-square: 0.916
  Adjusted R-square: 0.9088
  RMSE: 0.04187


scl3
Coefficients (with 95% confidence bounds):
       p00 =       155.9  (-110.6, 422.4)
       p10 =      -125.9  (-303.3, 51.62)
       p01 =       572.2  (-288.8, 1433)
       p20 =       35.25  (-6.612, 77.11)
       p11 =      -191.3  (-357.3, -25.25)
       p02 =      -331.1  (-1748, 1086)
       p30 =      -2.914  (-6.119, 0.2906)
       p21 =       16.57  (1.818, 31.32)
       p12 =       34.04  (-59.79, 127.9)
       p03 =       102.8  (-671.3, 877)

       155.9,-125.9,572.2,35.25,-191.3,-331.1,-2.914,16.57,34.04,102.8

Goodness of fit:
  SSE: 3.965e+05
  R-square: 0.2313
  Adjusted R-square: 0.1654
  RMSE: 61.45






    """

    # from pymatbridge import Matlab
    # mlab = Matlab('/Applications/MATLAB_R2014b.app/bin/matlab')
    # mlab.start()
    #
    # res = mlab.lsqcurvefit(sheet,param,posdntr,velodtr)
    # print res
    #
    # mlab.stop


    if fit_opt == 'snopt':
        optProb = Optimization('VAWTWake_Velo', obj_func)
        optProb.addObj('obj')

        optProb.addVarGroup('param', 140, 'c', lower=None, upper=None, value=param0)

        opt = SNOPT()
        opt.setOption('Scale option',2)
        if comp == 'mac':
            opt.setOption('Print file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptSheet/SNOPT_print'+str(argv[1])+str(argv[2])+str(argv[3])+str(argv[4])+str(argv[5])+'.out')
            opt.setOption('Summary file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptSheet/SNOPT_summary'+str(argv[1])+str(argv[2])+str(argv[3])+str(argv[4])+str(argv[5])+'.out')
        elif comp == 'fsl':
            opt.setOption('Print file','/fslhome/ebtingey/compute/VAWTWakeModel/CrossVal/SNOPT_print'+str(argv[1])+str(argv[2])+str(argv[3])+str(argv[4])+str(argv[5])+'.out')
            opt.setOption('Summary file','/fslhome/ebtingey/compute/VAWTWakeModel/CrossVal/SNOPT_summary'+str(argv[1])+str(argv[2])+str(argv[3])+str(argv[4])+str(argv[5])+'.out')
        result = opt(optProb, sens=None)

        res = result.xStar['param']

    elif fit_opt == 'scipy':
        res,_ = curve_fit(sheet,posdntr,velodtr,p0=param0,maxfev=2820000)


    coef0 = np.array([])
    coef1 = np.array([])
    coef2 = np.array([])
    coef3 = np.array([])
    coef4 = np.array([])
    coef5 = np.array([])
    coef6 = np.array([])
    coef7 = np.array([])
    coef8 = np.array([])
    coef9 = np.array([])
    coef10 = np.array([])
    coef11 = np.array([])
    coef12 = np.array([])
    coef13 = np.array([])

    k = 0
    for i in range(14):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',res[k])')
            k += 1


    if cv_on == 0:
        velc = _velcalc.sheet(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,coef10,coef11,coef12,coef13,sdv_gom,spr_gom)
        velf = _velcalc.sheet(xtts,ysts,posdnts[0,:],posdnts[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,coef10,coef11,coef12,coef13,sdv_gom,spr_gom)

        error_curve = np.sum((velc-velodtr)**2)
        error_cv = np.sum((velf-velodts)**2)

        print '\ncoef0 = np.array(',coef0.tolist(),')'
        print 'coef1 = np.array(',coef1.tolist(),')'
        print 'coef2 = np.array(',coef2.tolist(),')'
        print 'coef3 = np.array(',coef3.tolist(),')'
        print 'coef4 = np.array(',coef4.tolist(),')'
        print 'coef5 = np.array(',coef5.tolist(),')'
        print 'coef6 = np.array(',coef6.tolist(),')'
        print 'coef7 = np.array(',coef7.tolist(),')'
        print 'coef8 = np.array(',coef8.tolist(),')'
        print 'coef9 = np.array(',coef9.tolist(),')'
        print 'coef10 = np.array(',coef10.tolist(),')'
        print 'coef11 = np.array(',coef11.tolist(),')'
        print 'coef12 = np.array(',coef12.tolist(),')'
        print 'coef13 = np.array(',coef13.tolist(),')'

        print '\nparam = np.array([',res[0],',',res[1],',',res[2],',',res[3],',',res[4],',',res[5],',',res[6],',',res[7],',',res[8],',',res[9],',\t#men'
        print res[10],',',res[11],',',res[12],',',res[13],',',res[14],',',res[15],',',res[16],',',res[17],',',res[18],',',res[19],',\t#sdv1'
        print res[20],',',res[21],',',res[22],',',res[23],',',res[24],',',res[25],',',res[26],',',res[27],',',res[28],',',res[29],',\t#sdv2'
        print res[30],',',res[31],',',res[32],',',res[33],',',res[34],',',res[35],',',res[36],',',res[37],',',res[38],',',res[39],',\t#sdv3'
        print res[40],',',res[41],',',res[42],',',res[43],',',res[44],',',res[45],',',res[46],',',res[47],',',res[48],',',res[49],',\t#sdv4'
        print res[50],',',res[51],',',res[52],',',res[53],',',res[54],',',res[55],',',res[56],',',res[57],',',res[58],',',res[59],',\t#rat'
        print res[60],',',res[61],',',res[62],',',res[63],',',res[64],',',res[65],',',res[66],',',res[67],',',res[68],',',res[69],',\t#tns'
        print res[70],',',res[71],',',res[72],',',res[73],',',res[74],',',res[75],',',res[76],',',res[77],',',res[78],',',res[79],',\t#spr1'
        print res[80],',',res[81],',',res[82],',',res[83],',',res[84],',',res[85],',',res[86],',',res[87],',',res[88],',',res[89],',\t#spr2'
        print res[90],',',res[91],',',res[92],',',res[93],',',res[94],',',res[95],',',res[96],',',res[97],',',res[98],',',res[99],',\t#spr3'
        print res[100],',',res[101],',',res[102],',',res[103],',',res[104],',',res[105],',',res[106],',',res[107],',',res[108],',',res[109],',\t#spr4'
        print res[110],',',res[111],',',res[112],',',res[113],',',res[114],',',res[115],',',res[116],',',res[117],',',res[118],',',res[119],',\t#scl1'
        print res[120],',',res[121],',',res[122],',',res[123],',',res[124],',',res[125],',',res[126],',',res[127],',',res[128],',',res[129],',\t#scl2'
        print res[130],',',res[131],',',res[132],',',res[133],',',res[134],',',res[135],',',res[136],',',res[137],',',res[138],',',res[139],'])\t#scl3'

        print '\nCross Valiation Error:',error_cv
        print 'Curve Fit Error:',error_curve


    if cv_on == 1:
        velc = _velcalc.sheet(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,coef10,coef11,coef12,coef13,sdv_gom,spr_gom)

        error_curve = np.sum((velc-velodtr)**2)

        print '\ncoef0 = np.array(',coef0.tolist(),')'
        print 'coef1 = np.array(',coef1.tolist(),')'
        print 'coef2 = np.array(',coef2.tolist(),')'
        print 'coef3 = np.array(',coef3.tolist(),')'
        print 'coef4 = np.array(',coef4.tolist(),')'
        print 'coef5 = np.array(',coef5.tolist(),')'
        print 'coef6 = np.array(',coef6.tolist(),')'
        print 'coef7 = np.array(',coef7.tolist(),')'
        print 'coef8 = np.array(',coef8.tolist(),')'
        print 'coef9 = np.array(',coef9.tolist(),')'
        print 'coef10 = np.array(',coef10.tolist(),')'
        print 'coef11 = np.array(',coef11.tolist(),')'
        print 'coef12 = np.array(',coef12.tolist(),')'
        print 'coef13 = np.array(',coef13.tolist(),')'

        print '\nmen = np.array(',coef0.tolist(),')'
        print 'sdv1 = np.array(',coef1.tolist(),')'
        print 'sdv2 = np.array(',coef2.tolist(),')'
        print 'sdv3 = np.array(',coef3.tolist(),')'
        print 'sdv4 = np.array(',coef4.tolist(),')'
        print 'rat = np.array(',coef5.tolist(),')'
        print 'wdt = np.array(',coef6.tolist(),')'
        print 'spr1 = np.array(',coef7.tolist(),')'
        print 'spr2 = np.array(',coef8.tolist(),')'
        print 'spr3 = np.array(',coef9.tolist(),')'
        print 'spr4 = np.array(',coef10.tolist(),')'
        print 'scl1 = np.array(',coef11.tolist(),')'
        print 'scl2 = np.array(',coef12.tolist(),')'
        print 'scl3 = np.array(',coef13.tolist(),')'

        print '\nparam = np.array([',res[0],',',res[1],',',res[2],',',res[3],',',res[4],',',res[5],',',res[6],',',res[7],',',res[8],',',res[9],',\t#men'
        print res[10],',',res[11],',',res[12],',',res[13],',',res[14],',',res[15],',',res[16],',',res[17],',',res[18],',',res[19],',\t#sdv1'
        print res[20],',',res[21],',',res[22],',',res[23],',',res[24],',',res[25],',',res[26],',',res[27],',',res[28],',',res[29],',\t#sdv2'
        print res[30],',',res[31],',',res[32],',',res[33],',',res[34],',',res[35],',',res[36],',',res[37],',',res[38],',',res[39],',\t#sdv3'
        print res[40],',',res[41],',',res[42],',',res[43],',',res[44],',',res[45],',',res[46],',',res[47],',',res[48],',',res[49],',\t#sdv4'
        print res[50],',',res[51],',',res[52],',',res[53],',',res[54],',',res[55],',',res[56],',',res[57],',',res[58],',',res[59],',\t#rat'
        print res[60],',',res[61],',',res[62],',',res[63],',',res[64],',',res[65],',',res[66],',',res[67],',',res[68],',',res[69],',\t#tns'
        print res[70],',',res[71],',',res[72],',',res[73],',',res[74],',',res[75],',',res[76],',',res[77],',',res[78],',',res[79],',\t#spr1'
        print res[80],',',res[81],',',res[82],',',res[83],',',res[84],',',res[85],',',res[86],',',res[87],',',res[88],',',res[89],',\t#spr2'
        print res[90],',',res[91],',',res[92],',',res[93],',',res[94],',',res[95],',',res[96],',',res[97],',',res[98],',',res[99],',\t#spr3'
        print res[100],',',res[101],',',res[102],',',res[103],',',res[104],',',res[105],',',res[106],',',res[107],',',res[108],',',res[109],',\t#spr4'
        print res[110],',',res[111],',',res[112],',',res[113],',',res[114],',',res[115],',',res[116],',',res[117],',',res[118],',',res[119],',\t#scl1'
        print res[120],',',res[121],',',res[122],',',res[123],',',res[124],',',res[125],',',res[126],',',res[127],',',res[128],',',res[129],',\t#scl2'
        print res[130],',',res[131],',',res[132],',',res[133],',',res[134],',',res[135],',',res[136],',',res[137],',',res[138],',',res[139],'])\t#scl3'

        print '\nCurve Fit Error:',error_curve
