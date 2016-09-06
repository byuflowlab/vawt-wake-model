from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs
from sklearn.cross_validation import train_test_split
from scipy.optimize import curve_fit
from sys import argv
import _velcalcqme
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



def veldist(dn,lat,spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3):
    pow_v = pow1-pow2*dn**pow3
    exp_v = exp(-spr1*fabs(lat)**pow_v)
    quad_v = spr2*(lat-skw)**2 - 1.
    scl_v = scl3*scl2*scl1*exp(scl2*dn)*exp(-scl1*exp(scl2*dn))

    return exp_v*quad_v*scl_v + 1.



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


def sheet(pos,param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18,param19,param20,param21,param22,param23,param24,param25,param26,param27,param28,param29,param30,param31,param32,param33,param34,param35,param36,param37,param38,param39,param40,param41,param42,param43,param44,param45,param46,param47,param48,param49,param50,param51,param52,param53,param54,param55,param56,param57,param58,param59,param60,param61,param62,param63,param64,param65,param66,param67,param68,param69,param70,param71,param72,param73,param74,param75,param76,param77,param78,param79,param80,param81,param82,param83,param84,param85,param86,param87,param88,param89,param90):
    global xttr
    global ystr

    global spr1_ord
    global pow1_ord
    global pow2_ord
    global pow3_ord
    global spr2_ord
    global skw_ord
    global scl1_ord
    global scl2_ord
    global scl3_ord


    coef0 = np.array([])
    coef1 = np.array([])
    coef2 = np.array([])
    coef3 = np.array([])
    coef4 = np.array([])
    coef5 = np.array([])
    coef6 = np.array([])
    coef7 = np.array([])
    coef8 = np.array([])

    k = 1
    for i in range(9):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',param'+str(k)+')')
            k += 1

    posdn = pos[0,:]
    poslt = pos[1,:]

    # p00,p10,p01,p20,p11,p02,p30,p21,p12,p03
    if spr1_ord == 2:
        coef0[9] = 0.

    if pow1_ord == 1:
        coef1[3] = 0.
        coef1[6] = 0.
        coef1[7] = 0.
        coef1[8] = 0.
        coef1[9] = 0.
    elif pow1_ord == 2:
        coef1[6] = 0.
        coef1[7] = 0.
        coef1[8] = 0.
        coef1[9] = 0.

    if pow2_ord == 1:
        coef2[6] = 0.
        coef2[7] = 0.
        coef2[8] = 0.
        coef2[9] = 0.
    elif pow2_ord == 2:
        coef2[5] = 0.
        coef2[6] = 0.
        coef2[7] = 0.
        coef2[8] = 0.
        coef2[9] = 0.

    if pow3_ord == 1:
        coef3[9] = 0.

    if spr2_ord == 1:
        coef4[3] = 0.
        coef4[6] = 0.
        coef4[7] = 0.
        coef4[8] = 0.
        coef4[9] = 0.
    elif spr2_ord == 2:
        coef4[9] = 0.

    if skw_ord == 1:
        coef5[6] = 0.
    elif skw_ord == 2:
        coef5[9] = 0.

    if scl1_ord == 1:
        coef6[9] = 0.

    if scl2_ord == 1:
        coef7[9] = 0.

    if scl3_ord == 2:
        coef8[6] = 0.
        coef8[7] = 0.
        coef8[8] = 0.
        coef8[9] = 0.
    # p00,p10,p01,p20,p11,p02,p30,p21,p12,p03

    vel = _velcalcqme.sheet(xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8)
    # print vel

    ## Using Python instead of Fortran
    # vel = np.zeros_like(posdn)
    # for i in range(np.size(pos[0,:])):
    #     spr1 = overlay(xttr[i],ystr[i],coef0)
    #     pow1 = overlay(xttr[i],ystr[i],coef1)
    #     pow2 = overlay(xttr[i],ystr[i],coef2)
    #     pow3 = overlay(xttr[i],ystr[i],coef3)
    #     spr2 = overlay(xttr[i],ystr[i],coef4)
    #     skw = overlay(xttr[i],ystr[i],coef5)
    #     scl1 = overlay(xttr[i],ystr[i],coef6)
    #     scl2 = overlay(xttr[i],ystr[i],coef7)
    #     scl3 = overlay(xttr[i],ystr[i],coef8)
    #
    #
    #     vel[i] = veldist(posdn[i],poslt[i],spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3)
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

    k = 0
    for i in range(9):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',param['+str(k)+'])')
            k += 1

    posdn = posdntr_opt[0,:]
    poslt = posdntr_opt[1,:]

    vel = _velcalcqme.sheet(xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8)
    error = np.sum((vel-velodtr_opt)**2)

    ##Print
    # print error

    funcs['obj'] = error

    fail = False

    return funcs, fail

if __name__ == "__main__":

    # comp = 'mac'
    comp = 'fsl'

    global xttr
    global ystr
    global posdntr_opt
    global velodtr_opt

    global spr1_ord
    global pow1_ord
    global pow2_ord
    global pow3_ord
    global spr2_ord
    global skw_ord
    global scl1_ord
    global scl2_ord
    global scl3_ord

    read_data = int(argv[1])#4
    print 'Read Data:',read_data

    cv_on = int(argv[2])
    fitoptn = int(argv[3])
    if fitoptn == 0:
        fit_opt = 'scipy'
    elif fitoptn == 1:
        fit_opt = 'snopt'

    if cv_on == 0:
        print '********Cross Validating Wake Model********'
    elif cv_on == 1:
        print '********Curve Fitting Wake Model********'
    print '\nReading in '+str(read_data)+' data sets'

    print '\nUsing',fit_opt,'\n'

    spr1_ord = int(argv[4])
    pow1_ord = int(argv[5])
    pow2_ord = int(argv[6])
    pow3_ord = int(argv[7])
    spr2_ord = int(argv[8])
    skw_ord = int(argv[9])
    scl1_ord = int(argv[10])
    scl2_ord = int(argv[11])
    scl3_ord = int(argv[12])

    print 'spr1_ord =', spr1_ord
    print 'pow1_ord =', pow1_ord
    print 'pow2_ord =', pow2_ord
    print 'pow3_ord =', pow3_ord
    print 'spr2_ord =', spr2_ord
    print 'skw_ord =', skw_ord
    print 'scl1_ord =', scl1_ord
    print 'scl2_ord =', scl2_ord
    print 'scl3_ord =', scl3_ord,'\n'

    if comp == 'mac':
        opt_print = True
    elif comp == 'fsl':
        opt_print = False


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

            # print t[q],s[q]
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
    # p00,p10,p01,p20,p11,p02,p30,p21,p12,p03

    if spr1_ord == 3:
        spr10 = np.array([16.16,2.052,-61.49,-1.12,11.47,66.53,0.1071,-0.7078,-5.225,-24.49])
    elif spr1_ord == 1:
        spr10 = np.array([16.16,2.052,-61.49,-1.12,11.47,66.53,0.1071,-0.7078,-5.225,-24.49])
    elif spr1_ord == 2:
        spr10 = np.array([11.5,2.718,-34.83,-1.16,10.0,20.57,0.1052,-0.6225,-4.686,0.0])
    if pow1_ord == 3:
        pow10 = np.array([40.4,-18.3,-56.42,3.067,14.78,15.97,-0.14,-1.491,6.907,-4.616])
    elif pow1_ord == 1:
        pow10 = np.array([0.1634,1.32,14.67,0.0,-1.423,-6.478,0.0,0.0,0.0,0.0])
    elif pow1_ord == 2:
        pow10 = np.array([-0.261,1.511,15.01,-0.01954,-1.487,-6.58,0.0,0.0,0.0,0.0])
    if pow2_ord == 3:
        pow20 = np.array([-8.057,0.7212,95.71,0.6443,-30.43,-115.3,-0.07784,1.974,20.01,36.81])
    elif pow2_ord == 1:
        pow20 = np.array([0.2775,-0.08403,-0.9337,0.005547,0.1587,0.6667,0.0,0.0,0.0,0.0])
    elif pow2_ord == 2:
        pow20 = np.array([0.07601,-0.06402,0.007528,0.005877,0.08493,0.0,0.0,0.0,0.0,0.0])
    if pow3_ord == 3:
        pow30 = np.array([-7.561,6.675,2.328,-1.264,-2.159,4.598,0.0764,0.1978,-0.3486,-2.766])
    elif pow3_ord == 1:
        pow30 = np.array([-7.878,6.662,4.859,-1.258,-2.192,-0.3477,0.07585,0.1971,-0.3206,0.0])
    elif pow3_ord == 2:
        pow30 = np.array([-7.561,6.675,2.328,-1.264,-2.159,4.598,0.0764,0.1978,-0.3486,-2.766])
    if spr2_ord == 3:
        spr20 = np.array([5.244,-1.943,-6.702,0.2932,0.9833,4.555,-0.01451,-0.08345,0.005575,-1.745])
    elif spr2_ord == 1:
        spr20 = np.array([1.471,-0.1113,-2.713,0.0,0.08971,1.369,0.0,0.0,0.0,0.0])
    elif spr2_ord == 2:
        spr20 = np.array([4.519,-1.85,-3.583,0.3014,0.635,0.5372,-0.01591,-0.07089,0.1698,0.0])
    if skw_ord == 3:
        skw0 = np.array([-48.58,25.2,168.1,-4.004,-55.3,-143.9,0.2158,3.728,21.59,36.78])
    elif skw_ord == 1:
        skw0 = np.array([-35.13,13.69,169.7,-1.148,-53.99,-154.6,0.0,3.503,22.4,42.26])
    elif skw_ord == 2:
        skw0 = np.array([-38.41,22.41,122.4,-3.636,-51.22,-70.5,0.1959,3.671,19.07,0.0])
    if scl1_ord == 3:
        scl10 = np.array([0.7601,-0.159,0.2714,0.03194,-0.3049,-0.07724,-0.00221,0.03039,0.1698,-0.1953])
    elif scl1_ord == 1:
        scl10 = np.array([0.7559,-0.1726,0.4313,0.03529,-0.2998,-0.4329,-0.002441,0.02885,0.1777,0.0])
    elif scl1_ord == 2:
        scl10 = np.array([0.7601,-0.159,0.2714,0.03194,-0.3049,-0.07724,-0.00221,0.03039,0.1698,-0.1953])
    if scl2_ord == 3:
        scl20 = np.array([0.1466,-0.04902,-0.9488,0.004885,0.3807,0.9759,-0.00006812,-0.02005,-0.1374,-0.2565])
    elif scl2_ord == 1:
        scl20 = np.array([0.2543,0.09704,0.1465,-0.02001,0.03115,-0.002718,-0.0006896,-0.01884,-0.01955,0.0])
    elif scl2_ord == 2:
        scl20 = np.array([0.1466,-0.04902,-0.9488,0.004885,0.3807,0.9759,-0.00006812,-0.02005,-0.1374,-0.2565])
    if scl3_ord == 3:
        scl30 = np.array([-6.178,12.12,96.77,-1.478,-33.88,-94.01,0.05027,1.983,15.76,24.49])
    elif scl3_ord == 1:
        scl30 = np.array([-6.178,12.12,96.77,-1.478,-33.88,-94.01,0.05027,1.983,15.76,24.49])
    elif scl3_ord == 2:
        scl30 = np.array([13.46,-2.004,-4.436,0.007166,-0.6181,0.5103,0.0,0.0,0.0,0.0])
    # p00,p10,p01,p20,p11,p02,p30,p21,p12,p03

    param0 = np.array([])
    param0 = np.append(param0,spr10)
    param0 = np.append(param0,pow10)
    param0 = np.append(param0,pow20)
    param0 = np.append(param0,pow30)
    param0 = np.append(param0,spr20)
    param0 = np.append(param0,skw0)
    param0 = np.append(param0,scl10)
    param0 = np.append(param0,scl20)
    param0 = np.append(param0,scl30)


    """







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

        optProb.addVarGroup('param', 90, 'c', lower=None, upper=None, value=param0)

        opt = SNOPT()
        opt.setOption('Scale option',2)
        if comp == 'mac':
            opt.setOption('Print file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptSheet/SNOPT_print'+str(argv[1])+str(argv[2])+str(argv[3])+'ORD'+str(argv[4])+str(argv[5])+str(argv[6])+str(argv[7])+str(argv[8])+str(argv[9])+str(argv[10])+str(argv[11])+str(argv[12])+'.out')
            opt.setOption('Summary file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptSheet/SNOPT_summary'+str(argv[1])+str(argv[2])+str(argv[3])+'ORD'+str(argv[4])+str(argv[5])+str(argv[6])+str(argv[7])+str(argv[8])+str(argv[9])+str(argv[10])+str(argv[11])+str(argv[12])+'.out')
        elif comp == 'fsl':
            opt.setOption('Print file','/fslhome/ebtingey/compute/VAWTWakeModel/CrossVal/SNOPT_Files/SNOPT_print'+str(argv[1])+str(argv[2])+str(argv[3])+'ORD'+str(argv[4])+str(argv[5])+str(argv[6])+str(argv[7])+str(argv[8])+str(argv[9])+str(argv[10])+str(argv[11])+str(argv[12])+'.out')
            opt.setOption('Summary file','/fslhome/ebtingey/compute/VAWTWakeModel/CrossVal/SNOPT_Files/SNOPT_summary'+str(argv[1])+str(argv[2])+str(argv[3])+'ORD'+str(argv[4])+str(argv[5])+str(argv[6])+str(argv[7])+str(argv[8])+str(argv[9])+str(argv[10])+str(argv[11])+str(argv[12])+'.out')
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

    k = 0
    for i in range(9):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',res[k])')
            k += 1


    if cv_on == 0:
        velc = _velcalcqme.sheet(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8)
        velf = _velcalcqme.sheet(xtts,ysts,posdnts[0,:],posdnts[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8)

        error_curve = np.sum((velc-velodtr)**2)
        error_cv = np.sum((velf-velodts)**2)

        print '\nCross Validation Error:',error_cv
        print 'Curve Fit Error:',error_curve

        print '\ncoef0 = np.array(',coef0.tolist(),')'
        print 'coef1 = np.array(',coef1.tolist(),')'
        print 'coef2 = np.array(',coef2.tolist(),')'
        print 'coef3 = np.array(',coef3.tolist(),')'
        print 'coef4 = np.array(',coef4.tolist(),')'
        print 'coef5 = np.array(',coef5.tolist(),')'
        print 'coef6 = np.array(',coef6.tolist(),')'
        print 'coef7 = np.array(',coef7.tolist(),')'
        print 'coef8 = np.array(',coef8.tolist(),')'

        print '\nparam = np.array([',res[0],',',res[1],',',res[2],',',res[3],',',res[4],',',res[5],',',res[6],',',res[7],',',res[8],',',res[9],',\t#spr1'
        print res[10],',',res[11],',',res[12],',',res[13],',',res[14],',',res[15],',',res[16],',',res[17],',',res[18],',',res[19],',\t#pow1'
        print res[20],',',res[21],',',res[22],',',res[23],',',res[24],',',res[25],',',res[26],',',res[27],',',res[28],',',res[29],',\t#pow2'
        print res[30],',',res[31],',',res[32],',',res[33],',',res[34],',',res[35],',',res[36],',',res[37],',',res[38],',',res[39],',\t#pow3'
        print res[40],',',res[41],',',res[42],',',res[43],',',res[44],',',res[45],',',res[46],',',res[47],',',res[48],',',res[49],',\t#spr2'
        print res[50],',',res[51],',',res[52],',',res[53],',',res[54],',',res[55],',',res[56],',',res[57],',',res[58],',',res[59],',\t#skw'
        print res[60],',',res[61],',',res[62],',',res[63],',',res[64],',',res[65],',',res[66],',',res[67],',',res[68],',',res[69],',\t#scl1'
        print res[70],',',res[71],',',res[72],',',res[73],',',res[74],',',res[75],',',res[76],',',res[77],',',res[78],',',res[79],',\t#scl2'
        print res[80],',',res[81],',',res[82],',',res[83],',',res[84],',',res[85],',',res[86],',',res[87],',',res[88],',',res[89],'])\t#scl3'

    if cv_on == 1:
        velc = _velcalcqme.sheet(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,coef10,coef11,coef12,coef13,sdv_gom,spr_gom)

        error_curve = np.sum((velc-velodtr)**2)

        print '\nCurve Fit Error:',error_curve

        print '\ncoef0 = np.array(',coef0.tolist(),')'
        print 'coef1 = np.array(',coef1.tolist(),')'
        print 'coef2 = np.array(',coef2.tolist(),')'
        print 'coef3 = np.array(',coef3.tolist(),')'
        print 'coef4 = np.array(',coef4.tolist(),')'
        print 'coef5 = np.array(',coef5.tolist(),')'
        print 'coef6 = np.array(',coef6.tolist(),')'
        print 'coef7 = np.array(',coef7.tolist(),')'
        print 'coef8 = np.array(',coef8.tolist(),')'

        print '\nspr1 = np.array(',coef0.tolist(),')'
        print 'pow1 = np.array(',coef1.tolist(),')'
        print 'pow2 = np.array(',coef2.tolist(),')'
        print 'pow3 = np.array(',coef3.tolist(),')'
        print 'spr2 = np.array(',coef4.tolist(),')'
        print 'skw = np.array(',coef5.tolist(),')'
        print 'scl1 = np.array(',coef6.tolist(),')'
        print 'scl2 = np.array(',coef7.tolist(),')'
        print 'scl3 = np.array(',coef8.tolist(),')'

        print '\nparam = np.array([',res[0],',',res[1],',',res[2],',',res[3],',',res[4],',',res[5],',',res[6],',',res[7],',',res[8],',',res[9],',\t#spr1'
        print res[10],',',res[11],',',res[12],',',res[13],',',res[14],',',res[15],',',res[16],',',res[17],',',res[18],',',res[19],',\t#pow1'
        print res[20],',',res[21],',',res[22],',',res[23],',',res[24],',',res[25],',',res[26],',',res[27],',',res[28],',',res[29],',\t#pow2'
        print res[30],',',res[31],',',res[32],',',res[33],',',res[34],',',res[35],',',res[36],',',res[37],',',res[38],',',res[39],',\t#pow3'
        print res[40],',',res[41],',',res[42],',',res[43],',',res[44],',',res[45],',',res[46],',',res[47],',',res[48],',',res[49],',\t#spr2'
        print res[50],',',res[51],',',res[52],',',res[53],',',res[54],',',res[55],',',res[56],',',res[57],',',res[58],',',res[59],',\t#skw'
        print res[60],',',res[61],',',res[62],',',res[63],',',res[64],',',res[65],',',res[66],',',res[67],',',res[68],',',res[69],',\t#scl1'
        print res[70],',',res[71],',',res[72],',',res[73],',',res[74],',',res[75],',',res[76],',',res[77],',',res[78],',',res[79],',\t#scl2'
        print res[80],',',res[81],',',res[82],',',res[83],',',res[84],',',res[85],',',res[86],',',res[87],',',res[88],',',res[89],'])\t#scl3'
