from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs,log,sin,arctan,cosh
from sklearn.cross_validation import train_test_split
import matplotlib.pyplot as plt
import database_call as dbc
from sys import argv
# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'


def overlay(xt,ys,ordt,ords,coef):
    ntsr = np.size(ordt)
    nsol = np.size(ords)

    fun = 'f = 0.'
    k = 0
    for i in range(ntsr):
        for j in range(nsol):
            pow = '+'+str(coef[k])+'*xt**'+str(ntsr-1-i)+'*ys**'+str(nsol-1-i)
            fun = fun+pow
            k += 1

    exec(fun)
    return f



def veldist(dn,lat,men,sdv1,sdv2,sdv3,sdv4,rat,tns,spr1,spr2,spr3,spr4,scl1,scl2,scl3,sdv_gom,spr_gom):

    if sdv_gom == 0:
        sdv_v = sdv3*sdv2*sdv1*exp(sdv2*dn)*exp(-sdv1*exp(sdv2*dn))+sdv4
    elif sdv_gom == 1:
        sdv_v = sdv1

    if spr_gom == 0:
        spr_v = spr3*spr2*spr1*exp(spr2*dn)*exp(-spr1*exp(spr2*dn))+spr4
    elif spr_gom == 1:
        spr_v = 1.

    f1 = -1./(sdv_v*sqrt(2.*pi))*exp(-((lat/spr_v)-men)**2/(2.*sdv_v**2))*(1./(1.+exp(rat*fabs((lat/spr_v))-tns)))
    f2 = scl3*scl2*scl1*exp(scl2*dn)*exp(-scl1*exp(scl2*dn))

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


def obj_func(xdict):
    for i in range(80):
        name = str(i+1)
        exec('global posdn'+name+'tr')
        exec('global poslt'+name+'tr')
        exec('global velod'+name+'tr')
        exec('global xt'+name+'tr')
        exec('global ys'+name+'tr')
    global sdv_gom
    global spr_gom
    global ordt
    global ords

    menc = xdict['men']
    sdv1c = xdict['sdv1']
    sdv2c = xdict['sdv2']
    sdv3c = xdict['sdv3']
    sdv4c = xdict['sdv4']
    ratc = xdict['rat']
    tnsc = xdict['tns']
    spr1c = xdict['spr1']
    spr2c = xdict['spr2']
    spr3c = xdict['spr3']
    spr4c = xdict['spr4']
    scl1c = xdict['scl1']
    scl2c = xdict['scl2']
    scl3c = xdict['scl3']

    funcs = {}

    error = 0.
    for i in range(80):
        name = str(i+1)
        exec('mend = overlay(xt'+name+'tr,ys'+name+'tr,3,ords,menc)')
        exec('sdv1d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,sdv1c)')
        exec('sdv2d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,sdv2c)')
        exec('sdv3d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,sdv3c)')
        exec('sdv4d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,sdv4c)')
        exec('ratd = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,ratc)')
        exec('tnsd = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,tnsc)')
        exec('spr1d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,spr1c)')
        exec('spr2d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,spr2c)')
        exec('spr3d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,spr3c)')
        exec('spr4d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,spr4c)')
        exec('scl1d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,scl1c)')
        exec('scl2d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,scl2c)')
        exec('scl3d = overlay(xt'+name+'tr,ys'+name+'tr,ordt,ords,scl3c)')

        exec('for i in range(np.size(posdn'+name+'tr)):\n\tif posdn'+name+'tr[i] > 0.58:\n\t\tvel = veldist(posdn'+name+'tr[i],poslt'+name+'tr[i],mend,sdv1d,sdv2d,sdv3d,sdv4d,ratd,tnsd,spr1d,spr2d,spr3d,spr4d,scl1d,scl2d,scl3d,sdv_gom,spr_gom)\n\t\terror = error + (vel-velod'+name+'tr[i])**2')

    print error
    funcs['obj'] = error

    fail = False

    return funcs, fail

if __name__ == "__main__":

    comp = 'mac'
    # comp = 'fsl'

    global sdv_gom
    global spr_gom
    global ordt
    global ords

    if comp == 'mac':
        opt_print = True
        sdv_gom = int(argv[3])
        spr_gom = int(argv[4])
        print 'TSR Order:',int(argv[1])-1
        print 'Solidity Order:',int(argv[2])-1
        if int(argv[3]) == 0:
            print 'sdv_gom is True'
        elif int(argv[3]) == 1:
            print 'sdv_gom is False'
        if int(argv[4]) == 0:
            print 'spr_gom is True'
        elif int(argv[4]) == 1:
            print 'spr_gom is False'
    elif comp == 'fsl':
        opt_print = False
        sdv_gom = int(argv[3])
        spr_gom = int(argv[4])
        print 'TSR Order:',int(argv[1])-1
        print 'Solidity Order:',int(argv[2])-1
        if int(argv[3]) == 0:
            print 'sdv_gom is True'
        elif int(argv[3]) == 1:
            print 'sdv_gom is False'
        if int(argv[4]) == 0:
            print 'spr_gom is True'
        elif int(argv[4]) == 1:
            print 'spr_gom is False'
    read_data = 1

    for i in range(80):
        name = str(i+1)
        exec('global posdn'+name+'tr')
        exec('global poslt'+name+'tr')
        exec('global velod'+name+'tr')
        exec('global xt'+name+'tr')
        exec('global ys'+name+'tr')


    s1length = np.array([210,210,205,196,185,178,170,165,160,145,140,123,115,112,108,101,101,90,85,80,78,75,70])
    s2length = np.array([197,193,185,176,140,146,126,114,103,96,100,86,77,72,70,68,60,64,54,50,47,45,44])
    s3length = np.array([185,150,100,95,83,76,72,63,60,49,50,41,39,36,34,33,30,31,28,30,29,28,27])
    s4length = np.array([145,100,73,60,53,44,42,37,38,30,33,26,22,24,23,21,21,19,24,23,22,21,20])
    s5length = np.array([78,70,52,43,37,32,29,27,26,23,20,20,23,21,20,19,19,18,18,16,16,15,14])
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

    # Optimization
    optProb = Optimization('VAWTWake_Sheet', obj_func)
    optProb.addObj('obj')

    if comp == 'mac':
        ordt = int(argv[1])
        ords = int(argv[2])
    elif comp == 'fsl':
        ordt = int(argv[1])
        ords = int(argv[2])

    men0 = np.linspace(0.,0.,3*ords)
    sdv10 = np.linspace(0.,0.,ordt*ords)
    sdv20 = np.linspace(0.,0.,ordt*ords)
    sdv30 = np.linspace(0.,0.,ordt*ords)
    sdv40 = np.linspace(0.,0.,ordt*ords)
    rat0 = np.linspace(0.,0.,ordt*ords)
    tns0 = np.linspace(0.,0.,ordt*ords)
    spr10 = np.linspace(0.,0.,ordt*ords)
    spr20 = np.linspace(0.,0.,ordt*ords)
    spr30 = np.linspace(0.,0.,ordt*ords)
    spr40 = np.linspace(0.,0.,ordt*ords)
    scl10 = np.linspace(0.,0.,ordt*ords)
    scl20 = np.linspace(0.,0.,ordt*ords)
    scl30 = np.linspace(0.,0.,ordt*ords)

    men0[-1] = 0.
    sdv10[-1] = 0.1
    sdv20[-1] = 0.5
    sdv30[-1] = 10.
    sdv40[-1] = 0.5
    rat0[-1] = 10.
    tns0[-1] = 10.
    spr10[-1] = 0.1
    spr20[-1] = 0.5
    spr30[-1] = 10.
    spr40[-1] = 0.5
    scl10[-1] = 0.1
    scl20[-1] = 0.5
    scl30[-1] = 10.


    optProb.addVarGroup('men', 3*ords, 'c', lower=None, upper=None, value=men0)
    optProb.addVarGroup('sdv1', ordt*ords, 'c', lower=None, upper=None, value=sdv10)
    optProb.addVarGroup('sdv2', ordt*ords, 'c', lower=None, upper=None, value=sdv20)
    optProb.addVarGroup('sdv3', ordt*ords, 'c', lower=None, upper=None, value=sdv30)
    optProb.addVarGroup('sdv4', ordt*ords, 'c', lower=None, upper=None, value=sdv40)
    optProb.addVarGroup('rat', ordt*ords, 'c', lower=None, upper=None, value=rat0)
    optProb.addVarGroup('tns', ordt*ords, 'c', lower=None, upper=None, value=tns0)
    optProb.addVarGroup('spr1', ordt*ords, 'c', lower=None, upper=None, value=spr10)
    optProb.addVarGroup('spr2', ordt*ords, 'c', lower=None, upper=None, value=spr20)
    optProb.addVarGroup('spr3', ordt*ords, 'c', lower=None, upper=None, value=spr30)
    optProb.addVarGroup('spr4', ordt*ords, 'c', lower=None, upper=None, value=spr40)
    optProb.addVarGroup('scl1', ordt*ords, 'c', lower=None, upper=None, value=scl10)
    optProb.addVarGroup('scl2', ordt*ords, 'c', lower=None, upper=None, value=scl20)
    optProb.addVarGroup('scl3', ordt*ords, 'c', lower=None, upper=None, value=scl30)


    opt = SNOPT()
    opt.setOption('Scale option',0)
    if comp == 'mac':
        opt.setOption('Print file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptSheet/cv_cvx_'+argv[1]+argv[2]+argv[3]+argv[4]+'_SNOPT_print.out')
        opt.setOption('Summary file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptSheet/cv_cvx_SNOPT_'+argv[1]+argv[2]+argv[3]+argv[4]+'_summary.out')
    elif comp == 'fsl':
        opt.setOption('Print file','/fslhome/ebtingey/compute/VAWTWakeModel/OptSheet/cv_cvx_'+argv[1]+argv[2]+argv[3]+argv[4]+'_SNOPT_OptSheet_print.out')
        opt.setOption('Summary file','/fslhome/ebtingey/compute/VAWTWakeModel/OptSheet/cv_cvx_'+argv[1]+argv[2]+argv[3]+argv[4]+'_SNOPT_OptSheet_summary.out')
    res = opt(optProb, sens=None)
    if opt_print == True:
        print res

    error = res.fStar

    mencf = res.xStar['men']
    sdv1cf = res.xStar['sdv1']
    sdv2cf = res.xStar['sdv2']
    sdv3cf = res.xStar['sdv3']
    sdv4cf = res.xStar['sdv4']
    ratcf = res.xStar['rat']
    tnscf = res.xStar['tns']
    spr1cf = res.xStar['spr1']
    spr2cf = res.xStar['spr2']
    spr3cf = res.xStar['spr3']
    spr4cf = res.xStar['spr4']
    scl1cf = res.xStar['scl1']
    scl2cf = res.xStar['scl2']
    scl3cf = res.xStar['scl3']

    error_cv = 0.
    for i in range(35):
        name = str(i+1)
        exec('mend = overlay(xt'+name+'ts,ys'+name+'ts,3,ords,mencf)')
        exec('sdv1d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,sdv1cf)')
        exec('sdv2d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,sdv2cf)')
        exec('sdv3d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,sdv3cf)')
        exec('sdv4d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,sdv4cf)')
        exec('ratd = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,ratcf)')
        exec('tnsd = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,tnscf)')
        exec('spr1d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,spr1cf)')
        exec('spr2d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,spr2cf)')
        exec('spr3d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,spr3cf)')
        exec('spr4d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,spr4cf)')
        exec('scl1d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,scl1cf)')
        exec('scl2d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,scl2cf)')
        exec('scl3d = overlay(xt'+name+'ts,ys'+name+'ts,ordt,ords,scl3cf)')

        exec('for i in range(np.size(posdn'+name+'ts)):\n\tif posdn'+name+'ts[i] > 0.58:\n\t\tvel = veldist(posdn'+name+'ts[i],poslt'+name+'ts[i],mend,sdv1d,sdv2d,sdv3d,sdv4d,ratd,tnsd,spr1d,spr2d,spr3d,spr4d,scl1d,scl2d,scl3d,sdv_gom,spr_gom)\n\t\terror_cv = error_cv + (vel-velod'+name+'ts[i])**2')



    print '\n'
    print 'Total Error:',error,'\n'

    print 'men = np.array(',mencf.tolist(),')'
    print 'sdv1 = np.array(',sdv1cf.tolist(),')'
    print 'sdv2 = np.array(',sdv2cf.tolist(),')'
    print 'sdv3 = np.array(',sdv3cf.tolist(),')'
    print 'sdv4 = np.array(',sdv4cf.tolist(),')'
    print 'rat = np.array(',ratcf.tolist(),')'
    print 'tns = np.array(',tnscf.tolist(),')'
    print 'spr1 = np.array(',spr1cf.tolist(),')'
    print 'spr2 = np.array(',spr2cf.tolist(),')'
    print 'spr3 = np.array(',spr3cf.tolist(),')'
    print 'spr4 = np.array(',spr4cf.tolist(),')'
    print 'scl1 = np.array(',scl1cf.tolist(),')'
    print 'scl2 = np.array(',scl2cf.tolist(),')'
    print 'scl3 = np.array(',scl3cf.tolist(),')'

    print '\nCross Validation Error:',error_cv




