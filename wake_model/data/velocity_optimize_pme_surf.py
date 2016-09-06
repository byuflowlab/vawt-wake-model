from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs,log,sin,arctan,cosh
import matplotlib.pyplot as plt
import database_call as dbc
# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'


def veldist(dn,lat,spr1,pow1,pow2,pow3,spr2,skw,odr,scl1,scl2,scl3):
    pow = pow1-pow2*dn**pow3
    exp_v = exp(-spr1*fabs(lat)**pow)
    poly = spr2*fabs(lat-skw)**odr-1.
    scl_v = scl3*scl2*scl1*exp(scl2*dn)*exp(-scl1*exp(scl2*dn))

    return exp_v*poly*scl_v+1.
    # return exp_v*scl_v + 1.
    
def obj_func(xdict):
    global posdn
    global poslt
    global velod
    
    param = xdict['param']
    funcs = {}
    
    spr1 = param[0]
    pow1 = param[1]
    pow2 = param[2]
    pow3 = param[3]
    spr2 = param[4]
    skw = param[5]
    odr = param[6]
    scl1 = param[7]
    scl2 = param[8]
    scl3 = param[9]
    
    error = 0.

    for i in range(np.size(posdn)):
        if posdn[i] > 0.58:
            vel = veldist(posdn[i],poslt[i],spr1,pow1,pow2,pow3,spr2,skw,odr,scl1,scl2,scl3)
            error = error + (vel-velod[i])**2

    ##Print
    print error

    funcs['obj'] = error
    
    
    
    fail = False

    return funcs, fail

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


def starccm_read2(fdata,dia,windd,opt_print):
    for q in range(30):
        name = str(q+1)
        exec('pos'+name+'d = np.array([])')
        exec('velo'+name+'d = np.array([])')

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


        for i in range(30):
            name = str(i+1)

            #Ordering the data numerically by the position
            exec('pos'+name+', velo'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', velo'+name+'))))')
            #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
            exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)/dia\nvelo'+name+' = np.copy(velo'+name+'_0)/wind')
            #Deleting wall boundary data
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5. or pos'+name+'[j] < -5.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            #Deleting values greater than 1*wind
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif velo'+name+'[j] > 1. or fabs(pos'+name+'[j]) > 1.2:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            exec('pos'+name+'d = np.append(pos'+name+'d,pos'+name+')\nvelo'+name+'d = np.append(velo'+name+'d,velo'+name+')')


    # return pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d
    fac = 1.
    return pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d*fac,velo2d*fac,velo3d*fac,velo4d*fac,velo5d*fac,velo6d*fac,velo7d*fac,velo8d*fac,velo9d*fac,velo10d*fac,velo11d*fac,velo12d*fac,velo13d*fac,velo14d*fac,velo15d*fac,velo16d*fac,velo17d*fac,velo18d*fac,velo19d*fac,velo20d*fac,velo21d*fac,velo22d*fac,velo23d*fac,velo24d*fac,velo25d*fac,velo26d*fac,velo27d*fac,velo28d*fac,velo29d*fac,velo30d*fac


def fit(s,t,length,plot,comp,read_data,opt_print):
    global posdn
    global poslt
    global velod
    
    t2 = t+'.0'

    wfit = s+'_'+t2
    wfit2 = s+'_'+t2
    wfit3 = s+'_'+t2
    wfit4 = s+'_'+t2
    wfit5 = s+'_'+t2
    wfit6 = s+'_'+t2

    length2 = length
    length3 = length
    length4 = length
    length5 = length
    length6 = length
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
    elif comp == 'win':
        fdata = 'C://Users//TingeyPC//Documents//zStar-CCM//STAR-CCM//NACA0021//MoveForward//Velocity Sections//'+wfit+'.csv'
        fdata2 = 'C://Users//TingeyPC//Documents//zStar-CCM//STAR-CCM//NACA0021//MoveForward//CrossValidate//vel14//Velocity//'+wfit2+'.csv'
        fdata3 = 'C://Users//TingeyPC//Documents//zStar-CCM//STAR-CCM//NACA0021//MoveForward//CrossValidate//vel12//Velocity//'+wfit3+'.csv'
        fdata4 = 'C://Users//TingeyPC//Documents//zStar-CCM//STAR-CCM//NACA0021//MoveForward//CrossValidate//vel16//Velocity//'+wfit4+'.csv'
        fdata5 = 'C://Users//TingeyPC//Documents//zStar-CCM//STAR-CCM//NACA0021//MoveForward//CrossValidate//rot17//Velocity//'+wfit5+'.csv'
        fdata6 = 'C://Users//TingeyPC//Documents//zStar-CCM//STAR-CCM//NACA0021//MoveForward//CrossValidate//rot18//Velocity//'+wfit6+'.csv'



    if read_data ==1:
        posdn,poslt,velod = starccm_read(np.array([fdata]),dia,np.array([wind]),length,opt_print)
    if read_data ==2:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2]),dia,np.array([wind,wind2]),length,opt_print)
    if read_data ==3:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3]),dia,np.array([wind,wind2,wind3]),length,opt_print)
    if read_data ==4:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4]),dia,np.array([wind,wind2,wind3,wind4]),length,opt_print)
    if read_data ==5:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5]),dia,np.array([wind,wind2,wind3,wind4,wind5]),length,opt_print)
    if read_data ==6:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5,fdata6]),dia,np.array([wind,wind2,wind3,wind4,wind5,wind6]),length,opt_print)

    if plot == True:
        if read_data ==1:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata]),dia,np.array([wind]),opt_print)
        if read_data ==2:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2]),dia,np.array([wind,wind2]),opt_print)
        if read_data ==3:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3]),dia,np.array([wind,wind2,wind3]),opt_print)
        if read_data ==4:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3,fdata4]),dia,np.array([wind,wind2,wind3,wind4]),opt_print)
        if read_data ==5:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3,fdata4,fdata5]),dia,np.array([wind,wind2,wind3,wind4,wind5]),opt_print)
        if read_data ==6:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3,fdata4,fdata5,fdata6]),dia,np.array([wind,wind2,wind3,wind4,wind5,wind6]),opt_print)

        start = length/30.
        xd = np.linspace(start,length,30)/dia



## Optimization
    optProb = Optimization('VAWTWake_Velo', obj_func)
    optProb.addObj('obj')

    spr10 = 10.0
    pow10 = 5.0
    pow20 = 0.5
    pow30 = 1.0
    spr20 = 2.0
    skw0 = 0.0
    odr0 = 2.0
    scl10 = 0.5
    scl20 = 0.1
    scl30 = 20.0

    # spr10 = 213.8593169
    # pow10 = 10.39210953
    # pow20 = 2.086951239
    # pow30 = 0.035659319
    # spr20 = 0.007589688
    # skw0 = 10.63462155
    # odr0 = 2.0
    # scl10 = 0.537566448
    # scl20 = 0.041077603
    # scl30 = 56.74689143

    param0 = np.array([spr10,pow10,pow20,pow30,spr20,skw0,odr0,scl10,scl20,scl30])

    param_l = np.array([0.,0.,0.,0.,0.,None,0.,0.,0.,0.])
    param_u = np.array([None,None,None,None,None,None,None,1.,1.,None])

    nparam = np.size(param0)
    optProb.addVarGroup('param', nparam, 'c', lower=param_l, upper=param_u, value=param0)
    
    opt = SNOPT()
    opt.setOption('Scale option',2)
    if comp == 'mac':
        opt.setOption('Print file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptVel/SNOPT_print'+s+'_'+t+'.out')
        opt.setOption('Summary file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptVel/SNOPT_summary'+s+'_'+t+'.out')
    elif comp == 'fsl':
        opt.setOption('Print file','/fslhome/ebtingey/compute/VAWTWakeModel/OptVel/SNOPT_print'+s+'_'+t+'.out')
        opt.setOption('Summary file','/fslhome/ebtingey/compute/VAWTWakeModel/OptVel/SNOPT_summary'+s+'_'+t+'.out')
    elif comp == 'win':
        opt.setOption('Print file','C://Users//TingeyPC//Documents//FLOW Lab//VAWTWakeModel//wake_model//data//optVel//SNOPT_print'+s+'_'+t+'.out')
        opt.setOption('Summary file','C://Users//TingeyPC//Documents//FLOW Lab//VAWTWakeModel//wake_model//data//OptVel//SNOPT_summary'+s+'_'+t+'.out')
    res = opt(optProb, sens=None)
    if opt_print == True:
        print res
    
    pow = res.fStar
    paramf = res.xStar['param']
    if opt_print == True:
        print paramf[0]
        print paramf[1]
        print paramf[2]
        print paramf[3]
        print paramf[4]
        print paramf[5]
        print paramf[6]
        print paramf[7]
        print paramf[8]
        print paramf[9]


    spr1 = paramf[0]
    pow1 = paramf[1]
    pow2 = paramf[2]
    pow3 = paramf[3]
    spr2 = paramf[4]
    skw = paramf[5]
    odr = paramf[6]
    scl1 = paramf[7]
    scl2 = paramf[8]
    scl3 = paramf[9]

    paper = False

    if plot == True:
        if paper == True:
            for i in range(30):
                name = str(i+1)
                ind = str(i)
                plt.figure(1)
                ax1 = plt.subplot(5,6,i+1)
                color = 'bo'
                color2 = 'r-'
                fs = 15
                lab = 'CFD'
                lab2 = 'Trend'
                tex = '$x/D$ = '+str("{0:.2f}".format(x[i]/dia))
                exec('xfit = np.linspace(min(pos'+name+'/dia)-1.,max(pos'+name+'/dia)+1.,500)')
                if i == 5:
                    exec('xfit = np.linspace(min(pos'+name+'d)-1.,max(pos'+name+'d)+1.,500)')
                    exec('plt.plot(velo'+name+'d,pos'+name+'d,color,label=lab)')
                    skw_v,spr_v,scl_v,rat_v,spr_v = paramfit(xd[i],skw,spr,scl,rat,spr)
                    plt.plot(veldist(xfit,skw_v,spr_v,scl_v,rat_v,spr_v),xfit,'r-',linewidth=2,label=lab2)
                    plt.xlim(0.,1.5)
                    # plt.ylim(-4.,4.)
                    plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
                else:
                    exec('xfit = np.linspace(min(pos'+name+'d)-1.,max(pos'+name+'d)+1.,500)')
                    exec('plt.plot(velo'+name+'d,pos'+name+'d,color)')
                    skw_v,spr_v,scl_v,rat_v,spr_v = paramfit(xd[i],skw,spr,scl,rat,spr)
                    plt.plot(veldist(xfit,skw_v,spr_v,scl_v,rat_v,spr_v),xfit,'r-',linewidth=2)
                    plt.xlim(0.,1.5)
                    # plt.ylim(-4.,4.)
                plt.text(0.3,0.8,tex,fontsize=fs)
                if i <= 23:
                    plt.setp(ax1.get_xticklabels(), visible=False)
                else:
                    plt.xlabel('$y/D$',fontsize=fs)
                    plt.xticks(fontsize=fs)
                if i == 0 or i == 6 or i == 12 or i == 18 or i ==24:
                    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
                    plt.yticks(fontsize=fs)
                else:
                    plt.setp(ax1.get_yticklabels(), visible=False)


        elif paper == False:
            for i in range(30):
                name = str(i+1)
                plt.figure(1)
                plt.subplot(5,6,i+1)
                color = 'bo'
                exec('xfit = np.linspace(min(pos'+name+'d)-1.,max(pos'+name+'d)+1.,500)')
                exec('plt.plot(velo'+name+'d,pos'+name+'d,color)')
                plt.plot(veldist(xd[i],xfit,spr1,pow1,pow2,pow3,spr2,skw,odr,scl1,scl2,scl3),xfit,'r-',linewidth=2)
                plt.xlim(0.,1.5)
                # plt.ylim(-4.,4.)
                # plt.legend(loc=1)
                plt.xlabel('Normalized Velocity')
                plt.ylabel('$y/D$')
    
    return spr1,pow1,pow2,pow3,spr2,skw,odr,scl1,scl2,scl3

## Main File
if __name__ == "__main__":
    cv = False

# 1.5   1.75	2	    2.25	2.5	    2.75	3	3.25	3.5	    3.75	4	4.25	4.5	    4.75	5	5.25	5.5	    5.75	6	6.25	6.5	    6.75	7
# 210	210	    205	    196	    185	    178	    170	165	    160	    145	    140	123	    115	    112	    108	101	    101	    90	    85	80	    78	    75	    70
# 197	193	    185	    176	    140	    146	    126	114	    103	    96	    100	86	    77	    72	    70	68	    60	    64	    54	50	    47	    45	    44
# 185	150	    100	    95	    83	    76	    72	63	    60	    49	    50	41	    39	    36	    34	33	    30	    31	    28	30	    29	    28  	27
# 145	100	    73	    60	    53	    44	    42	37	    38	    30	    33	26	    22	    24	    23	21	    21	    19	    24	23	    22	    21	    20
# 78	70	    52	    43	    37	    32	    29	27	    26	    23	    20	20	    23	    21	    20	19	    19	    18	    18	16	    16	    15	    14



    s = 's4'
    t = '150'
    length = 145.

    dia = 6.
    
    comp = 'mac'
    # comp = 'fsl'
    # comp = 'win'

    spr1,pow1,pow2,pow3,spr2,skw,odr,scl1,scl2,scl3 = fit(s,t,length,True,comp,4,True)
    
    print '\n'
    print 'spr1 =',spr1
    print 'pow1 =',pow1
    print 'pow2 =',pow2
    print 'pow3 =',pow3
    print 'spr2 =',spr2
    print 'skw =',skw
    print 'odr =',odr
    print 'scl1 =',scl1
    print 'scl2 =',scl2
    print 'scl3 =',scl3

plt.show()