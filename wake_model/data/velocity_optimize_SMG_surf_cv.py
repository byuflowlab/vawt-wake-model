from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs,log,sin,arctan,cosh
import matplotlib.pyplot as plt
from sklearn.cross_validation import train_test_split
import database_call as dbc
# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'


def veldist(dn,lat,men,sdv1,sdv2,sdv3,sdv4,rat1,bow4,spr1,bow1,bow2,bow3,scl1,scl2,scl3):

    bow = bow3*bow2*bow1*exp(bow2*dn)*exp(-bow1*exp(bow2*dn))#+bow4
    # bow = 1.#bow1*dn**2+bow2*dn+bow3

    sdv_v = sdv3*sdv2*sdv1*exp(sdv2*dn)*exp(-sdv1*exp(sdv2*dn)) + sdv4
    # sdv_v = sdv1
    rat_v = rat1#rat1*dn + rat2
    # spr_v = spr2#spr1*dn + spr2
    # spr_v = spr1*dn**2 + spr2*dn + spr3
    spr_v = spr1#spr3*spr2*spr1*exp(spr2*dn)*exp(spr1)*exp(-spr1*exp(spr2*dn)) + spr4

    f1 = -1./(sdv_v*sqrt(2.*pi))*exp(-((lat/bow)-men)**2/(2.*sdv_v**2))*(1./(1.+exp(rat_v*fabs((lat/bow))-spr_v)))
    f2 = scl3*scl2*scl1*exp(scl2*dn)*exp(scl1)*exp(-scl1*exp(scl2*dn))

    return f1*f2 + 1.

    
def obj_func(xdict):
    global posdntr
    global poslttr
    global velodtr
    
    param = xdict['param']
    funcs = {}
    
    men = param[0]
    sdv1 = param[1]
    sdv2 = param[2]
    sdv3 = param[3]
    sdv4 = param[4]
    rat1 = param[5]
    rat2 = param[6]
    spr1 = param[7]
    spr2 = param[8]
    spr3 = param[9]
    spr4 = param[10]
    scl1 = param[11]
    scl2 = param[12]
    scl3 = param[13]
    
    error = 0.

    for i in range(np.size(posdntr)):
        if posdntr[i] > 0.58:
            vel = veldist(posdntr[i],poslttr[i],men,sdv1,sdv2,sdv3,sdv4,rat1,rat2,spr1,spr2,spr3,spr4,scl1,scl2,scl3)
            error = error + (vel-velodtr[i])**2

    ##Print
    # print error

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


    return pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d


def fit(s,t,length,plot,comp,read_data,opt_print):
    global posdntr
    global poslttr
    global velodtr
    
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

    cvtest = 0.3
    posdntr,posdnts,poslttr,posltts,velodtr,velodts = train_test_split(posdn,poslt,velod,test_size=cvtest)

## Optimization
    optProb = Optimization('VAWTWake_Velo', obj_func)
    optProb.addObj('obj')

    men0 = 0.
    sdv10 = 0.5
    sdv20 = 0.1
    sdv30 = 10.
    sdv40 = 0.5
    rat0 = 10.
    spr0 = 10.
    bow1 = 0.5
    bow2 = 0.1
    bow3 = 20.
    bow4 = 1.
    # bow1 = -1.
    # bow2 = 1.
    # bow3 = 1.
    # bow4 = 1.
    
    param0 = np.array([men0,sdv10,sdv20,sdv30,sdv40,rat0,bow4,spr0,bow1,bow2,bow3,0.5,0.1,20.])

    param_l = np.array([None,0.,0.,0.,0.,None,0.,None,0.,0.,0.,0.])
    param_u = np.array([None,None,None,None,None,0.,None,0.,None,1.,1.,50.])

    param_l = np.array([None,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
    param_u = np.array([None,10.,1.,50.,None,None,None,None,1.,1.,50.,1.,1.,None])

    # param_l = np.array([None,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
    # param_u = np.array([None,None,None,None,None,None,None,None,None,None,None,1.,1.,None])

    # param_l = np.array([None,0.,0.,0.,0.,0.,None,0.,None,None,None,0.,0.,0.])
    # param_u = np.array([None,1.,1.,50.,None,None,None,None,None,None,None,1.,1.,50.])
    
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
        print paramf[10]
        print paramf[11]
        print paramf[12]
        print paramf[13]

    men = paramf[0]
    sdv1 = paramf[1]
    sdv2 = paramf[2]
    sdv3 = paramf[3]
    sdv4 = paramf[4]
    rat1 = paramf[5]
    rat2 = paramf[6]
    spr1 = paramf[7]
    spr2 = paramf[8]
    spr3 = paramf[9]
    spr4 = paramf[10]
    scl1 = paramf[11]
    scl2 = paramf[12]
    scl3 = paramf[13]

    cv_error = 0.
    for i in range(np.size(posdnts)):
        if posdnts[i] > 0.58:
            vel = veldist(posdnts[i],posltts[i],men,sdv1,sdv2,sdv3,sdv4,rat1,rat2,spr1,spr2,spr3,spr4,scl1,scl2,scl3)
            cv_error = cv_error + (vel-velodts[i])**2
    
    # men = np.array([paramf[0],paramf[1],paramf[2]])
    # spr = np.array([paramf[3],paramf[4],paramf[5],paramf[6]])
    # scl = np.array([paramf[7],paramf[8],paramf[9]])
    # rat = np.array([paramf[10],paramf[11]])
    # spr = np.array([paramf[12],paramf[13]])
    #
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
                    men_v,spr_v,scl_v,rat_v,spr_v = paramfit(xd[i],men,spr,scl,rat,spr)
                    plt.plot(veldist(xfit,men_v,spr_v,scl_v,rat_v,spr_v),xfit,'r-',linewidth=2,label=lab2)
                    plt.xlim(0.,1.5)
                    # plt.ylim(-4.,4.)
                    plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
                else:
                    exec('xfit = np.linspace(min(pos'+name+'d)-1.,max(pos'+name+'d)+1.,500)')
                    exec('plt.plot(velo'+name+'d,pos'+name+'d,color)')
                    men_v,spr_v,scl_v,rat_v,spr_v = paramfit(xd[i],men,spr,scl,rat,spr)
                    plt.plot(veldist(xfit,men_v,spr_v,scl_v,rat_v,spr_v),xfit,'r-',linewidth=2)
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
                plt.plot(veldist(xd[i],xfit,men,sdv1,sdv2,sdv3,sdv4,rat1,rat2,spr1,spr2,spr3,spr4,scl1,scl2,scl3),xfit,'r-',linewidth=2)
                plt.xlim(0.,1.5)
                # plt.ylim(-4.,4.)
                # plt.legend(loc=1)
                plt.xlabel('Normalized Velocity')
                plt.ylabel('$y/D$')
    
    return men,sdv1,sdv2,sdv3,sdv4,rat1,rat2,spr1,spr2,spr3,spr4,scl1,scl2,scl3,cv_error

## Main File
if __name__ == "__main__":
    cv = False

    s = 's4'
    t = '600'
    length = 24.

    dia = 6.
    
    comp = 'mac'
    # comp = 'fsl'

    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error0 = fit('s1','150',210.,False,comp,4,False)
    print 0,cv_error0
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error1 = fit('s1','325',165.,False,comp,4,False)
    print 1,cv_error1
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error2 = fit('s1','500',108.,False,comp,4,False)
    print 2,cv_error2
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error3 = fit('s2','200',185.,False,comp,4,False)
    print 3,cv_error3
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error4 = fit('s2','375',96.,False,comp,4,False)
    print 4,cv_error4
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error5 = fit('s2','550',60.,False,comp,4,False)
    print 5,cv_error5
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error6 = fit('s3','250',83.,False,comp,4,False)
    print 6,cv_error6
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error7 = fit('s3','425',41.,False,comp,4,False)
    print 7,cv_error7
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error8 = fit('s3','600',28.,False,comp,4,False)
    print 8,cv_error8
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error9 = fit('s4','300',42.,False,comp,4,False)
    print 9,cv_error9
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error10 = fit('s4','475',24.,False,comp,4,False)
    print 10,cv_error10
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error11 = fit('s4','650',22.,False,comp,4,False)
    print 11,cv_error11
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error12 = fit('s5','350',26.,False,comp,4,False)
    print 12,cv_error12
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error13 = fit('s5','525',19.,False,comp,4,False)
    print 13,cv_error13
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,cv_error14 = fit('s5','700',14.,False,comp,4,False)
    print 14,cv_error14
    
    # print '\n'
    # print 'men =',men
    # print 'sdv1 =',sdv1
    # print 'sdv2 =',sdv2
    # print 'sdv3 =',sdv3
    # print 'sdv4 =',sdv4
    # print 'rat1 =',rat1
    # print 'rat2 =',rat2
    # print 'spr1 =',spr1
    # print 'spr2 =',spr2
    # print 'spr3 =',spr3
    # print 'spr4 =',spr4
    # print 'scl1 =',scl1
    # print 'scl2 =',scl2
    # print 'scl3 =',scl3

    cv_errortot = np.average([cv_error0, cv_error1, cv_error2, cv_error3, cv_error4, cv_error5, cv_error6, cv_error7, cv_error8, cv_error9, cv_error10, cv_error11, cv_error12, cv_error13,cv_error14])


    print 'error:',cv_errortot

plt.show()