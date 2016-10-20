from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs,log,sin,arctan,cosh
from scipy.special import erf
import matplotlib.pyplot as plt
# import database_call as dbc
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


def veldist(dn,lat,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3):
    loc = loc1*dn**2 + loc2*dn + loc3
    spr = spr1*dn + spr2
    skw = skw1*dn + skw2
    scl = scl1/(1.0 + exp(scl2*(dn - scl3)))

    gam1 = scl*skw/2.0*exp(skw/2.0*(2.0*loc+skw*spr**2.0-2.0*lat))*(1.0-erf((loc + skw*spr**2.0 - lat)/(sqrt(2.0)*spr)))

    loc = -loc
    spr = -spr
    skw = -skw
    scl = -scl

    gam2 = scl*skw/2.0*exp(skw/2.0*(2.0*loc+skw*spr**2.0-2.0*lat))*(1.0-erf((loc + skw*spr**2.0 - lat)/(sqrt(2.0)*spr)))

    return gam1-gam2



def obj_func(xdict):
    global posdn
    global poslt
    global velod

    param = xdict['param']
    funcs = {}

    loc1 = param[0]
    loc2 = param[1]
    loc3 = param[2]
    spr1 = param[3]
    spr2 = param[4]
    skw1 = param[5]
    skw2 = param[6]
    scl1 = param[7]
    scl2 = param[8]
    scl3 = param[9]

    error = 0.

    for i in range(np.size(posdn)):
        if posdn[i] > 0.58:
            vel = veldist(posdn[i],poslt[i],loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3)
            error = error + (vel-velod[i])**2

    ##Print
    print error

    funcs['obj'] = error



    fail = False

    return funcs, fail

def starccm_read(fdata,dia,rotd,length,opt_print):
    start = length/30.
    lendat =  np.linspace(start,length,30)/dia

    for j in range(np.size(fdata)):
        for k in range(30):
            name = str(k+1)
            exec('pos'+name+' = np.array([])')
            exec('velo'+name+' = np.array([])')
        rot = rotd[j]

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
            exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)/dia\nvelo'+name+' = np.copy(velo'+name+'_0)/rot')
            #Deleting wall boundary data
            # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5. or pos'+name+'[j] < -5.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            #Deleting values greater than 1*wind
            # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif velo'+name+'[j] > 1. or fabs(pos'+name+'[j]) > 1.2:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            exec('lensize = np.size(pos'+name+')')
            exec('posdn = np.append(posdn,np.ones(lensize)*lendat['+ind+'])\nposlt = np.append(poslt,pos'+name+')\nvelod = np.append(velod,velo'+name+')')
        #Deleting values inside turbine region (and shortly after)
        exec('delvec = np.array([])\nfor j in range(np.size(posdn)):\n\tif posdn[j] < 0.58 or poslt[j] > 1.08 or poslt[j] < -1.08:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tposdn = np.delete(posdn,delvec[j])\n\tposlt = np.delete(poslt,delvec[j])\n\tvelod = np.delete(velod,delvec[j])')


    return posdn,poslt,velod


def starccm_read2(fdata,dia,rotd,opt_print):
    for q in range(30):
        name = str(q+1)
        exec('pos'+name+'d = np.array([])')
        exec('velo'+name+'d = np.array([])')

    for j in range(np.size(fdata)):
        for k in range(30):
            name = str(k+1)
            exec('pos'+name+' = np.array([])')
            exec('velo'+name+' = np.array([])')
        rot = rotd[j]

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
            exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)/dia\nvelo'+name+' = np.copy(velo'+name+'_0)/rot')
            #Deleting wall boundary data
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 1.08 or pos'+name+'[j] < -1.08:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            #Deleting values greater than 1*wind
            # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif velo'+name+'[j] > 1. or fabs(pos'+name+'[j]) > 1.2:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            exec('pos'+name+'d = np.append(pos'+name+'d,pos'+name+')\nvelo'+name+'d = np.append(velo'+name+'d,velo'+name+')')



    # return pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d
    fac = 1.
    return pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d*fac,velo2d*fac,velo3d*fac,velo4d*fac,velo5d*fac,velo6d*fac,velo7d*fac,velo8d*fac,velo9d*fac,velo10d*fac,velo11d*fac,velo12d*fac,velo13d*fac,velo14d*fac,velo15d*fac,velo16d*fac,velo17d*fac,velo18d*fac,velo19d*fac,velo20d*fac,velo21d*fac,velo22d*fac,velo23d*fac,velo24d*fac,velo25d*fac,velo26d*fac,velo27d*fac,velo28d*fac,velo29d*fac,velo30d*fac


def fit(s,t,length,plot,comp,read_data,opt_print):
    global posdn
    global poslt
    global velod

    t2 = t+'.0'

    wfit = s+'_'+t
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
        fdata = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/Vorticity Sections/'+wfit+'.csv'
        fdata2 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel14/Vorticity/'+wfit2+'t.csv'
        fdata3 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel12/Vorticity/'+wfit3+'t.csv'
        fdata4 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel16/Vorticity/'+wfit4+'t.csv'
        fdata5 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot17/Vorticity/'+wfit5+'t.csv'
        fdata6 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot18/Vorticity/'+wfit6+'t.csv'
    elif comp == 'fsl':
        fdata = '/fslhome/ebtingey/compute/moveForward/Vorticity/'+wfit+'.csv'
        fdata2 = '/fslhome/ebtingey/compute/moveForward/vel14/Vorticity/'+wfit2+'t.csv'
        fdata3 = '/fslhome/ebtingey/compute/moveForward/vel12/Vorticity/'+wfit3+'t.csv'
        fdata4 = '/fslhome/ebtingey/compute/moveForward/vel16/Vorticity/'+wfit4+'t.csv'
        fdata5 = '/fslhome/ebtingey/compute/moveForward/rot17/Vorticity/'+wfit5+'.csv'
        fdata6 = '/fslhome/ebtingey/compute/moveForward/rot18/Vorticity/'+wfit6+'.csv'



    if read_data ==1:
        posdn,poslt,velod = starccm_read(np.array([fdata]),dia,np.array([rot]),length,opt_print)
    if read_data ==2:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2]),dia,np.array([rot,rot2]),length,opt_print)
    if read_data ==3:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3]),dia,np.array([rot,rot2,rot3]),length,opt_print)
    if read_data ==4:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4]),dia,np.array([rot,rot2,rot3,rot4]),length,opt_print)
    if read_data ==5:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5]),dia,np.array([rot,rot2,rot3,rot4,rot5]),length,opt_print)
    if read_data ==6:
        posdn,poslt,velod = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5,fdata6]),dia,np.array([rot,rot2,rot3,rot4,rot5,rot6]),length,opt_print)

    if plot == True:
        if read_data ==1:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata]),dia,np.array([rot]),opt_print)
        if read_data ==2:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2]),dia,np.array([rot,rot2]),opt_print)
        if read_data ==3:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3]),dia,np.array([rot,rot2,rot3]),opt_print)
        if read_data ==4:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3,fdata4]),dia,np.array([rot,rot2,rot3,rot4]),opt_print)
        if read_data ==5:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3,fdata4,fdata5]),dia,np.array([rot,rot2,rot3,rot4,rot5]),opt_print)
        if read_data ==6:
            pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read2(np.array([fdata,fdata2,fdata3,fdata4,fdata5,fdata6]),dia,np.array([rot,rot2,rot3,rot4,rot5,rot6]),opt_print)

        start = length/30.
        xd = np.linspace(start,length,30)/dia



## Optimization
    optProb = Optimization('VAWTWake_Velo', obj_func)
    optProb.addObj('obj')

    coef0 = np.array( [0.0025703809856661534, -0.0007386258659065129, 0.004595508188667984, 0.000380123563204793, -0.0005090098755683027, 0.005744581813281894, -4.103393770815313e-05, -0.0014146918534486358, -0.013975958482495927, 0.0] )
    coef1 = np.array( [-0.5047504670963536, 0.23477391362058556, 0.8414256436198028, -0.04252528528617351, -0.06962875967504166, -0.6566907653208429, 0.002839318332370807, 0.00571803958194812, 0.0070744372783060295, 0.22805286438890995] )
    coef2 = np.array( [0.2878345841026334, 0.11512552658662782, 0.7303949879914625, -0.007035517839387948, -0.18284850673545897, -0.5241921153256568, -0.0003704899921255296, 0.010972527139685873, 0.04380801537377295, 0.1724129349605399] )
    coef3 = np.array( [0.08234816067475287, -0.03530687906626052, -0.3662863944976986, 0.003240141344532779, 0.12172015102204112, 0.2993048183466721, 0.0, -0.009253185586804007, -0.057469126406649716, -0.07257633583877886] )
    coef4 = np.array( [-0.07083579909945328, 0.016182024377569406, 0.1985436342461859, 0.0017738254727425816, -0.09111094817943823, -0.06561408122153217, -0.0005115133402638633, 0.009434288536679505, 0.022392136905926813, 0.0] )
    coef5 = np.array( [-1.6712830849073221, 1.5625053380692426, -6.180392756736983, -0.20407668040293722, -4.6476103643607685, 29.380064536220306, 0.0, 0.7502978877582536, -0.16358232641365608, -19.937609244085568] )
    coef6 = np.array( [-3.423561091777921, -9.228795430171687, 86.95722105482042, 2.772872601988039, -11.968168333741515, -150.61261090270446, -0.24715316589674527, 0.5283723108899993, 4.537286811245538, 82.50581844010263] )
    coef7 = np.array( [-0.19815381951708524, 0.08438758133540872, 1.2650146439483734, -0.007606115512168328, -0.2747023984740461, -0.8844640101378567, 0.0, 0.01870057580949183, 0.0699898278743648, 0.2794360008051127] )
    coef8 = np.array( [2.3932787625531815, -2.020874419612962, -8.938221963838357, 0.576323845480877, 2.8782448498416944, 16.598492450314534, -0.04746016700352029, -0.197101203594028, -1.3860007472886064, -8.289767128060362] )
    coef9 = np.array( [104.40501489600803, -29.942999569370276, -174.42008279158216, 3.708514822202037, 25.14336546356742, 132.35546551746415, -0.16479555172343271, -1.351556690339512, -6.721810844025761, -40.39565289044579] )

    if s == 's1':
        solidity = 0.15
    elif s == 's2':
        solidity = 0.25
    elif s == 's3':
        solidity = 0.5
    elif s == 's4':
        solidity = 0.75
    elif s == 's5':
        solidity = 1.0

    tsr = float(t)/100.


    loc10 = overlay(tsr,solidity,coef0)
    loc20 = overlay(tsr,solidity,coef1)
    loc30 = overlay(tsr,solidity,coef2)
    spr10 = overlay(tsr,solidity,coef3)
    spr20 = overlay(tsr,solidity,coef4)
    skw10 = overlay(tsr,solidity,coef5)
    skw20 = overlay(tsr,solidity,coef6)
    scl10 = overlay(tsr,solidity,coef7)
    scl20 = overlay(tsr,solidity,coef8)
    scl30 = overlay(tsr,solidity,coef9)

    param0 = np.array([loc10,loc20,loc30,spr10,spr20,skw10,skw20,scl10,scl20,scl30])

    param_l = np.array([None,0.,0.,None,None,None,None,0.,0.,0.])
    param_u = np.array([0.,None,None,0.,0.,None,0.,None,None,None])

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


    loc1 = paramf[0]
    loc2 = paramf[1]
    loc3 = paramf[2]
    spr1 = paramf[3]
    spr2 = paramf[4]
    skw1 = paramf[5]
    skw2 = paramf[6]
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
                    skw1_v,spr_v,scl_v,rat_v,spr_v = paramfit(xd[i],skw1,spr,scl,rat,spr)
                    plt.plot(veldist(xfit,skw1_v,spr_v,scl_v,rat_v,spr_v),xfit,'r-',linewidth=2,label=lab2)
                    plt.xlim(0.,1.5)
                    # plt.ylim(-4.,4.)
                    plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
                else:
                    exec('xfit = np.linspace(min(pos'+name+'d)-1.,max(pos'+name+'d)+1.,500)')
                    exec('plt.plot(velo'+name+'d,pos'+name+'d,color)')
                    skw1_v,spr_v,scl_v,rat_v,spr_v = paramfit(xd[i],skw1,spr,scl,rat,spr)
                    plt.plot(veldist(xfit,skw1_v,spr_v,scl_v,rat_v,spr_v),xfit,'r-',linewidth=2)
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
                plt.plot(veldist(xd[i],xfit,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3),xfit,'r-',linewidth=2)
                # plt.xlim(0.,1.5)
                plt.ylim(-1.08,1.08)
                # plt.legend(loc=1)
                plt.xlabel('Normalized Velocity')
                plt.ylabel('$y/D$')

    return loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

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
    t = '600'
    length = 24.

    dia = 6.

    comp = 'mac'
    # comp = 'fsl'
    # comp = 'win'

    loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3 = fit(s,t,length,True,comp,4,True)

    print '\n'
    print 'loc1 =',loc1
    print 'loc2 =',loc2
    print 'loc3 =',loc3
    print 'spr1 =',spr1
    print 'spr2 =',spr2
    print 'skw1 =',skw1
    print 'skw2 =',skw2
    print 'scl1 =',scl1
    print 'scl2 =',scl2
    print 'scl3 =',scl3

plt.show()
