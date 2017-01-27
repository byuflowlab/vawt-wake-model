from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs,log,sin,arctan,cosh
import matplotlib.pyplot as plt
import database_call as dbc
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'


def veldist(x,men,spr,scl,rat,tns):
    g = 1./(1 + exp(rat*fabs(x)-tns))
    y = (-scl/(spr*sqrt(2.*pi))*exp(-(x+men)**2/(2.*spr**2)))*g + 1.

    return y
    
def paramfit(xd,men,spr,scl,rat,tns):
    men_v = men[0]*xd**2 + men[1]*xd + men[2]
    # men_v = men[0]*xd + men[1]
    # if men_v > 0.5:
    #     men_v = 0.5
    # elif men_v < -0.5:
    #     men_v = -0.5
    
    # spr_v = spr[2]/(spr[1]*sqrt(2.*pi))*exp(-(xd-spr[0])**2/(2.*spr[1]**2))
    # spr_v = spr[0]*xd**2 + spr[1]*xd + spr[2]
    spr_v = spr[2]*spr[1]*spr[0]*exp(spr[1]*xd)*exp(spr[0])*exp(-spr[0]*exp(spr[1]*xd)) + spr[3] #gom
    # if spr_v < 0.35:
    #     spr_v = 0.35
    # elif spr_v > 4.:
    #     spr_v = 4.
    
    # scl_v = scl[2]/(scl[1]*sqrt(2.*pi))*exp(-(xd-scl[0])**2/(2.*scl[1]**2))
    # scl_v = scl[2]/(scl[1]*sqrt(2.*pi))*(exp(-(xd-scl[0])**2/(2.*scl[1]**2))+exp(-(xd+scl[0])**2/(2.*scl[1]**2))) #fold
    scl_v = scl[2]*scl[1]*scl[0]*exp(scl[1]*xd)*exp(scl[0])*exp(-scl[0]*exp(scl[1]*xd)) #gom
    # scl_v = scl[0]*xd**2 + scl[1]*xd + scl[2]
    
    rat_v = rat[0]*xd + rat[1]
    # if rat_v < 0.:
    #     rat_v = 0.
    
    tns_v = tns[0]*xd + tns[1]
    # if tns_v < 0.:
    #     tns_v = 0.
        
    return men_v,spr_v,scl_v,rat_v,tns_v
    
    
def obj_func(xdict):
    global xd
    
    global pos1d
    global pos2d
    global pos3d
    global pos4d
    global pos5d
    global pos6d
    global pos7d
    global pos8d
    global pos9d
    global pos10d
    global pos11d
    global pos12d
    global pos13d
    global pos14d
    global pos15d
    global pos16d
    global pos17d
    global pos18d
    global pos19d
    global pos20d
    global pos21d
    global pos22d
    global pos23d
    global pos24d
    global pos25d
    global pos26d
    global pos27d
    global pos28d
    global pos29d
    global pos30d
    
    global velo1d
    global velo2d
    global velo3d
    global velo4d
    global velo5d
    global velo6d
    global velo7d
    global velo8d
    global velo9d
    global velo10d
    global velo11d
    global velo12d
    global velo13d
    global velo14d
    global velo15d
    global velo16d
    global velo17d
    global velo18d
    global velo19d
    global velo20d
    global velo21d
    global velo22d
    global velo23d
    global velo24d
    global velo25d
    global velo26d
    global velo27d
    global velo28d
    global velo29d
    global velo30d
    
    param = xdict['param']
    funcs = {}
    
    men1 = param[0]
    men2 = param[1]
    men3 = param[2]
    spr1 = param[3]
    spr2 = param[4]
    spr3 = param[5]
    spr4 = param[6]
    scl1 = param[7]
    scl2 = param[8]
    scl3 = param[9]
    rat1 = param[10]
    rat2 = param[11]
    tns1 = param[12]
    tns2 = param[13]
    
    men = np.array([men1,men2,men3])
    spr = np.array([spr1,spr2,spr3,spr4])
    scl = np.array([scl1,scl2,scl3])
    rat = np.array([rat1,rat2])
    tns = np.array([tns1,tns2])
    
    error = 0.
        
    if xd[0] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[0],men,spr,scl,rat,tns)
        for j in range(np.size(pos1d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos1d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos1d[j])-tns_v))) + 1.
            error = error + (vel-velo1d[j])**2
    if xd[1] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[1],men,spr,scl,rat,tns)
        for j in range(np.size(pos2d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos2d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos2d[j])-tns_v))) + 1.
            error = error + (vel-velo2d[j])**2
    if xd[2] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[2],men,spr,scl,rat,tns)
        for j in range(np.size(pos3d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos3d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos3d[j])-tns_v))) + 1.
            error = error + (vel-velo3d[j])**2
    if xd[3] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[3],men,spr,scl,rat,tns)
        for j in range(np.size(pos4d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos4d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos4d[j])-tns_v))) + 1.
            error = error + (vel-velo4d[j])**2
    if xd[4] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[4],men,spr,scl,rat,tns)
        for j in range(np.size(pos5d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos5d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos5d[j])-tns_v))) + 1.
            error = error + (vel-velo5d[j])**2
    if xd[5] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[5],men,spr,scl,rat,tns)
        for j in range(np.size(pos6d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos6d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos6d[j])-tns_v))) + 1.
            error = error + (vel-velo6d[j])**2
    if xd[6] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[6],men,spr,scl,rat,tns)
        for j in range(np.size(pos7d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos7d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos7d[j])-tns_v))) + 1.
            error = error + (vel-velo7d[j])**2
    if xd[7] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[7],men,spr,scl,rat,tns)
        for j in range(np.size(pos8d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos8d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos8d[j])-tns_v))) + 1.
            error = error + (vel-velo8d[j])**2
    if xd[8] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[8],men,spr,scl,rat,tns)
        for j in range(np.size(pos9d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos9d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos9d[j])-tns_v))) + 1.
            error = error + (vel-velo9d[j])**2
    if xd[9] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[9],men,spr,scl,rat,tns)
        for j in range(np.size(pos10d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos10d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos10d[j])-tns_v))) + 1.
            error = error + (vel-velo10d[j])**2
    if xd[10] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[10],men,spr,scl,rat,tns)
        for j in range(np.size(pos11d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos11d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos11d[j])-tns_v))) + 1.
            error = error + (vel-velo11d[j])**2
    if xd[11] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[11],men,spr,scl,rat,tns)
        for j in range(np.size(pos12d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos12d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos12d[j])-tns_v))) + 1.
            error = error + (vel-velo12d[j])**2
    if xd[12] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[12],men,spr,scl,rat,tns)
        for j in range(np.size(pos13d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos13d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos13d[j])-tns_v))) + 1.
            error = error + (vel-velo13d[j])**2
    if xd[13] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[13],men,spr,scl,rat,tns)
        for j in range(np.size(pos14d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos14d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos14d[j])-tns_v))) + 1.
            error = error + (vel-velo14d[j])**2
    if xd[14] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[14],men,spr,scl,rat,tns)
        for j in range(np.size(pos15d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos15d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos15d[j])-tns_v))) + 1.
            error = error + (vel-velo15d[j])**2
    if xd[15] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[15],men,spr,scl,rat,tns)
        for j in range(np.size(pos16d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos16d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos16d[j])-tns_v))) + 1.
            error = error + (vel-velo16d[j])**2
    if xd[16] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[16],men,spr,scl,rat,tns)
        for j in range(np.size(pos17d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos17d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos17d[j])-tns_v))) + 1.
            error = error + (vel-velo17d[j])**2
    if xd[17] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[17],men,spr,scl,rat,tns)
        for j in range(np.size(pos18d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos18d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos18d[j])-tns_v))) + 1.
            error = error + (vel-velo18d[j])**2
    if xd[18] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[18],men,spr,scl,rat,tns)
        for j in range(np.size(pos19d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos19d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos19d[j])-tns_v))) + 1.
            error = error + (vel-velo19d[j])**2
    if xd[19] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[19],men,spr,scl,rat,tns)
        for j in range(np.size(pos20d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos20d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos20d[j])-tns_v))) + 1.
            error = error + (vel-velo20d[j])**2
    if xd[20] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[20],men,spr,scl,rat,tns)
        for j in range(np.size(pos21d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos21d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos21d[j])-tns_v))) + 1.
            error = error + (vel-velo21d[j])**2
    if xd[21] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[21],men,spr,scl,rat,tns)
        for j in range(np.size(pos22d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos22d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos22d[j])-tns_v))) + 1.
            error = error + (vel-velo22d[j])**2
    if xd[22] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[22],men,spr,scl,rat,tns)
        for j in range(np.size(pos23d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos23d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos23d[j])-tns_v))) + 1.
            error = error + (vel-velo23d[j])**2
    if xd[23] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[23],men,spr,scl,rat,tns)
        for j in range(np.size(pos24d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos24d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos24d[j])-tns_v))) + 1.
            error = error + (vel-velo24d[j])**2
    if xd[24] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[24],men,spr,scl,rat,tns)
        for j in range(np.size(pos25d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos25d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos25d[j])-tns_v))) + 1.
            error = error + (vel-velo25d[j])**2
    if xd[25] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[25],men,spr,scl,rat,tns)
        for j in range(np.size(pos26d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos26d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos26d[j])-tns_v))) + 1.
            error = error + (vel-velo26d[j])**2
    if xd[26] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[26],men,spr,scl,rat,tns)
        for j in range(np.size(pos27d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos27d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos27d[j])-tns_v))) + 1.
            error = error + (vel-velo27d[j])**2
    if xd[27] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[27],men,spr,scl,rat,tns)
        for j in range(np.size(pos28d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos28d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos28d[j])-tns_v))) + 1.
            error = error + (vel-velo28d[j])**2
    if xd[28] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[28],men,spr,scl,rat,tns)
        for j in range(np.size(pos29d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos29d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos29d[j])-tns_v))) + 1.
            error = error + (vel-velo29d[j])**2
    if xd[29] > 0.58:
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[29],men,spr,scl,rat,tns)
        for j in range(np.size(pos30d)):
            vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(pos30d[j]+men_v)**2/(2.*spr_v**2)))*(1./(1. + exp(rat_v*fabs(pos30d[j])-tns_v))) + 1.
            error = error + (vel-velo30d[j])**2
    #Print
    print error

    funcs['obj'] = error
    
    
    
    fail = False

    return funcs, fail

def starccm_read(fdata,dia,windd,opt_print):
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
            # #Deleting values greater than 1*wind
            # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif velo'+name+'[j] > 1. or fabs(pos'+name+'[j]) > 1.2:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            exec('pos'+name+'d = np.append(pos'+name+'d,pos'+name+')\nvelo'+name+'d = np.append(velo'+name+'d,velo'+name+')')




    

    


    return pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d


def fit(s,t,length,plot,comp,read_data,opt_print):
    global xd
    
    global pos1d
    global pos2d
    global pos3d
    global pos4d
    global pos5d
    global pos6d
    global pos7d
    global pos8d
    global pos9d
    global pos10d
    global pos11d
    global pos12d
    global pos13d
    global pos14d
    global pos15d
    global pos16d
    global pos17d
    global pos18d
    global pos19d
    global pos20d
    global pos21d
    global pos22d
    global pos23d
    global pos24d
    global pos25d
    global pos26d
    global pos27d
    global pos28d
    global pos29d
    global pos30d
    
    global velo1d
    global velo2d
    global velo3d
    global velo4d
    global velo5d
    global velo6d
    global velo7d
    global velo8d
    global velo9d
    global velo10d
    global velo11d
    global velo12d
    global velo13d
    global velo14d
    global velo15d
    global velo16d
    global velo17d
    global velo18d
    global velo19d
    global velo20d
    global velo21d
    global velo22d
    global velo23d
    global velo24d
    global velo25d
    global velo26d
    global velo27d
    global velo28d
    global velo29d
    global velo30d
    
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
        pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read(np.array([fdata]),dia,np.array([wind]),opt_print)
    if read_data ==2:
        pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read(np.array([fdata,fdata2]),dia,np.array([wind,wind2]),opt_print)
    if read_data ==3:
        pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read(np.array([fdata,fdata2,fdata3]),dia,np.array([wind,wind2,wind3]),opt_print)
    if read_data ==4:
        pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read(np.array([fdata,fdata2,fdata3,fdata4]),dia,np.array([wind,wind2,wind3,wind4]),opt_print)
    if read_data ==5:
        pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5]),dia,np.array([wind,wind2,wind3,wind4,wind5]),opt_print)
    if read_data ==6:
        pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read(np.array([fdata,fdata2,fdata3,fdata4,fdata5,fdata6]),dia,np.array([wind,wind2,wind3,wind4,wind5,wind6]),opt_print)
    
    start = length/30.
    xd = np.linspace(start,length,30)/dia
    
    
    
## Optimization
    optProb = Optimization('VAWTWake_Velo', obj_func)
    optProb.addObj('obj')
    
    param0 = np.array([2.91638655e-04,  -1.70286993e-03 ,  2.38051673e-02 , -7.65610623e-01,6.40509205e-02  , 6.99046413e-01,   7.83484187e-01  , 4.55408268e-01, 1.18716383e-01  , 2.05484572e+01 , -2.67741935e+00  , 4.43022575e+01,-2.10925147e+00  , 3.30400554e+01])
    
    # param0 = np.array([-0.000168794,0.005413905,-0.151372907,0.194333896,-29.52514866,-5.555943879,0.350540531,106.116331,0.000304188,17.4218187,30.06169144,47.3929253,32.97190468,33.4188042])
    # param0 = np.array([-0.430957043,0.202977737,0.20361798,-0.818862214,-1.568705068,0.753568492,0.255557818,0.162385094,3.7137047,19.67995783,31.32756058,77.24259087,13.72492175,48.22628446])
    
    # param0 = np.array([-0.0002347,0.0094285,-0.2262108,-0.0847296,-0.7634520,5.0331194,0.3964091,6.0135974,0.0052853,19.9466507,-2.0322239,47.3381111,-1.1426171,27.5907093])
    param_l = np.array([-1.,-1,-1.,-1.,-1,-1.,0.,0.,0.,None,0.,None,0.]) 
    param_u = np.array([1.,1.,1.,1.,1.,1.,None,None,None,0.,None,0.,None])
    param_l = np.array([None,None,None,None,None,None,None,None,None,None,None,0.,None,0.]) 
    param_u = np.array([None,None,None,None,None,None,None,None,None,None,0.,None,0.,None])
    
    nparam = np.size(param0)
    optProb.addVarGroup('param', nparam, 'c', lower=param_l, upper=param_u, value=param0)
    
    opt = SNOPT()
    opt.setOption('Scale option',2)
    if comp == 'mac':
        print 'mac'
        # opt.setOption('Print file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptVel/SNOPT_print'+s+'_'+t+'.out')
        # opt.setOption('Summary file','/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/OptVel/SNOPT_summary'+s+'_'+t+'.out')
    elif comp == 'fsl':
        opt.setOption('Print file','/fslhome/ebtingey/compute/VAWTWakeModel/OptVel/SNOPT_print'+s+'_'+t+'.out')
        opt.setOption('Summary file','/fslhome/ebtingey/compute/VAWTWakeModel/OptVel/SNOPT_summary'+s+'_'+t+'.out')
    # res = opt(optProb, sens=None)
    # if opt_print == True:
    #     print res
    #
    # pow = res.fStar
    # paramf = res.xStar['param']
    # if opt_print == True:
    #     print paramf[0]
    #     print paramf[1]
    #     print paramf[2]
    #     print paramf[3]
    #     print paramf[4]
    #     print paramf[5]
    #     print paramf[6]
    #     print paramf[7]
    #     print paramf[8]
    #     print paramf[9]
    #     print paramf[10]
    #     print paramf[11]
    #     print paramf[12]
    #     print paramf[13]
    #
    # men = np.array([paramf[0],paramf[1],paramf[2]])
    # spr = np.array([paramf[3],paramf[4],paramf[5],paramf[6]])
    # scl = np.array([paramf[7],paramf[8],paramf[9]])
    # rat = np.array([paramf[10],paramf[11]])
    # tns = np.array([paramf[12],paramf[13]])

    men = np.array( [0.0002809298385469616, -0.002152904786611757, 0.02421458254508763] )
    spr = np.array( [-0.7619069186757337, 0.054654574491074574, 0.693209512634655, 0.7589702331063731] )
    scl = np.array( [0.43439424298020607, 0.1208447537121723, 20.53828795693433] )
    rat = np.array( [-2.4841293302971597, 44.32673189170779] )
    tns = np.array( [-1.9925689536368407, 33.01510013478536] )
    
    
    
    if plot == True:
        # for i in range(30):
        #     name = str(i+1)
        #     plt.figure(1)
        #     plt.subplot(5,6,i+1)
        #     color = 'bo'
        #     exec('xfit = np.linspace(min(pos'+name+'d)-1.,max(pos'+name+'d)+1.,500)')
        #     exec('plt.plot(velo'+name+'d,pos'+name+'d,color)')
        #     men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[i],men,spr,scl,rat,tns)
        #     plt.plot(veldist(xfit,men_v,spr_v,scl_v,rat_v,tns_v),xfit,'r-',linewidth=2)
        #     plt.xlim(0.,1.5)
        #     # plt.ylim(-4.,4.)
        #     # plt.legend(loc=1)
        #     plt.xlabel('Normalized Velocity')
        #     plt.ylabel('$y/D$')

        fs = 25
        fig2 = plt.figure(2,figsize=(9,7))
        fig2.subplots_adjust(left=.14,bottom=.13,right=.95)
        color = 'bo'
        xfit = np.linspace(min(pos9d)-1.,max(pos9d)+1.,500)
        plt.plot(velo9d,pos9d,'b.')
        # plt.subplot(1,2,1)
        # plt.plot(velo5d/velf,pos5d/dia,'b.')
        men_v,spr_v,scl_v,rat_v,tns_v = paramfit(xd[8],men,spr,scl,rat,tns)
        plt.plot(veldist(xfit,men_v,spr_v,scl_v,rat_v,tns_v),xfit,'r-',linewidth=2)
        plt.xlim(0,1.4)
        plt.ylim(-5,5)
        plt.xlabel(r'$u/U_\infty$',fontsize=fs)
        plt.ylabel(r'$y/D$',fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
    
    return men,spr,scl,rat,tns

## Main File
if __name__ == "__main__":
    cv = False

    s = 's2'
    t = '400'
    length = 100.

    dia = 6.
    
    comp = 'mac'
    # comp = 'fsl'

    men,spr,scl,rat,tns = fit(s,t,length,True,'mac',4,True)
    
    print '\n'
    print 'men = np.array(',men.tolist(),')'
    print 'spr = np.array(',spr.tolist(),')'
    print 'scl = np.array(',scl.tolist(),')'
    print 'rat = np.array(',rat.tolist(),')'
    print 'tns = np.array(',tns.tolist(),')'
    
    plt.show()