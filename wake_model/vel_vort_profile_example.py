import csv
import numpy as np
from numpy import exp,sqrt
from scipy.special import erf
import matplotlib.pyplot as plt

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


def vortdist(dn,lat,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3):
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



def starccm_read(fdata,dia,rotd,length,opt_print):
    start = length/30.
    lendat =  np.linspace(start,length,30)/dia

    for j in range(np.size(fdata)):
        for k in range(30):
            name = str(k+1)
            exec('pos'+name+' = np.array([])')
            exec('vort'+name+' = np.array([])')
        rot = rotd

        f = open(fdata)

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
                vort1 = np.append(vort1,float(row[1]))
            if row[3] != 'null' and i != 0:
                vort2 = np.append(vort2,float(row[3]))
            if row[5] != 'null' and i != 0:
                vort3 = np.append(vort3,float(row[5]))
            if row[7] != 'null' and i != 0:
                vort4 = np.append(vort4,float(row[7]))
            if row[9] != 'null' and i != 0:
                vort5 = np.append(vort5,float(row[9]))
            if row[11] != 'null' and i != 0:
                vort6 = np.append(vort6,float(row[11]))
            if row[13] != 'null' and i != 0:
                vort7 = np.append(vort7,float(row[13]))
            if row[15] != 'null' and i != 0:
                vort8 = np.append(vort8,float(row[15]))
            if row[17] != 'null' and i != 0:
                vort9 = np.append(vort9,float(row[17]))
            if row[19] != 'null' and i != 0:
                vort10 = np.append(vort10,float(row[19]))
            if row[21] != 'null' and i != 0:
                vort11 = np.append(vort11,float(row[21]))
            if row[23] != 'null' and i != 0:
                vort12 = np.append(vort12,float(row[23]))
            if row[25] != 'null' and i != 0:
                vort13 = np.append(vort13,float(row[25]))
            if row[27] != 'null' and i != 0:
                vort14 = np.append(vort14,float(row[27]))
            if row[29] != 'null' and i != 0:
                vort15 = np.append(vort15,float(row[29]))
            if row[31] != 'null' and i != 0:
                vort16 = np.append(vort16,float(row[31]))
            if row[33] != 'null' and i != 0:
                vort17 = np.append(vort17,float(row[33]))
            if row[35] != 'null' and i != 0:
                vort18 = np.append(vort18,float(row[35]))
            if row[37] != 'null' and i != 0:
                vort19 = np.append(vort19,float(row[37]))
            if row[39] != 'null' and i != 0:
                vort20 = np.append(vort20,float(row[39]))
            if row[41] != 'null' and i != 0:
                vort21 = np.append(vort21,float(row[41]))
            if row[43] != 'null' and i != 0:
                vort22 = np.append(vort22,float(row[43]))
            if row[45] != 'null' and i != 0:
                vort23 = np.append(vort23,float(row[45]))
            if row[47] != 'null' and i != 0:
                vort24 = np.append(vort24,float(row[47]))
            if row[49] != 'null' and i != 0:
                vort25 = np.append(vort25,float(row[49]))
            if row[51] != 'null' and i != 0:
                vort26 = np.append(vort26,float(row[51]))
            if row[53] != 'null' and i != 0:
                vort27 = np.append(vort27,float(row[53]))
            if row[55] != 'null' and i != 0:
                vort28 = np.append(vort28,float(row[55]))
            if row[57] != 'null' and i != 0:
                vort29 = np.append(vort29,float(row[57]))
            if row[59] != 'null' and i != 0:
                vort30 = np.append(vort30,float(row[59]))
            i += 1

        f.close()

        if opt_print == True:
            print 'Imported Data Set',j+1
        # print 'DATA IMPORTED'

        posdn = np.array([])
        poslt = np.array([])
        vortd = np.array([])
        for i in range(30):
            name = str(i+1)
            ind = str(i)

            #Ordering the data numerically by the position
            exec('pos'+name+', vort'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', vort'+name+'))))')
            #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
            exec('pos'+name+'_0 = np.array([])\nvort'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvort'+name+'_0 = np.append(vort'+name+'_0,vort'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvort'+name+' = np.copy(vort'+name+'_0)')
            #Deleting wall boundary data
            # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5. or pos'+name+'[j] < -5.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+' = np.delete(vort'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
            #Deleting values greater than 1*wind
            # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif vort'+name+'[j] > 1. or fabs(pos'+name+'[j]) > 1.2:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+' = np.delete(vort'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')


    return pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9,pos10,pos11,pos12,pos13,pos14,pos15,pos16,pos17,pos18,pos19,pos20,pos21,pos22,pos23,pos24,pos25,pos26,pos27,pos28,pos29,pos30,vort1,vort2,vort3,vort4,vort5,vort6,vort7,vort8,vort9,vort10,vort11,vort12,vort13,vort14,vort15,vort16,vort17,vort18,vort19,vort20,vort21,vort22,vort23,vort24,vort25,vort26,vort27,vort28,vort29,vort30



fdata = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/Vorticity Sections/s2_400.csv'
fdata2 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/Velocity Sections/s2_400.0.csv'

velf = 15.
tsr = 4.
sol = 0.25
dia = 6.
rot = tsr*velf/(dia/2.)
length = 100.

dn = length/30.*5.


pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9,pos10,pos11,pos12,pos13,pos14,pos15,pos16,pos17,pos18,pos19,pos20,pos21,pos22,pos23,pos24,pos25,pos26,pos27,pos28,pos29,pos30,vort1,vort2,vort3,vort4,vort5,vort6,vort7,vort8,vort9,vort10,vort11,vort12,vort13,vort14,vort15,vort16,vort17,vort18,vort19,vort20,vort21,vort22,vort23,vort24,vort25,vort26,vort27,vort28,vort29,vort30 = starccm_read(fdata,dia,rot,length,False)

#TSR: 4.0, Solidity: 0.25; optimized at point
loc1 = -0.00457201560733
loc2 = 0.058092890556
loc3 = 0.641889970022
spr1 = -0.00567160127755
spr2 = -0.0186909979419
skw1 = 0.119497879971
skw2 = -8.142755325
scl1 = 0.103363575966
scl2 = 1.15582186488
scl3 = 14.5968426352

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

loc1 = overlay(tsr,sol,coef0)
loc2 = overlay(tsr,sol,coef1)
loc3 = overlay(tsr,sol,coef2)
spr1 = overlay(tsr,sol,coef3)
spr2 = overlay(tsr,sol,coef4)
skw1 = overlay(tsr,sol,coef5)
skw2 = overlay(tsr,sol,coef6)
scl1 = overlay(tsr,sol,coef7)
scl2 = overlay(tsr,sol,coef8)
scl3 = overlay(tsr,sol,coef9)

pos1d,pos2d,pos3d,pos4d,pos5d,pos6d,pos7d,pos8d,pos9d,pos10d,pos11d,pos12d,pos13d,pos14d,pos15d,pos16d,pos17d,pos18d,pos19d,pos20d,pos21d,pos22d,pos23d,pos24d,pos25d,pos26d,pos27d,pos28d,pos29d,pos30d,velo1d,velo2d,velo3d,velo4d,velo5d,velo6d,velo7d,velo8d,velo9d,velo10d,velo11d,velo12d,velo13d,velo14d,velo15d,velo16d,velo17d,velo18d,velo19d,velo20d,velo21d,velo22d,velo23d,velo24d,velo25d,velo26d,velo27d,velo28d,velo29d,velo30d = starccm_read(fdata2,dia,rot,length,False)

fs = 15

plt.figure(1)
plt.plot(vort5/rot,pos5/dia,'b.',label='CFD Data')
plt.plot(vortdist(dn/dia,pos5/dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3),pos5/dia,'r-',linewidth=2,label='EMG Fitting')
plt.ylim(-1.5,1.5)
plt.legend(loc=1)
plt.xlabel(r'$\gamma/\Omega$',fontsize=fs)
plt.ylabel(r'$y/D$',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)

plt.figure(2)
plt.plot(velo5d/velf,pos5d/dia,'b.')
plt.xlim(0,1.4)
plt.ylim(-5,5)
plt.xlabel(r'$u/U_\infty$',fontsize=fs)
plt.ylabel(r'$y/D$',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)

plt.show()