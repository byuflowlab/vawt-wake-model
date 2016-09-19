import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs,log
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,leastsq,minimize
from scipy.stats import beta, norm
from scipy.special import erf,i0
from scipy.integrate import simps
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'


def pdf(x):
    return 1/sqrt(2*pi) * exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/sqrt(2))) / 2

def skew(x,e,w,a,c):
# def skew(x,b,a,c,d):
    # w = w*sqrt(c)
    # t = (x - e)/w
    # return c*(2./w)*pdf(t)*cdf(a*t)
    
    return c*a/2.0*exp(a/2.0*(2*e+a*w**2-2*x))*(1-erf((e + a*w**2 -x)/(sqrt(2)*w)))
    
    # return c*(a*b*exp(b*x)*exp(a)*exp(-a*exp(b*x)))
    # return c*(b*exp(-b*x)*exp(-a*exp(-b*x))*(1+a*(1-exp(-b*x))))
    # y = np.zeros_like(x)
    # for i in range(np.size(x)):
    #     if x[i] >= 0.:
    #         y[i]= d*((c/b)*(x[i]/b)**(c-1)*exp(-(x[i]/b)**c))
    #         # y[i] = d*((c/b)*(x[i]/b)**(c-1))/(1 + (x[i]/b)**c)**2
    #         # y[i] = d*((sqrt(b/(2*pi*x[i])**3))*exp(-b*(x[i]-c)**2/(2*c**2*x[i])))
    #         # y[i] = d*((x[i]/b**2)*exp(-(x[i]**2+c**2)/(2*b**2)*i0(c*x[i]/b**2)))
    #     else:
    #         y[i] = 0.
    # 
    # return y

def error(params):
    y = skew(x, *params)
    # print params, np.sum((y - vort)**2)
    return np.sum((y - vort)**2)


def residuals(p,x,y,i):
    # if p[3] < 0.:
    #     penalty = penalty + 1000.
    # # if fabs(p[0]) < 0:
    # #     penalty = penalty + 1000.
    # # if fabs(p[3]-fabs(int)) > 0.5:
    # #     penalty = penalty + 1000.
    # if fabs(p[2]) < 0.:
    #     penalty = penalty + 1000.
    return (y - skew(x,p[0],p[1],p[2],i))
    # return np.sum((y - skew(x,p[0],p[1],p[2],p[3]))**2)

def residualfun(x):
    return x[5] - skew(x[4],x[0],x[1],x[2],x[3])
    
def expfit(x,lam,horz,lat):
    # return exp(lam*(x-horz))+lat
    return horz*x+lat
    
def piecewise(x,a,b,c,d,e,horza,scal):
    # global decayu0
    # global decayl0
    #quadratic fit/folded normal distribution
    
    y = np.zeros_like(x)
    
    # for i in range(np.size(x)):
    #     if x[i] > horz:
    #         y[i] = c*(1/(b*sqrt(2*pi)))*(exp(-(scal*(x[i]-horz)-a)**2/(2*b**2))+exp(-(scal*(x[i]-horz)+a)**2/(2*b**2)))
    #     else:
    #         zero = c*(1/(b*sqrt(2*pi)))*(exp(-(-a)**2/(2*b**2))+exp(-(a)**2/(2*b**2)))
    #         y[i] = d*(x[i]-horz)**2+e*(x[i]-horz)+zero
    # return y
    # return horz*x+scal
    
    # for i in range(np.size(x)):
    #     if x[i] < horza:
    #         y[i] = b*x[i] + c
    #     else:
    #         zero = b*horza + c
    #         horz = -(log(zero)/d - horza)
    #         y[i] = exp(d*x[i]-horz)
    #         
    # return y
    
    return a/(1 + exp(b*(x - c)))
    
def piecewiseint(x,a,b,c,d,e,horza,scal):
    #quadratic fit/folded normal distribution
    
    y = np.zeros_like(x)
    
    # for i in range(np.size(x)):
    #     if x[i] > horz:
    #         y[i] = c*(1/(b*sqrt(2*pi)))*(exp(-(scal*(x[i]-horz)-a)**2/(2*b**2))+exp(-(scal*(x[i]-horz)+a)**2/(2*b**2)))
    #     else:
    #         zero = c*(1/(b*sqrt(2*pi)))*(exp(-(-a)**2/(2*b**2))+exp(-(a)**2/(2*b**2)))
    #         y[i] = d*(x[i]-horz)**2+e*(x[i]-horz)+zero
    # return y
    # return horz*x+scal
    
    # for i in range(np.size(x)):
    #     if x[i] < horza:
    #         y[i] = b*x[i] + c
    #     else:
    #         zero = b*horza + c
    #         horz = -(log(zero)/d - horza)
    #         y[i] = exp(d*x[i]-horz)
    #         
    # return y
    
    return a/(1 + exp(b*(x - c)))
    
    


def poly(x,a,b,c,d):
    # return a*x**3+b*x**2+c*x+d
    # return b*x**2+c*x+d
    return c*x+d
    # return a/(1. + exp(b*(x-c))) + 1.
    
def sigmoid(x,max,decay,infl):
    return max/(1 + exp(decay*(x - infl)))
    
def starccm_read(fdata,length,rad,wind,rot,tsr,s):
    integrate = True
    ploti = True
    
    dia = rad*2.
    
    for i in range(30):
        name = str(i+1)
        exec('pos'+name+' = np.array([])')
        exec('vort'+name+' = np.array([])')
        
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
    
    for i in range(30):
        name = str(i+1)
        
        #Ordering the data numerically by the position
        exec('pos'+name+', vort'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', vort'+name+'))))')
        #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
        exec('pos'+name+'_0 = np.array([])\nvort'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvort'+name+'_0 = np.append(vort'+name+'_0,vort'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvort'+name+' = np.copy(vort'+name+'_0)')
    
    
    start = length/30.
    x = np.linspace(start,length,30)
    
    # fst = 15
    # plt.figure(53)
    # plt.subplot(1,4,1)
    # plt.plot(vort2/rot,pos2/dia,'bo-')
    # plt.xlabel('Norm. Vorticity',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.title('$x/D$ = 1.11')
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(-1,1)
    # plt.ylim(-1.2,1.2)
    # 
    # plt.subplot(1,4,2)
    # plt.plot(vort15/rot,pos15/dia,'bo-')
    # plt.xlabel('Norm. Vorticity',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.title('$x/D$ = 8.33',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(-1,1)
    # plt.ylim(-1.2,1.2)
    # 
    # plt.subplot(1,4,3)
    # plt.plot(vort21/rot,pos21/dia,'bo-')
    # plt.xlabel('Norm. Vorticity',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.title('$x/D$ = 11.67',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(-1,1)
    # plt.ylim(-1.2,1.2)
    # 
    # plt.subplot(1,4,4)
    # plt.plot(vort24/rot,pos24/dia,'bo-')
    # plt.xlabel('Norm. Vorticity',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.title('$x/D$ = 13.33',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(-1,1)
    # plt.ylim(-1.2,1.2)
    # 
    # plt.show()
    

## Integrating over the vorticity profiles
    if integrate == True:
        intu = np.array([])
        intl = np.array([])
        
        #Making copies to integrate over
        for i in range(30):
            name = str(i+1)
            
            exec('pos'+name+'intu = np.copy(pos'+name+')')
            exec('pos'+name+'intl = np.copy(pos'+name+')')
            exec('vort'+name+'intu = np.copy(vort'+name+')')
            exec('vort'+name+'intl = np.copy(vort'+name+')')
        
        for i in range(30):
            name = str(i+1)
            ind = str(i)
            #Deleting the wrong half of the vorticity data in each set
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'intu)):\n\tif vort'+name+'intu[j] > 0. or pos'+name+'intu[j] < -1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'intu = np.delete(vort'+name+'intu,delvec[j])\n\tpos'+name+'intu = np.delete(pos'+name+'intu,delvec[j])')
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'intl)):\n\tif vort'+name+'intl[j] < 0. or pos'+name+'intl[j] > 1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'intl = np.delete(vort'+name+'intl,delvec[j])\n\tpos'+name+'intl = np.delete(pos'+name+'intl,delvec[j])')
            #Normalizing the data by diameter and rotation
            exec('pos'+name+'intu = (1./dia)*pos'+name+'intu')
            exec('pos'+name+'intl = (1./dia)*pos'+name+'intl')
            exec('vort'+name+'intu = (1./rot)*vort'+name+'intu')
            exec('vort'+name+'intl = (1./rot)*vort'+name+'intl')
            #Integrating
            exec('intu0 = np.trapz(vort'+name+'intu, x=pos'+name+'intu)')
            exec('intl0 = np.trapz(vort'+name+'intl, x=pos'+name+'intl)')
            intu = np.append(intu,intu0)
            intl = np.append(intl,intl0)
    
        # plt.figure(8)
        # plt.plot(x/dia,intu,'bo')
        # plt.plot(x/dia,intl,'ro')
        # plt.show()
    else:
        intu = 0
        intl = 0
    


## Making Copies for Picking Values
    for i in range(30):
        name = str(i+1)
        
        exec('pos'+name+'gup = np.copy(pos'+name+')')
        exec('pos'+name+'glow = np.copy(pos'+name+')')
        exec('vort'+name+'gup = np.copy(vort'+name+')')
        exec('vort'+name+'glow = np.copy(vort'+name+')')
    
## Chosing vortex sheet locations and strengths
    region = 0.3
    fix = 0.55 # Taking into account short wakes
    
    # Finding maximum vortex sheet
    for i in range(30):
        ind = str(i)
        name = str(i+1)
        exec('if x['+ind+'] < rad+fix:\n\twhile True:\n\t\tmin'+name+' = np.argmin(vort'+name+')\n\t\tmax'+name+' = np.argmax(vort'+name+')\n\t\tif pos'+name+'[min'+name+'] < rad+region:\n\t\t\tvort'+name+' = np.delete(vort'+name+',min'+name+')\n\t\t\tpos'+name+' = np.delete(pos'+name+',min'+name+')\n\t\telse:\n\t\t\tif pos'+name+'[max'+name+'] > -rad-region:\n\t\t\t\tvort'+name+' = np.delete(vort'+name+',max'+name+')\n\t\t\t\tpos'+name+' = np.delete(pos'+name+',max'+name+')\n\t\t\telse:\n\t\t\t\tbreak\nelse:\n\tmin'+name+' = np.argmin(vort'+name+')\n\tmax'+name+' = np.argmax(vort'+name+')')
        # print i, 'Max'
    
    # Naming maximum values
    for i in range(30):
        name = str(i+1)
        exec('vort'+name+'u = vort'+name+'[min'+name+']')
        exec('vort'+name+'l = vort'+name+'[max'+name+']')
        exec('pos'+name+'u = pos'+name+'[min'+name+']')
        exec('pos'+name+'l = pos'+name+'[max'+name+']')
        
    
## Compiling max values    
    gamma_up = np.array([vort1[min1],vort2[min2],vort3[min3],vort4[min4],vort5[min5],vort6[min6],vort7[min7],vort8[min8],vort9[min9],vort10[min10],vort11[min11],vort12[min12],vort13[min13],vort14[min14],vort15[min15],vort16[min16],vort17[min17],vort18[min18],vort19[min19],vort20[min20],vort21[min21],vort22[min22],vort23[min23],vort24[min24],vort25[min25],vort26[min26],vort27[min27],vort28[min28],vort29[min29],vort30[min30]])
    gamma_low = np.array([vort1[max1],vort2[max2],vort3[max3],vort4[max4],vort5[max5],vort6[max6],vort7[max7],vort8[max8],vort9[max9],vort10[max10],vort11[max11],vort12[max12],vort13[max13],vort14[max14],vort15[max15],vort16[max16],vort17[max17],vort18[max18],vort19[max19],vort20[max20],vort21[max21],vort22[max22],vort23[max23],vort24[max24],vort25[max25],vort26[max26],vort27[max27],vort28[max28],vort29[max29],vort30[max30]])
    lat_up = np.array([pos1[min1],pos2[min2],pos3[min3],pos4[min4],pos5[min5],pos6[min6],pos7[min7],pos8[min8],pos9[min9],pos10[min10],pos11[min11],pos12[min12],pos13[min13],pos14[min14],pos15[min15],pos16[min16],pos17[min17],pos18[min18],pos19[min19],pos20[min20],pos21[min21],pos22[min22],pos23[min23],pos24[min24],pos25[min25],pos26[min26],pos27[min27],pos28[min28],pos29[min29],pos30[min30]])
    lat_low = np.array([pos1[max1],pos2[max2],pos3[max3],pos4[max4],pos5[max5],pos6[max6],pos7[max7],pos8[max8],pos9[max9],pos10[max10],pos11[max11],pos12[max12],pos13[max13],pos14[max14],pos15[max15],pos16[max16],pos17[max17],pos18[max18],pos19[max19],pos20[max20],pos21[max21],pos22[max22],pos23[max23],pos24[max24],pos25[max25],pos26[max26],pos27[max27],pos28[max28],pos29[max29],pos30[max30]])

## SKEW DISTRIBUTION FITTING
## Creating arrays for fitting
    regionfix = (rad)/dia+.01
    fix = 1. # Taking into account short wakes (.55)
    vortmax = 0.1
    sreg = -0.3 #eliminating double vortex on s1 profiles
    svortmax = 0.03
    
    eu = np.zeros(30)
    wu = np.zeros(30)
    au = np.zeros(30)
    cu = np.zeros(30)
    el = np.zeros(30)
    wl = np.zeros(30)
    al = np.zeros(30)
    cl = np.zeros(30)
    
    latmax1 = max(np.fabs(lat_up))
    latmax2 = max(np.fabs(lat_low))
    latmax = max(latmax1,latmax2)/dia
    
    for i in range(30):
        name = str(i+1)
        ind = str(i)
        exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'gup)):\n\tif vort'+name+'gup[j] > 0. or pos'+name+'gup[j] > 6.5 or pos'+name+'gup[j] < -1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'gup = np.delete(vort'+name+'gup,delvec[j])\n\tpos'+name+'gup = np.delete(pos'+name+'gup,delvec[j])')
        exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'glow)):\n\tif vort'+name+'glow[j] < 0. or pos'+name+'glow[j] < -6.5 or pos'+name+'glow[j] > 1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'glow = np.delete(vort'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
        # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'gup)):\n\tif vort'+name+'gup[j] > 0. or pos'+name+'gup[j] < -1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'gup = np.delete(vort'+name+'gup,delvec[j])\n\tpos'+name+'gup = np.delete(pos'+name+'gup,delvec[j])')
        # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'glow)):\n\tif vort'+name+'glow[j] < 0. or pos'+name+'glow[j] > 1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'glow = np.delete(vort'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
        exec('pos'+name+'gup = (-1./dia)*pos'+name+'gup')
        exec('pos'+name+'glow = (1./dia)*pos'+name+'glow')
        exec('vort'+name+'gup = (-1./rot)*vort'+name+'gup')
        exec('vort'+name+'glow = (1./rot)*vort'+name+'glow')
        exec('delvec = np.array([])\nif x['+ind+'] < rad+fix:\n\tfor j in range(np.size(vort'+name+'gup)):\n\t\tif pos'+name+'gup[j] > -regionfix and vort'+name+'gup[j] > vortmax:\n\t\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'gup = np.delete(vort'+name+'gup,delvec[j])\n\tpos'+name+'gup = np.delete(pos'+name+'gup,delvec[j])')
        exec('delvec = np.array([])\nif x['+ind+'] < rad+fix:\n\tfor j in range(np.size(vort'+name+'glow)):\n\t\tif pos'+name+'glow[j] > -regionfix and vort'+name+'glow[j] > vortmax:\n\t\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'glow = np.delete(vort'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
        # exec('delvec = np.array([])\nif x['+ind+'] < rad+fix:\n\tfor j in range(np.size(vort'+name+'glow)):\n\t\tif pos'+name+'glow[j] < regionfix and vort'+name+'glow[j] > vortmax:\n\t\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'glow = np.delete(vort'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
        if s == 's1' and tsr <= 3.5:
            exec('delvec = np.array([])\nfor j in range(np.size(vort'+name+'glow)):\n\tif pos'+name+'glow[j] > sreg and vort'+name+'glow[j] > svortmax:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'glow = np.delete(vort'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
        if s == 's2' and tsr <= 3.:
            exec('delvec = np.array([])\nfor j in range(np.size(vort'+name+'glow)):\n\tif pos'+name+'glow[j] > sreg and vort'+name+'glow[j] > svortmax:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvort'+name+'glow = np.delete(vort'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')

## Skewed Distribution Fitting
        diap = 1.
        rotp = 1.
        gam = np.fabs(gamma_up[i]/rot)
        mpos = lat_up[i]
        # e = 0.6
        # w = 0.5
        # a = 5.
        # c = 0.1
        
        # e = 0.8
        exec('e = -pos'+name+'u/dia')
        w = 0.05
        a = 15.
        c = -0.2
        
        x0 = np.array([e,w,a,c])
        bounds = ((0.0,1.0), (-0.2, 0.0), (-60.0, 0.0), (-0.2, 0.0))
        # meth = 'L-BFGS-B'
        meth = 'SLSQP'
        # meth = 'Nelder-Mead'
        dis = 'disp'
        tol = 'xtol'
        
        ft = 1e-8
        xt = 1e-8
        mf = int(1e6)
        
        # exec('betaparam,_ = leastsq(residuals,np.array([e,w,a]),args=(pos'+name+'gup,vort'+name+'gup,fabs(intu)),ftol=ft, xtol=xt,maxfev=mf)\neu['+ind+'] = betaparam[0]\nwu['+ind+'] = betaparam[1]\nau['+ind+'] = betaparam[2]\ncu['+ind+'] = betaparam[3]')
        exec('betaparam,_ = leastsq(residuals,np.array([e,w,a]),args=(pos'+name+'gup,vort'+name+'gup,-intu['+ind+']),ftol=ft, xtol=xt,maxfev=mf)\neu['+ind+'] = -betaparam[0]\nwu['+ind+'] = -betaparam[1]\nau['+ind+'] = -betaparam[2]')
        # exec('betaparam,_ = curve_fit(skew,pos'+name+'gup,vort'+name+'gup,p0=(e,w,a,c))\neu['+ind+'] = betaparam[0]\nwu['+ind+'] = betaparam[1]\nau['+ind+'] = betaparam[2]\ncu['+ind+'] = betaparam[3]')
        # exec('betaparam,_ = curve_fit(skew,pos'+name+'gup,vort'+name+'gup,p0=(e,w,a,c))\neu['+ind+'] = betaparam[0]\nwu['+ind+'] = betaparam[1]\nau['+ind+'] = betaparam[2]\ncu['+ind+'] = betaparam[3]')
        # exec('res = minimize(residuals,x0,args=(pos'+name+'gup,vort'+name+'gup),method=meth,bounds=bounds,options={tol: 1e-8, dis: False})\neu['+ind+'] = res.x[0]\nwu['+ind+'] = res.x[1]\nau['+ind+'] = res.x[2]\ncu['+ind+'] = res.x[3]')
        # exec('res = minimize(residuals,x0,args=(pos'+name+'gup,vort'+name+'gup),method=meth,options={tol: 1e-8, dis: False})\neu['+ind+'] = res.x[0]\nwu['+ind+'] = res.x[1]\nau['+ind+'] = res.x[2]\ncu['+ind+'] = res.x[3]')

        # exec('print eu['+ind+'],wu['+ind+'],au['+ind+'],-intu['+ind+']')
        ceu = np.array( [-0.051533741732590066, 0.603645756517617, 3.404058219228679] )
        cwu = np.array( [0.0012367836421159789, -0.016292716663412185, 0.08228677443522553, 0.6891188133055813] )
        cau = np.array( [0.006775814669580242, -0.09990512632914271, 0.6437647406588762, -5.547388443086547] )
        ccu = np.array( [-0.001095215841509046, 0.025087453250091336, -0.17923773313904467, 0.6219114941017065, 9.724403752075677] )
        # Plotting
        if ploti == True:
            plt.figure(1)
            ax1 = plt.subplot(5,6,i+1)
            color = 'bo'
            color2 = 'r-'
            color3 = 'c.'
            color4 = 'g.'
            fs = 15
            lab = 'CFD'
            lab2 = 'Trend'
            tex = '$x/D$ = '+str("{0:.2f}".format(x[i]/dia))
            exec('xfit = np.linspace(min(pos'+name+'gup),max(pos'+name+'gup),500)')
            if i == 5 and wind == 15.:
                exec('plt.plot(-pos'+name+'gup/diap,vort'+name+'gup/rotp,color,label=lab)')
                exec('plt.plot(-xfit,skew(xfit,-eu['+ind+'],-wu['+ind+'],-au['+ind+'],-intu['+ind+'])/rotp,color2,linewidth=2,label= lab2)')
                plt.xlim(0.,1.2)
                plt.ylim(0.,1.)
                plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
            else:
                exec('plt.plot(-pos'+name+'gup/diap,vort'+name+'gup/rotp,color)')
                exec('plt.plot(-xfit,skew(xfit,-eu['+ind+'],-wu['+ind+'],-au['+ind+'],-intu['+ind+'])/rotp,color2,linewidth=2)')
                plt.xlim(0.,1.2)
                plt.ylim(0.,1.)
            if wind == 15.:
                plt.text(0.3,0.8,tex,fontsize=fs)
            if i <= 23:
                plt.setp(ax1.get_xticklabels(), visible=False)
            else:
                plt.xlabel('$y/D$',fontsize=fs)
                plt.xticks(fontsize=fs)
            if i == 0 or i == 6 or i == 12 or i == 18 or i ==24:
                plt.ylabel('Norm. Vorticity',fontsize=fs)
                plt.yticks(fontsize=fs)
            else:
                plt.setp(ax1.get_yticklabels(), visible=False)
            
            # if i==1:
            #     plt.figure(52)
            #     if wind == 15.:
            #         exec('plt.plot(-pos'+name+'gup/diap,vort'+name+'gup/rotp,color,label=lab)')
            #         exec('plt.plot(-xfit,skew(xfit,-eu['+ind+'],-wu['+ind+'],-au['+ind+'],-intu['+ind+'])/rotp,color2,linewidth=2,label= lab2)')
            #     
            #     else:
            #         exec('plt.plot(-pos'+name+'gup/diap,vort'+name+'gup/rotp,color)')
            #         exec('plt.plot(-xfit,skew(xfit,-eu['+ind+'],-wu['+ind+'],-au['+ind+'],-intu['+ind+'])/rotp,color2,linewidth=2)')
            #         
            #         plt.xlim(0.,1.2)
            #         plt.ylim(0.,1.)
            #         plt.text(0.1,0.9,tex,fontsize=fs+2)
            #         plt.xlabel('$y/D$',fontsize=fs)
            #         plt.xticks(fontsize=fs)
            #         plt.ylabel('Norm. Vorticity',fontsize=fs)
            #         plt.yticks(fontsize=fs)
            #         plt.legend(loc=1,fontsize=fs)
                
            
            
            
            
            
            
            
            
        
        gam = np.fabs(gamma_low[i]/rotp)
        mpos = lat_low[i]
        # e = -0.6
        # w = -0.5
        # a = -5.
        # c = 0.1
        
        # e = -0.8
        exec('e = pos'+name+'l/dia')
        w = 0.05
        a = 15.
        c = 0.2
        
        x0 = np.array([e,w,a,c])
        bounds = ((-1.0,0.0), (0.0, 0.2), (0.0, 60.0), (0.0, 0.2))
        
        
        # exec('betaparam,_ = leastsq(residuals,np.array([e,w,a]),args=(pos'+name+'glow,vort'+name+'glow,fabs(intl)),ftol=ft, xtol=xt,maxfev=mf)\nel['+ind+'] = betaparam[0]\nwl['+ind+'] = betaparam[1]\nal['+ind+'] = betaparam[2]\ncl['+ind+'] = betaparam[3]')
        exec('betaparam,_ = leastsq(residuals,np.array([e,w,a]),args=(pos'+name+'glow,vort'+name+'glow,intl['+ind+']),ftol=ft, xtol=xt,maxfev=mf)\nel['+ind+'] = betaparam[0]\nwl['+ind+'] = betaparam[1]\nal['+ind+'] = betaparam[2]')
        # exec('betaparam,_ = curve_fit(skew,pos'+name+'glow,vort'+name+'glow,p0=(e,w,a,c))\nel['+ind+'] = betaparam[0]\nwl['+ind+'] = betaparam[1]\nal['+ind+'] = betaparam[2]\ncl['+ind+'] = betaparam[3]')
        # exec('betaparam,_ = curve_fit(skew,pos'+name+'glow,vort'+name+'glow)\nel['+ind+'] = betaparam[0]\nwl['+ind+'] = betaparam[1]\nal['+ind+'] = betaparam[2]\ncl['+ind+'] = betaparam[3]')
        # exec('betaparam,_ = curve_fit(skew,pos'+name+'glow,vort'+name+'glow,p0=(e,w,a,c))\nel['+ind+'] = betaparam[0]\nwl['+ind+'] = betaparam[1]\nal['+ind+'] = betaparam[2]\ncl['+ind+'] = betaparam[3]')
        # exec('res = minimize(residuals,x0,args=(pos'+name+'glow,vort'+name+'glow),method=meth,bounds=bounds,options={tol: 1e-8, dis: False})\nel['+ind+'] = res.x[0]\nwl['+ind+'] = res.x[1]\nal['+ind+'] = res.x[2]\ncl['+ind+'] = res.x[3]')
        # exec('res = minimize(residuals,x0,args=(pos'+name+'glow,vort'+name+'glow),method=meth,options={tol: 1e-8, dis: False})\nel['+ind+'] = res.x[0]\nwl['+ind+'] = res.x[1]\nal['+ind+'] = res.x[2]\ncl['+ind+'] = res.x[3]')

        # exec('print el['+ind+'],wl['+ind+'],al['+ind+'],intl['+ind+']')
        ced = np.array( [0.038792890464154184, -0.4886209519396058, -3.635922359449307] )
        cwd = np.array( [-0.00036771375242720776, 0.013902999931376486, -0.07867012563750761, 0.9456061873427617] )
        cad = np.array( [-0.006352533336562251, 0.13097620069131152, -1.0268938055673231, 6.401473482471596] )
        ccd = np.array( [-0.0010705538708980264, 0.02408395748957303, -0.16376834249765918, 0.4692981648345209, 10.545454098450099] )
        # Plotting
        if ploti == True:
            plt.figure(2)
            plt.subplot(5,6,i+1)
            color = 'bo'
            color2 = 'r-'
            color3 = 'c.'
            color4 = 'g.'
            lab = str(x[i]/dia)+'D'
            exec('xfit = np.linspace(min(pos'+name+'glow),max(pos'+name+'glow),500)')
            exec('plt.plot(pos'+name+'glow/diap,vort'+name+'glow/rotp,color,label=lab)')
            exec('plt.plot(xfit,skew(xfit,el['+ind+'],wl['+ind+'],al['+ind+'],intl['+ind+'])/rotp,color2,linewidth=2)')
            plt.ylim(0.,1.)
            # plt.legend(loc=1)
            
            
            
            
            
            plt.xlabel('Lateral Distance')
            plt.ylabel('Vorticity')

            
        print i+1,'EMG  Wind =',wind,'Rot =',rot,s,'TSR =',tsr
    # plt.show()
    
    return x,eu,wu,au,el,wl,al,intu,intl


def fit(s,t,length,comp):
    read_data = 4
    
    plott = False
    plotac = False
    plottrend = True
    
    t2 = t+'.0t'
    
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
        fdata2 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel14/Vorticity/'+wfit2+'.csv'
        fdata3 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel12/Vorticity/'+wfit3+'.csv'
        fdata4 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel16/Vorticity/'+wfit4+'.csv'
        fdata5 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot17/Vorticity/'+wfit5+'.csv'
        fdata6 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot18/Vorticity/'+wfit6+'.csv'
    elif comp == 'fsl':
        fdata = '/fslhome/ebtingey/compute/moveForward/Vorticity/'+wfit+'.csv'
        fdata2 = '/fslhome/ebtingey/compute/moveForward/vel14/Vorticity/'+wfit2+'t.csv'
        fdata3 = '/fslhome/ebtingey/compute/moveForward/vel12/Vorticity/'+wfit3+'t.csv'
        fdata4 = '/fslhome/ebtingey/compute/moveForward/vel16/Vorticity/'+wfit4+'t.csv'
        fdata5 = '/fslhome/ebtingey/compute/moveForward/rot17/Vorticity/'+wfit5+'.csv'
        fdata6 = '/fslhome/ebtingey/compute/moveForward/rot18/Vorticity/'+wfit6+'.csv'

    
    
    if read_data >=1:
        xd,eud,wud,aud,eld,wld,ald,intud,intld = starccm_read(fdata,length,rad,wind,rot,tsr,s)
    if read_data >=2:
        xd2,eud2,wud2,aud2,eld2,wld2,ald2,intud2,intld2 = starccm_read(fdata2,length2,rad,wind2,rot2,tsr,s)
    if read_data >=3:
        xd3,eud3,wud3,aud3,eld3,wld3,ald3,intud3,intld3 = starccm_read(fdata3,length3,rad,wind3,rot3,tsr,s)
    if read_data >=4:
        xd4,eud4,wud4,aud4,eld4,wld4,ald4,intud4,intld4 = starccm_read(fdata4,length4,rad,wind4,rot4,tsr,s)
    if read_data >=5:
        xd5,eud5,wud5,aud5,eld5,wld5,ald5,intud5,intld5 = starccm_read(fdata5,length5,rad,wind5,rot5,tsr,s)
    if read_data >=6:
        xd6,eud6,wud6,aud6,eld6,wld6,ald6,intud6,intld6 = starccm_read(fdata6,length6,rad,wind6,rot6,tsr,s)
    
    x = np.zeros(30*read_data)
    eu = np.zeros(30*read_data)
    el = np.zeros(30*read_data)
    wu = np.zeros(30*read_data)
    wl = np.zeros(30*read_data)
    au = np.zeros(30*read_data)
    al = np.zeros(30*read_data)
    # cu = np.zeros(30*read_data)
    # cl = np.zeros(30*read_data)
    intu = np.zeros(30*read_data)
    intl = np.zeros(30*read_data)
    
    
    for i in range(30):
        if read_data >=1:
            x[read_data*i] = xd[i]
            eu[read_data*i] = eud[i]
            el[read_data*i] = eld[i]
            wu[read_data*i] = wud[i]
            wl[read_data*i] = wld[i]
            au[read_data*i] = aud[i]
            al[read_data*i] = ald[i]
            # cu[read_data*i] = cud[i]
            # cl[read_data*i] = cld[i]
            intu[read_data*i] = intud[i]
            intl[read_data*i] = intld[i]
        if read_data >=2:
            x[read_data*i+1] = xd2[i]
            eu[read_data*i+1] = eud2[i]
            el[read_data*i+1] = eld2[i]
            wu[read_data*i+1] = wud2[i]
            wl[read_data*i+1] = wld2[i]
            au[read_data*i+1] = aud2[i]
            al[read_data*i+1] = ald2[i]
            # cu[read_data*i+1] = cud2[i]
            # cl[read_data*i+1] = cld2[i]
            intu[read_data*i+1] = intud2[i]
            intl[read_data*i+1] = intld2[i]
        if read_data >=3:
            x[read_data*i+2] = xd3[i]
            eu[read_data*i+2] = eud3[i]
            el[read_data*i+2] = eld3[i]
            wu[read_data*i+2] = wud3[i]
            wl[read_data*i+2] = wld3[i]
            au[read_data*i+2] = aud3[i]
            al[read_data*i+2] = ald3[i]
            # cu[read_data*i+2] = cud3[i]
            # cl[read_data*i+2] = cld3[i]
            intu[read_data*i+2] = intud3[i]
            intl[read_data*i+2] = intld3[i]
        if read_data >=4:
            x[read_data*i+3] = xd4[i]
            eu[read_data*i+3] = eud4[i]
            el[read_data*i+3] = eld4[i]
            wu[read_data*i+3] = wud4[i]
            wl[read_data*i+3] = wld4[i]
            au[read_data*i+3] = aud4[i]
            al[read_data*i+3] = ald4[i]
            # cu[read_data*i+3] = cud4[i]
            # cl[read_data*i+3] = cld4[i]
            intu[read_data*i+3] = intud4[i]
            intl[read_data*i+3] = intld4[i]
        if read_data >=5:
            x[read_data*i+4] = xd5[i]
            eu[read_data*i+4] = eud5[i]
            el[read_data*i+4] = eld5[i]
            wu[read_data*i+4] = wud5[i]
            wl[read_data*i+4] = wld5[i]
            au[read_data*i+4] = aud5[i]
            al[read_data*i+4] = ald5[i]
            # cu[read_data*i+4] = cud5[i]
            # cl[read_data*i+4] = cld5[i]
            intu[read_data*i+4] = intud5[i]
            intl[read_data*i+4] = intld5[i]
        if read_data >=6:
            x[read_data*i+5] = xd6[i]
            eu[read_data*i+5] = eud6[i]
            el[read_data*i+5] = eld6[i]
            wu[read_data*i+5] = wud6[i]
            wl[read_data*i+5] = wld6[i]
            au[read_data*i+5] = aud6[i]
            al[read_data*i+5] = ald6[i]
            # cu[read_data*i+5] = cud6[i]
            # cl[read_data*i+5] = cld6[i]
            intu[read_data*i+5] = intud6[i]
            intl[read_data*i+5] = intld6[i]
    
    x = x/dia
 
## Deleting Values (outlier)    
    xeu = np.copy(x)
    xel = np.copy(x)
    xwu = np.copy(x)
    xwl = np.copy(x)
    xau = np.copy(x)
    xal = np.copy(x)
    # xcu = np.copy(x)
    # xcl = np.copy(x)
    xintu = np.copy(x)
    xintl = np.copy(x)
    
    xeup = np.copy(xeu)
    xelp = np.copy(xel)
    xwup = np.copy(xwu)
    xwlp = np.copy(xwl)
    xaup = np.copy(xau)
    xalp = np.copy(xal)
    eup = np.copy(eu)
    elp = np.copy(el)
    wup = np.copy(wu)
    wlp = np.copy(wl)
    aup = np.copy(au)
    alp = np.copy(al)
    
    deleu = np.array([])
    delel = np.array([])
    delwu = np.array([])
    delwl = np.array([])
    delau = np.array([])
    delal = np.array([])
    delcu = np.array([])
    delcl = np.array([])
    delintu = np.array([])
    delintl = np.array([])
    for i in range(np.size(x)):
        if fabs(eu[i]) > 1.:
            deleu = np.append(deleu,i)
        if fabs(el[i]) > 1.:
            delel = np.append(delel,i)
        if fabs(wu[i]) > 1.:
            delwu = np.append(delwu,i)
        if fabs(wl[i]) > 1.:
            delwl = np.append(delwl,i)
        if fabs(au[i]) > 20.:
            delau = np.append(delau,i)
        if fabs(al[i]) > 20.:
            delal = np.append(delal,i)
        # if fabs(au[i]) < 10. and xau[i] < 1.:
        #     delau = np.append(delau,i)
        # if fabs(cu[i]) > 1.:
        #     delcu = np.append(delcu,i)
        # if fabs(cl[i]) > 1.:
        #     delcl = np.append(delcl,i)
        if fabs(intu[i]) > 0.3:
            delintu = np.append(delintu,i)
        if fabs(intl[i]) > 0.3:
            delintl = np.append(delintl,i)
    deleu = deleu[::-1]
    delel = delel[::-1]
    delwu = delwu[::-1]
    delwl = delwl[::-1]
    delau = delau[::-1]
    delal = delal[::-1]
    # delcu = delcu[::-1]
    # delcl = delcl[::-1]
    delintu = delintu[::-1]
    delintl = delintl[::-1]
    for i in range(np.size(deleu)):
        eu = np.delete(eu,deleu[i])
        xeu = np.delete(xeu,deleu[i])
    for i in range(np.size(delel)):
        el = np.delete(el,delel[i])
        xel = np.delete(xel,delel[i])
    for i in range(np.size(delwu)):
        wu = np.delete(wu,delwu[i])
        xwu = np.delete(xwu,delwu[i])
    for i in range(np.size(delwl)):
        wl = np.delete(wl,delwl[i])
        xwl = np.delete(xwl,delwl[i])
    for i in range(np.size(delau)):
        au = np.delete(au,delau[i])
        xau = np.delete(xau,delau[i])
    for i in range(np.size(delal)):
        al = np.delete(al,delal[i])
        xal = np.delete(xal,delal[i])
    # for i in range(np.size(delcu)):
    #     cu = np.delete(cu,delcu[i])
    #     xcu = np.delete(xcu,delcu[i])
    # for i in range(np.size(delcl)):
    #     cl = np.delete(cl,delcl[i])
    #     xcl = np.delete(xcl,delcl[i])
    for i in range(np.size(delintu)):
        intu = np.delete(intu,delintu[i])
        xintu = np.delete(xintu,delintu[i])
    for i in range(np.size(delintl)):
        intl = np.delete(intl,delintl[i])
        xintl = np.delete(xintl,delintl[i])
    
    
    if plott == True:
        diap = 1.
        plt.figure(4)
        plt.subplot(2,2,1)
        plt.plot(xeu/diap,eu/diap,'bo')
        plt.plot(xel/diap,el/diap,'ro')
        plt.title(r'Location ($\xi$)')
        # plt.ylim(-1,1)
        
        plt.subplot(2,2,2)
        plt.plot(xwu/diap,wu,'bo')
        plt.plot(xwl/diap,wl,'ro')
        plt.title(r'Scale ($\omega$)')
        # plt.ylim(-0.05,0.15)
        
        plt.subplot(2,2,3)
        plt.plot(xau/diap,au,'bo')
        plt.plot(xal/diap,al,'ro')
        plt.title(r'Skew ($\alpha$)')
        # plt.ylim(-5,5)
        
        # plt.subplot(2,2,4)
        # plt.plot(xcu/diap,cu,'bo')
        # plt.plot(xcl/diap,cl,'ro')
        # plt.title('Resize ($c$)')
        # plt.show()
    # print xal,al
## Deletinge Values (turbulent)
    devrange = 3.0 #0.4
    devrangew = 0.05
    devrangea = 5.
    averange = 0.07
    skewtest = 2.8
    avechange = 0.04
    
    #Spread Variation
    xtest = xwu[0]
    wutest = np.array([])
    kdel = np.array([])
    k = 0
    for i in range(np.size(xwu)):
        xtest0 = xwu[i]
        wut0 = wu[i]
        if xtest == xtest0:
            wutest = np.append(wutest,wut0)
            kdel = np.append(kdel,i)
        else:
            xtest = xtest0
            wut = wut0
            testave = fabs(np.average(wutest))
            # testrange = fabs(max(wutest)-min(wutest))
            # if testrange > devrangew and xwu[i] > 0.5*(length/6):
            #     wudel0 = xwu[kdel[0]]
            #     break
            # else:
            #     wutest = np.array([])
            #     kdel = np.array([])
            #     wutest = np.append(wutest,wut0)
            #     kdel = np.append(kdel,i)
            #     wudel0 = 1e8
                
            if k == 0:
                testave0 = testave
                k += 1
            else:
                if testave - testave0 > avechange and xwu[i] > 0.5*(length/6):
                    wudel0 = xwu[kdel[0]]
                    break
                else:
                    wutest = np.array([])
                    kdel = np.array([])
                    wutest = np.append(wutest,wut0)
                    kdel = np.append(kdel,i)
                    wudel0 = 1e8
                    k += 1
                
    xtest = xwl[0]
    wltest = np.array([])
    kdel = np.array([])
    for i in range(np.size(xwl)):
        xtest0 = xwl[i]
        wlt0 = wl[i]
        if xtest == xtest0:
            wltest = np.append(wltest,wlt0)
            kdel = np.append(kdel,i)
        else:
            xtest = xtest0
            wlt = wlt0
            testave = fabs(np.average(wltest))
            # testrange = fabs(max(wltest)-min(wltest))
            # if testrange > devrangew and xwl[i] > 0.5*(length/6):
            #     wldel0 = xwl[kdel[0]]
            #     break
            # else:
            #     wltest = np.array([])
            #     kdel = np.array([])
            #     wltest = np.append(wltest,wlt0)
            #     kdel = np.append(kdel,i)
            #     wldel0 = 1e8
            
            if k == 0:
                testave0 = testave
                k += 1
            else:
                if testave - testave0 > avechange and xwl[i] > 0.5*(length/6):
                    wldel0 = xwu[kdel[0]]
                    break
                else:
                    wltest = np.array([])
                    kdel = np.array([])
                    wltest = np.append(wltest,wlt0)
                    kdel = np.append(kdel,i)
                    wldel0 = 1e8
                    k += 1
    
    #Skew Variation
    xtest = xau[0]
    autest = np.array([])
    kdel = np.array([])
    for i in range(np.size(xau)):
        xtest0 = xau[i]
        aut0 = au[i]
        if xtest == xtest0:
            autest = np.append(autest,aut0)
            kdel = np.append(kdel,i)
        else:
            xtest = xtest0
            aut = aut0
            testrange = fabs(max(autest)-min(autest))
            if testrange > devrangea and xau[i] > 0.5*(length/6):
                audel0 = xau[kdel[0]]
                break
            else:
                autest = np.array([])
                kdel = np.array([])
                autest = np.append(autest,aut0)
                kdel = np.append(kdel,i)
                audel0 = 1e8
        
    xtest = xal[0]
    altest = np.array([])
    kdel = np.array([])
    for i in range(np.size(xal)):
        xtest0 = xal[i]
        alt0 = al[i]
        if xtest == xtest0:
            altest = np.append(altest,alt0)
            kdel = np.append(kdel,i)
        else:
            xtest = xtest0
            alt = alt0
            testrange = fabs(max(altest)-min(altest))
            if testrange > devrangea and xal[i] > 0.5*(length/6):
                aldel0 = xal[kdel[0]]
                break
            else:
                altest = np.array([])
                kdel = np.array([])
                altest = np.append(altest,alt0)
                kdel = np.append(kdel,i)
                aldel0 = 1e8
    
    
    if wudel0 == 1e8 or wldel0 == 1e8:
        udel = (audel0+aldel0)/2.#(audel+aldel)/2
        ldel = (audel0+aldel0)/2.#(audel+aldel)/2
    elif audel0 == 1e8 or aldel0 == 1e8:
        udel = (wudel0+wldel0)/2.#(audel+aldel)/2
        ldel = (wudel0+wldel0)/2.#(audel+aldel)/2
    else:
        udel = (wudel0+wldel0+audel0+aldel0)/4.#(audel+aldel)/2
        ldel = (wudel0+wldel0+audel0+aldel0)/4.#(audel+aldel)/2
    # print udel,ldel
    
    
    return udel,ldel


## Main File
if __name__ == "__main__":
    
    s = 's2'
    t = '400'
    length = 100.

    comp = 'mac'
    # comp = 'fsl'

    s1length = np.array([210,210,205,196,185,178,170,165,160,145,140,123,115,112,108,101,101,90,85,80,78,75,70])*1.
    s2length = np.array([197,193,185,176,140,146,126,114,103,96,100,86,77,72,70,68,60,64,54,50,47,45,44])*1.
    s3length = np.array([185,150,100,95,83,76,72,63,60,49,50,41,39,36,34,33,30,31,28,30,29,28,27])*1.
    s4length = np.array([145,100,73,60,53,44,42,37,38,30,33,26,22,24,23,21,21,19,24,23,22,21,20])*1.
    s5length = np.array([78,70,52,43,37,32,29,27,26,23,20,20,23,21,20,19,19,18,18,16,16,15,14])*1.

    slength = np.array([])
    slength = np.append(slength,s1length)
    slength = np.append(slength,s2length)
    slength = np.append(slength,s3length)
    slength = np.append(slength,s4length)
    slength = np.append(slength,s5length)
    slength = slength*1.


    solidity = np.array(['s1','s2','s3','s4','s5'])
    tsr = np.linspace(150,700,23)
    
    dia = 6.


    sturb = np.zeros(115)

    k = 0
    for i in range(5):
        for j in range(23):
            s = solidity[i]
            t = str(int(tsr[j]))
            length = slength[k]

            udel,_ = fit(s,t,length,comp)
            sturb[k] = udel

            k += 1

    print 'sturb = np.array(',sturb.tolist(),')'