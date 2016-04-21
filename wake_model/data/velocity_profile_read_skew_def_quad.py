import csv
import numpy as np
from numpy import pi,sqrt,exp,fabs,log,sin,arctan,cosh
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,leastsq,minimize
from pyoptsparse import Optimization, SNOPT, pyOpt_solution
from scipy.stats import beta, norm
from scipy.special import erf,i0
from scipy.integrate import simps
from matplotlib import rcParams
import database_call as dbc
# rcParams['font.family'] = 'Times New Roman'


def pdf(x):
    return 1/sqrt(2*pi) * exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/sqrt(2))) / 2

def veldist(x,e,w,a,c,b,d,f):
    # y = np.zeros_like(x)
    # if b < 1.5:
    #     b = 1.5
    
    # y = w*fabs(x)**b + d
    # for i in range(np.size(y)):
    #     if y[i] > 1.:
    #         y[i] = 1.
    

    # g = -a/(sqrt(2.*pi))*exp(-(x)**2/(2.*d**2)) + 1.
    g = 1./(1 + exp(b*fabs(x)-a))
    # g = 1.
    
    # y = -b*(2./w)*pdf(t)*cdf(d*t) + 1.
    # y = -b*pdf(x) + 1.
    
    y = (-d/(w*sqrt(2.*pi))*exp(-(x+e)**2/(2.*w**2)))*g + 1.

    return y


def error(params):
    y = veldist(x, *params)
    # print params, np.sum((y - velo)**2)
    return np.sum((y - velo)**2)


def residuals(p,x,y):
    # if p[3] < 0.:
    #     penalty = penalty + 1000.
    # # if fabs(p[0]) < 0:
    # #     penalty = penalty + 1000.
    # # if fabs(p[3]-fabs(int)) > 0.5:
    # #     penalty = penalty + 1000.
    # if fabs(p[2]) < 0.:
    #     penalty = penalty + 1000.
    return (y - veldist(x,p[0],p[1],p[2],p[3],p[4],p[5],p[6]))**4
    # return np.sum((y - veldist(x,p[0],p[1],p[2],p[3]))**2)

def normgauss(x,e,w,c):
    # return c/(w*sqrt(2.*pi))*exp(-(x-e)**2/(2.*w**2))
    
    # return c*(e**2/(e+1.)*(x+1.)*exp(-e*x)) #lindley
    # return c/w*exp((e-x)/w)*exp(-exp((e-x)/w)) #gumbel
    
    
    # return c/(w*sqrt(2.*pi))*(exp(-(x-e)**2/(2.*w**2))+exp(-(x+e)**2/(2.*w**2))) #folded-normal
    return c*w*e*exp(w*x)*exp(e)*exp(-e*exp(w*x)) #Gompertz
    
    
    
    # t = (x - e)/w
    # return c*(2./w)*pdf(t)*cdf(a*t)

def normgaussw(x,e,w,c,d):
    # return c/(w*sqrt(2.*pi))*exp(-(x-e)**2/(2.*w**2))
    
    # return c*(e**2/(e+1.)*(x+1.)*exp(-e*x)) #lindley
    # return c/w*exp((e-x)/w)*exp(-exp((e-x)/w)) #gumbel
    
    
    # return c/(w*sqrt(2.*pi))*(exp(-(x-e)**2/(2.*w**2))+exp(-(x+e)**2/(2.*w**2))) #folded-normal
    return c*w*e*exp(w*x)*exp(e)*exp(-e*exp(w*x)) + d #Gompertz

def residualfun(x):
    return x[5] - veldist(x[4],x[0],x[1],x[2],x[3])
    
def residualnorm(p,x,y):
    return (y - normgauss(x,p[0],p[1],p[2]))**2

def obj_func(xdict):
    global yopt
    global xopt
    p = xdict['param']
    funcs = {}
    # yopt = kwargs['yopt']
    # xopt = kwargs['xopt']
    
    funcs['obj'] = max((yopt - veldist(xopt,p[0],p[1],p[2],p[3],p[4],p[5],p[6]))**4)
    
    fail = False

    return funcs, fail
    
def obj_funcng(xdict):
    global xdata
    global ydata
    p = xdict['param']
    funcs = {}
    # yopt = kwargs['yopt']
    # xopt = kwargs['xopt']
    
    funcs['obj'] = max((ydata - normgauss(xdata,p[0],p[1],p[2]))**2)
    
    fail = False

    return funcs, fail


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

def starccm_read(fdata,length,rad,wind,rot,tsr,s,paper,ploti):
    global yopt
    global xopt
    integrate = False
    # ploti = True

    dia = rad*2.

    for i in range(61):
        name = str(i+1)
        exec('pos'+name+' = np.array([])')
        exec('velo'+name+' = np.array([])')

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
    # print 'DATA IMPORTED'


    for i in range(30):
        name = str(i+1)

        #Ordering the data numerically by the position
        exec('pos'+name+', velo'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', velo'+name+'))))')
        #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
        exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvelo'+name+' = np.copy(velo'+name+'_0)')
        #Deleting wall boundary data
        exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5.*dia or pos'+name+'[j] < -5.*dia:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
        
        exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif velo'+name+'[j] > 1.*wind or fabs(pos'+name+'[j]) > 1.2*dia:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')




    start = length/30.
    x = np.linspace(start,length,30)

    #Deleting turbine region data
    # wr = np.linspace(-rad,rad,21)
    # k = 0
    # for i in range(30,51):
    #     name = str(i+1)
    #     exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] < sqrt(wr[k]**2+3.3**2) and pos'+name+'[j] > -sqrt(wr[k]**2+3.3**2):\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
    #     k += 1

    # print 'pos = np.array(',pos2.tolist(),')'
    # print 'velo = np.array(',velo2.tolist(),')'
    # print "a;lskdjf;lasjfas;dlkfjas;ldkfja;sldkjfa;sldkjfa;sldkjfa;sldkjfa;slkdfja;slkdjfa;slkdfja;sldkfja;lskfja"

    # fwrite = 'C:\\Users\\TingeyPC\\Documents\\zStar-CCM\\NACA0021\\MoveForward\\posvel.csv'
    # with open(fwrite,'w') as fp:
    #     a = csv.writer(fp)
    #
    #     data = pos2
    #     data = np.vstack([data,velo2])
    #
    #     a.writerows(data)



    # fst = 15
    # plt.figure(53)
    # # plt.subplot(1,4,1)
    # plt.plot(velo5/wind,pos5/dia,'b.')
    # plt.xlabel(r'$u/U_{\infty}$',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # # plt.title('$x/D$ = 1.11')
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(0,1.4)
    # plt.ylim(-5,5)
    #
    # plt.subplot(1,4,2)
    # plt.plot(velo45/wind,pos45/dia,'bo-')
    # plt.xlabel('Norm. velocity',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.title('$x/D$ = 8.33',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(0,2)
    # plt.ylim(-4,4)
    #
    # plt.subplot(1,4,3)
    # plt.plot(velo52/wind,pos52/dia,'bo-')
    # plt.xlabel('Norm. velocity',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.title('$x/D$ = 11.67',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(0,2)
    # plt.ylim(-4,4)
    #
    # plt.subplot(1,4,4)
    # plt.plot(velo56/wind,pos56/dia,'bo-')
    # plt.xlabel('Norm. velocity',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.title('$x/D$ = 13.33',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(0,2)
    # plt.ylim(-4,4)
    #
    # plt.show()




## Creating arrays for fitting
    regionfix = (rad)-0.98
    fix = 0.55 # Taking into account short wakes (.55)
    velomax = 0.1
    sreg = -0.3 #eliminating double veloex on s1 profiles
    svelomax = 0.03

    eu = np.zeros(30)
    wu = np.zeros(30)
    au = np.zeros(30)
    cu = np.zeros(30)
    bu = np.zeros(30)
    du = np.zeros(30)
    fu = np.zeros(30)
    # el = np.zeros(30)
    # wl = np.zeros(30)
    # al = np.zeros(30)
    # cl = np.zeros(30)

    # latmax1 = max(np.fabs(lat_up))
    # latmax2 = max(np.fabs(lat_low))
    # latmax = max(latmax1,latmax2)/dia

    for i in range(30):
        name = str(i+1)
        ind = str(i)
    #     exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'gup)):\n\tif velo'+name+'gup[j] > 0. or pos'+name+'gup[j] > 6.5 or pos'+name+'gup[j] < -1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'gup = np.delete(velo'+name+'gup,delvec[j])\n\tpos'+name+'gup = np.delete(pos'+name+'gup,delvec[j])')
    #     exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'glow)):\n\tif velo'+name+'glow[j] < 0. or pos'+name+'glow[j] < -6.5 or pos'+name+'glow[j] > 1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'glow = np.delete(velo'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
    #     # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'gup)):\n\tif velo'+name+'gup[j] > 0. or pos'+name+'gup[j] < -1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'gup = np.delete(velo'+name+'gup,delvec[j])\n\tpos'+name+'gup = np.delete(pos'+name+'gup,delvec[j])')
    #     # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+'glow)):\n\tif velo'+name+'glow[j] < 0. or pos'+name+'glow[j] > 1.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'glow = np.delete(velo'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
    #     exec('pos'+name+'gup = (-1./dia)*pos'+name+'gup')
    #     exec('pos'+name+'glow = (1./dia)*pos'+name+'glow')
    #     exec('velo'+name+'gup = (-1./rot)*velo'+name+'gup')
    #     exec('velo'+name+'glow = (1./rot)*velo'+name+'glow')
    
    
    
    
        # exec('delvec = np.array([])\nif x['+ind+'] < rad+fix:\n\tfor j in range(np.size(pos'+name+')):\n\t\tif pos'+name+'[j] > -regionfix and pos'+name+'[j] <regionfix:\n\t\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
    
    
    
    
    
    
    #     exec('delvec = np.array([])\nif x['+ind+'] < rad+fix:\n\tfor j in range(np.size(velo'+name+'glow)):\n\t\tif pos'+name+'glow[j] > -regionfix and velo'+name+'glow[j] > velomax:\n\t\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'glow = np.delete(velo'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
    #     # exec('delvec = np.array([])\nif x['+ind+'] < rad+fix:\n\tfor j in range(np.size(velo'+name+'glow)):\n\t\tif pos'+name+'glow[j] < regionfix and velo'+name+'glow[j] > velomax:\n\t\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'glow = np.delete(velo'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
    #     if s == 's1' and tsr <= 3.5:
    #         exec('delvec = np.array([])\nfor j in range(np.size(velo'+name+'glow)):\n\tif pos'+name+'glow[j] > sreg and velo'+name+'glow[j] > svelomax:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'glow = np.delete(velo'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')
    #     if s == 's2' and tsr <= 3.:
    #         exec('delvec = np.array([])\nfor j in range(np.size(velo'+name+'glow)):\n\tif pos'+name+'glow[j] > sreg and velo'+name+'glow[j] > svelomax:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+'glow = np.delete(velo'+name+'glow,delvec[j])\n\tpos'+name+'glow = np.delete(pos'+name+'glow,delvec[j])')

## Velocity Distribution Fitting
        diap = dia
        rotp = 1.
        # gam = np.fabs(gamma_up[i]/rot)
        # mpos = lat_up[i]
        
        pmen = 0. # 0.
        pspr = 2. # 2.
        pscl = 2. # 2.
        prat = 2. # 2.
        ptns = 1. # 1.
        
        c = .8
        f = 1.
        # if 
        
        
        

        # x0 = np.array([e,w,a,c])
        bounds = ((0.0,1.0), (-0.2, 0.0), (-60.0, 0.0), (-0.2, 0.0))
        # meth = 'L-BFGS-B'
        meth = 'SLSQP'
        # meth = 'Nelder-Mead'
        dis = 'disp'
        tol = 'xtol'

        ft = 1e-8
        xt = 1e-8
        mf = int(1e6)

        exec('betaparam,_ = leastsq(residuals,np.array([pmen,pspr,ptns,c,prat,pscl,f]),args=(pos'+name+'/dia,velo'+name+'/wind),ftol=ft, xtol=xt,maxfev=mf)\neu['+ind+'] = betaparam[0]\nwu['+ind+'] = betaparam[1]\nau['+ind+'] = betaparam[2]\ncu['+ind+'] = betaparam[3]\nbu['+ind+'] = betaparam[4]\ndu['+ind+'] = betaparam[5]\nfu['+ind+'] = betaparam[6]')
        # exec('betaparam = np.polyfit(pos'+name+'/dia,velo'+name+'/wind,4)')#\nwu['+ind+'] = betaparam[0]\nau['+ind+'] = betaparam[1]\ndu['+ind+'] = betaparam[2]')
        
        #SNOPT#######
        # params0 = np.array([pmen,pspr,ptns,c,prat,pscl,f])
        # exec('yoptt = velo'+name+'/wind')
        # exec('xoptt = pos'+name+'/dia')
        # yopt = yoptt
        # xopt = xoptt
        # optProb = Optimization('Vel_Fit', obj_func)
        # optProb.addObj('obj')
        # 
        # 
        # optProb.addVarGroup('param', 7, 'c', lower=0., upper=None, value=params0)
        # 
        # opt = SNOPT()#yopt=yopt, xopt=xopt)
        # opt.setOption('Scale option',0)
        # res = opt(optProb, sens=None)#, yopt=yopt, xopt=xopt)
        # # print res
        # 
        # paramf = res.xStar['param']
        # 
        # exec('eu['+ind+'] = paramf[0]\nwu['+ind+'] = paramf[1]\nau['+ind+'] = paramf[2]\ncu['+ind+'] = paramf[3]\nbu['+ind+'] = paramf[4]\ndu['+ind+'] = paramf[5]\nfu['+ind+'] = paramf[6]')
        
        
        
        
        
        
        
        
        
        # exec('betaparam,_ = leastsq(residuals,np.array([e,w,a,c]),args=(pos'+name+',velo'+name+'),ftol=ft, xtol=xt,maxfev=mf)\neu['+ind+'] = betaparam[0]\nwu['+ind+'] = betaparam[1]\nau['+ind+'] = betaparam[2]\ncu['+ind+'] = betaparam[3]')
        # exec('betaparam,_ = curve_fit(veldist,pos'+name+'/dia,velo'+name+'/wind,p0=(e,w,a,c,b,d,f))\neu['+ind+'] = betaparam[0]\nwu['+ind+'] = betaparam[1]\nau['+ind+'] = betaparam[2]\ncu['+ind+'] = betaparam[3]\nbu['+ind+'] = betaparam[4]\ndu['+ind+'] = betaparam[5]\nfu['+ind+'] = betaparam[6]')
        # exec('betaparam,_ = curve_fit(veldist,pos'+name+'gup,velo'+name+'gup,p0=(e,w,a,c))\neu['+ind+'] = betaparam[0]\nwu['+ind+'] = betaparam[1]\nau['+ind+'] = betaparam[2]\ncu['+ind+'] = betaparam[3]')
        # exec('res = minimize(residuals,x0,args=(pos'+name+'gup,velo'+name+'gup),method=meth,bounds=bounds,options={tol: 1e-8, dis: False})\neu['+ind+'] = res.x[0]\nwu['+ind+'] = res.x[1]\nau['+ind+'] = res.x[2]\ncu['+ind+'] = res.x[3]')
        # exec('res = minimize(residuals,x0,args=(pos'+name+'gup,velo'+name+'gup),method=meth,options={tol: 1e-8, dis: False})\neu['+ind+'] = res.x[0]\nwu['+ind+'] = res.x[1]\nau['+ind+'] = res.x[2]\ncu['+ind+'] = res.x[3]')

        # exec('print eu['+ind+'],wu['+ind+'],au['+ind+'],-intu['+ind+']')
        # ceu = np.array( [-0.051533741732590066, 0.603645756517617, 3.404058219228679] )
        # cwu = np.array( [0.0012367836421159789, -0.016292716663412185, 0.08228677443522553, 0.6891188133055813] )
        # cau = np.array( [0.006775814669580242, -0.09990512632914271, 0.6437647406588762, -5.547388443086547] )
        # ccu = np.array( [-0.001095215841509046, 0.025087453250091336, -0.17923773313904467, 0.6219114941017065, 9.724403752075677] )
        # Plotting
        
        # men = np.array( [-0.0007448786610163438, 0.011700465818493566, -0.005332505770684337] )
        # spr = np.array( [6.462355161093711, 7.079901300173991, 12.102886237210939] )
        # scl = np.array( [8.509272717226171, 7.023483471068396, 27.707846411384697] )
        # rat = np.array( [-2.107186196351149, 44.93845180541949] )
        # tns = np.array( [-1.4660542829002265, 30.936653231840257] )
        # 
        # men = np.array( [-0.00059737414699399, 0.009890587474506057, -0.0016721254639608882] )
        # spr = np.array( [6.470586725452819, 8.205138442690668, 13.656983372616262] )
        # scl = np.array( [6.639808894608728, 5.477607580787858, 21.13678312202297] )
        # rat = np.array( [-2.0794010451530873, 44.798557035611] )
        # tns = np.array( [-1.43164706657537, 30.761785195818447] )
        
        #s2_400
        # men = np.array( [-0.00059737414699399, 0.009890587474506057, -0.0016721254639608882] )
        # spr = np.array( [-0.005652314031253564, 0.06923002880544946, 0.526304136118912] )
        # scl = np.array( [6.639808894608728, 5.477607580787858, 21.13678312202297] )
        # rat = np.array( [-2.0794010451530873, 44.798557035611] )
        # tns = np.array( [-1.43164706657537, 30.761785195818447] )
        # 
        # men = np.array( [-0.0006344223751663201, 0.01055675755786011, -0.004073212523707764] )
        # spr = np.array( [-0.005187125854670714, 0.06397918461247416, 0.543874357807372] )
        # scl = np.array( [6.667328694868336, 5.617498827673229, 21.520026361522778] )
        # rat = np.array( [-2.129054494312758, 45.17191461412915] )
        # tns = np.array( [-1.5569348878268718, 31.913143231782648] )
        
        # #s1_300
        # men = np.array( [4.0348056629563466e-05, -0.0005838894472673057, -0.12319864474746939] )
        # spr = np.array( [-0.00014193456712206312, 0.004286965394049379, 0.4333125819471815] )
        # scl = np.array( [12.041023363327602, 17.35527083827463, 18.614541065902408] )
        # rat = np.array( [-0.7860582227038313, 39.02603416198663] )
        # tns = np.array( [-0.4503431986903151, 21.43333223988026] )
        # 
        # #s4_600
        # men = np.array( [-0.003958655969610645, -0.012497079414312735, 0.0015000042748494803] )
        # spr = np.array( [0.06719808524021777, -0.2665046830015697, 0.5869144851566086] )
        # scl = np.array( [1.2729325247012917, 0.7531051165932696, 15.479368036459785] )
        # rat = np.array( [-9.613284996801559, 47.15828538608972] )
        # tns = np.array( [-9.26983482585302, 39.87384954967598] )
        
        # s3_150
        # men = np.array( [0.00039827881005907734, -0.00904675311786268, -0.10910433324181813] )
        # spr = np.array( [0.0052653322222484915, -0.19787486142936342, 2.154796622458931] )
        # scl = np.array( [8.831676552790272, 5.463124302900828, 16.562600343794855] )
        # rat = np.array( [0.18543157941280936, 6.293602107581222] )
        # tns = np.array( [0.2555826800561036, 1.831038353492003] )
        
        if s == 's1':
            solidity = 0.15
        elif s == 's2':
            solidity = 0.25
        elif s == 's3':
            solidity = 0.5
        elif s == 's4':
            solidity = 0.75
        elif s == 's5':
            solidity = 1.
        
        men,spr,scl,rat,tns = dbc.velocity(tsr,solidity)
        # men = np.array( [0.00039827881005907734, -0.00904675311786268, -0.10910433324181813] )
        # spr = np.array( [0.005711491642366975, -0.20557616791507677, 2.088706809453395] )
        # scl = np.array( [8.871204327482078, 6.4750000000000005, 19.218163360943127] )
        # rat = np.array( [0.18543157941280936, 6.293602107581222] )
        # tns = np.array( [0.2555826800561036, 1.831038353492003] )
        
        # men = np.array( [0.00039827881005907734, -0.00904675311786268, -0.10910433324181813] )
        # spr = np.array( [-0.0002761554016230923, 0.009356520540766187, 0.4079297314900556] )
        # scl = np.array( [8.911515879192555, 7.40031677482411, 21.695418281284383] )
        # rat = np.array( [0.18543157941280936, 6.293602107581222] )
        # tns = np.array( [0.2555826800561036, 1.831038353492003] )
        # men = np.array( [0.00040389238585919425, 0.0058920789865672795, -0.6245463883539139] )
        # spr = np.array( [-0.0016450772484362373, 0.02455176915483611, 1.026988702523543] )
        # scl = np.array( [10.944444445983667, 7.88, 25.0] )
        # rat = np.array( [-0.39583161297691366, 26.03967341597475] )
        # tns = np.array( [-0.1563482845263192, 12.248022357162302] )
        xd = x/dia
        men_v = men[0]*xd[i]**2 + men[1]*xd[i] + men[2]
        if men_v > 0.5:
            men_v = 0.5
        elif men_v < -0.5:
            men_v = -0.5
        
        # spr_v = spr[2]/(spr[1]*sqrt(2.*pi))*exp(-(xd[i]-spr[0])**2/(2.*spr[1]**2))
        spr_v = spr[0]*xd[i]**2 + spr[1]*xd[i] + spr[2]
        if spr_v < 0.35:
            spr_v = 0.35
        elif spr_v > 4.:
            spr_v = 4.
        
        # scl_v = scl[2]/(scl[1]*sqrt(2.*pi))*exp(-(xd[i]-scl[0])**2/(2.*scl[1]**2))
        # scl_v = scl[2]/(scl[1]*sqrt(2.*pi))*(exp(-(xd[i]-scl[0])**2/(2.*scl[1]**2))+exp(-(xd[i]+scl[0])**2/(2.*scl[1]**2))) #fold
        scl_v = scl[2]*scl[1]*scl[0]*exp(scl[1]*xd[i])*exp(scl[0])*exp(-scl[0]*exp(scl[1]*xd[i])) #gom
        
        rat_v = rat[0]*xd[i] + rat[1]
        if rat_v < 0.:
            rat_v = 0.
        
        tns_v = tns[0]*xd[i] + tns[1]
        if tns_v < 0.:
            tns_v = 0.
        
        if ploti == True:
            if paper == False:
                plt.figure(1)
                plt.subplot(5,6,i+1)
                color = 'bo'
                color2 = 'r-'
                color3 = 'c.'
                color4 = 'g.'
                lab = str(x[i]/dia)+'D'
                exec('xfit = np.linspace(min(pos'+name+'/dia)-1.,max(pos'+name+'/dia)+1.,500)')
                exec('plt.plot(velo'+name+'/wind,pos'+name+'/diap,color,label=lab)')
                exec('plt.plot(veldist(xfit,eu['+ind+'],wu['+ind+'],au['+ind+'],cu['+ind+'],bu['+ind+'],du['+ind+'],fu['+ind+']),xfit,color2,linewidth=2)')
                plt.plot(veldist(xfit,men_v,spr_v,tns_v,1.,rat_v,scl_v,1.),xfit,'c',linewidth=1)
                # exec('plt.plot(np.polyval(betaparam,xfit),xfit,color2,linewidth=2)')
                plt.xlim(0.,1.5)
                # plt.ylim(-4.,4.)
                # plt.legend(loc=1)
                plt.xlabel('Normalized Velocity')
                plt.ylabel('$y/D$')

            elif paper == True:
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
                exec('xfit = np.linspace(min(pos'+name+'/dia)-1.,max(pos'+name+'/dia)+1.,500)')
                if i == 5 and wind == 15.:
                    exec('plt.plot(velo'+name+'/wind,pos'+name+'/diap,color,label=lab)')
                    exec('plt.plot(veldist(xfit,eu['+ind+'],wu['+ind+'],au['+ind+'],cu['+ind+'],bu['+ind+'],du['+ind+'],fu['+ind+']),xfit,color2,linewidth=2,label= lab2)')
                    plt.xlim(0.,1.5)
                    # plt.ylim(-4.,4.)
                    plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
                else:
                    exec('plt.plot(velo'+name+'/wind,pos'+name+'/diap,color)')
                    exec('plt.plot(veldist(xfit,eu['+ind+'],wu['+ind+'],au['+ind+'],cu['+ind+'],bu['+ind+'],du['+ind+'],fu['+ind+']),xfit,color2,linewidth=2)')
                    plt.xlim(0.,1.5)
                    # plt.ylim(-4.,4.)
                if wind == 15.:
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



            
        print i+1,'Velo: Wind =',wind,'Rot =',rot,s,'TSR =',tsr

    # fst = 15
    # plt.figure(55)
    # plt.subplot(1,2,1)
    # color = 'bo'
    # color2 = 'r-'
    # color3 = 'c.'
    # color4 = 'g.'
    # lab = str(x[i]/dia)+'D'
    # xfit = np.linspace(min(pos5/dia),max(pos5/dia),500)
    # plt.plot(velo5/wind,pos5/diap,color)
    # plt.plot(veldist(xfit,eu[4],wu[4],au[4],cu[4],bu[4],du[4],fu[4]),xfit,color2,linewidth=2)
    # plt.xlabel(r'$u/U_{\infty}$',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(0,1.4)
    # plt.ylim(-5,5)
    # plt.subplot(1,4,3)
    # color = 'bo'
    # color2 = 'r-'
    # color3 = 'c.'
    # color4 = 'g.'
    # lab = str(x[i]/dia)+'D'
    # xfit = np.linspace(min(pos5/dia),max(pos5/dia),500)
    # plt.plot(velo5/wind,pos5/diap,color)
    # plt.plot(veldist(xfit,eu[4],wu[4],au[4],cu[4],bu[4],du[4],fu[4]),xfit,color2,linewidth=2)
    # plt.xlabel(r'$u/U_{\infty}$',fontsize=fst)
    # plt.ylabel('$y/D$',fontsize=fst)
    # plt.xticks(fontsize=fst)
    # plt.yticks(fontsize=fst)
    # plt.xlim(1,1.2)
    # plt.ylim(0,5)

    # plt.show()

    

        # print i+1,'veldist  Wind =',wind,'Rot =',rot,s,'TSR =',tsr
    # plt.show()

    xw = np.linspace(-rad,rad,21)
    xh = np.linspace(-12.*rad,-12*3./10.,10)

    # return x,xw,xh,eu,wu,au,cu,bu,du,euw,wuw,auw,cuw,buw,duw,euh,wuh,auh,cuh,buh,duh
    return x,eu,wu,au,cu,bu,du


def fit(s,t,length,plottrend,ploti):
    global xdata
    global ydata
    
    read_data = 4

    plott = False
    plotac = False
    # plottrend = True
    paper = False

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

    # fdata = '/Users/ning1/Documents/Flow Lab/STAR-CCM+/NACA0021/MoveForward/test.csv'
    fdata = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/Velocity Sections/'+wfit+'.csv'
    fdata2 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel14/Velocity/'+wfit2+'.csv'
    fdata3 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel12/Velocity/'+wfit3+'.csv'
    fdata4 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/vel16/Velocity/'+wfit4+'.csv'
    fdata5 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot17/Velocity/'+wfit5+'.csv'
    fdata6 = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/CrossValidate/rot18/Velocity/'+wfit6+'.csv'



    if read_data >=1:
        xd,eud,wud,aud,cud,bud,dud = starccm_read(fdata,length,rad,wind,rot,tsr,s,paper,ploti)
    if read_data >=2:
        xd2,eud2,wud2,aud2,cud2,bud2,dud2 = starccm_read(fdata2,length2,rad,wind2,rot2,tsr,s,paper,ploti)
    if read_data >=3:
        xd3,eud3,wud3,aud3,cud3,bud3,dud3 = starccm_read(fdata3,length3,rad,wind3,rot3,tsr,s,paper,ploti)
    if read_data >=4:
        xd4,eud4,wud4,aud4,cud4,bud4,dud4 = starccm_read(fdata4,length4,rad,wind4,rot4,tsr,s,paper,ploti)
    if read_data >=5:
        xd5,eud5,wud5,aud5,cud5,bud5,dud5 = starccm_read(fdata5,length5,rad,wind5,rot5,tsr,s,paper,ploti)
    if read_data >=6:
        xd6,eud6,wud6,aud6,cud6,bud6,dud6 = starccm_read(fdata6,length6,rad,wind6,rot6,tsr,s,paper,ploti)

    x = np.zeros(30*read_data)
    eu = np.zeros(30*read_data)
    wu = np.zeros(30*read_data)
    au = np.zeros(30*read_data)
    cu = np.zeros(30*read_data)
    bu = np.zeros(30*read_data)
    du = np.zeros(30*read_data)



    for i in range(30):
        if read_data >=1:
            x[read_data*i] = xd[i]
            eu[read_data*i] = eud[i]
            wu[read_data*i] = wud[i]
            au[read_data*i] = aud[i]
            cu[read_data*i] = cud[i]
            bu[read_data*i] = bud[i]
            du[read_data*i] = dud[i]
        if read_data >=2:
            x[read_data*i+1] = xd2[i]
            eu[read_data*i+1] = eud2[i]
            wu[read_data*i+1] = wud2[i]
            au[read_data*i+1] = aud2[i]
            cu[read_data*i+1] = cud2[i]
            bu[read_data*i+1] = bud2[i]
            du[read_data*i+1] = dud2[i]
        if read_data >=3:
            x[read_data*i+2] = xd3[i]
            eu[read_data*i+2] = eud3[i]
            wu[read_data*i+2] = wud3[i]
            au[read_data*i+2] = aud3[i]
            cu[read_data*i+2] = cud3[i]
            bu[read_data*i+2] = bud3[i]
            du[read_data*i+2] = dud3[i]
        if read_data >=4:
            x[read_data*i+3] = xd4[i]
            eu[read_data*i+3] = eud4[i]
            wu[read_data*i+3] = wud4[i]
            au[read_data*i+3] = aud4[i]
            cu[read_data*i+3] = cud4[i]
            bu[read_data*i+3] = bud4[i]
            du[read_data*i+3] = dud4[i]
        if read_data >=5:
            x[read_data*i+4] = xd5[i]
            eu[read_data*i+4] = eud5[i]
            wu[read_data*i+4] = wud5[i]
            au[read_data*i+4] = aud5[i]
            cu[read_data*i+4] = cud5[i]
            bu[read_data*i+4] = bud5[i]
            du[read_data*i+4] = dud5[i]
        if read_data >=6:
            x[read_data*i+5] = xd6[i]
            eu[read_data*i+5] = eud6[i]
            wu[read_data*i+5] = wud6[i]
            au[read_data*i+5] = aud6[i]
            cu[read_data*i+5] = cud6[i]
            bu[read_data*i+5] = bud6[i]
            du[read_data*i+5] = dud6[i]



    x = x/dia
    if paper == True:
        fsp = 20
    elif paper == False:
        fsp = 15
    if plottrend == True:
        if paper == True:
            # plt.figure(10)
            plt.figure(11)
            # plt.subplot(3,2,1)
            plt.plot(x,eu,'bo')
            plt.xlabel('$x/D$',fontsize=fsp)
            plt.ylabel('Parameter Value',fontsize=fsp)
            plt.xticks(fontsize=fsp)
            plt.yticks(fontsize=fsp)
            # plt.title('Mean')
            plt.ylim(-0.1,0.1)
            plt.figure(12)
            # plt.subplot(3,2,2)
            plt.plot(x,wu,'bo')
            plt.xlabel('$x/D$',fontsize=fsp)
            plt.ylabel('Parameter Value',fontsize=fsp)
            plt.xticks(fontsize=fsp)
            plt.yticks(fontsize=fsp)
            # plt.title('St Dev')
            plt.ylim(0,1.5)
            plt.figure(13)
            # plt.subplot(3,2,3)
            plt.plot(x,du,'bo')
            plt.xlabel('$x/D$',fontsize=fsp)
            plt.ylabel('Parameter Value',fontsize=fsp)
            plt.xticks(fontsize=fsp)
            plt.yticks(fontsize=fsp)
            # plt.title('Scale')
            plt.ylim(0,2)
            plt.figure(14)
            # plt.subplot(3,2,4)
            plt.plot(x,bu,'bo')
            plt.xlabel('$x/D$',fontsize=fsp)
            plt.ylabel('Parameter Value',fontsize=fsp)
            plt.xticks(fontsize=fsp)
            plt.yticks(fontsize=fsp)
            # plt.title('Rate')
            plt.ylim(0,50)
            plt.figure(15)
            # plt.subplot(3,2,5)
            plt.plot(x,au,'bo')
            plt.xlabel('$x/D$',fontsize=fsp)
            plt.ylabel('Parameter Value',fontsize=fsp)
            plt.xticks(fontsize=fsp)
            plt.yticks(fontsize=fsp)
            # plt.title('Spread')
            plt.ylim(0,40)
            # plt.ylim(-4,4)
            # plt.subplot(2,3,4)
            # plt.plot(x,cu,'bo')
            # plt.title('Sig. Translation')
            # plt.subplot(2,3,5)
            # plt.plot(x,bu,'bo')
            # plt.title('Exp. Rate')
            # plt.subplot(2,3,3)
            # plt.plot(x,du,'bo')
            # plt.title('Sig. Rate')
            # plt.ylim(-4,4)
        elif paper == False:
            plt.figure(10)
            plt.subplot(3,2,1)
            plt.plot(x,eu,'bo')
            plt.xlabel('$x/D$')
            plt.ylabel('Parameter Value')
            plt.title('Mean')
            plt.ylim(-0.1,0.1)
            plt.subplot(3,2,2)
            plt.plot(x,wu,'bo')
            plt.xlabel('$x/D$')
            plt.ylabel('Parameter Value')
            plt.title('St Dev')
            plt.ylim(0,1.5)
            plt.subplot(3,2,3)
            plt.plot(x,du,'bo')
            plt.xlabel('$x/D$')
            plt.ylabel('Parameter Value')
            plt.title('Scale')
            plt.ylim(0,2)
            plt.subplot(3,2,4)
            plt.plot(x,bu,'bo')
            plt.xlabel('$x/D$')
            plt.ylabel('Parameter Value')
            plt.title('Rate')
            plt.ylim(0,50)
            plt.subplot(3,2,5)
            plt.plot(x,au,'bo')
            plt.xlabel('$x/D$')
            plt.ylabel('Parameter Value')
            plt.title('Spread')
            plt.ylim(0,40)
            # plt.ylim(-4,4)
            # plt.subplot(2,3,4)
            # plt.plot(x,cu,'bo')
            # plt.title('Sig. Translation')
            # plt.subplot(2,3,5)
            # plt.plot(x,bu,'bo')
            # plt.title('Exp. Rate')
            # plt.subplot(2,3,3)
            # plt.plot(x,du,'bo')
            # plt.title('Sig. Rate')
            # plt.ylim(-4,4)


    # plt.show()



## Deleting Values (outlier)
    xeu = np.copy(x)
    xwu = np.copy(x)
    xdu = np.copy(x)
    xbu = np.copy(x)
    xau = np.copy(x)
    dele = np.array([])
    delw = np.array([])
    deld = np.array([])
    delb = np.array([])
    dela = np.array([])
    lendel = 9./10.
    radfix = 0.5+0.08
    for i in range(np.size(x)):
        if fabs(eu[i]) > 0.5 or xeu[i] > lendel*length/dia or xeu[i] < radfix:
            dele = np.append(dele,i)
        if fabs(wu[i]) > 2. or bu[i] < 0. or au[i] < 0. or wu[i] < 0. or xwu[i] > lendel*length/dia or xwu[i] < radfix:
            delw = np.append(delw,i)
        if fabs(du[i]) > 12. or bu[i] < 0. or au[i] < 0. or du[i] < 0. or xdu[i] > lendel*length/dia or xdu[i] < radfix:
            deld = np.append(deld,i)
        if fabs(bu[i]) > 50. or bu[i] < 0. or au[i] < 0. or xbu[i] > lendel*length/dia or xbu[i] < radfix:
            delb = np.append(delb,i)
        if fabs(au[i]) > 70. or bu[i] < 0. or au[i] < 0. or xau[i] > lendel*length/dia or xau[i] < radfix:
            dela = np.append(dela,i)
    dele = dele[::-1]
    delw = delw[::-1]
    deld = deld[::-1]
    delb = delb[::-1]
    dela = dela[::-1]
    
    for i in range(np.size(dele)):
        eu = np.delete(eu,dele[i])
        xeu = np.delete(xeu,dele[i])
    for i in range(np.size(delw)):
        wu = np.delete(wu,delw[i])
        xwu = np.delete(xwu,delw[i])
    for i in range(np.size(deld)):
        du = np.delete(du,deld[i])
        xdu = np.delete(xdu,deld[i])
    for i in range(np.size(delb)):
        bu = np.delete(bu,delb[i])
        xbu = np.delete(xbu,delb[i])
    for i in range(np.size(dela)):
        au = np.delete(au,dela[i])
        xau = np.delete(xau,dela[i])
    
    # Deleting values that the sigmoid curve does not fit correctly
    # delp = np.array([])
    # for i in range(np.size(xbu)):
    #     if or bu[i] < 0. or au[i] < 0.:
    #         delp = np.append(delb,i)
    # delp = delp[::-1]
    # for i in range(np.size(delp)):
    #     wu = np.delete(wu,delp[i])
    #     xwu = np.delete(xwu,delp[i])
    #     du = np.delete(du,delp[i])
    #     xdu = np.delete(xdu,delp[i])
    #     bu = np.delete(bu,delp[i])
    #     xbu = np.delete(xbu,delp[i])
    #     au = np.delete(au,delp[i])
    #     xau = np.delete(xau,delp[i])
        
    #
    #
    # if plott == True:
    #     diap = 1.
    #     plt.figure(4)
    #     plt.subplot(2,2,1)
    #     plt.plot(xeu/diap,eu/diap,'bo')
    #     plt.plot(xel/diap,el/diap,'ro')
    #     plt.title(r'Location ($\xi$)')
    #     plt.ylim(-1,1)
    #
    #     plt.subplot(2,2,2)
    #     plt.plot(xwu/diap,wu,'bo')
    #     plt.plot(xwl/diap,wl,'ro')
    #     plt.title(r'Scale ($\omega$)')
    #     plt.ylim(-0.05,0.15)
    #
    #     plt.subplot(2,2,3)
    #     plt.plot(xau/diap,au,'bo')
    #     plt.plot(xal/diap,al,'ro')
    #     plt.title(r'veldist ($\alpha$)')
    #     plt.ylim(-5,5)
    #
    #     plt.subplot(2,2,4)
    #     plt.plot(xcu/diap,cu,'bo')
    #     plt.plot(xcl/diap,cl,'ro')
    #     plt.title('Resize ($c$)')
    #     plt.show()
    # print xal,al
## Deletinge Values (turbulent)
    # devrange = 3.0 #0.4
    # devrangew = 0.05
    # devrangea = 5.
    # averange = 0.07
    # veldisttest = 2.8
    # avechange = 0.04
    #
    # #Spread Variation
    # xtest = xwu[0]
    # wutest = np.array([])
    # kdel = np.array([])
    # k = 0
    # for i in range(np.size(xwu)):
    #     xtest0 = xwu[i]
    #     wut0 = wu[i]
    #     if xtest == xtest0:
    #         wutest = np.append(wutest,wut0)
    #         kdel = np.append(kdel,i)
    #     else:
    #         xtest = xtest0
    #         wut = wut0
    #         testave = fabs(np.average(wutest))
    #         # testrange = fabs(max(wutest)-min(wutest))
    #         # if testrange > devrangew and xwu[i] > 0.5*(length/6):
    #         #     wudel0 = xwu[kdel[0]]
    #         #     break
    #         # else:
    #         #     wutest = np.array([])
    #         #     kdel = np.array([])
    #         #     wutest = np.append(wutest,wut0)
    #         #     kdel = np.append(kdel,i)
    #         #     wudel0 = 1e8
    #
    #         if k == 0:
    #             testave0 = testave
    #             k += 1
    #         else:
    #             if testave - testave0 > avechange and xwu[i] > 0.5*(length/6):
    #                 wudel0 = xwu[kdel[0]]
    #                 break
    #             else:
    #                 wutest = np.array([])
    #                 kdel = np.array([])
    #                 wutest = np.append(wutest,wut0)
    #                 kdel = np.append(kdel,i)
    #                 wudel0 = 1e8
    #                 k += 1
    #
    # xtest = xwl[0]
    # wltest = np.array([])
    # kdel = np.array([])
    # for i in range(np.size(xwl)):
    #     xtest0 = xwl[i]
    #     wlt0 = wl[i]
    #     if xtest == xtest0:
    #         wltest = np.append(wltest,wlt0)
    #         kdel = np.append(kdel,i)
    #     else:
    #         xtest = xtest0
    #         wlt = wlt0
    #         testave = fabs(np.average(wltest))
    #         # testrange = fabs(max(wltest)-min(wltest))
    #         # if testrange > devrangew and xwl[i] > 0.5*(length/6):
    #         #     wldel0 = xwl[kdel[0]]
    #         #     break
    #         # else:
    #         #     wltest = np.array([])
    #         #     kdel = np.array([])
    #         #     wltest = np.append(wltest,wlt0)
    #         #     kdel = np.append(kdel,i)
    #         #     wldel0 = 1e8
    #
    #         if k == 0:
    #             testave0 = testave
    #             k += 1
    #         else:
    #             if testave - testave0 > avechange and xwl[i] > 0.5*(length/6):
    #                 wldel0 = xwu[kdel[0]]
    #                 break
    #             else:
    #                 wltest = np.array([])
    #                 kdel = np.array([])
    #                 wltest = np.append(wltest,wlt0)
    #                 kdel = np.append(kdel,i)
    #                 wldel0 = 1e8
    #                 k += 1
    #
    # #veldist Variation
    # xtest = xau[0]
    # autest = np.array([])
    # kdel = np.array([])
    # for i in range(np.size(xau)):
    #     xtest0 = xau[i]
    #     aut0 = au[i]
    #     if xtest == xtest0:
    #         autest = np.append(autest,aut0)
    #         kdel = np.append(kdel,i)
    #     else:
    #         xtest = xtest0
    #         aut = aut0
    #         testrange = fabs(max(autest)-min(autest))
    #         if testrange > devrangea and xau[i] > 0.5*(length/6):
    #             audel0 = xau[kdel[0]]
    #             break
    #         else:
    #             autest = np.array([])
    #             kdel = np.array([])
    #             autest = np.append(autest,aut0)
    #             kdel = np.append(kdel,i)
    #             audel0 = 1e8
    #
    # xtest = xal[0]
    # altest = np.array([])
    # kdel = np.array([])
    # for i in range(np.size(xal)):
    #     xtest0 = xal[i]
    #     alt0 = al[i]
    #     if xtest == xtest0:
    #         altest = np.append(altest,alt0)
    #         kdel = np.append(kdel,i)
    #     else:
    #         xtest = xtest0
    #         alt = alt0
    #         testrange = fabs(max(altest)-min(altest))
    #         if testrange > devrangea and xal[i] > 0.5*(length/6):
    #             aldel0 = xal[kdel[0]]
    #             break
    #         else:
    #             altest = np.array([])
    #             kdel = np.array([])
    #             altest = np.append(altest,alt0)
    #             kdel = np.append(kdel,i)
    #             aldel0 = 1e8
    #
    # # xtest = xau[0]
    # # autest = np.array([])
    # # kdel = np.array([])
    # # k = 0
    # # for i in range(np.size(xau)):
    # #     xtest0 = xau[i]
    # #     aut0 = au[i]
    # #     if xtest == xtest0:
    # #         autest = np.append(autest,aut0)
    # #         kdel = np.append(kdel,i)
    # #     else:
    # #         xtest = xtest0
    # #         aut = aut0
    # #         testave = fabs(np.average(autest))
    # #         if testave < veldisttest:
    # #             audel0 = xau[kdel[0]]
    # #             break
    # #         else:
    # #             autest = np.array([])
    # #             kdel = np.array([])
    # #             autest = np.append(autest,aut0)
    # #             kdel = np.append(kdel,i)
    # #             testave0 = testave
    # #             k += 1
    # #             audel0 = 1e8
    # #
    # #         # if k == 0:
    # #         #     testave0 = testave
    # #         #     k += 1
    # #         # else:
    # #         #     if xau[i] > 0.3*(length/6):
    # #         #         if testave < testave0+averange and testave > testave0-averange:
    # #         #             audel0 = xau[kdel[0]]
    # #         #             break
    # #         #         else:
    # #         #             autest = np.array([])
    # #         #             kdel = np.array([])
    # #         #             autest = np.append(autest,aut0)
    # #         #             kdel = np.append(kdel,i)
    # #         #             testave0 = testave
    # #         #             k += 1
    # #         #     else:
    # #         #         autest = np.array([])
    # #         #         kdel = np.array([])
    # #         #         autest = np.append(autest,aut0)
    # #         #         kdel = np.append(kdel,i)
    # #         #         testave0 = testave
    # #         #         k += 1
    # #
    # # xtest = xal[0]
    # # altest = np.array([])
    # # kdel = np.array([])
    # # k = 0
    # # for i in range(np.size(xal)):
    # #     xtest0 = xal[i]
    # #     alt0 = al[i]
    # #     if xtest == xtest0:
    # #         altest = np.append(altest,alt0)
    # #         kdel = np.append(kdel,i)
    # #     else:
    # #         xtest = xtest0
    # #         alt = alt0
    # #         testave = fabs(np.average(altest))
    # #         if testave < veldisttest:
    # #             aldel0 = xal[kdel[0]]
    # #             break
    # #         else:
    # #             altest = np.array([])
    # #             kdel = np.array([])
    # #             altest = np.append(altest,alt0)
    # #             kdel = np.append(kdel,i)
    # #             testave0 = testave
    # #             k += 1
    # #             aldel0 = 1e8
    # #
    # #         # if k == 0:
    # #         #     testave0 = testave
    # #         #     k += 1
    # #         # else:
    # #         #     if xal[i] > 0.3*(length/6):
    # #         #         if testave < testave0+averange and testave > testave0-averange:
    # #         #             aldel0 = xal[kdel[0]]
    # #         #             break
    # #         #         else:
    # #         #             altest = np.array([])
    # #         #             kdel = np.array([])
    # #         #             altest = np.append(altest,alt0)
    # #         #             kdel = np.append(kdel,i)
    # #         #             testave0 = testave
    # #         #             k += 1
    # #         #     else:
    # #         #         altest = np.array([])
    # #         #         kdel = np.array([])
    # #         #         altest = np.append(altest,alt0)
    # #         #         kdel = np.append(kdel,i)
    # #         #         testave0 = testave
    # #         #         k += 1
    #
    #
    # if wudel0 == 1e8 or wldel0 == 1e8:
    #     udel = (audel0+aldel0)/2.#(audel+aldel)/2
    #     ldel = (audel0+aldel0)/2.#(audel+aldel)/2
    # elif audel0 == 1e8 or aldel0 == 1e8:
    #     udel = (wudel0+wldel0)/2.#(audel+aldel)/2
    #     ldel = (wudel0+wldel0)/2.#(audel+aldel)/2
    # else:
    #     udel = (wudel0+wldel0+audel0+aldel0)/4.#(audel+aldel)/2
    #     ldel = (wudel0+wldel0+audel0+aldel0)/4.#(audel+aldel)/2
    # # print udel,ldel
    #
    # deleu = np.array([])
    # delel = np.array([])
    # delwu = np.array([])
    # delwl = np.array([])
    # delau = np.array([])
    # delal = np.array([])
    # # delcu = np.array([])
    # # delcl = np.array([])
    # for i in range(np.size(xeu)):
    #     if xeu[i] >= udel:
    #         deleu = np.append(deleu,i)
    # for i in range(np.size(xel)):
    #     if xel[i] >= ldel:
    #         delel = np.append(delel,i)
    # for i in range(np.size(xwu)):
    #     if xwu[i] >= udel:
    #         delwu = np.append(delwu,i)
    # for i in range(np.size(xwl)):
    #     if xwl[i] >= ldel:
    #         delwl = np.append(delwl,i)
    # for i in range(np.size(xau)):
    #     if xau[i] >= udel:
    #         delau = np.append(delau,i)
    # for i in range(np.size(xal)):
    #     if xal[i] >= ldel:
    #         delal = np.append(delal,i)
    # # for i in range(np.size(xcu)):
    # #     if xcu[i] >= udel:
    # #         delcu = np.append(delcu,i)
    # # for i in range(np.size(xcl)):
    # #     if xcl[i] >= ldel:
    # #         delcl = np.append(delcl,i)
    # deleu = deleu[::-1]
    # delel = delel[::-1]
    # delwu = delwu[::-1]
    # delwl = delwl[::-1]
    # delau = delau[::-1]
    # delal = delal[::-1]
    # # delcu = delcu[::-1]
    # # delcl = delcl[::-1]
    # for i in range(np.size(deleu)):
    #     eu = np.delete(eu,deleu[i])
    #     xeu = np.delete(xeu,deleu[i])
    # for i in range(np.size(delel)):
    #     el = np.delete(el,delel[i])
    #     xel = np.delete(xel,delel[i])
    # for i in range(np.size(delwu)):
    #     wu = np.delete(wu,delwu[i])
    #     xwu = np.delete(xwu,delwu[i])
    # for i in range(np.size(delwl)):
    #     wl = np.delete(wl,delwl[i])
    #     xwl = np.delete(xwl,delwl[i])
    # for i in range(np.size(delau)):
    #     au = np.delete(au,delau[i])
    #     xau = np.delete(xau,delau[i])
    # for i in range(np.size(delal)):
    #     al = np.delete(al,delal[i])
    #     xal = np.delete(xal,delal[i])
    # # for i in range(np.size(delcu)):
    # #     cu = np.delete(cu,delcu[i])
    # #     xcu = np.delete(xcu,delcu[i])
    # # for i in range(np.size(delcl)):
    # #     cl = np.delete(cl,delcl[i])
    # # xcl = np.delete(xcl,delcl[i])
    #
    #
    #
    #
    # if plott == True:
    #     diap = 1.
    #     plt.figure(4)
    #     plt.subplot(2,2,1)
    #     plt.plot(xeu/diap,eu/diap,'bo')
    #     plt.plot(xel/diap,el/diap,'ro')
    #     plt.title(r'Location ($\xi$)')
    #     # plt.ylim(-1,1)
    #
    #     plt.subplot(2,2,2)
    #     plt.plot(xwu/diap,wu,'bo')
    #     plt.plot(xwl/diap,wl,'ro')
    #     plt.title(r'Scale ($\omega$)')
    #     # plt.ylim(-0.05,0.15)
    #
    #     plt.subplot(2,2,3)
    #     plt.plot(xau/diap,au,'bo')
    #     plt.plot(xal/diap,al,'ro')
    #     plt.title(r'veldist ($\alpha$)')
    #     # plt.ylim(-5,5)
    #
    #     # plt.subplot(2,2,4)
    #     # plt.plot(xcu/diap,cu,'bo')
    #     # plt.plot(xcl/diap,cl,'ro')
    #     # plt.title('Resize ($c$)')
    #
    #     plt.figure(9)
    #     plt.plot(xintu/diap,intu,'bo')
    #     plt.plot(xintl/diap,intl,'ro')
    #     plt.title('Integral')
    #     plt.show()
    #
    #




## Curve Fitting        
    ft = 1e-8
    xt = 1e-8
    mf = int(1e6)
    
    # efit,_ = curve_fit(normgauss,xeu,eu)
    efit = np.polyfit(xeu,eu,2)
    # efit,_ = leastsq(residualnorm,np.array([3.,1.,0.5]),args=(xeu,eu),ftol=ft,xtol=xt,maxfev=mf)
    
    wfit = np.polyfit(xwu,wu,2)
    # wfit,_ = curve_fit(sigmoid,xwu,wu,p0=(2.,1.,0.75*xdu[-1]))
    # wfit,_ = curve_fit(normgaussw,xwu,wu,p0=(3.,1.,0.5))
    # wfit,_ = leastsq(residualnorm,np.array([3.,1.,0.5]),args=(xwu,wu),ftol=ft,xtol=xt,maxfev=mf)
    
    # dfit = np.polyfit(xdu,du,2)
    # dfit,_ = curve_fit(sigmoid,xdu,du,p0=(2.,1.,0.75*xdu[-1]))
    try:
        # dfit,_ = curve_fit(normgauss,xdu,du,p0=(3.,1.,0.5))
        dfit,_ = curve_fit(normgauss,xdu,du,p0=(0.05,1.,1.))
    except:
        dfit = np.array([10,10,10])
    # dfit,_ = leastsq(residualnorm,np.array([3.,1.,0.5]),args=(xdu,du),ftol=ft,xtol=xt,maxfev=mf)
    # #SNOPT#######
    # xdata = xdu
    # ydata = du
    # dfit0 = np.array([3.,3.,15.])
    # dfitl = np.array([1.,0.04*length,0.])
    # dfitu = np.array([length/dia,300.,25.])
    # optProb = Optimization('Vel_Fit', obj_funcng)
    # optProb.addObj('obj')
    # optProb.addVarGroup('param', 3, 'c', lower=dfitl, upper=dfitu, value=dfit0)
    # opt = SNOPT()
    # opt.setOption('Scale option',0)
    # res = opt(optProb, sens=None)
    # # print res
    # dfit = res.xStar['param']
    
    bfit = np.polyfit(xbu,bu,1)
    afit = np.polyfit(xau,au,1)
    
    xfite = np.linspace(xeu[0],xeu[-1],100)
    xfitw = np.linspace(xwu[0],xwu[-1],100)
    xfitd = np.linspace(xdu[0],xdu[-1],100)
    xfitb = np.linspace(xbu[0],xbu[-1],100)
    xfita = np.linspace(xau[0],xau[-1],100)
    
    efit = np.array( [0.00011143447157452207, -6.548322990064805e-05, 0.018529320355285157] )
    wfit = np.array( [-0.7652190255937096, 0.06107150050623949, 0.6982471030158889, 0.7652436757392653] )
    dfit = np.array( [0.46007275370246564, 0.11572807674332522, 20.54895467382088] )
    bfit = np.array( [-2.6889951822073024, 44.30363297368563] )
    afit = np.array( [-2.1215370317567306, 33.035244743684075] )
    
    efit = np.array( [0.0002809298385469616, -0.002152904786611757, 0.02421458254508763] )
    wfit = np.array( [-0.7619069186757337, 0.054654574491074574, 0.693209512634655, 0.7589702331063731] )
    dfit = np.array( [0.43439424298020607, 0.1208447537121723, 20.53828795693433] )
    bfit = np.array( [-2.4841293302971597, 44.32673189170779] )
    afit = np.array( [-1.9925689536368407, 33.01510013478536] )
    
    if plottrend == True:
        # plt.figure(10)
        # plt.subplot(5,1,1)
        # plt.plot(xeu,eu,'bo')
        # plt.title('Mean')
        # plt.ylim(-0.5,0.5)
        # plt.subplot(5,1,2)
        # plt.plot(xwu,wu,'bo')
        # plt.title('St Dev')
        # plt.ylim(0,10)
        # plt.subplot(5,1,3)
        # plt.plot(xdu,du,'bo')
        # plt.title('Scale')
        # plt.ylim(0,10)
        # plt.subplot(5,1,4)
        # plt.plot(xbu,bu,'bo')
        # plt.title('Rate')
        # plt.ylim(0,50)
        # plt.subplot(5,1,5)
        # plt.plot(xau,au,'bo')
        # plt.title('Trans')
        # plt.ylim(0,40)
        
        # men = np.array( [-0.0006344223751663201, 0.01055675755786011, -0.004073212523707764] )
        # spr = np.array( [-0.005187125854670714, 0.06397918461247416, 0.543874357807372] )
        # dfit = np.array( [6.667328694868336, 5.617498827673229, 21.520026361522778] )
        # rat = np.array( [-2.129054494312758, 45.17191461412915] )
        # tns = np.array( [-1.5569348878268718, 31.913143231782648] )
        
        if paper == True:
            plt.figure(11)
            # plt.subplot(3,2,1)
            # plt.plot(xfite,normgauss(xfite,efit[0],efit[1],efit[2]),'r')
            plt.plot(xfite,np.polyval(efit,xfite),'r',linewidth=2)
            plt.figure(12)
            # plt.subplot(3,2,2)
            # plt.plot(xfitw,np.polyval(wfit,xfitw),'r',linewidth=2)
            # plt.plot(xfitw,sigmoid(xfitw,wfit[0],wfit[1],wfit[2]),'r')
            plt.plot(xfitw,normgaussw(xfitw,wfit[0],wfit[1],wfit[2],wfit[3]),'r',linewidth=2)
            plt.figure(13)
            # plt.subplot(3,2,3)
            # plt.plot(xfitd,np.polyval(dfit,xfitd),'r')
            # plt.plot(xfitd,sigmoid(xfitd,dfit[0],dfit[1],dfit[2]),'r')
            plt.plot(xfitd,normgauss(xfitd,dfit[0],dfit[1],dfit[2]),'r',linewidth=2)
            plt.figure(14)
            # plt.subplot(3,2,4)
            plt.plot(xfitb,np.polyval(bfit,xfitb),'r',linewidth=2)
            plt.figure(15)
            # plt.subplot(3,2,5)
            plt.plot(xfita,np.polyval(afit,xfita),'r',linewidth=2)
        elif paper == False:
            plt.figure(10)
            plt.subplot(3,2,1)
            # plt.plot(xfite,normgauss(xfite,efit[0],efit[1],efit[2]),'r')
            plt.plot(xfite,np.polyval(efit,xfite),'r',linewidth=2)
            plt.subplot(3,2,2)
            # plt.plot(xfitw,np.polyval(wfit,xfitw),'r',linewidth=2)
            # plt.plot(xfitw,sigmoid(xfitw,wfit[0],wfit[1],wfit[2]),'r')
            plt.plot(xfitw,normgaussw(xfitw,wfit[0],wfit[1],wfit[2],wfit[3]),'r',linewidth=2)
            plt.subplot(3,2,3)
            # plt.plot(xfitd,np.polyval(dfit,xfitd),'r')
            # plt.plot(xfitd,sigmoid(xfitd,dfit[0],dfit[1],dfit[2]),'r')
            plt.plot(xfitd,normgauss(xfitd,dfit[0],dfit[1],dfit[2]),'r',linewidth=2)
            plt.subplot(3,2,4)
            plt.plot(xfitb,np.polyval(bfit,xfitb),'r',linewidth=2)
            plt.subplot(3,2,5)
            plt.plot(xfita,np.polyval(afit,xfita),'r',linewidth=2)
    
    # print wfit,dfit,bfit,afit
    if plottrend == True:
        print '\n'
        print efit[0]
        print efit[1]
        print efit[2]
        print wfit[0]
        print wfit[1]
        print wfit[2]
        print dfit[0]
        print dfit[1]
        print dfit[2]
        print bfit[0]
        print bfit[1]
        print afit[0]
        print afit[1]
        print '\n'
    
    plt.show()


    return x,efit,wfit,dfit,bfit,afit#,xeu,xwu,xau,xcu,eu,wu,au,cu

## Main File
if __name__ == "__main__":
    cv = False

    s = 's2'
    t = '400'
    length = 100.

    dia = 6.

    x,efit,wfit,dfit,bfit,afit = fit(s,t,length,True,True)

## Output parameters
    print '\n'
    print 'men = np.array(',efit.tolist(),')'
    print 'spr = np.array(',wfit.tolist(),')'
    print 'scl = np.array(',dfit.tolist(),')'
    print 'rat = np.array(',bfit.tolist(),')'
    print 'tns = np.array(',afit.tolist(),')'
    
    
    # # print '\nx = np.array(',x.tolist(),')'
    # # print '\neu = np.array(',eu.tolist(),')'
    # # print 'wu = np.array(',wu.tolist(),')'
    # # print 'au = np.array(',au.tolist(),')'
    # # print 'cu = np.array(',cu.tolist(),')'
    # # print 'el = np.array(',el.tolist(),')'
    # # print 'wl = np.array(',wl.tolist(),')'
    # # print 'al = np.array(',al.tolist(),')'
    # # print 'cl = np.array(',cl.tolist(),')\n'
    # # print 'ceu = np.array(',ceu.tolist(),')'
    # # print 'cwu = np.array(',cwu.tolist(),')'
    # print 'cau = np.array(',cau.tolist(),')'
    # # print 'ccu = np.array(',ccu.tolist(),')'
    # # print 'ced = np.array(',cel.tolist(),')'
    # # print 'cwd = np.array(',cwl.tolist(),')'
    # print 'cad = np.array(',cal.tolist(),')'
    # # print 'ccd = np.array(',ccl.tolist(),')'

    # Order 3 (2 for location)
    # print '\neu = '+str(ceu[0])+'_dp*sd**2 + '+str(ceu[1])+'_dp*sd + '+str(ceu[2])+'_dp'
    # print 'wu = '+str(cwu[0])+'_dp*exp('+str(cwu[1])+'_dp*sd) + '+str(cwu[2])+'_dp'
    # print 'au = '+str(cau[0])+'_dp*sd**3 + '+str(cau[1])+'_dp*sd**2 + '+str(cau[2])+'_dp*sd + '+str(cau[3])+'_dp'
    # # print 'cu = '+str(ccu[0])+'_dp*sd**3 + '+str(ccu[1])+'_dp*sd**2 + '+str(ccu[2])+'_dp*sd + '+str(ccu[3])+'_dp'
    # print 'ed = '+str(cel[0])+'_dp*sd**2 + '+str(cel[1])+'_dp*sd + '+str(cel[2])+'_dp'
    # print 'wd = '+str(cwl[0])+'_dp*exp('+str(cwl[1])+'_dp*sd) + '+str(cwl[2])+'_dp'
    # print 'ad = '+str(cal[0])+'_dp*sd**3 + '+str(cal[1])+'_dp*sd**2 + '+str(cal[2])+'_dp*sd + '+str(cal[3])+'_dp'
    # # print 'cd = '+str(ccl[0])+'_dp*sd**3 + '+str(ccl[1])+'_dp*sd**2 + '+str(ccl[2])+'_dp*sd + '+str(ccl[3])+'_dp'
    #
    # # Order 4 (2 for location)
    # # print '\neu = '+str(ceu[0])+'_dp*sd**2 + '+str(ceu[1])+'_dp*sd + '+str(ceu[2])+'_dp'
    # # print 'wu = '+str(cwu[0])+'_dp*exp('+str(cwu[1])+'_dp*sd) + '+str(cwu[2])+'_dp'
    # # print 'au = '+str(cau[0])+'_dp*sd**4 + '+str(cau[1])+'_dp*sd**3 + '+str(cau[2])+'_dp*sd**2 &\n+ '+str(cau[3])+'_dp*sd + '+str(cau[4])+'_dp'
    # print 'cu = '+str(ccu[0])+'_dp*sd**4 + '+str(ccu[1])+'_dp*sd**3 + '+str(ccu[2])+'_dp*sd**2 &\n+ '+str(ccu[3])+'_dp*sd + '+str(ccu[4])+'_dp'
    # # print 'ed = '+str(cel[0])+'_dp*sd**2 + '+str(cel[1])+'_dp*sd + '+str(cel[2])+'_dp'
    # # print 'wd = '+str(cwl[0])+'_dp*exp('+str(cwl[1])+'_dp*sd) + '+str(cwl[2])+'_dp'
    # # print 'ad = '+str(cal[0])+'_dp*sd**4 + '+str(cal[1])+'_dp*sd**3 + '+str(cal[2])+'_dp*sd**2 &\n+ '+str(cal[3])+'_dp*sd + '+str(cal[4])+'_dp'
    # print 'cd = '+str(ccl[0])+'_dp*sd**4 + '+str(ccl[1])+'_dp*sd**3 + '+str(ccl[2])+'_dp*sd**2 &\n+ '+str(ccl[3])+'_dp*sd + '+str(ccl[4])+'_dp'

    # print '\n*****************AVERAGE VALUES*****************\n'

    # print 'e = np.array([',e1,',',e2,',',e3,'])'
    # print 'w = np.array([',w1,',',w2,',',w3,'])'
    # print 'a = np.array([',a1,',',a2,',',a3,',',a4,'])'
    # # print 'c = np.array([',c1,',',c2,',',c3,'])'
    # # print 'aed = np.array([',-e1,',',-e2,',',-e3,'])'
    # # print 'awd = np.array([',w1,',',w2,',',w3,'])'
    # # print 'aau = np.array([',-a1,',',-a2,',',-a3,',',-a4,'])'
    # # print 'acd = np.array([',c1,',',c2,',',c3,',',c4,',',c5,'])'
    # print 'i = np.array([',i1,',',i2,',',i3,'])'
    # # print 'aid = np.array([',-i1,',',i2,',',i3,'])'

    # Order 3 (2 for location)
    # print '\neu = '+str(e1)+'_dp*sd**2 + '+str(e2)+'_dp*sd + '+str(e3)+'_dp'
    # print 'wu = '+str(w1)+'_dp*exp('+str(w2)+'_dp*sd) + '+str(w3)+'_dp'
    # print 'au = '+str(a1)+'_dp*sd**3 + '+str(a2)+'_dp*sd**2 + '+str(a3)+'_dp*sd + '+str(a4)+'_dp'
    # # print 'cu = '+str(c1)+'_dp*sd**3 + '+str(c2)+'_dp*sd**2 + '+str(c3)+'_dp*sd + '+str(c4)+'_dp'
    # print 'ed = '+str(-e1)+'_dp*sd**2 + '+str(-e2)+'_dp*sd + '+str(-e3)+'_dp'
    # print 'wd = '+str(w1)+'_dp*exp('+str(w2)+'_dp*sd) + '+str(w3)+'_dp'
    # print 'ad = '+str(-a1)+'_dp*sd**3 + '+str(-a2)+'_dp*sd**2 + '+str(-a3)+'_dp*sd + '+str(-a4)+'_dp'
    # # print 'cd = '+str(c1)+'_dp*sd**3 + '+str(c2)+'_dp*sd**2 + '+str(c3)+'_dp*sd + '+str(c4)+'_dp'
    #
    # # Order 4 (2 for location)
    # # print '\neu = '+str(e1)+'_dp*sd**2 + '+str(e2)+'_dp*sd + '+str(e3)+'_dp'
    # # print 'wu = '+str(w1)+'_dp*exp('+str(w2)+'_dp*sd) + '+str(w3)+'_dp'
    # # print 'au = '+str(a1)+'_dp*sd**4 + '+str(a2)+'_dp*sd**3 + '+str(a3)+'_dp*sd**2 &\n+ '+str(a4)+'_dp*sd + '+str(a5)+'_dp'
    # print 'cu = '+str(c1)+'_dp*sd**4 + '+str(c2)+'_dp*sd**3 + '+str(c3)+'_dp*sd**2 &\n+ '+str(c4)+'_dp*sd + '+str(c5)+'_dp'
    # # print 'ed = '+str(-e1)+'_dp*sd**2 + '+str(-e2)+'_dp*sd + '+str(-e3)+'_dp'
    # # print 'wd = '+str(w1)+'_dp*exp('+str(w2)+'_dp*sd) + '+str(w3)+'_dp'
    # # print 'ad = '+str(-a1)+'_dp*sd**4 + '+str(-a2)+'_dp*sd**3 + '+str(-a3)+'_dp*sd**2 &\n+ '+str(-a4)+'_dp*sd + '+str(-a5)+'_dp'
    # print 'cd = '+str(c1)+'_dp*sd**4 + '+str(c2)+'_dp*sd**3 + '+str(c3)+'_dp*sd**2 &\n+ '+str(c4)+'_dp*sd + '+str(c5)+'_dp'

    # # Plotting parameter trends
    # xfit = np.linspace(0.,max(x),500)
    # feu = np.polyval(ceu,xfit)
    # fwu = expfit(xfit,cwu[0],cwu[1],cwu[2])
    # fau = np.polyval(cau,xfit)
    # # fcu = np.polyval(ccu,xfit)
    # fcu = piecewise(xfit,ccu[0],ccu[1],ccu[2],ccu[3],ccu[4],ccu[5],ccu[6])
    # fel = np.polyval(cel,xfit)
    # fwl = expfit(xfit,cwl[0],cwl[1],cwl[2])
    # fal = np.polyval(cal,xfit)
    # # fcl = np.polyval(ccl,xfit)
    # fcl = piecewise(xfit,ccl[0],ccl[1],ccl[2],ccl[3],ccl[4],ccl[5],ccl[6])
    #
    # geu = np.polyval([e1,e2,e3],xfit)
    # gwu = expfit(xfit,w1,w2,w3)
    # gau = np.polyval([a1,a2,a3,a4],xfit)
    # # gcu = np.polyval([c1,c2,c3,c4,c5],xfit)
    # gcu = piecewise(xfit,c1,c2,c3,c4,c5,c6,c7)
    # gel = np.polyval([-e1,-e2,-e3],xfit)
    # gwl = expfit(xfit,w1,w2,w3)
    # gal = np.polyval([-a1,-a2,-a3,-a4],xfit)
    # # gcl = np.polyval([c1,c2,c3,c4,c5],xfit)
    # gcl = piecewise(xfit,c1,c2,c3,c4,c5,c6,c7)
    #
    # fs = 12
    # diap = 1.
    #
    # plt.figure(3)
    # plt.subplot(2,2,1)
    # plt.plot(xeu/diap,eu/diap,'bo')
    # plt.plot(xel/diap,el/diap,'ro')
    # plt.plot(xfit/diap,feu/diap,'c-')
    # plt.plot(xfit/diap,fel/diap,'m-')
    # plt.plot(xfit/diap,geu/diap,'g-')
    # plt.plot(xfit/diap,gel/diap,'g-')
    # plt.title(r'Location ($\xi$)')
    # # plt.ylim(-1,1)
    #
    # plt.subplot(2,2,2)
    # plt.plot(xwu/diap,wu,'bo')
    # plt.plot(xwl/diap,wl,'ro')
    # plt.plot(xfit/diap,fwu,'c-')
    # plt.plot(xfit/diap,fwl,'m-')
    # plt.plot(xfit/diap,gwu,'g-')
    # plt.plot(xfit/diap,gwl,'g-')
    # plt.title(r'Spread ($\omega$)')
    # # plt.ylim(-0.05,0.15)
    #
    # plt.subplot(2,2,3)
    # plt.plot(xau/diap,au,'bo')
    # plt.plot(xal/diap,al,'ro')
    # plt.plot(xfit/diap,fau,'c-')
    # plt.plot(xfit/diap,fal,'m-')
    # plt.plot(xfit/diap,gau,'g-')
    # plt.plot(xfit/diap,gal,'g-')
    # plt.title(r'veldist ($\alpha$)')
    # # plt.ylim(-5,5)
    #
    # plt.subplot(2,2,4)
    # plt.plot(xcu/diap,cu,'bo')
    # plt.plot(xcl/diap,cl,'ro')
    # plt.plot(xfit/diap,fcu,'c-')
    # plt.plot(xfit/diap,fcl,'m-')
    # plt.plot(xfit/diap,gcu,'g-')
    # plt.plot(xfit/diap,gcl,'g-')
    # plt.title('Scale ($c$)')


## Cross Validation (K Fold)
    if cv == True:
        from sklearn.cross_validation import train_test_split
        def compute_error(x,y,p):
            yfit = np.polyval(p,x)
            return sqrt(np.mean((y-yfit)**2))

        test = 0.1

        euxtrain,euxtest,euytrain,euytest = train_test_split(xeu,eu,test_size = test)
        edxtrain,edxtest,edytrain,edytest = train_test_split(xel,el,test_size = test)
        wuxtrain,wuxtest,wuytrain,wuytest = train_test_split(xwu,wu,test_size = test)
        wdxtrain,wdxtest,wdytrain,wdytest = train_test_split(xwl,wl,test_size = test)
        auxtrain,auxtest,auytrain,auytest = train_test_split(xau,au,test_size = test)
        adxtrain,adxtest,adytrain,adytest = train_test_split(xal,al,test_size = test)
        # cuxtrain,cuxtest,cuytrain,cuytest = train_test_split(xcu,cu,test_size = test)
        # cdxtrain,cdxtest,cdytrain,cdytest = train_test_split(xcl,cl,test_size = test)

        degrees = np.arange(8)
        eutrain_error = np.zeros(len(degrees))
        euvalidation_error = np.zeros(len(degrees))
        edtrain_error = np.zeros(len(degrees))
        edvalidation_error = np.zeros(len(degrees))
        wutrain_error = np.zeros(len(degrees))
        wuvalidation_error = np.zeros(len(degrees))
        wdtrain_error = np.zeros(len(degrees))
        wdvalidation_error = np.zeros(len(degrees))
        autrain_error = np.zeros(len(degrees))
        auvalidation_error = np.zeros(len(degrees))
        adtrain_error = np.zeros(len(degrees))
        advalidation_error = np.zeros(len(degrees))
        # cutrain_error = np.zeros(len(degrees))
        # cuvalidation_error = np.zeros(len(degrees))
        # cdtrain_error = np.zeros(len(degrees))
        # cdvalidation_error = np.zeros(len(degrees))


        for i,d in enumerate(degrees):
            peu = np.polyfit(euxtrain,euytrain,d)
            eutrain_error[i] = compute_error(euxtrain,euytrain,peu)
            euvalidation_error[i] = compute_error(euxtest,euytest,peu)
            ped = np.polyfit(edxtrain,edytrain,d)
            edtrain_error[i] = compute_error(edxtrain,edytrain,ped)
            edvalidation_error[i] = compute_error(edxtest,edytest,ped)
            pwu = np.polyfit(wuxtrain,wuytrain,d)
            wutrain_error[i] = compute_error(wuxtrain,wuytrain,pwu)
            wuvalidation_error[i] = compute_error(wuxtest,wuytest,pwu)
            pwd = np.polyfit(wdxtrain,wdytrain,d)
            wdtrain_error[i] = compute_error(wdxtrain,wdytrain,pwd)
            wdvalidation_error[i] = compute_error(wdxtest,wdytest,pwd)
            pau = np.polyfit(auxtrain,auytrain,d)
            autrain_error[i] = compute_error(auxtrain,auytrain,pau)
            auvalidation_error[i] = compute_error(auxtest,auytest,pau)
            pad = np.polyfit(adxtrain,adytrain,d)
            adtrain_error[i] = compute_error(adxtrain,adytrain,pad)
            advalidation_error[i] = compute_error(adxtest,adytest,pad)
            # pcu = np.polyfit(cuxtrain,cuytrain,d)
            # cutrain_error[i] = compute_error(cuxtrain,cuytrain,pcu)
            # cuvalidation_error[i] = compute_error(cuxtest,cuytest,pcu)
            # pcd = np.polyfit(cdxtrain,cdytrain,d)
            # cdtrain_error[i] = compute_error(cdxtrain,cdytrain,pcd)
            # cdvalidation_error[i] = compute_error(cdxtest,cdytest,pcd)

        plt.figure(4)
        plt.subplot(2,4,1)
        plt.plot(degrees,eutrain_error,'b',label='Training Error')
        plt.plot(degrees,euvalidation_error,'r',label='Cross Validation Error')
        plt.xlabel('Order of the Polynomial Fit')
        plt.ylabel('Error Value')
        plt.title('Location (up)')

        plt.subplot(2,4,2)
        plt.plot(degrees,edtrain_error,'b',label='Training Error')
        plt.plot(degrees,edvalidation_error,'r',label='Cross Validation Error')
        plt.xlabel('Order of the Polynomial Fit')
        plt.ylabel('Error Value')
        plt.title('Location (down)')


        plt.subplot(2,4,3)
        plt.plot(degrees,wutrain_error,'b',label='Training Error')
        plt.plot(degrees,wuvalidation_error,'r',label='Cross Validation Error')
        plt.xlabel('Order of the Polynomial Fit')
        plt.ylabel('Error Value')
        plt.title('Spread (up)')

        plt.subplot(2,4,4)
        plt.plot(degrees,wdtrain_error,'b',label='Training Error')
        plt.plot(degrees,wdvalidation_error,'r',label='Cross Validation Error')
        plt.xlabel('Order of the Polynomial Fit')
        plt.ylabel('Error Value')
        plt.title('Spread (down)')
        plt.legend(loc="upper left", bbox_to_anchor=(1,1))

        plt.subplot(2,4,5)
        plt.plot(degrees,autrain_error,'b',label='Training Error')
        plt.plot(degrees,auvalidation_error,'r',label='Cross Validation Error')
        plt.xlabel('Order of the Polynomial Fit')
        plt.ylabel('Error Value')
        plt.title('veldist (up)')

        plt.subplot(2,4,6)
        plt.plot(degrees,adtrain_error,'b',label='Training Error')
        plt.plot(degrees,advalidation_error,'r',label='Cross Validation Error')
        plt.xlabel('Order of the Polynomial Fit')
        plt.ylabel('Error Value')
        plt.title('veldist (down)')

        # plt.subplot(2,4,7)
        # plt.plot(degrees,cutrain_error,'b',label='Training Error')
        # plt.plot(degrees,cuvalidation_error,'r',label='Cross Validation Error')
        # plt.xlabel('Order of the Polynomial Fit')
        # plt.ylabel('Error Value')
        # plt.title('Scale (up)')
        #
        # plt.subplot(2,4,8)
        # plt.plot(degrees,cdtrain_error,'b',label='Training Error')
        # plt.plot(degrees,cdvalidation_error,'r',label='Cross Validation Error')
        # plt.xlabel('Order of the Polynomial Fit')
        # plt.ylabel('Error Value')
        # plt.title('Scale (down)')


    plt.show()