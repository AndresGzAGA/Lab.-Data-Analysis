# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math as m

def data_matrix(file_name):
    """
    From a file_name which contains a table of numbers separated with tabulator,
    returns a matrix with those numbers, not ignoring the first row. 
    """
    infile=open(file_name,"r")
    #infile.readline() #skip the first line
    matrix=[]
    for line in infile:
        num_row=[]
        str_row=line.split()
        for n in str_row:
            #n=float(n)
            num_row.append(n)
        matrix.append(num_row)
    return matrix
    infile.close()

def data_dictionary(file_name,location_list):
    """
    Creates a dictionary from file_name, using location_list with the syntax:
    location_list=[first term(list with: term name(str), first column(int), last column(int)), second term(...),...]
    location_list does not includes the keys, they should always be in the first column.
    Example of what it returns: {5.0: {'first': [7.99, 8.42, 7.96],'second': [8.11, 7.87, 7.93],'third': [7.9]},7.0: {'first': [15.89, 15.76, 15.7], 'second': [15.7, 16.1, 16.01], 'third': [15.76]},9.0: {'first': [21.73, 21.38, 21.0], 'second': [21.37, 21.65, 21.39], 'third': [21.84]},11.0: {'first': [27.08, 27.32, 27.14], 'second': [27.4, 27.52, 27.36], 'third': [27.21]},13.0: {'first': [32.24, 31.92, 31.85], 'second': [31.88, 32.4, 32.57], 'third': [31.95]},15.0: {'first': [37.06, 37.19, 36.65], 'second': [36.51, 36.0, 36.04], 'third': [36.46]},17.0: {'first': [41.32, 40.26, 41.6], 'second': [41.69, 40.88, 41.39], 'third': [41.64]}}
    """
    infile=open(file_name,"r")
    infile.readline() #skip the first line
    dic={}
    for line in infile:
        num_row={} #Dictionary for each line, which will go in dic with the main key.
        str_row=line.split() #Raw list with all the numbers in the working line.
        for term in location_list:
            term_list=[] #List for each different term that will go in num_row.
            for i in range(term[1],term[2]+1):
                n=float(str_row[i])
                term_list.append(n)
            num_row[term[0]]=term_list
        n=float(str_row[0]) #Main key
        dic[n]=num_row
    return dic
    
def write(file_name,something,name):
    """
    Writes something (string) on the last column with the name of the column.
    """
    D=data_matrix(file_name)
    D[0].append(name)
    for i in range(1,len(D)):
        D[i].append(something)
    outfile=open(file_name,"w")
    for i in D:
        j=i[0]
        outfile.write(j)
        for j in i[1:]:
            outfile.write("\t" + j)
        outfile.write('\n')
    outfile.close()
    
    
def plot_points(file_name,ci=1,cf=1):
    """
    Plots all of the data from the column ci to the column cf, against the first key column.
    """
    location_list=[["y",ci,cf]]
    d=data_dictionary(file_name,location_list)
    for i in range(0,cf-ci+1):
        yi=[]
        for k in d.keys():
            yi.append(d[k]["y"][i])
        plt.plot(d.keys(),yi,'b+')
        plt.hold('on')

def avrg(file_name,ci=1,cf=1):
    """
    Computes the average of the data from column ci to column cf,
    returns a dictionary with the keys from the first column.
    """
    location_list=[["y",ci,cf]]
    d=data_dictionary(file_name,location_list)
    av={}
    for k in d.keys():
        a=0
        for i in range(0,cf-ci+1):
            a+=d[k]["y"][i]
        a=a/(cf-ci+1)
        av[k]=a
    return av
    
def multi_avrg(file_names,ci_l,cf_l):
    """
    file_names: List with the file names.
    ci_l: List with the initial columns
    cf_l: List with the final columns.
    Computes the average of the data of the file_names, from their respective
    column ci to column cf, returns a dictionary with the keys from the first column. 
    """
    av={}
    N=0
    for i in range(0,len(file_names)):
        file_name=file_names[i]
        ci=ci_l[i]
        cf=cf_l[i]
        N+=cf-ci+1
        location_list=[["y",ci,cf]]
        d=data_dictionary(file_name,location_list)
        for k in d.keys():
            a=0
            for j in range(0,cf-ci+1):
                a+=d[k]["y"][j]
            if i==0:
                av[k]=a
            else:
                av[k]+=a
    for k in d.keys():
        av[k]=av[k]/N
    return av
       
def std(file_name,ci=1,cf=1):
    """
    Calculates the sample standard deviation of the data of each row from ci to cf.
    Returns a dictionary with the keys from the first column.
    """
    av=avrg(file_name,ci,cf)
    location_list=[["y",ci,cf]]
    d=data_dictionary(file_name,location_list)
    dev={}
    for k in d.keys():
        devk=0
        for i in range(0,cf-ci+1):
            devk+=(d[k]["y"][i]-av[k])**2
        devk=devk/(cf-ci)
        devk=m.sqrt(devk)
        dev[k]=devk
    return dev
    
def multi_std(file_names,ci_l,cf_l):
    """
    file_names: List with the file names.
    ci_l: List with the initial columns
    cf_l: List with the final columns.
    Computes the sample standard deviation of the data of the file_names,
    from their respective column ci to column cf,
    returns a dictionary with the keys from the first column.
    """
    av=multi_avrg(file_names,ci_l,cf_l)
    dev={}
    N=0
    for i in range(0,len(file_names)):
        file_name=file_names[i]
        ci=ci_l[i]
        cf=cf_l[i]
        N+=cf-ci+1
        location_list=[["y",ci,cf]]
        d=data_dictionary(file_name,location_list)
        for k in d.keys():
            devk=0
            for j in range(0,cf-ci+1):
                devk+=(d[k]["y"][j]-av[k])**2
            if i==0:
                dev[k]=devk
            else:
                dev[k]+=devk
    for k in d.keys():
        dev[k]=dev[k]/(N-1)
        dev[k]=m.sqrt(dev[k])
    return dev

    
def dist_plot(file_name,ci=1,cf=1,xerr=2,yerr_type="Standard deviation",yerr=3):
    """
    Plots the average of the distribution from ci of cf with its
    asociated standard deviation if err_type="Standard deviation";
    or, if err_type="Measurement" with its mesurement error in the err column,
    plotting against the keys from column 0 with their measurement error.
    """
    avd=avrg(file_name,ci,cf)
    avl=[]
    yerrl=[]
    xerrl=[]
    if yerr_type=="Standard deviation":
        yerrd=std(file_name,ci,cf)
        xerrd=data_dictionary(file_name,[["xerr",xerr,xerr]])
        yerrdkeys=yerrd.keys()
        yerrdkeys.sort()
        for k in yerrdkeys:
            avl.append(avd[k])
            yerrl.append(yerrd[k])
            xerrl.append(xerrd[k]["xerr"][0])
    elif yerr_type=="Measurement":
        errd=data_dictionary(file_name,[["xerr",xerr,xerr],["yerr",yerr,yerr]])
        errdkeys=errd.keys()
        errdkeys.sort()
        for k in errdkeys:
            avl.append(avd[k])
            yerrl.append(errd[k]["yerr"][0]/m.sqrt(cf-ci+1))
            xerrl.append(errd[k]["xerr"][0])

    x=avd.keys()
    x.sort()
    plt.errorbar(x,avl,xerr=xerrl,yerr=yerrl,ls=':',label='Data')
    plt.show()

def err_a(o, do):
    return np.sqrt(978**2*(np.cos(o)**2)*do**2/4) #NOTA: Checar 9.78 en cm o m, o (theta) tiene que estar en radianes!!!

def modelx(x, a):
    return a * x**2

def err_modelx(x, dx, a, da):
    return np.sqrt(4*a**2*x**2*dx**2+x**4*da**2)

def r(o):
    return 29.7+12.2*np.cos(o)

def dr(o,do):
    return np.sqrt(0.05**2+np.cos(o)**2*0.05**2+12.2**2*np.sin(o)**2*do**2)

def X(Ro,o):
    return 125.2+(Ro-15.3)*np.tan(o)+r(o)/np.cos(o)
    
def err_X(o, do, D, dD, dr, dP, dRo):
    return np.sqrt(dP**2+(dRo**2+dD**2)*(np.tan(o)**2)+dr**2*(np.cos(o)**(-2))+do**2*((do-D)*(np.cos(o)**(-2))+r(o)/np.cos(o)*np.tan(o))**2)

def p_A(Ro,o):
    return 978*m.sin(o)/(2*X(Ro,o))

def p_dA(Ro,o,do,D,dD,dr,dP,dRo):
    return np.sqrt(err_a(o,do)**2/(X(Ro,o)**2)+(978*m.sin(o)/2)**2*err_X(o,do,D,dD,dr,dP,dRo)**2/(X(Ro,o)**4))

def modelR(x, Ro, A):
    return 15.3/(1-A*np.power(x,2))+Ro-15.3
    
#def err_modelR(x, dx, Ro, A, dRo, dA):
#    return np.sqrt(dRo**2/np.power((1-A*np.power(x,2)),2)+(Ro**2*np.power(x,4)*dA**2+4*Ro**2*A**2*np.power(x,2)*dx**2)/np.power(1-A*np.power(x,2),4))

def err_modelR(x, dx, Ro, A, dRo, dA):
    D=20.3
    return np.sqrt(A**2*np.power(x,4)*0.05**2/np.power((1-A*np.power(x,2)),2)+dRo**2+(D**2*np.power(x,4)*dA**2+4*A**2*D**2*np.power(x,2)*dx**2)/np.power(1-A*np.power(x,2),4))

def fitR(file_name, Po, Po_err, ci=1, cf=1, xerr=2, yerr_type="Measurement", yerr=3):
    """
    Fits the function func to the average of the data from ci to cf, 
    with the associated error for each row if err_type is "Measurement";
    if err_type is "Standard deviation", then it calculates the standard
    deviation of the data and uses it as the error.
    Requires a pre-defined func.
    Po: List of the initial guess of the parameters.
    Po_err: List of the respective measurement error of the parameters.
    """
    dic_ydata=avrg(file_name,ci,cf)
    ydata=[]
    
    dicydatakeys=dic_ydata.keys()
    dicydatakeys.sort()
    for k in dicydatakeys:
        ydata.append(dic_ydata[k])
    xdata=dicydatakeys
    sigma=[]
    if yerr_type=="Standard deviation":
        dic_sigma=std(file_name,ci,cf)
        dicsigmakeys=dic_sigma.keys()
        dicsigmakeys.sort()
        for k in dicsigmakeys:
            sigma.append(dic_sigma[k])
    elif yerr_type=="Measurement":
        dic_sigma=data_dictionary(file_name,[["yerr",yerr,yerr]])
        dicsigmakeys=dic_sigma.keys()
        dicsigmakeys.sort()
        for k in dicsigmakeys:
            sigma.append(dic_sigma[k]["yerr"][0])
    else:
        print "Invalid err_type"
    
    dist_plot(file_name,ci,cf,xerr,yerr_type,yerr)
    
    xdataE=[]
    E=(xdata[-1]-xdata[0])/len(xdata)
    xdataE.append(xdata[0]-E)
    for i in xdata:
        xdataE.append(i)
    xdataE.append(xdata[-1]+E)
    
    popt, pcov = curve_fit(modelR, xdata, ydata, Po, sigma)
    plt.plot(xdataE, modelR(xdataE, *popt), 'r-', label='fit')
    plt.fill_between(xdataE,modelR(xdataE, *popt)-err_modelR(xdataE, 0.016, popt[0], popt[1], pcov[0][0], pcov[1][1]),modelR(xdataE, *popt)+err_modelR(xdataE, 0.016, popt[0], popt[1], pcov[0][0], pcov[1][1]),alpha=0.5,edgecolor='#CC4F1B', facecolor='#FF9848')
    
    plt.plot(xdataE, modelR(xdataE, *Po), 'g-', label='Model')
    plt.fill_between(xdataE, modelR(xdataE,*Po)-err_modelR(xdataE,0.016,Po[0],Po[1],Po_err[0],Po_err[1]), modelR(xdataE, *Po)+err_modelR(xdataE,0.016,Po[0],Po[1],Po_err[0],Po_err[1]),alpha=0.5, edgecolor='#3F7F4C', facecolor='#7EFF99',linewidth=0)

    plt.xlabel('Tiempo [s]')
    plt.ylabel('R [cm]')
    plt.legend()
    plt.title('R(t) - Vista en perspectiva a 21 grados')
    plt.show()
    
    print popt
    print pcov


def P_file(file_name):
    """
    Escribe en file_name todos los parámetros esperados.
    """
    alphal=["alpha[deg]"]
    for a in range(3,24,3):
        alphal.append(a)
    
    dal=["da[deg]"]
    for a in range(3,24,3):
        dal.append(1.5)
    
    thetal=["theta[rad]"]
    for a in alphal[1:]:
        thetal.append(np.deg2rad(a))
    
    dol=["do[rad]"]
    for a in dal[1:]:
        dol.append(a*3.141593/180)
    
    r=["r[cm]"]
    for i in range(1,len(thetal)):
        r.append(29.7+12.2*m.cos(thetal[i]))
        
    dr=["dr[cm]"]
    for i in range(1,len(thetal)):
        dr.append(m.sqrt(0.05**2+m.cos(thetal[i])**2*0.05**2+12.2**2*m.sin(thetal[i])**2*dol[i]**2))
    
    Ro=["Ro[cm]",8.48,-9.63,-0.89,-3.88,1.44,-1.55,-2.93]
    
    dRo=["dRo[cm]",0.125,0.125,0.125,0.125,0.125,0.125,0.125]
    
    Xl=["X[cm]"]
    for i in range(1,len(thetal)):
        Xl.append(X(Ro[i],thetal[i]))
        
    dX=["dX[cm]"]
    for i in range(1,len(thetal)):
        dX.append(err_X(thetal[i], dol[i], 15.3, r[i], dr[i], 0.05, dRo[i]))
    
    A=["A[s^-2]"]
    for i in range(1,len(thetal)):
        A.append(p_A(Ro[i],thetal[i]))
    
    dA=["dA[s^-2]"]
    for i in range(1,len(thetal)):
        dA.append(p_dA(Ro[i],thetal[i],dol[i],15.3,0.05,dr[i],0.05,0.125))
    
    D=[]
    D.append([alphal[0],dal[0],thetal[0],dol[0],r[0],dr[0],Ro[0],dRo[0],A[0],dA[0],Xl[0],dX[0]])
    for k in range(1,len(thetal)):
        row=[]
        row=[str(alphal[k]),str(dal[k]),str(round(thetal[k],3)),str(round(dol[k],4)),str(round(r[k],2)),str(round(dr[k],2)),str(round(Ro[k],2)),str(round(dRo[k],2)),str(round(A[k],2)),str(round(dA[k],2)),str(round(Xl[k],2)),str(round(dX[k],2))]
        D.append(row)
    
    outfile=open(file_name,"w")
    for i in D:
        j=i[0]
        outfile.write(j)
        for j in i[1:]:
            outfile.write("\t" + j)
        outfile.write('\n')
    outfile.close()