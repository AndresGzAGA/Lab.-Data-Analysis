import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import sqrt

def data_matrix(file_name):
    """
    From a file_name which contains a table of numbers separated with tabulator,
    returns a matrix with those numbers, ignoring the first row. 
    """
    infile=open(file_name,"r")
    infile.readline() #skip the first line
    matrix=[]
    for line in infile:
        num_row=[]
        str_row=line.split()
        for n in str_row:
            n=float(n)
            num_row.append(n)
        matrix.append(num_row)
    return matrix

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
        devk=sqrt(devk)
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
        dev[k]=sqrt(dev[k])
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
        for k in yerrd.keys():
            avl.append(avd[k])
            yerrl.append(yerrd[k])
            xerrl.append(xerrd[k]["xerr"][0])
    elif yerr_type=="Measurement":
        errd=data_dictionary(file_name,[["xerr",xerr,xerr],["yerr",yerr,yerr]])
        for k in errd.keys():
            avl.append(avd[k])
            yerrl.append(errd[k]["yerr"][0]/sqrt(cf-ci+1))
            xerrl.append(errd[k]["xerr"][0])

    plt.errorbar(avd.keys(),avl,xerr=xerrl,yerr=yerrl,ls=':',label='Data')
    plt.show()

def err_a(r1,r2,dr):
    return np.sqrt(8*r1**2*dr**2/(981*r2**4)+8*r1**4*dr**2/(981*r2**6))
    
def err_b(r1,r2,ho,dr,dh):
    return np.sqrt((8*r1**2*ho/(981*r2**4)+8*r1**4*ho/(981*r2**6))*dr**2+r1**4*dh**2/(2*981*ho*r2**4))

def model(x, a, b):
    return a * np.sqrt(x) + b

def err_model(x,a,da,db,dh):
    return np.sqrt(np.sqrt(x)**2*da**2+(a**2/(2*np.sqrt(x)))*dh**2+db**2)

def fit(file_name, Po, ci=1, cf=1, xerr=2, yerr_type="Measurement", yerr=3):
    """
    Fits the function func to the average of the data from ci to cf, 
    with the associated error for each row if err_type is "Measurement";
    if err_type is "Standard deviation", then it calculates the standard
    deviation of the data and uses it as the error.
    Requires a pre-defined func.
    Po: List of the initial guess of the parameters.
    Po_err: List of the respective measurement error of the parameters
    """
    dic_ydata=avrg(file_name,ci,cf)
    ydata=[]
    for k in dic_ydata.keys():
        ydata.append(dic_ydata[k])
    xdata=dic_ydata.keys()
    sigma=[]
    if yerr_type=="Standard deviation":
        dic_sigma=std(file_name,ci,cf)
        for k in dic_sigma.keys():
            sigma.append(dic_sigma[k])
    elif yerr_type=="Measurement":
        dic_sigma=data_dictionary(file_name,[["yerr",yerr,yerr]])
        for k in dic_sigma.keys():
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
    
    popt, pcov = curve_fit(model, xdata, ydata, Po, sigma)
    plt.plot(xdataE, model(xdataE, *popt), 'r-', label='fit')
    plt.fill_between(xdataE,model(xdataE, *popt)-err_model(xdataE,popt[0],pcov[0][0],pcov[1][1],0.4),model(xdataE, *popt)+err_model(xdataE,popt[0],pcov[0][0],pcov[1][1],0.4),alpha=0.5,edgecolor='#CC4F1B', facecolor='#FF9848')
    
    plt.plot(xdataE, model(xdataE, *Po), 'g-', label='Model')
    plt.fill_between(xdataE, model(xdataE,*Po)-err_model(xdataE,Po[0],2.44,4.55,0.4), model(xdataE, *Po)+err_model(xdataE,Po[0],2.44,4.55,0.4),alpha=0.5, edgecolor='#3F7F4C', facecolor='#7EFF99',linewidth=0)

    plt.xlabel('Altura inicial [cm]')
    plt.ylabel('Tiempo [s]')
    plt.legend()
    plt.title('Promedios de todas las medidas')
    plt.show()
    
    print popt
    print pcov
