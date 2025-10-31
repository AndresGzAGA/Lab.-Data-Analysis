import numpy as np
import pandas as pd

def data_frame(file_name, latex_name):
    """
    From a file_name which contains a table of numbers separated with tabulator,
    returns a DataFrame with those numbers, not ignoring the first row. 
    """
    infile=open(file_name,"r")
    #infile.readline() #skip the first line
    matrix=[]
    T=True
    for line in infile:
        num_row=[]
        str_row=line.split()
        for n in str_row:
            #if len(matrix)!=0:
                #n=float(n)
            num_row.append(n)
        if T==True: #Guardar los nombres de las columnas aparte
            name_col=line.split()
        else:
            matrix.append(num_row)
        T=False
    matrix=pd.DataFrame(data=matrix, columns=name_col)
    infile.close()
    
    with open(latex_name,'w') as tf:
        tf.write(matrix.to_latex())
        
    
    
#d=data_frame(filename)
#d.to_latex()
    
    
    
    
    
    