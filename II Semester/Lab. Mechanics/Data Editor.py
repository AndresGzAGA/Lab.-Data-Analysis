import plotly.plotly as py
py.sign_in('andresgz', 'wKB8HNViDkHsCLhTIn0S')
from plotly.tools import FigureFactory as ff

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
            #if len(matrix)!=0:
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
    
def writes(file_name,something,column,name):
    """
    Writes something (string) on the column (int) with the name (string) of the column.
    """
    D=data_matrix(file_name)
    D[0].insert(column,name)
    for i in range(1,len(D)):
        D[i].insert(column,something)
    outfile=open(file_name,"w")
    for i in D:
        j=i[0]
        outfile.write(j)
        for j in i[1:]:
            outfile.write("\t" + j)
        outfile.write('\n')
    outfile.close()

def writel(file_name,l,column):
    """
    Writes the elements of the list l on the column. The list should include the name of the column
    """
    D=data_matrix(file_name)
    for i in range(0,len(D)):
        D[i].insert(column,l[i])
    outfile=open(file_name,"w")
    for i in D:
        j=i[0]
        outfile.write(j)
        for j in i[1:]:
            outfile.write("\t" + j)
        outfile.write('\n')
    outfile.close()


def delete (file_name,column):
    """
    Deletes the column (integer which represents its position).
    """
    D=data_matrix(file_name)
    for i in range(0,len(D)):
        del D[i][column]
    outfile=open(file_name,"w")
    for i in D:
        j=i[0]
        outfile.write(j)
        for j in i[1:]:
            outfile.write("\t" + j)
        outfile.write('\n')
    outfile.close()
    
def table(file_name,name):
    """
    Makes a table from a file_name and saves it in https://plot.ly/organize/andresgz:home with its name (string).
    """
    matrix=data_matrix(file_name)
    table = ff.create_table(matrix)
    py.iplot(table, filename=name)


"""
data_matrix = [['Country', 'Year', 'Population'],
               ['United States', 2000, 282200000],
               ['Canada', 2000, 27790000],
               ['United States', 2005, 295500000],
               ['Canada', 2005, 32310000],
               ['United States', 2010, 309000000],
               ['Canada', 2010, 34000000]]

table = ff.create_table(data_matrix)
py.iplot(table, filename='simple_table')
"""

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    