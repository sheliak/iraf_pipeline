#! /usr/bin/env python
import ebf
import pdb
import numpy
import sys
import time
import os
import glob


def write(filename,data,keys=None,nsize=None):
    if type(data)==dict:
        if keys == None:
            keys=data.keys()
        for key in keys:
            if type(data[key])!=numpy.ndarray:
                raise RuntimeError('keys must be numpy array')
        nsize1=data[keys[0]].size
        for key in keys:
            if data[key].size!=nsize1:
                raise RuntimeError('dict items must be of same size')
        if nsize == None:
            nsize=nsize1
        else:
            if nsize>nsize1:
                raise RuntimeError('input size too large')

        with open(filename,'w') as fp:
            fp.write(str(keys)[1:-1].replace("'",'')+'\n')        
            for i in range(nsize):
                s=",".join([str(data[key][i]) for key in keys])        
                fp.write(s+'\n')        
    elif type(data)==numpy.ndarray:
        if data.dtype.char == 'V':            
            with open(filename,'w') as fp:
                s=str(data.dtype.names)
                s=s[1:-1].replace("'",'')+'\n'
                fp.write(s)
                for temp in data:
                    s=str(temp)
                    s=s[1:-1].replace("'",'')+'\n'
                    fp.write(s)
    print 'Written file:',filename



def isfloat(mystr):
    try:
        float(mystr)
        return True
    except ValueError:
        return False

def isint(mystr):
    try:
        int(mystr)
        return True
    except ValueError:
        return False

def cformat2type_width(myformat):
    myformat=myformat.strip(' \t\n\r\"\'')
    myformatl=myformat.split('%')
    if len(myformatl[0]) == 0:
        myformatl=myformatl[1:]
    datatypel=[]
    widthl=[]
    for x in myformatl:
        x=x.lstrip(' +-0#')
        width=x.count(' ')
        x=x.strip()
#        temp=x[-1]
#        width=width+int(x[0:-1].split('.')[0])
        temp=x[0]
        width=width+int(x[1:].split('.')[0])
        if temp == 'f':
            datatype='float64'
        elif temp == 'F':
            datatype='float64'
        elif temp == 'g':
            datatype='float64'
        elif temp == 'G':
            datatype='float64'
        elif temp == 'e':
            datatype='float64'
        elif temp == 'E':
            datatype='float64'
        elif temp == 'i':
            datatype='int64'
        elif temp == 'I':
            datatype='int64'
        elif temp == 'd':
            datatype='int64'
        elif temp == 'u':
            datatype='uint64'
        elif temp == 'A':
            datatype='S'
        elif temp == 'S':
            datatype='S'
        elif temp == 's':
            datatype='S'
        else:
            print temp
            datatype=''
            raise RuntimeError('Unknown Data type format')
        datatypel.append(datatype)
        widthl.append(width)
    print datatypel    
    print widthl    
    return [datatypel,widthl]


def get_datatypes_fast(f,datastart,delimiter,cols):
    datatypes=[]
    f.seek(0)
    i=-1
    for line in f:
        i=i+1
        dataline=line.strip(' \t\n\r')
        if (line[0] != '#')and(i >= datastart)and(len(dataline)>0):
            values=dataline.split(delimiter)
            if len(datatypes) == 0:
                for value in values:
                    datatypes.append('')
            status=0
            j=0
            for value in values:
                value1=value.strip(' \t\n\r\'\"')                
                if len(datatypes[j]) == 0:
                    if len(value1) == 0:
                        status=1
                    else:    
                        if isint(value1):
                            datatypes[j]='int64'
                        elif isfloat(value1):
                            datatypes[j]='float64'
                        else:
                            datatypes[j]='S'
                j=j+1
            if status == 0:
                break
    return datatypes

def str2typednum(x):
    if x.isdigit():
        return numpy.int64(x)
    elif isfloat(x):
        return numpy.float64(x)
    else:
        return numpy.array(x.strip())

def get_datatypes(f,datastart,delimiter,cols):
    status=0
    nfloat=numpy.zeros(cols,dtype=numpy.int64)
    ndigit=numpy.zeros(cols,dtype=numpy.int64)
    nvalid=numpy.zeros(cols,dtype=numpy.int64)
    f.seek(0)
    i=-1
    for line in f:
        i=i+1
        dataline=line.strip(' \t\n\r')
        if (line[0] != '#')and(i >= datastart)and(len(dataline)>0):
            values=dataline.split(delimiter)
            if len(values) < cols:
#                print 'In line no:',i,' cols expected=',cols,' actual cols=',len(values)
                while len(values) != cols:
                    values.append('null')
                # print line
                # raise RuntimeError('All rows do not have same number of colums')
            for j in range(cols):
                value1=values[j].strip(' \t\n\r\'\"').lower()
                if (value1!='null')and(value1!=''):
                    nvalid[j]=nvalid[j]+1
                    if isint(value1):
                        ndigit[j]=ndigit[j]+1
                    elif isfloat(value1):
                        nfloat[j]=nfloat[j]+1
    
    datatypes=[]
    for i in range(cols):
        if ndigit[i] == nvalid[i]:
            datatypes.append('int64')
        elif (ndigit[i]+nfloat[i]) == nvalid[i]:
            datatypes.append('float64')
        else:
            datatypes.append('S')
    return datatypes


def ascii_init(filename):
    f=open(filename,'r')
    params={}
    fields=[]
    datatypes=[]
    widths=[]
    units=[]
    ncols=0

#    print 'reading header'
    # read from header, i.e., commented portion
    for line in f:
        line1=line.strip(' \t\n\r')
        if (line[0] != '#')and(len(line1)>0):
            break
        line=line.strip(' \t\n\r#').split('=')
        quant=line[0].strip()
        if len(line) > 1:
            if quant == 'format':
                [datatypes,widths]=cformat2type_width(line[1])
            else:
                values1=line[1].strip(" []").split(',')
                values=[]
                for value in values1:
                    values.append(value.strip(' \t\n\r\'\"'))
                if quant == 'fields':
                    fields=values
                elif quant == 'datatypes':
                    datatypes=values
                elif quant == 'units':
                    units=values
                elif quant == 'params':
                    params={}
                    for temp in values:
                        key,value=temp.split(':')
                        params[key.strip().lower()]=str2typednum(value)
                elif quant == 'widths':
                    for value in values:
                        widths.append(int(value))


    if len(fields) > 0:
        fields1=[]
        i=0
        for field in fields:
            if len(field) == 0:
                fields1.append('col'+str(i))
                i=i+1
            else:
                fields1.append(field)
        fields=fields1
            

#    print 'reading delimiter'
    # get delimiter
    f.seek(0)
    for line in f:
        dataline=line.strip(' \t\n\r')
        if (line[0] != '#')and(len(dataline)>0):
            if ',' in dataline:
                delimiter=','
            elif '|' in dataline:
                delimiter='|'
            elif ';' in dataline:
                delimiter=';'
            else:
                delimiter=None
            break

#    print 'delimiter=',delimiter

#    print 'reading datastart'
    # get datastart and if csv format then also fields
    f.seek(0)
    i=-1
    status=0
    datastart=-1
    for line in f:
        i=i+1
        dataline=line.strip(' \t\n\r')
        if (line[0] != '#')and(datastart == -1)and(len(dataline)>0):
            values=dataline.split(delimiter)
            for value in values:
                value1=value.strip(' \t\n\r\'\"')
                if value1.isdigit() or isfloat(value1):
                    datastart=i
                    break
            if (datastart == -1)and(len(fields)==0):
                for value in values:
                    fields.append(value.strip(' \t\n\r\'\"'))
            else:
                datastart=i
                ncols=len(values)
                break


#    print 'reading datatype'
    # get datatype

    if len(datatypes) == 0:
        datatypes=get_datatypes(f,datastart,delimiter,ncols)

    f.close()

    # initialize to default values  and check length with ncols
    if len(datatypes) != ncols:
        print len(datatypes),ncols,'delimiter=',delimiter
        raise RuntimeError("len(datatypes) != ncols")

    if len(fields) == 0:
        for i in range(0,ncols):
            fields.append('col'+str(i))            
    else:
        if len(fields) != ncols:
            print len(fields),ncols,'delimiter=',delimiter
            raise RuntimeError("len(fields) != ncols")
        
    if len(units) == 0:
        for i in range(0,ncols):
            units.append('')            
    else:
        if len(units) != ncols:
            print len(units),ncols,'delimiter=',delimiter
            raise RuntimeError("len(units) != ncols")

    if len(widths) == 0:
        for i in range(0,ncols):
            widths.append(0)
    else:
        if len(widths) != ncols:
            print len(widths),ncols,'delimiter=',delimiter
            raise RuntimeError("len(widths) != ncols")

    return [fields,datatypes,widths,units,delimiter,datastart,params]

def read_from_schema(schemafile):
    f=open(schemafile,'r')
    mystr=f.read()
    f.close()
    mystr=mystr.replace('\n','')
    items=mystr.partition('(')[2].partition(')')[0].split(',')
    fields=[]
    datatypes=[]
    widths=[]
    units=[]
    for item in items:
        item=item.strip()
        if len(item) > 0:
            word=item.split()
            fields.append(word[0].strip(' \'\"'))
            datatypes.append(word[1].strip('\'\"'))
            if len(word) > 2:
                widths.append(int(word[2].strip('\'\"')))
            else:
                widths.append(0)            
            if len(word) > 3:
                units.append(word[3].strip('\'\"'))
            else:
                units.append('')
            if len(word) > 4:
                raise RuntimeError('Max 4 space separated word allowed. Enclose units in quotes if needed')

    return [fields,datatypes,widths,units]



def str2data(data,fields,datatypes,struct,params):    
#    print fields
#    print datatypes
    verbose=0
    dt=[]
    for j,field in enumerate(fields):
        if j == 0:
            nsize=len(data[field])
        else:
            if len(data[field]) != nsize:
                print field,len(data[field]),nsize
                print 'size of fields are not same',
                raise ValueError('')
        dt.append((field,datatypes[j]))

    s=set(zip(*dt)[0])
    for key,value in params.iteritems():
        if key not in s:
            dt.append((key,value.dtype))
        

#    data1=numpy.zeros(nsize,dtype=dt)
    if type(struct) == numpy.dtype:
        if len(dt) != len(struct):
            print 'Supplied size of dtype does not match with data on file'
        for temp in dt:
            if temp[0] not in struct.names:
                print struct.names
                print 'colname=',temp[0]
                print 'Supplied colname in dtype does not match with data on file:'            
        data1=numpy.zeros(nsize,dtype=struct)
    elif struct == True:
        data1=numpy.zeros(nsize,dtype=dt)
    else: 
        data1={}

    for key,value in params.iteritems():
        if key not in s:
            data1[key]=numpy.zeros(nsize,dtype=value.dtype)+value
        
        
    for j,field in enumerate(fields):                
        if verbose == 1:
            print "{0:20s} {1:10s} {2:d}".format(field,datatypes[j],len(data[field]))
        missing=''
        if datatypes[j][0] != 'S':
            if datatypes[j][0] == 'i':
                missing=numpy.iinfo(datatypes[j]).max
            elif datatypes[j][0] == 'u':
                missing=numpy.iinfo(datatypes[j]).max
            else:
                missing='inf'
        temp=numpy.array(data[field],dtype='S')
        ind=numpy.where((temp == 'NULL') |(temp == 'null') | (temp=='\N'))[0]
        if len(ind) > 0:
            if verbose == 1:
                print 'missing elements found using default missing value=',missing, len(data[field])
                print ind,len(ind)
            for i in ind:
                data[field][i]=missing
            
        try:        
            data1[field]=numpy.array(data[field],dtype=datatypes[j])
        except ValueError:
#            print 'ValueError',field,datatypes[j]
            try:        
                data1[field]=numpy.array(data[field],dtype='float64')
            except ValueError:
                for x in data[field]:
                    try:
                        temp=numpy.array(x,dtype=datatypes[j])
                    except ValueError:
                        print x,' Is it NULL?',(x == 'NULL')
                        raise ValueError('')
        except TypeError:
            print datatypes[j],data1[field]
            raise TypeError(' datattype='+datatypes[j])

    
    return data1



def read(filename,struct=False,schema=None,maxcols=None,params=False,sortby=None):
    t1=time.time()    
    [fields,datatypes,widths,units,delimiter,datastart,params1]=ascii_init(filename)        
    if schema != None:
        [fields,datatypes,widths,units]=read_from_schema(schema)
    # print 'fields=',fields
    # print 'datatypes=',datatypes
    # print 'widths=',widths
    # print 'units=',units
    # print 'delimiter=',delimiter
    # print 'datastart=',datastart
    fields=[temp.lower() for temp in fields]
    
    
    data={}
    for field in fields:
        data[field]=[]

    f=open(filename,'r')
    if sum(widths) == 0:
        i=0
        di=0
        for line in f:
            if (line[0] != '#')and(i >= datastart)and(len(line.strip(' \t\r\n'))>0):
                j=0
                for word in line.strip().strip('\r').split(delimiter):
                    word=word.strip()
                    if len(word) == 0:
                        word='NULL'
                    if j<len(fields):
                        data[fields[j]].append(word) 
                    j=j+1
                if (j>0)&(j < len(fields)):
                    # print line
                    # raise RuntimeError('not enough columns in line')
                    for l in range(j,len(fields)):
                        data[fields[l]].append('NULL') 
                di=di+1
            i=i+1
    else:
        start=[]
        end=[]
        temp=0
        for x in widths:
            start.append(temp)
            temp=temp+x
            end.append(temp)
            
        i=0
        di=0
        for line in f:
            if (line[0] != '#')and(i >= datastart)and(len(line.strip(' \t\r\n'))>0):
                for j in range(0,len(widths)):
                    data[fields[j]].append(line[start[j]:end[j]])                     
                di=di+1
            i=i+1

    f.close()        
    if maxcols != None:
        fields=fields[0:maxcols]
#    print fields
    if params==False:
        params1={}
        
    data=str2data(data,fields,datatypes,struct,params1)
    if sortby!=None:
        ind=numpy.argsort(data[sortby])
        if type(data) == dict:
            data1=data
            data={}
            for key in data1.keys():
                data[key]=data1[key][ind]
        else:
            data=data[ind]
    return data
    
        
 
def read_concat(pathname,params=False,sortby=None):
    dtype=None
    data1=[]
    struct=True
    for filename in glob.glob(pathname):
        data=read(filename,params=params,struct=struct,sortby=sortby)
        if type(struct) != numpy.dtype:
            struct=data.dtype
        data1.append(data)
    return numpy.concatenate(data1)
