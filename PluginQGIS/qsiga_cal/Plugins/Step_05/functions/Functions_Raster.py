# -*- coding: utf-8 -*-
"""
Module Raster
-------------

Author
------
Rendón-Álvarez, J. P., Tue Apr 25 13:02:46 2017.
E-mails: jprendona@unal.edu.co, pablo.rendon@udea.edu.co

Description
-----------
This module contains functions to read and write raster files

License
-------
This work is licensed under the Creative Commons Attribution-NonComercial-
ShareAlike 4.0 International License. To view a copy o this license, visit
http://creativecommons.org/licenses/by-nc-sa/4.0

Notices
-------
This program is distributed in the hope that it will be used, but without any
warranty. No author or distributor accepts responsability to anyone for
consequences of using them or for whether they serve any particular purpose
or work at all, unless he says so in writing. Everyone is granted permission
to copy, modify an distribute this program, but only under the condition that
this notice and above authority and license notices remain intact
"""

import numpy as np
import struct
from osgeo import gdal
from netCDF4 import Dataset
from rasterio.warp import calculate_default_transform, reproject, Resampling
import copy as copy
import rasterio
from rasterio.features import rasterize
import geopandas as gpd
import os
import pandas as pd

def ReadRaster(path):

    u"""
    It loads a raster after identifiying the file extension

    Parameters
    ----------
    path : text
        Raster file direction in the hard disk

    Returns
    -------
    rst : dictionary
        Raster dictionary with the next attribute naming convention:
            - ncols : Number of columns
            - nrows : Number of rows
            - xll   : X-axis coordinate of the lower left corner
            - yll   : Y-axis coordinate of the lower left corner
            - xur   : X-axis coordinate of the upper right corner
            - yur   : Y-axis coordinate of the upper right corner
            - clsz  : Cell size
            - nodt  : Missing data value
            - mtrx  : Data matrix
    """

    # It identifies the file extension and reads the raster
    ext=path.lower().split('.')[-1]
    if ext=='asc':
        rst=ReadTifRaster(path)
    elif ext=='tif' or ext=='tiff':
        rst=ReadTifRaster(path)
    else:
        print('File format not recognized')
        rst=BuildRaster(0.0,0.0,0.0,0.0,np.array(0))
    return rst

def WriteRaster(path,rst,prcsn=gdal.GDT_Float32):

    u"""
    It writes a raster after identifiying the file extension

    Parameters
    ----------
    path : text
        Raster file direction in the hard disk

    rst : dictionary
        Raster dictionary with the next attribute naming convention:
            - ncols : Number of columns
            - nrows : Number of rows
            - xll   : X-axis coordinate of the lower left corner
            - yll   : Y-axis coordinate of the lower left corner
            - xur   : X-axis coordinate of the upper right corner
            - yur   : Y-axis coordinate of the upper right corner
            - clsz  : Cell size
            - nodt  : Missing data value
            - mtrx  : Data matrix
    """

    # It identifies the file extension and writes the raster
    ext=path.lower().split('.')[-1]
    if ext=='asc':
        rst=WriteAsciiRaster(path,rst)
    elif ext=='tif' or ext=='tiff':
        rst=WriteTifRaster(path,rst,prcsn)
    else:
        print('File format not recognized')

def ReadAsciiRaster(path):

    u"""
    It loads an Arc/ASCII raster file from the hard disk

    Parameters
    ----------
    path : text
        Raster file direction in the hard disk

    Returns
    -------
    rst : dictionary
        Raster dictionary with the next attribute naming convention:
            - ncols : Number of columns
            - nrows : Number of rows
            - xll   : X-axis coordinate of the lower left corner
            - yll   : Y-axis coordinate of the lower left corner
            - xur   : X-axis coordinate of the upper right corner
            - yur   : Y-axis coordinate of the upper right corner
            - clsz  : Cell size
            - nodt  : Missing data value
            - mtrx  : Data matrix
    """

    #It reads the header information as text
    asc_file=open(path,'r')
    asc_line=asc_file.readline()
    temp,ncols=str.split(asc_line)
    asc_line=asc_file.readline()
    temp,nrows=str.split(asc_line)
    asc_line=asc_file.readline()
    temp,xll=str.split(asc_line)
    asc_line=asc_file.readline()
    temp,yll=str.split(asc_line)
    asc_line=asc_file.readline()
    temp,clsz=str.split(asc_line)
    asc_line=asc_file.readline()
    temp,nodt=str.split(asc_line)
    asc_file.close()

    #It converts header information from text to numbers
    ncols=int(ncols)
    nrows=int(nrows)
    xll=float(xll)
    yll=float(yll)
    clsz=float(clsz)
    nodt=float(nodt)

    #It reads the matrix information, builds the raster and returns
    mtrx=np.loadtxt(path,skiprows=6)

    #It builds the Raste dictionary and returns
    rst=BuildRaster(xll,yll,clsz,nodt,mtrx)
    ChangeNoData(rst,-9999.0)
    return rst

def ReadNetCDFRaster(path,vrbl):

    u"""
    It loads a 2D netCDF raster file from the hard disk

    Parameters
    ----------
    path : text
        Raster file direction in the hard disk
    vrbl : text
        Variable to be loaded

    Returns
    -------
    rst : dictionary
        Raster dictionary with the next attribute naming convention:
            - ncols : Number of columns
            - nrows : Number of rows
            - xll   : X-axis coordinate of the lower left corner
            - yll   : Y-axis coordinate of the lower left corner
            - xur   : X-axis coordinate of the upper right corner
            - yur   : Y-axis coordinate of the upper right corner
            - clsz  : Cell size
            - nodt  : Missing data value
            - mtrx  : Data matrix
    """

    # It reads the georeferenciation properties
    dataset=Dataset(path)
    dims=[dim for dim in dataset.dimensions]
    xdim=np.array(dataset.variables[dims[0]])
    ydim=np.array(dataset.variables[dims[1]])
    dx=np.average(xdim[1:]-xdim[0:np.size(xdim)-1])
    dy=np.average(ydim[1:]-ydim[0:np.size(ydim)-1])
    clsz=0.5*(dx+dy)
    xll=np.min(xdim)-0.5*dx
    yll=np.min(ydim)-0.5*dy
    nodt=dataset.variables[vrbl]._FillValue

    # It loads the netCDF Dataset
    mtrx=np.fliplr(np.array(dataset.variables[vrbl])).T
    mtrx[mtrx==np.NaN]=nodt
    dataset.close()

    # It builds the Raster dictionary and returns
    rst=BuildRaster(xll,yll,clsz,nodt,mtrx)
    ChangeNoData(rst,-9999.0)
    return rst

def ReadTifRaster(path):

    u"""
    It loads a GeoTiff raster file from the hard disk

    Parameters
    ----------
    path : text
        Raster file direction in the hard disk
    vrbl : text
        Variable to be loaded

    Returns
    -------
    rst : dictionary
        Raster dictionary with the next attribute naming convention:
            - ncols : Number of columns
            - nrows : Number of rows
            - xll   : X-axis coordinate of the lower left corner
            - yll   : Y-axis coordinate of the lower left corner
            - xur   : X-axis coordinate of the upper right corner
            - yur   : Y-axis coordinate of the upper right corner
            - clsz  : Cell size
            - nodt  : Missing data value
            - mtrx  : Data matrix
    """

    # It reads the georeferenciation properties and the data matrix
    tif=gdal.Open(path)
    mtrx=tif.GetRasterBand(1).ReadAsArray().astype(float)
    georef=tif.GetGeoTransform()
    if georef[1]<0:
        xll=georef[0]+georef[1]*mtrx.shape[1]#np.size(mtrx,1)
    else:
        xll=georef[0]
    if georef[5]<0:
        yll=georef[3]+georef[5]*mtrx.shape[0]#np.size(mtrx,0)
    else:
        yll=georef[3]
    clsz=0.5*(np.abs(georef[1])+np.abs(georef[5]))
    nodt=tif.GetRasterBand(1).GetNoDataValue()
    tif=None

    # It builds the Raster dictionary and returns
    rst=BuildRaster(xll,yll,clsz,nodt,mtrx)
    ChangeNoData(rst,-9999.0)
    return rst

def WriteAsciiRaster(path,rst):
    
    u"""
    It saves a Raster dictionary in the hard disk using the Arc/ASCII format

    Parameters
    ----------
    path: text
        Direction in the hard disk in which the Arc/ASCII file will be written
    rst: dictionary
        AsciiRaster dictionary to export
    """

    #It builds the header information string
    hdr='ncols        '+str(rst['ncols'])\
    +'\nnrows        '+str(rst['nrows'])\
    +'\nxllcorner    '+str(rst['xll'])\
    +'\nyllcorner    '+str(rst['yll'])\
    +'\ncellsize     '+str(rst['clsz'])\
    +'\nnodata_value '+str(rst['nodt'])

    #It saves data in the hard disk and finishes
    file=open(path,'w')
    file.write(hdr)
    for row in rst['mtrx']:
        line='\n'+' '.join(row.astype('str'))
        file.write(line)
    file.close()

def WriteTifRaster(path,rst,prcsn=gdal.GDT_Float32):

    u"""
    It saves a Raster dictionary in the hard disk using the GeoTiff format

    Parameters
    ----------
    path: text
        Direction in the hard disk in which the GeoTiff file will be written
    rst: dictionary
        AsciiRaster dictionary to export
    prcsn: optional, integer
        Data precision. By default gdal.GDT_Float32, but any other data type
        can be specified
    """

    # It creates the driver and the file
    driver=gdal.GetDriverByName("GTiff")
    tif=driver.Create(path,rst['ncols'],rst['nrows'],1,prcsn)

    # It sets the metadata, band data and finishes
    tif.SetGeoTransform((rst['xll'],rst['clsz'],0.0,rst['yur'],0.0,-rst['clsz']))
    band=tif.GetRasterBand(1)
    band.SetNoDataValue(rst['nodt'])
    band.WriteArray(rst['mtrx'],0,0)
    band.FlushCache()
    tif=None;band=None

def BuildRaster(xll,yll,clsz,nodt,mtrx):

    u"""
    It builds a Raster dictionary when all of its attributes are specified

    Parameters
    ----------
    xll : float
        X-axis coordinate of the lower left corner
    yll : float
        Y-axis coordinate of the lower left corner
    clsz : float
        Cell size
    nodt : number
        Missing data value
    mtrx : float
        Data array

    Returns
    -------
    rst : dictionary
        A dictionary with the raster information
    """

    #It calculates aditional necessary information
    nrows=np.size(mtrx,0)
    ncols=np.size(mtrx,1)
    nclls=nrows*ncols
    xur=xll+ncols*clsz
    yur=yll+nrows*clsz

    #It assambles the Raster dictinary and returns
    rst={'ncols':ncols,'nrows':nrows,'nclls':nclls,'xll':xll,'yll':yll,\
    'xur':xur,'yur':yur,'clsz':clsz,'nodt':nodt,'mtrx':mtrx}
    return rst

def ChangeNoData(rst,new_nodt):

    u"""
    It replaces the no data value in the data array for an specified new no data
    value

    Parameters
    ----------
    rst : dictionary
        Raster dictionary
    new_nodt : number
        New no data value
    """
    rst['mtrx'][np.where(rst['mtrx']==rst['nodt'])]=new_nodt
    rst['nodt']=new_nodt

def ClipRaster(rst,limit):

    u"""
    It clips a raster within a square window limit

    Parameters
    ----------
    rst : dictionary
        Raster dictionary
    limit : flat
        Vector of square window limit coordinates with the next order
        [xmin,ymin,xmax,ymax]
    """

    # It defines the row-col index limits
    col_min=max(0,np.int((limit[0]-rst['xll'])/rst['clsz']))
    col_max=min(np.int((limit[1]-rst['xll'])/rst['clsz'])+1,rst['ncols'])
    row_min=max(0,np.int((rst['yur']-limit[3])/rst['clsz']))
    row_max=min(np.int((rst['yur']-limit[2])/rst['clsz'])+1,rst['nrows'])

    # It clips the data matrix and defines the new georeferenciation
    mtrx=rst['mtrx'][row_min:row_max,col_min:col_max]
    xll=rst['xll']+col_min*rst['clsz']
    yll=rst['yur']-row_max*rst['clsz']

    # It builds the clipped raster and returns
    clip_rst=BuildRaster(xll,yll,rst['clsz'],rst['nodt'],mtrx)
    return clip_rst

def NullsInBounds(rst):

    u"""
    It sets no data values in the bounds of raster

    Parameters
    ----------
    rst : dictionary
        Raster dictionary
    """

    # It sets no data values in the bounds
    rst['mtrx'][0,:]  = rst['nodt']
    rst['mtrx'][-1,:] = rst['nodt']
    rst['mtrx'][:,0]  = rst['nodt']
    rst['mtrx'][:,-1] = rst['nodt']

def XYtoij(X,Y,raster):
    """
    Transform X Y coordinates to i j indices

    Parameters
    ----------
    X:float or list of floats
        X coordinate

    Y: float or list of floats
        Y coordinate

    raster: raster object
            Base raster to obtain matrix indices

    Returns
    -------
    i: int
        Row index

    j: int
        Colum index

    """
    i=np.floor(raster['nrows']-((Y-raster['yll'])/raster['clsz'])).astype(int)
    j=(np.floor((X-raster['xll'])/raster['clsz'])).astype(int)     
    return i,j

def IJtoxy(i,j,raster):
    """
    Transform X Y coordinates to i j indices

    Parameters
    ----------
    i: int or list of int
        Row index

    j: int or list of int
        Colum index

    raster: raster object
            Base raster to obtain matrix indices

    Returns
    -------
    X:float
        X coordinate

    Y: float
        Y coordinate
    """
    Y=((raster['nrows']-i)*(raster['clsz']))+raster['yll']-raster['clsz']/2
    X=raster['xll']+j*raster['clsz']+raster['clsz']/2
    return X,Y

def DrainageMap(dir_rst,FlowVectorA,AuxReajuste,topic,path=''):

    u"""
    Crea mapas en red

    Entradas
    ----------
    dir_rst: Ráster
             Mapa de dirección de flujo
    FlowVector :  DataFrame
                  Con atributos de la cuenca, las columnas iniciales son:
                  ['ID'] Id celda
                  ['IDD] Id celda destino
                  ['Area_Acum'] Area acumulada
                  ['Lon_NoaAcum'] Longitud no acumulada
                  ['Elevacion'] Cota del DEM
                  ['X'] Coordenada este
                  ['Y'] Coordenada norte
                  ['i'] índice fila de la matriz ráster
                  ['j'] índice columna de la matriz ráster
                  ['Pendiente'] pendiente
                  ['Tipo'] tipo celda [0: ladera; 1:cauce; 2: drenaje efímero 3:embalse, 4: Pie de presa]
                  ['Nodos'] Nodos [0:es cauce; 1:cabecera; 2:hidrológico; 3:topográfico]
                  ['Elev_corr'] Cota corregida cauce
                  ['ID_tramo'] Número de tramo al que corresponde cada celda que
                  sea cauce
                  ['Beechie'] Clasificación de Beechie para celda tipo cauce
                  ['Flores'] Clasificación de Flores para celda tipo cauce

    topic: str
            Tipo de mapa que se va a crear:
            'Tipo':Dibuja red de drenaje
            'Red': Grafica los nodos en la red de drenaje y la red
            'ID_tramo': Grafica los tramos de la red
            'Beachie': Grafica la clasificacion segun flores
            'Flores': Grafica la clasificacion segun Beechie

    path: string
            Donde se va a almacenar el mapa creado

    Salidas
    -------
    Map : ASC
            Devuelve un mapa en formato ASC

    """
    #Para guardar en raster
    mapT1     = np.zeros((dir_rst['nrows'],dir_rst['ncols']))-9999.0
    isnetwork = np.where(np.logical_or(FlowVectorA['Tipo']==1.0, FlowVectorA['Tipo']>3.0))[0]
    i         = FlowVectorA.loc[isnetwork,'i'].astype(int)
    j         = FlowVectorA.loc[isnetwork,'j'].astype(int)
    mapT1[i,j]= FlowVectorA.loc[isnetwork,topic]
    map_rst_1 ={'ncols':dir_rst['ncols'],'nrows':dir_rst['nrows'],'nclls':dir_rst['nclls'],'xll':dir_rst['xll'],'yll':dir_rst['yll'],'xur':dir_rst['xur'],'yur':dir_rst['yur'],'clsz':dir_rst['clsz'],'nodt':-9999,'mtrx':mapT1,'prj':''}
    if path is not '':
        WriteRaster(path,map_rst_1)
    return map_rst_1

def SampleRastToRast(rasterA,rasterB):
    """
    Turns values from rasterA into shape and resolution of rasterB

    Parameters
    ----------

    rasterA:
        raster object with data to sample.

    rasterB:
        raster object with desired spatial resolution and domain.

    Returns
    -------

    raster:
        raster object with sampled data on the desired resolution and domain.
    """
    xll=rasterB['xll']
    yll=rasterB['yll']
    clsz=rasterB['clsz']
    nodt=rasterB['nodt']
    mtrx=rasterB['mtrx'] * 0. + nodt
    rasterA['mtrx'][np.where(rasterA['mtrx']==rasterA['nodt'])]=nodt
    iB=np.argwhere(np.isreal(mtrx))[:,0]
    jB=np.argwhere(np.isreal(mtrx))[:,1]
    xB,yB=IJtoxy(iB,jB,rasterB)
    iBA,jBA=XYtoij(xB,yB,rasterA)
    iA=np.argwhere(np.isreal(rasterA['mtrx']))[:,0]
    jA=np.argwhere(np.isreal(rasterA['mtrx']))[:,1]

    corrBA=np.where((iBA>=iA.min())*(iBA<=iA.max())*(jBA>=jA.min())*(jBA<=jA.max()))[0]
    mtrx[iB[corrBA],jB[corrBA]]=rasterA['mtrx'][iBA[corrBA],jBA[corrBA]]

    raster=BuildRaster(xll,yll,clsz,nodt,mtrx)
    return raster

def MergeRst(path_rst1='', path_rst2='', raster1='',raster2=''):
    if path_rst1 is not '':
            raster1=ReadRaster(path_rst1)
            raster2=ReadRaster(path_rst2)
    merged=copy.copy(raster1)

    yur=max(raster1['yur'],raster2['yur'])
    xur=max(raster1['xur'],raster2['xur'])
    yll=min(raster1['yll'],raster2['yll'])
    xll=min(raster1['xll'],raster2['xll'])

    nrowup=max(np.ceil((yur-raster1['yur'])/raster1['clsz']),0)
    nrowdown=max(np.ceil((raster1['yll']-yll)/raster1['clsz']),0)
    ncolr=max(np.ceil((xur-raster1['xur'])/raster1['clsz']),0)
    ncoll=max(np.ceil((raster1['xll']-xll)/raster1['clsz']),0)

    merged['ncols']=int(merged['ncols']+ncolr+ncoll)
    merged['nrows']=int(merged['nrows']+nrowup+nrowdown)

    merged['yur']=merged['yur']+merged['clsz']*nrowup
    merged['xur']=merged['xur']+merged['clsz']*ncolr
    merged['yll']=merged['yll']-merged['clsz']*nrowdown
    merged['xll']=merged['xll']-merged['clsz']*ncoll

    merged['mtrx']=np.insert(merged['mtrx'],list(np.zeros((nrowup)).astype(int)),raster1['nodt'],axis=0)
    merged['mtrx']=np.insert(merged['mtrx'],list(np.ones((nrowdown)).astype(int)*np.shape(merged['mtrx'])[0]),raster1['nodt'],axis=0)
    merged['mtrx']=np.insert(merged['mtrx'],list(np.zeros((ncoll)).astype(int)),raster1['nodt'],axis=1)
    merged['mtrx']=np.insert(merged['mtrx'],list(np.ones((ncolr)).astype(int)*np.shape(merged['mtrx'])[1]),raster1['nodt'],axis=1)

    raster2=SampleRastToRast(raster2,merged)
    merged['mtrx'][np.where(raster2['mtrx']!=raster2['nodt'])]=raster2['mtrx'][np.where(raster2['mtrx']!=raster2['nodt'])]
    #WriteRaster(generalpth+'merged.tif',merged)
    return merged

def MultiMerge(lst_rst_pth, consecutive_val=False):
    if len(lst_rst_pth)==1:
        r_1=ReadRaster(lst_rst_pth[0])
        merged=copy.copy(r_1)
        return merged
    else:
        r_1=ReadRaster(lst_rst_pth[0])
        merged=copy.copy(r_1)
        for i in lst_rst_pth[1:]:
            r_2=ReadRaster(i)
            if consecutive_val:
                max1=np.nanmax(merged['mtrx'])
                r_2['mtrx'][np.where(r_2['mtrx']!=r_2['nodt'])]=r_2['mtrx'][np.where(r_2['mtrx']!=r_2['nodt'])] + max1 + 1
            merged=MergeRst(raster1=merged,raster2=r_2)
        return merged

def raster2sgabr(path):
    
    u"""
    Convierte un archivo ráster (con extensión .asc o .tif) en un archivo .sgabr compatible con el 
    modelo SIGA.
    
    Parámetros
    ----------
    path : texto
        Dirección del archivo ráster de entrada en el disco duro.
    """
    
    # Lectura del archivo ráster
    rst = ReadRaster(path)
    rst['mtrx'] = rst['mtrx'].astype('f')

    # Empaque de los metadatos en una variable binaria
    bin_mtdt = struct.pack('f'*6,rst['ncols'],rst['nrows'],rst['xll'],rst['yll'],rst['clsz'],rst['nodt'])
    # Empaque de la matriz de datos en una variable binaria
    bin_mtrx = rst['mtrx'].tobytes()

    # Escritura de los metadatos en un archivo .sgabr homónimo al ráster de entrada
    path_out = path[:-len(path.split('.')[-1])]+'sgabr'
    new_rst = open(path_out,'wb')
    new_rst.write(bin_mtdt)
    new_rst.close()
    # Escritura de la matriz de datos en el mismo archivo .sgabr
    new_rst = open(path_out,'ab')
    new_rst.write(bin_mtrx)
    new_rst.close()
    
    del bin_mtdt, bin_mtrx
    return rst

def sgabr2raster(path):
    
    u"""
    Convierte un archivo binario .sgabr en un archivo ráster (con extensión .asc o .tif).
    
    Parámetros
    ----------
    path : texto
        Dirección del archivo binario de entrada en el disco duro.
    """
    
    # Lectura del archivo .sgabr
    bin_file = open(path,'rb')
    bin_content = bin_file.read()
    bin_file.close()

    # Creación del diccionario que almacena los datos necesarios para construir el ráster
    rst = {}

    # Rastreo de la información entre los datos binarios
    (rst['ncols'],rst['nrows'],rst['xll'],rst['yll'],rst['clsz'],rst['nodt'])=struct.unpack('ffffff',bin_content[0:24])
    rst['ncols'] = int(rst['ncols'])
    rst['nrows'] = int(rst['nrows'])
    rst['nclls'] = rst['ncols']*rst['nrows']
    rst['xur'] = rst['xll']+rst['ncols']*rst['clsz']
    rst['yur'] = rst['yll']+rst['nrows']*rst['clsz']
    rst['mtrx'] = np.reshape(struct.unpack('f'*int(rst['nclls']),bin_content[24:]),(rst['nrows'],rst['ncols']))
    del bin_content

    # Exportación del archivo ráster
    path_out = path[:-len(path.split('.')[-1])-1]+'_conv.tif'
    WriteRaster(path_out,rst)
    
    return rst

def PolyponToRaster(PathShp,PathOutput,PathRaster='',Value='',shape=(1000,1000)):

    # Leer archivo shapefile
    df = gpd.read_file(PathShp)

    # Número de filas y columnas
    if PathRaster != '':
        # Leer Raster
        DEM = rasterio.open(PathRaster)

        # Leer tamaño del raster
        shape = DEM.shape

        # Leer transformación
        transform = DEM.transform
    else:
        # Asignar transformación
        transform = rasterio.transform.from_bounds(*df['geometry'].total_bounds, *shape)

    if PathRaster == '':
        dtype = rasterio.uint16
        # Convertir shp a raster
        rasterize_rivernet = rasterize(
            [(geo, 1) for geo in df['geometry']],
            out_shape=shape,
            transform=transform,
            fill=0,
            all_touched=True,
            dtype=dtype)
    else:
        if np.issubdtype(df[Value], np.integer):
            dtype = np.int32
        elif np.issubdtype(df[Value].max(), np.floating):
            dtype = np.float32

        # Convertir shp a raster
        rasterize_rivernet = rasterize(
            [(geo, va) for geo, va in zip(df['geometry'],df[Value])],
            out_shape=shape,
            transform=transform,
            fill=0,
            all_touched=True,
            dtype=dtype)

    # Guardar raster
    with rasterio.open(
        PathOutput, 'w',
        driver='GTiff',
        dtype=dtype,
        count=1,
        nodata=-99,
        width=shape[1],
        height=shape[0],
        transform=transform
    ) as dst:
        dst.write(rasterize_rivernet, indexes=1)

def resample_to_match3(input_raster, reference_raster, output_raster):
    with rasterio.open(reference_raster) as ref:
        transform, width, height = calculate_default_transform(
            ref.crs, ref.crs, ref.width, ref.height, *ref.bounds)
        kwargs = ref.meta.copy()
        kwargs.update({
            'crs': ref.crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(input_raster) as src:
            with rasterio.open(output_raster, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=ref.crs,
                        resampling=Resampling.bilinear)
def resample_to_match(input_raster, reference_raster, output_raster, integer_output=False):
    with rasterio.open(reference_raster) as ref:
        transform, width, height = calculate_default_transform(
            ref.crs, ref.crs, ref.width, ref.height, *ref.bounds)
        kwargs = ref.meta.copy()
        kwargs.update({
            'crs': ref.crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        # Modificar el tipo de datos si integer_output es True
        if integer_output:
            kwargs['dtype'] = 'int16'  # Tipo entero de 16 bits, ajustable según necesidad
            kwargs['nodata'] = -32768  # Definir un valor nodata si es necesario para enteros
            resampling_method = Resampling.nearest  # Vecino más cercano para mantener enteros
        else:
            resampling_method = Resampling.bilinear  # Interpolación bilineal para datos continuos

        with rasterio.open(input_raster) as src:
            with rasterio.open(output_raster, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=ref.crs,
                        resampling=resampling_method)

def Erodability(PathK, PathDEM, PathSand, PathSilt, PathClay, PathOMC):
    """
    Está función estima la erodabilidad del suelo a partir de
    la textura del suelo y el contenido de materia orgánica.
    Los valores se toman del trabajo de:
    Wischmeier, Johnson y Cross 1971 (reportado en Roose, 1996).
    La tabla utilizada se puede encontrar en el siguiente link:
    http://releases.naturalcapitalproject.org/invest-userguide/latest/en/sdr.html#soil-erodibility-k

    Nota: Los datos de la tabla están multiplicados por 0.1317
    para convertir de EE. UU. a (ton*ha*hr)/(ha*MJ*mm)

    Parameters
    ----------
    PathK: text
        Ruta del raster de salida de erodabilidad del suelo en (ton*ha*hr)/(ha*MJ*mm)
    PathSand: text
        Ruta del raster del contenido de arenas con valores de [0 - 1]
    PathSilt: text
        Ruta del raster del contenido de limos con valores de [0 - 1]
    PathClay: text
        Ruta del raster del contenido de arcillas con valores de [0 - 1]
    PathOMC: text
        Ruta del raster del contenido de materia orgánica con valores de [0 - 1]
    """

    # Raster Sand
    DEM = ReadRaster(PathDEM)

    # Raster Sand
    sand = SampleRastToRast(ReadRaster(PathSand),DEM)

    # Raster Silt
    silt = SampleRastToRast(ReadRaster(PathSilt),DEM)

    # Raster Clay
    clay = SampleRastToRast(ReadRaster(PathClay),DEM)

    # OM
    OMC = SampleRastToRast(ReadRaster(PathOMC),DEM)

    # Crear raster vacio
    K = DEM
    K['mtrx'] = K['mtrx']*0

    # ------------------------------------------------------------------------------------------------------------------
    # Sand
    # ------------------------------------------------------------------------------------------------------------------
    id          = (silt['mtrx'] + (1.5*clay['mtrx'])) < 0.15
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.004
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)]  = 0.001

    # ------------------------------------------------------------------------------------------------------------------
    # Loamy sand
    # ------------------------------------------------------------------------------------------------------------------
    id          = ((silt['mtrx'] + (1.5*clay['mtrx'])) >= 0.15) & ((silt['mtrx'] + (2*clay['mtrx'])) < 0.3)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.007
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.005

    # ------------------------------------------------------------------------------------------------------------------
    # Sandy loam
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.07) & (clay['mtrx']<= 0.2) & (sand['mtrx'] > 0.52) & ((silt['mtrx'] + (2*clay['mtrx'])) >= 0.3)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.018
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.016

    id          = (clay['mtrx']< 0.07) & (silt['mtrx'] < 0.5) & ((silt['mtrx'] + (2*clay['mtrx'])) >= 0.3)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.018
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.016

    # ------------------------------------------------------------------------------------------------------------------
    # loam
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.07) & (clay['mtrx']<= 0.27) & (silt['mtrx'] >= 0.28) & (silt['mtrx'] < 0.5) & (sand['mtrx'] <= 0.52)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.045
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.034

    # ------------------------------------------------------------------------------------------------------------------
    # silt Loam
    # ------------------------------------------------------------------------------------------------------------------
    id          = ((silt['mtrx']>=0.5) & (clay['mtrx']>= 0.12) & (clay['mtrx']<0.27)) | ((silt['mtrx']>=0.5) & (silt['mtrx']<0.8) & (clay['mtrx']<0.12))
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.054
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.049

    # ------------------------------------------------------------------------------------------------------------------
    # silt
    # ------------------------------------------------------------------------------------------------------------------
    id          = (silt['mtrx'] >= 0.8) & (clay['mtrx']< 0.12)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.054
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.049

    # ------------------------------------------------------------------------------------------------------------------
    # Sandy Clay Loam
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.2) & (clay['mtrx']< 0.35) & (silt['mtrx'] < 0.28) & (sand['mtrx'] > 0.45)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.026
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.026

    # ------------------------------------------------------------------------------------------------------------------
    # Clay Loam
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.27) & (clay['mtrx']< 0.4) & (sand['mtrx'] > 0.2) & (sand['mtrx'] <= 0.45)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.043
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.037

    # ------------------------------------------------------------------------------------------------------------------
    # Silty Clay Loam
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.27) & (clay['mtrx']< 0.4) & (sand['mtrx'] <= 0.2)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.046
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.046

    # ------------------------------------------------------------------------------------------------------------------
    # Sandy Clay
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.35) & (sand['mtrx'] >= 0.45)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.026
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.026

    # ------------------------------------------------------------------------------------------------------------------
    # Silty Clay
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.4) & (silt['mtrx'] >= 0.4)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.036
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.034

    # ------------------------------------------------------------------------------------------------------------------
    # Clay
    # ------------------------------------------------------------------------------------------------------------------
    id          = (clay['mtrx']>= 0.4) & (sand['mtrx'] <= 0.45) & (silt['mtrx'] < 0.4)
    id_1        = OMC['mtrx'] <= 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.032
    id_1        = OMC['mtrx'] > 0.02
    K['mtrx'][np.logical_and(id,id_1)] = 0.028

    # Guardar raster
    WriteTifRaster(PathK, K)


def shp_to_csv_muestreador(path_out, shp_name):
    gdf = gpd.read_file(shp_name)
    # Convertir a DataFrame de pandas
    if 'Cod_P' in gdf.columns:
        # Pivoting the DataFrame
        df = pd.DataFrame(gdf)
        pivot_df = df.pivot(index='Cod_P', columns='Cod_T', values=['X', 'Y'])
        pivot_df.columns = [f'{x1}{x2}' for x1, x2 in pivot_df.columns]
        # Adding the 'Nombre' column
        pivot_df['Nombre'] = df.drop_duplicates('Cod_P')['Nombre'].values
        # Resetting the index
        pivot_df.reset_index(drop=True, inplace=True)
        # Renaming the columns to match the required format
        pivot_df.columns = ['X1', 'X2', 'Y1', 'Y2', 'Nombre']
        # Reordering the columns
        pivot_df = pivot_df[['X1', 'Y1', 'X2', 'Y2', 'Nombre']]
        pivot_df.to_csv(os.path.join(path_out, '04_TramosControl.csv'), index=False)
    else:
        df = gdf[gdf.columns[1:-1]]
        df.columns = ['coord_x_or_aju', 'coord_y_or_aju', 'Nombre', 'Estadistico']
        df.to_csv(os.path.join(path_out, '02_PuntosControl.csv'), index=False)