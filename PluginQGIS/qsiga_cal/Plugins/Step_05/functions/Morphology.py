# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:49:07 2017

@author: Gotta
"""
import os.path

#==============================================================================
#--Se importan paquetes requeridos
#==============================================================================
import numpy as np
from .Functions_Raster import *
import copy as copy
from scipy import spatial
from scipy import signal
from scipy import cluster
import pickle as cPickle
from osgeo import gdal
import pandas as pd
import geopandas as gpd
from shapely import wkt
from shapely.geometry import Point, LineString
from pyproj import Proj,transform
#==============================================================================
# Para imprimir
#==============================================================================
def print_message(msg):
    u"""
    Función para imprimir
    """
    print(msg)
#==============================================================================
#--Trazar cuenca
#==============================================================================

def TraceBasin(X,Y,dir_rst, return_raster=False, return_mask=False,XY_ini=''):
    #
    u"""
    Función para trazar una cuenca

    Entradas
    ----------
    X : float_like
        coordenada este de cierre de la cuenca
    Y : float_like
        coordenada norte de cierre de la cuenca
    dir_rst : RasterObject_like
        raster de dirección de flujo
    return_raster : Boolean_like (opcional)
        'True':  TraceBasin devuelve dos rasters (ID y IDD). 
        'False': TraceBasin devuelve a flow vector DataFrame and a list with matrix indices
    return_mask : Boolean_like (opcional)
        'True' TraceBasin devuelve máscara de cuenca con 1 donde está la cuenca y -9999 donde no
    XY_ini: array (opcional)
        array con las coordenadas este y norte de la(s) cabecera(s)


    Salidas
    -------
    rID : RasterObject_like
        Ráster con ID de las celdas que se encuentran dentro de la cuenca
    rIDD: RasterObject_like
        Ráster con ID de la celda a la que drena 
    fvec: DataFrame
        DataFrame con ID y IDD
    ijs: ndarray
        Lista de indices de la matriz ID 

    """
    rdir = copy.copy(dir_rst['mtrx'])
    # Borrar fronteras
    rdir[0,:]  = dir_rst['nodt']
    rdir[:,0]  = dir_rst['nodt']
    rdir[-1,:] = dir_rst['nodt']
    rdir[:,-1] = dir_rst['nodt']
    # Definir algunas variables locales
    cont1  = 0
    cont2  = 0
    sumijs = [[-1,-1,-1,0,0,0,1,1,1],[-1,0,1,-1,0,1,-1,0,1]]
    
    # Array con las direcciones de flujo
    dirHS = np.array([3,2,1,6,0,4,9,8,7])
    fvec  = [[0],[0]]
    
    i,j   = XYtoij(X,Y,raster=dir_rst)

    ijs   = np.array([[i],[j]],dtype=np.int64)
    
    # Se encuentra la respectiva fila y columna en matriz ráster para las coordenadas de cabecera
    if len(XY_ini) > 0 :
        rr = np.array([XYtoij(XY_ini[k,0],XY_ini[k,1],raster=dir_rst) for k in range(len(XY_ini))])
        
    # Se seleccionan las celdas dentro de la cuenca
    while cont2 >= cont1:
        # Coordenadas de celda a evaluar, inicia con el punto de cierre
        i = ijs[0,cont1]
        j = ijs[1,cont1]
        # Se saltan las celdas que son cabecera
        if len(XY_ini)>0 :
            while np.any(np.all(rr == [i,j], axis=1)):
                print(i,j) 
                cont1 = cont1 + 1
                i     = ijs[0,cont1]
                j     = ijs[1,cont1]  
        # Para evaluar qué celdas drenan a la celda actual
        newijs = np.array([(i+l,j+k) for l,k in zip(sumijs[0],sumijs[1])]).T
        aux1   = rdir[newijs[0,:],newijs[1,:]]
        sidren = np.equal(aux1,dirHS)
        # Se pone un ID diferente a cada celda
        newid  = (sidren.cumsum()+cont2)*sidren
        newijs = newijs[:,sidren]
        fvec[0].extend(newid[sidren])
        fvec[1].extend(np.zeros(newid[sidren].shape, dtype=np.int64)+cont1)
        cont2  = fvec[0][-1] # Número de ID máximo asignado
        ijs    = np.append(ijs,newijs,axis=1)
        cont1  = cont1+1     # Número de celda en la que va evaluando los drenajes 
    
    # Se cambia orden
    fvec       = np.array(fvec).T
    fvec       = cont2-fvec
    fvec       = fvec[::-1,:]
    ijs        = list(ijs[:,::-1])
    if return_raster == True or return_mask == True:
        rID  = np.ones((np.shape(rdir)))*dir_rst['nodt']
        rIDD = np.ones((np.shape(rdir)))*dir_rst['nodt']
        rID[ijs[0],ijs[1]]  = fvec[:,0]
        rIDD[ijs[0],ijs[1]] = fvec[:,1]
        # Crear raster donde se almacenarán los nuevos ID
        rID  = BuildRaster(dir_rst['xll'],dir_rst['yll'],dir_rst['clsz'],dir_rst['nodt'],rID)
        rIDD = BuildRaster(dir_rst['xll'],dir_rst['yll'],dir_rst['clsz'],dir_rst['nodt'],rIDD)
        # Se crea máscara de la cuenca
        mask                 = rID.copy()
        mask['mtrx']         = np.copy(rID['mtrx'])
        mask['mtrx'][mask['mtrx']>=0] = 1
        if return_mask == True and return_raster == False:
            return mask
        elif return_mask == False and return_raster == True:
            return rID, rIDD
        elif return_mask == True and return_raster == True:
            return mask, rID, rIDD
    else:
        fvec_print        = pd.DataFrame() 
        fvec_print['ID']  = fvec[:,0]
        fvec_print['IDD'] = fvec[:,1]
        return fvec_print, ijs

#==============================================================================
#--Crear vector de flujo con los atributos de la cuenca
#==============================================================================
def FlowVector(dem_rst,dir_rst,areas_rst,bsnmask_rst,Uarea,xyA_ini=''):

    u"""
    Construye  DataFrame que contiene información para todas las celdas dentro de la cuenca


    Entradas
    ----------
    dem_rst: RasterObject_like
             raster del modelo de elevación digital

    dir_rst: RasterObject_like
             raster de dirección de flujo

    areas_rst: RasterObject_like
               raster de área de flujo acumulada [km^2]

    bsnmask_rst: RasterObject_like
                 mácara de cuenca

    Uarea: float_like
           Umbral área para definición de la red de drenaje

    Salidas
    -------
    FlowVector : DataFrame
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

    AuxReajuste: numpy array
                Array con lista de índices para obtener mapas ráster

    """

    # Se llenan los bordes de los rásters con faltantes
    NullsInBounds(dem_rst)
    NullsInBounds(dir_rst)
    NullsInBounds(areas_rst)
    NullsInBounds(bsnmask_rst)

    # Se crean las matrices con los datos necesarios
    RS_dir=np.copy(dir_rst['mtrx'])
    RS_area=np.copy(areas_rst['mtrx'])
    RS_dem=np.copy(dem_rst['mtrx'])
    RS_bsnmask=np.copy(bsnmask_rst['mtrx'])

    AuxReajuste=np.where(np.reshape(RS_bsnmask,RS_bsnmask.size)>=0.0)

    # Se reemplazan los ceros del ráster de direcciones con faltantes
    RS_dir[RS_dir==0]=-9999
    
    # Se crea el Dataframe de flujo usando únicamente las celdas dentro de la cuenca
    celdas=np.where(RS_bsnmask!=-9999.)

    n             = celdas[0].shape[0]            # Número de celdas en la cuenca
    RS_ID         = np.copy(bsnmask_rst['mtrx']) 
    RS_ID[celdas] = np.arange(n,dtype=np.float32) # Se pone un ID para cada celda dentro de la cuenca

    FlowVectorA   = pd.DataFrame()                # Se crea DataFrame

    # Aux para longitud no acumulada
    aux_long = np.array([np.sqrt(2),1.,np.sqrt(2),1.,0.,1.,np.sqrt(2),1.,np.sqrt(2)])
    aux_row  = np.array([1,1,1,0,0,0,-1,-1,-1])
    aux_col  = np.array([-1,0,1,-1,0,1,-1,0,1])

    # Se crean los vectores de filas y columnas con información del índice de celda a la que drena la celda respectiva 
    des_row  = celdas[0] + aux_row[RS_dir[celdas].astype(int)-1]
    des_col  = celdas[1] + aux_col[RS_dir[celdas].astype(int)-1]


    #Se inicia la construcción de DataFrame con atributos de la cuenca
    FlowVectorA['ID']  = np.arange(n,dtype=np.float32) #ID de la celda
    FlowVectorA['IDD'] = RS_ID[des_row,des_col] #IDD: ID de celda a la que drena

    FlowVectorA['Area_Acum']          = RS_area[celdas] #Area
    FlowVectorA['Lon_NoAcum']         = aux_long[RS_dir[celdas].astype(int)-1]*dem_rst['clsz'] #Long
    FlowVectorA['Elevacion']          = RS_dem[celdas] #cota inicial
    FlowVectorA['X'],FlowVectorA['Y'] = IJtoxy(celdas[0],celdas[1],dem_rst) # X, y
    FlowVectorA['i']                  = celdas[0] # i
    FlowVectorA['j']                  = celdas[1] # j
    FlowVectorA['Pendiente']          = (RS_dem[celdas]-RS_dem[des_row,des_col])/FlowVectorA['Lon_NoAcum'] #SlopeMax

    # Se corrigen las áreas en cuencas incompletas
    if len(xyA_ini) >0:
        FlowVectorA=AreasFit(FlowVectorA,xyA_ini,dir_rst,areas_rst)
        print(u'Se corrigieron las áreas de la cuenca incompleta')

    # Se ordena de acuerdo al área de drenaje acumulada 
    FlowVectorA = FlowVectorA.sort_values(by=['Area_Acum'],ignore_index=True)

    # Se cambia el index de la columna ID, así la celda con ID 0 es la que menor área acumulada tiene 
    orig_index         = FlowVectorA['ID'].values.argsort()
    search             = np.searchsorted(np.arange(0,n),FlowVectorA['IDD'].values)
    search[search<0]   = 0; search[search>=n] = 0
    new_index          = orig_index[search]
    outlets            = FlowVectorA['ID'].iloc[new_index].values != FlowVectorA['IDD'].values
    FlowVectorA['ID']  = np.arange(0,n)
    new_index[outlets] = FlowVectorA.loc[outlets,'ID'].values
    FlowVectorA['IDD'] = new_index

    # En la salida de la cuenca, la celda se drena a sí misma: IDD=ID
    FlowVectorA['IDD'].iloc[-1] = len(FlowVectorA['IDD']) -1
    
    #define red de drenaje con umbral de area
    FlowVectorA['Tipo']                                  = np.zeros((len(FlowVectorA),1))
    FlowVectorA.loc[FlowVectorA['Area_Acum'] >= Uarea,'Tipo'] = 1.0  # 0:ladera, 1:cauce
    FlowVectorA['Nodos']                                   = np.zeros((len(FlowVectorA),1))
    FlowVectorA.loc[FlowVectorA['Tipo'] != 1.0,'Nodos']    = -9999.0 # columna [11]=0 donde es red
    
    return (FlowVectorA,AuxReajuste)


#==============================================================================
#--Corrección de pendientes negativas
#==============================================================================

def SlopeFit(FlowVectorA):
    u"""
    Corrige pendientes negativas
    
    Entradas
    ----------
    FlowVector : DataFrame
                 Con atributos de la cuenca:
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

        
    Salidas
    -------
    FlowVector : DataFrame
                Columna 'Pendiente' corregida
    """
    # Seleccionar los pixeles con pendiente negativa
    aux1 = FlowVectorA.index[FlowVectorA['Pendiente']<=0.]
    for i in aux1[::-1]:
        # poner la pendiente de la celda de destino
        FlowVectorA.loc[i,'Pendiente'] = FlowVectorA.loc[FlowVectorA.loc[i,'IDD'],'Pendiente']
    # Si todavía existen pendientes negativas se reemplaza la pendiente por un valor bajo
    aux1 = FlowVectorA.index[FlowVectorA['Pendiente']<=0.]
    for i in aux1:
        FlowVectorA['Pendiente'][i] = 0.001
    
    return FlowVectorA


#==============================================================================
#--Trazar cabeceras desde Shapefile
#==============================================================================
def TraceRiverheads(FlowVectorA,rivhead_pth,bsnmask_rst):
    u"""
    Agregar redes de drenaje basadas en un shapefile predefinido
   
    Entradas
    ----------
    FlowVector : DataFrame
                 Con atributos de la cuenca:
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

    
    rivhead_pth : str
                  Directorio donde está el shapefile con las cabeceras
    
    bsnmask_rst : RasterObject_like
                  máscara donde se encuentra la cuenca
        
    Salidas
    -------
    FlowVector : DataFrame
        Correcciones en columnas 'Tipo' y 'Nodos'
    """
    RS_bsnmask = np.copy(bsnmask_rst['mtrx'])
    
    # Crea un ID para cada celda en bsnmask 
    celdas = np.where(RS_bsnmask!=-9999.)
    
    n             = celdas[0].shape[0]
    RS_id         = np.copy(bsnmask_rst['mtrx'])
    RS_id[celdas] = np.arange(n,dtype=np.float32)
    
    # Lee cabeceras desde .shp y define ID
    shape       = gpd.read_file(rivhead_pth)
    shape['x']  = shape.geometry.x
    shape['y']  = shape.geometry.y
    shape['i']  = np.nan
    shape['j']  = np.nan
    shape['ID'] = np.nan
    
    imax = np.shape(RS_bsnmask)[0]
    jmax = np.shape(RS_bsnmask)[1]
    for nn in range(len(shape)):
        x = shape.x[nn]
        y = shape.y[nn]
        i,j = XYtoij(X=x, Y=y, raster=bsnmask_rst)
        if i>=0 and i<imax:
            if j>=0 and j<jmax:
                if RS_bsnmask[i,j] >= 0:
                    shape.loc[nn,'i']  = i
                    shape.loc[nn,'j']  = j
                    shape.loc[nn,'ID'] = FlowVectorA['ID'][(FlowVectorA['i']==i)&(FlowVectorA['j']==j)].iloc[0]
    
    shape = shape.dropna()
    
    # Selección de las cabeceras
    cab = shape.ID.astype(int).values.T
    cab.sort(axis=0)
    
    
    for p in range(len(cab)):
        celda=int(cab[p])
    
        # Traza drenaje desde las cabeceras hasta las salidas y define las celdas como cauce en FlowVector
        while FlowVectorA['ID'][celda] != FlowVectorA['IDD'][celda]:
            if FlowVectorA['Tipo'][celda] == 1.0:
                break
            else:
                FlowVectorA.loc[celda,'Tipo']  = 1.0
                FlowVectorA.loc[celda,'Nodos'] = 0.0
                celda = FlowVectorA.loc[celda,'IDD'].astype(int)
    return (FlowVectorA)

#==============================================================================
#--Cabeceras
#==============================================================================
def IdCab(FlowVectorA):

    u"""
    Identifica cabeceras de la red de drenaje definidas con el umbral de área

    Entradas
    ----------
    FlowVector : DataFrame
                 Con atributos de la cuenca

    Salidas
    -------
    FlowVector : DataFrame
                 Definición de cabeceras en columna 'Nodos', estas tienen un índice de 1
    """
    # ID de celdas definidas como red de drenaje
    iddren = FlowVectorA.index[FlowVectorA['Tipo']==1.0]
    # selecciona las celdas de cauce que no son drenadas por ninguna otra
    cabeceras = np.setdiff1d(FlowVectorA.loc[iddren,'ID'],FlowVectorA.loc[iddren,'IDD']).astype(int)
    FlowVectorA.loc[cabeceras,'Nodos'] = 1.0
    return (FlowVectorA)

#==============================================================================
#--Nodos Hidrológicos
#==============================================================================
def HydrologicNodes(FlowVectorA):

    u"""
    Identificar nodos de confluencia hidrológica

    Entradas
    ----------
    FlowVector : DataFrame
                 Con atributos de la cuenca:
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

    Salidas
    -------
    FlowVector : DataFrame
                 Llena columna 'Nodos' de FlowVector con 2 donde hay confluencias
    """
    # id de celdas que hacen parte de la red de drenaje
    iddren = FlowVectorA.index[FlowVectorA['Tipo']==1.0]
    
    # id de las celdas a las que drenan celdas de red,
    # y número de celdas tipo red que drenan a ellas.
    aux1,aux2 = np.unique(FlowVectorA.loc[iddren,'IDD'], return_counts=True)
    
    # id de las celdas a las que drenan más de 1 celda de red
    aux3 = aux1[aux2>1].astype(int)
    
    # marca dichas celdas como confluencia o nodo hidrológico
    FlowVectorA.loc[aux3,'Nodos'] = 2.0
    
    return (FlowVectorA)

#==============================================================================
#--Umbral de longitud
#==============================================================================
def UmbralLong(FlowVectorA,Ulong):

    u"""
    Corrección de la red de drenaje, ignorando corrientes de primer orden
    (Strahler-Horton) con longitudes menor que Ulen.

    Entradas
    ----------
    FlowVector : DataFrame
                 Con atributos de las celdas de la cuenca

    Ulong : Float
            Umbral de longitud

    Salidas
    -------
    FlowVector :  DataFrame
                  Cambia columna 'Tipo' creando nuevas celdas de ladera 
                  y cambia columna 'Nodos' removiendo los nodos hidrológicos correspondientes

    """
    # id de celdas marcadas como cabeceras
    idcabe = FlowVectorA.index[FlowVectorA['Nodos']==1.0]
    
    for p in idcabe:
        # marca la cabecera como celda de inicio y crea un array para almacenar
        # el id de cada celda del tramo, la longitud acumulada, y si es nodo
        celda = int(p)
        Tramo = np.array([[p,FlowVectorA.loc[p,'Lon_NoAcum'],FlowVectorA.loc[p,'Nodos']]])#; d=0

        
        while Tramo[-1,2] != 2.0:
            # si todavía no se ha encontrado un nodo hidrológico
            celda = int(FlowVectorA.loc[celda,'IDD'])# Calcula celda siguiente
            
            # adiciona las propiedades de la nueva celda
            nueva_celda = np.array([[ celda , Tramo[-1,1] + FlowVectorA.loc[celda,'Lon_NoAcum'] , FlowVectorA.loc[celda,'Nodos']]])
            Tramo=np.append(Tramo, nueva_celda,axis=0)
            
        if Tramo[-1,1]<Ulong:
            # marca todas las celdas del tramo como ladera
            FlowVectorA.loc[Tramo[:-1,0].astype(int),'Tipo'] = 0.0
            # quita todos los nodos en las celdas del tramo
            FlowVectorA.loc[Tramo[:-1,0].astype(int),'Nodos']  = -9999.0
            
    return (FlowVectorA)

#==============================================================================
#--Umbral de longitud por percentil
#==============================================================================
def PercUmbralLong(FlowVectorA,bsnmask_rst,Ulong=100.0,percent=75.0,sbn_pth=''):
    u"""
    Corrección de la red de drenaje, ignorando corrientes de primer orden
    (Strahler-Horton) con longitud menor que un percentil de todas las
    longitudes de las corrientes de primer orden definido por el usuario,
    manteniendo las corrientes que intersectan algún área de interés definida
    en el ráster.

    Entradas
    ----------
    FlowVector : DataFrame
                 Con atributos de las celdas de la cuenca
    bsnmask_rst : Raster
                  máscara donde está la cuenca
    Ulong : Float
            Umbral de longitud
    percent : Float
            Percentil umbral para las longitudes de primer orden
    sbn_pth : str
            Directorio del ráster de áreas de interés, estas están definidas por 1        

    Salidas
    -------
    FlowVector :  DataFrame
                Cambio a la columna 'Tipo': se crean nuevas celdas de ladera, 
                y cambios en la columna 'Nodos': se remueven las correspondientes celdas de cauce
    """
   
    if sbn_pth != '':
        sbn_rst=ReadRaster(sbn_pth)
        RS_sbn=sbn_rst['mtrx']
    else:
        RS_sbn=bsnmask_rst['mtrx']*0.0
    
    # id de celdas marcadas como cabeceras
    idcabe      = FlowVectorA.index[FlowVectorA['Nodos']==1.0]
    
    lenghts     = []
    for p in idcabe:
        # marca la cabecera como celda de inicio y crea un array para almacenar
        # el id de cada celda del tramo, la longitud acumulada, y si es nodo
        celda = int(p)
        Tramo = np.array([[p,FlowVectorA.loc[p,'Lon_NoAcum'],FlowVectorA.loc[p,'Nodos']]])#; d=0
        sbn_check = 0
        while Tramo[-1,2] != 2.0:
            # se verifica que el tramo no cruce una zona SBN
            i = int(FlowVectorA.loc[p,'i'])
            j = int(FlowVectorA.loc[p,'j'])
            if RS_sbn[i,j] > 0.0:
                sbn_check = 1
    
            # si todavía no se ha encontrado un nodo hidrológico
            celda = int(FlowVectorA.loc[celda,'IDD'])# Calcula celda siguiente
            
            # adiciona las propiedades de la nueva celda
            nueva_celda = np.array([[ celda , Tramo[-1,1] + FlowVectorA.loc[celda,'Lon_NoAcum'],FlowVectorA.loc[celda,'Nodos']]])
            Tramo=np.append(Tramo, nueva_celda,axis=0)
                  
        if percent is None:
            if sbn_check == 0:
                if Tramo[-1,1]<Ulong:
                    # marca todas las celdas del tramo como ladera
                    FlowVectorA.loc[Tramo[:-1,0].astype(int),'Tipo'] = 0.0
                    # quita todos los nodos en las celdas del tramo
                    FlowVectorA.loc[Tramo[:-1,0].astype(int),'Nodos']  = -9999.0

                    
        if len(Tramo)<3:
            # marca todas las celdas del tramo como ladera
            FlowVectorA.loc[Tramo[:-1,0].astype(int),'Tipo'] = 0.0
            # quita todos los nodos en las celdas del tramo
            FlowVectorA.loc[Tramo[:-1,0].astype(int),'Nodos']  = -9999.0
                        
        # se agregan todas las longiudes de tramo orden 1
        lenghts.append(Tramo[-1,1])
    
    if percent is not None:
        lenghts = np.array(lenghts)
        Uperc = np.percentile(lenghts,percent)
        
        for p in idcabe:
            # marca la cabecera como celda de inicio y crea un array para almacenar
            # el id de cada celda del tramo, la longitud acumulada, y si es nodo
            celda     = int(p)
            Tramo     = np.array([[p,FlowVectorA.loc[p,'Lon_NoAcum'],FlowVectorA.loc[p,'Nodos']]])
            sbn_check = 0
            while Tramo[-1,2] != 2.0:
                # se verifica que el tramo no cruce una zona SBN
                i = int(FlowVectorA.loc[p,'i'])
                j = int(FlowVectorA.loc[p,'j'])
                if RS_sbn[i,j] > 0.0:
                    sbn_check = 1
        
                # si todavía no se ha encontrado un nodo hidrológico
                celda = int(FlowVectorA.loc[celda,'IDD'])# Calcula celda siguiente
                
                # adiciona las propiedades de la nueva celda
                nueva_celda = np.array([[ celda , Tramo[-1,1] + FlowVectorA.loc[celda,'Lon_NoAcum'] , FlowVectorA.loc[celda,'Nodos']]])
                Tramo       = np.append(Tramo, nueva_celda,axis=0)
    
            if sbn_check == 0:
                if Tramo[-1,1]<Uperc:
                    # marca todas las celdas del tramo como ladera
                    FlowVectorA.loc[Tramo[:-1,0].astype(int),'Tipo'] = 0.0
                    # quita todos los nodos en las celdas del tramo
                    FlowVectorA.loc[Tramo[:-1,0].astype(int),'Nodos']  = -9999.0
            
            if len(Tramo)<3:
                # marca todas las celdas del tramo como ladera
                FlowVectorA.loc[Tramo[:-1,0].astype(int),'Tipo'] = 0.0
                # quita todos los nodos en las celdas del tramo
                FlowVectorA.loc[Tramo[:-1,0].astype(int),'Nodos']  = -9999.0
    
    return (FlowVectorA)

#==============================================================================
#--Nodos donde drenajes efímeros se vuelven permanentes
#==============================================================================
def KnEphtoPerm(FlowVectorA,Uarea,Uareaper):
    u"""
    Define nodos topográficos, donde el drenaje efímero se convierte en permanente

    Entradas
    ----------
    FlowVector : DataFrame
                 con atributos de la cuenca

    Uareaper: float
                Umbral de área que define drenaje permanente

    Uarea: float
                Umbral de área que define drenaje efímero

    Salidas
    -------
    FlowVector :  DataFrame
                con actualización en columna 'Nodos'

    """
    Aux  = copy.copy(FlowVectorA['Tipo'])
    Aux2 = copy.copy(FlowVectorA['Nodos'])
    Aux[np.logical_and(FlowVectorA['Area_Acum']>Uarea,FlowVectorA['Area_Acum']<Uareaper)] = 2
    #FlowVectorA[np.where(FlowVectorA[:,10]==1)[0],11]=0 # columna [11]=0 donde es red
    iddren = np.where(Aux==1)[0]
    Aux2[np.setdiff1d(FlowVectorA.loc[iddren,'ID'],FlowVectorA.loc[iddren,'IDD']).astype(int)] = 3
    FlowVectorA.loc[np.logical_and(Aux2==3,FlowVectorA['Nodos']==0),'Nodos']    = 3
    return FlowVectorA

#==============================================================================
#-- Nodos Topográficos
#==============================================================================
def KnickNodes(FlowVectorA,NoCllTram,ventana,peso_media,tol,output_pth=''):

    u"""
    Suaviza el perfil de flujo y calcula los nodos topográficos

    Entradas
    ----------
    FlowVector : DataFrame con los atributos de cuenca

    Salidas
    -------
    FlowVector :  numpy array
                [12] New corrected elevation in the drainage network
                [11] New Topographic nodes

    """
    FlowVectorA['Elev_corr'] = -9999.0
    
    # Parametros para el lienzo de gráficos
    if output_pth != '':
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        rcParams['font.family']='sans-serif'
        rcParams['font.sans-serif']=['Arial']
        rcParams['font.size']=16
        plt.ioff()

    if output_pth != '':
        fig=plt.figure(facecolor='w',edgecolor='w',figsize=(10,10))
    
    # Identificacion de cabeceras
    cabeceras = FlowVectorA.loc[FlowVectorA['Nodos']==1,'ID']
    
    for p,cabecera in enumerate(cabeceras.astype(int),1):
        # Iniciar arreglo para guardar propiedades de tramos
        tramo=[]
        abscisa=0.0
        celda=cabecera
        
        # Traza una corriente desde la cabecera hasta la salida de la cuenca
        while FlowVectorA.loc[celda,'ID'] != FlowVectorA.loc[celda,'IDD']:
            tramo.append([celda,abscisa,FlowVectorA.loc[celda,'Elevacion']]) #Celda inicio/log acum/cota
            abscisa = abscisa + FlowVectorA.loc[celda,'Lon_NoAcum']
            celda   = FlowVectorA.loc[celda,'IDD'].astype(int)
        tramo.append([celda,abscisa,FlowVectorA.loc[celda,'Elevacion']])
        tramo=np.array(tramo)

        #- Extensión del perfil hacia aguas arriba y aguas abajo
        if np.mod(ventana,2) == 0:
            # El valor de la ventana debe ser impar
            ventana = ventana + 1
        media_ventana = int((ventana - 1 ) / 2)# la mitad de la ventana
        
        # Arreglo para guardar la información de elevación
        cotas = np.zeros(np.size(tramo, 0) + 2 * media_ventana)
        
        # Pendiente e intercepto al inicio del tramo
        rec_pen1 = (tramo[1,2] - tramo[0,2]) / (tramo[1,1] - tramo[0,1])
        rec_int1 = tramo[0,2] - rec_pen1 * tramo[0,1]
        
        # Pendiente e intercepto al final del tramo
        rec_pen2 = (tramo[-1,2] - tramo[-2,2]) / (tramo[-1,1]-tramo[-2,1])
        rec_int2 = tramo[-1,2] - rec_pen2 * tramo[-1,1]
        
        # Extender el perfil adelante en el final del tramo y atrás al 
        # inicio del tramo, en un tamaño igual a media ventana
        for i in range(0,media_ventana):
            cotas[media_ventana - i - 1] = rec_pen1 * (tramo[0,1] - (i + 1) * tramo[1,1]) + rec_int1
            cotas[-(i + 1)] = rec_pen2 * (tramo[-1,1] + (media_ventana - i) * (tramo[-1,1] - tramo[-2,1])) + rec_int2
        cotas[media_ventana : - media_ventana] = tramo[:,2]

        # Filtro de suavizado
        filtro_media = np.convolve(cotas, np.ones(ventana)/ventana, mode='valid')
        filtro_mediana = signal.medfilt(cotas, ventana)[media_ventana:-media_ventana]
        
        # guardar cota suavizada en arreglo de tramo
        tramo      = np.insert(tramo, 3, np.NaN, axis=1)
        tramo[:,3] = peso_media*filtro_media + (1.0 - peso_media) * filtro_mediana

        # Filtro monotono decreciente
        AuxCota = tramo[0,3]# cota de inicio de tramo
        for i, cota in enumerate(tramo[1:,3], 1):
            if cota < AuxCota:
                # si la pendiente es negativa: OK, cambie la altura de referencia
                AuxCota = cota
            else:
                # si la pendiente es positiva: NO, borre la cota para recalcularla
                tramo[i,3] = -9999.0
        
        if tramo[-1,3] == -9999.0:
            # Recalcule elevación de la celda final si tiene pendiente positiva
            tramo[-1,3] = AuxCota - 0.001
        
        # Hace interpolación lineal para reclacular elevaciones de 
        # zonas con pendientes positivas
        IdTramoNan = np.where(tramo[:,3] == -9999.0)[0]
        for tnan in IdTramoNan:
            down          = np.where(~(tramo[tnan:,3] == -9999))[0][0] + tnan
            pend_local    = (tramo[tnan - 1, 3] - tramo[down,3]) / (tramo[down,1] - tramo[tnan - 1, 1])
            tramo[tnan,3] = tramo[down,3] + ((tramo[down, 1] - tramo[tnan, 1]) * pend_local)

        # Reemplaza el perfil corregido en el vector de flujo
        FlowVectorA.loc[tramo[:,0].astype(int),'Elev_corr']=tramo[:,3]

        # Encuentra los nodos topográficos
        if np.size(tramo,0)>=5: # Para evaluar la segunda derivada se necesitan 5 puntos

            # Extrae propiedades del tramo para graficar posteriormente
            Ltramo = np.round(tramo[-1,1]/1000.0,3)
            dH     = np.round(tramo[0,3]-tramo[-1,3],3)
            Xini   = int(FlowVectorA.loc[int(tramo[0,0]),'X']); Yini=int(FlowVectorA.loc[int(tramo[0,0]),'Y'])
            Xfin   = int(FlowVectorA.loc[int(tramo[-1,0]),'X']); Yfin=int(FlowVectorA.loc[int(tramo[-1,0]),'Y'])

            # Normaliza el perfil de flujo
            tramo[:,1]=tramo[:,1]/tramo[:,1].max()
            cota_min=tramo[:,3].min()
            cota_ran=tramo[:,3].max()-cota_min
            tramo[:,3]=(tramo[:,3]-cota_min)/cota_ran
            tramo[:,2]=(tramo[:,2]-cota_min)/cota_ran
#            tol=tol*np.std(np.abs(tramo[:,2]-(filtro_media-cota_min)/cota_ran))
#            tol=tol*np.std(np.abs(tramo[:,2]-filtro_media))

            # Filtro de banda digital
            nuevo_tramo=[]
            nuevo_tramo.append([tramo[0,1],tramo[0,3],0])
            for i,cota in enumerate(tramo[1:,3],1):
                mask=np.arange(nuevo_tramo[-1][2],i+1)
                # Recta entre el primer y el último punto del tramo
                rec_pen=(cota-nuevo_tramo[-1][1])/(tramo[i,1]-nuevo_tramo[-1][0])
                rec_int=nuevo_tramo[-1][1]-rec_pen*nuevo_tramo[-1][0]
                # Recta de regresión en mínimos cuadrados
#                rec_pen,rec_int=np.linalg.lstsq(np.vstack([tramo[mask,1],np.ones(len(tramo[mask,1]))]).T,tramo[mask,3])[0]
                lim_inf=rec_pen*tramo[mask,1]+rec_int-tol
                lim_sup=rec_pen*tramo[mask,1]+rec_int+tol
                if ~np.all(np.logical_and(tramo[mask,3]>=lim_inf,tramo[mask,3]<=lim_sup)):
                    nuevo_tramo.append([tramo[i,1],tramo[i,3],i])
            if nuevo_tramo[-1][2]!=i:
                nuevo_tramo.append([tramo[-1,1],tramo[-1,3],i])
            nuevo_tramo=np.array(nuevo_tramo)

            # Interpolacion lineal del perfil filtrado
            tramo=np.insert(tramo,4,tramo[:,3],axis=1)
            for i in range(1,nuevo_tramo.shape[0]):
                rec_pen=(nuevo_tramo[i-1,1]-nuevo_tramo[i,1])/(nuevo_tramo[i-1,0]-nuevo_tramo[i,0])
                rec_int=nuevo_tramo[i-1,1]-rec_pen*nuevo_tramo[i-1,0]
                tramo[int(nuevo_tramo[i-1,2]):int(nuevo_tramo[i,2]),3]=\
                      rec_pen*tramo[int(nuevo_tramo[i-1,2]):int(nuevo_tramo[i,2]),1]+rec_int

            # Calcula la primera y la segunda derivada del perfil corregido
            dx_adela=np.roll(tramo[:,1],-1)-tramo[:,1]
            dx_atras=tramo[:,1]-np.roll(tramo[:,1],1)
            alfa=dx_adela/(dx_adela+dx_atras)
            prim_deri=(1.0-alfa)*(tramo[:,3]-np.roll(tramo[:,3],1))/dx_atras+\
                           alfa*(np.roll(tramo[:,3],-1)-tramo[:,3])/dx_adela
            prim_deri[0]=(tramo[1,3]-tramo[0,3])/(tramo[1,1]-tramo[0,1])
            prim_deri[-1]=(tramo[-1,3]-tramo[-2,3])/(tramo[-1,1]-tramo[-2,1])
            segu_deri=(1.0-alfa)*(prim_deri-np.roll(prim_deri,1))/dx_atras+\
                           alfa*(np.roll(prim_deri,-1)-prim_deri)/dx_adela
            segu_deri[0]=(prim_deri[1]-prim_deri[0])/(tramo[1,1]-tramo[0,1])
            segu_deri[-1]=(prim_deri[-1]-prim_deri[-2])/(tramo[-1,1]-tramo[-2,1])
#            curvatura=segu_deri/(1+segu_deri**2.0)**(3/2)
            terce_deri=(1.0-alfa)*(segu_deri-np.roll(segu_deri,1))/dx_atras+\
                           alfa*(np.roll(segu_deri,-1)-segu_deri)/dx_adela
            terce_deri[0]=(segu_deri[1]-segu_deri[0])/(tramo[1,1]-tramo[0,1])
            terce_deri[-1]=(segu_deri[-1]-segu_deri[-2])/(tramo[-1,1]-tramo[-2,1])
            cambio_signo=terce_deri*np.roll(terce_deri,-1)
            cambio_signo[-1]=1.0
            segu_deri[np.where(cambio_signo>0.0)[0]]=0.0
            segu_deri[0]=0.0

            # Encuentra los puntos de inflexion en el perfil de flujo
            lim_inf=-50.0#segu_deri.mean()-2.0*segu_deri.std()
            lim_sup=50.0#segu_deri.mean()+2.0*segu_deri.std()
            puntos_inflex=np.where(np.logical_or(segu_deri>lim_sup,segu_deri<lim_inf))[0]

            # Depura y asigna los nodos topográficos al arreglo de flujo
            for i in puntos_inflex.astype(int):
                j_ar = max(i-ventana,0)
                j_ab = min(i+ventana,np.size(tramo,0)-1)+1
                if np.max(FlowVectorA.loc[tramo[j_ar:j_ab,0].astype(int),'Nodos'])<=0.0:
                    FlowVectorA.loc[tramo[i,0].astype(int),'Nodos']=3.0

            # Gráfico
            if output_pth != '':
                ax1=fig.add_subplot(211)
                ax1.set_xlabel(u'Abscisa [adim.]'); ax1.set_ylabel(u'Cota [adim.]')
                ax1.plot(tramo[:,1],tramo[:,2],'b',lw=0.8)
                ax1.plot(tramo[:,1],tramo[:,3],'g',lw=0.8)
                ax1.plot(tramo[:,1],tramo[:,4],'m',lw=0.8)
                ax2=fig.add_subplot(212)
                ax1.text(0.65,0.95,u'$L_{Tramo}$ = '+str(Ltramo)+' km',fontsize=14)
                ax1.text(0.65,0.85,u'$\Delta$H = '+str(dH)+' m ',fontsize=14)
                ax1.text(0.65,0.75,u'$(x,y)_{Ini}$ = ('+str(Xini)+', '+str(Yini)+')',fontsize=14)
                ax1.text(0.65,0.65,u'$(x,y)_{Fin}$ = ('+str(Xfin)+', '+str(Yfin)+')',fontsize=14)
                ax2.plot(tramo[:,1],segu_deri,'b',lw=0.8)
                ax2.set_xlabel(u'Abscisa [adim.]'); ax2.set_ylabel(u'Segunda derivada')
                for celda in puntos_inflex:
                    x=tramo[celda,1]; y=segu_deri[celda]
                    ax1.axvline(x=x,c='g',linewidth=0.5,zorder=0, clip_on=False)
    #                ax2.axvline(x=x,c='g',linewidth=0.5, zorder=0,clip_on=False)
                    ax2.plot(x,y,'ro',markersize=0.5)
                ax2.axhline(y=segu_deri.mean(),c='g',linewidth=0.5,zorder=0,clip_on=False)
                ax2.axhline(y=lim_inf,c='g',linewidth=0.5,zorder=0, clip_on=False)
                ax2.axhline(y=lim_sup,c='g',linewidth=0.5,zorder=0, clip_on=False)
                num_perfil=str(p)
                for cero in range(len(str(np.size(cabeceras,0)))-len(num_perfil)):
                    num_perfil='0'+num_perfil
                fig.savefig(os.path.join(output_pth,'perfil_'+num_perfil+'.png'),dpi=150)
                plt.clf()
    if output_pth != '':
        plt.close(fig)
    return FlowVectorA

#==============================================================================
#--Define celdas de drenaje efímero
#==============================================================================

def Carcavas(FlowVectorA,Uarea,Uareaper):
    u"""
    Define celdas de drenaje efímero

    Entradas
    ----------
    FlowVector : DataFrame
                 Con atributos de cuenca

    Uareaper:    float
                 Umbral de Área que define drenajes permanentes

    Uarea:       float
                 Umbral de Área que define drenajes efímeros
    Salidas
    -------
    FlowVector :  DataFrame
                  Se actualizan drenajes efímeros en columna 'Tipo'

    """
    FlowVectorA.loc[np.logical_and(FlowVectorA['Tipo']==1,np.logical_and(FlowVectorA['Area_Acum']>Uarea,FlowVectorA['Area_Acum']<Uareaper)),'Tipo'] = 2.0
    return FlowVectorA

#==============================================================================
#--Get Link properties
#==============================================================================
def fillDic(x,dic):

    u"""
    Función para llenar diccionario usando una tupla

    Entradas
    ----------
    x: Tuple [1] ->key [0]->valor
    dic: Diccionario relacionado con la tupla

    Salidas
    -------
    aux: valor del diccionario correspondiente a la llave
    """
    auxdic = dic[x[1]]=x[0]
    return auxdic

def GetReaches(FlowVectorA):

    u"""
    Enumerate each link after having the consolidated hydrologyc and topographic nodes

    Entradas
    ----------
    FlowVector : DataFrame
                 con atributos de la cuenca

    Salidas
    -------
    FlowVector :  DataFrame
                 con cambios en la columna 'Pendiente'
                 con una nueva columna 'ID_tramo' con el ID del tramo al que pertenece la celda
    Tramos:     DataFrame
                con las propiedades de los tramos
                [0] Link ID
                [1] Link length
                [2] Link Slope
                [3] Link mean Basin Area
                [4] Id of destination link
                [5] initial cell Id
                [6] end cell Id
                [7] destination link initial cell id
                [8] corresponding Hillslope Area

    """
    FlowVectorA['ID_tramo'] = -9999.0
    
    #- crear variables para almacenar resultados
    Ntram       = 0
    TramosA     = []
    cabeceras   = FlowVectorA.loc[FlowVectorA['Nodos']==1.,'ID']
    
    for p,cabecera in enumerate(cabeceras.astype(int),1):
        tramo   = []
        abscisa = 0.0
        celda   = cabecera
        #- se recorre todo el tramo hasta la salida de la cuenca
        while FlowVectorA.loc[celda,'ID'] != FlowVectorA.loc[celda,'IDD']:
            abscisa = abscisa + FlowVectorA.loc[celda,'Lon_NoAcum']#- acumula longitudes
            tramo.append([celda,abscisa,FlowVectorA.loc[celda,'Elevacion']]) #- Celda inicio/long acum/cota
            celda   = FlowVectorA.loc[celda,'IDD'].astype(int) #- siguiente celda
        #- para la última celda
        abscisa = abscisa + FlowVectorA.loc[celda,'Lon_NoAcum']
        tramo.append([celda,abscisa,FlowVectorA.loc[celda,'Elevacion']])
        
        Tramo   = np.array(tramo)
        i       = 0
        # se recorren todas las celdas del tramo hasta la salida de la cuenca
        # o hasta encontrar una a la que se le haya asignado ya id de tramo
        while i<len(Tramo[:,0])-1 and FlowVectorA.loc[Tramo[i,0].astype(int),'ID_tramo'] == -9999:
            #- Marcar celda con ID de tramo
            Idcelli = Tramo[i,0].astype(int)
            celli   = i
            FlowVectorA.loc[Idcelli, 'ID_tramo'] = Ntram
            i += 1 #- siguiente celda
            
            #- recorre tramo entre nodos
            #- mientras la celda siguiente no sea nodo y no sea la última
            while FlowVectorA.loc[Tramo[i,0].astype(int),'Nodos']==0 and i<len(Tramo[:,0])-1:
                #- marcar la celda con id de tramo
                FlowVectorA.loc[Tramo[i,0].astype(int),'ID_tramo'] = Ntram
                i += 1 #- siguiente celda
            
            #- si la celda es la última
            if i == len(Tramo[:,0])-1:
                #- marcar con id de tramo
                FlowVectorA.loc[Tramo[i,0].astype(int),'ID_tramo'] = Ntram
            
            #- celda final  del tramo
            Idcellf = Tramo[i - 1,0].astype(int)
            cellf   = i
            #- celda inicial del tramo de destino
            Idcellides = int(FlowVectorA.loc[Idcellf,'IDD'])
            
            #- cálculo del área de ladera aferente
            if FlowVectorA.loc[Idcelli,'Nodos'] == 1:
                # si el inicio es cabecera
                Aladera = FlowVectorA.loc[Idcellf,'Area_Acum']
            else:
                # si es un tramo intermedio
                Aladera = FlowVectorA.loc[Idcellf,'Area_Acum'] - FlowVectorA.loc[Idcelli,'Area_Acum']
            # Longitud acumulada
            Longitud = Tramo[cellf,1] - Tramo[celli,1]
            
            #- calcular pendiente
            if Idcelli == Idcellf:
                #- si es un tramo de una sola celda
                if Idcelli != Idcellides:
                    #- si no es el último tramo
                    s = (FlowVectorA.loc[Idcelli,'Elev_corr']-FlowVectorA.loc[Idcellides,'Elev_corr'])/Longitud
                elif Idcelli == Idcellides:
                    s = 0.001
            else:
                #- para cualquier otro caso
                s = (FlowVectorA.loc[Idcelli,'Elev_corr']-FlowVectorA.loc[Idcellf,'Elev_corr'])/Longitud

            #@@@@@ JPRA: Noté que considerar pendientes con magnitudes inferiores provoca
            #@@@@@ resultados indeseados en las estimaciones de geometría hidráulica, como
            #@@@@@ calados de banca muy elevados. Recursivamente fijé el valor mínimo de
            #@@@@@ pendiente en 0.001. Cualquier sugerencia es bienvenida.
            if s <= 0.001:
                s = 0.001

            Acuenca = (FlowVectorA.loc[Idcellf,'Area_Acum'] + FlowVectorA.loc[Idcelli,'Area_Acum'])/2
            TramosA.append([Ntram,Longitud,s,Acuenca,-9999.0,Idcelli,Idcellf,Idcellides,Aladera])
            Ntram += 1
    TramosA = np.array(TramosA)

    #Cambiar pendiente de las celdas cauce
    mascara = FlowVectorA['ID_tramo'] != -9999.0
    filas   = FlowVectorA.loc[mascara,'ID_tramo'].astype(int)
    if True in mascara:
        FlowVectorA.loc[mascara,'Pendiente'] = TramosA[filas,2]
#    aux=zip(TramosA[:,2],TramosA[:,0])
#    dic={}
#    map(lambda x:fillDic(x,dic),aux)
#    FlowVectorA[np.where(FlowVectorA[:,13]!=-9999.0)[0],9]=np.array(map(lambda x: dic[x],FlowVectorA[np.where(FlowVectorA[:,13]!=-9999)[0],13]))
    if True in mascara:
        TramosA[:,4] = FlowVectorA.loc[TramosA[:,7].astype(int),'ID_tramo']
        
    TramosA = pd.DataFrame(TramosA,columns=['ID','Longitud','Pendiente','Area_cuenca_prom','IDD','ID_i','ID_f','ID_des','Area_ladera'])

    return FlowVectorA,TramosA

#==============================================================================
# Guardar vectores de tramos a shapefile
#==============================================================================
def Reaches2Shp(flowvec,reachesvec,output_pth=''):
    u"""
    Creación de shapefile con los tramos

    Entradas
    ----------
    flowvec :  DataFrame
               con atributos de la cuenca
    reachesvec: DataFrame
                con información de los tramos
    percent : Float
            Umbral de percentil de las longitudes de tramos de primer ordén
    output_pth : str
            Directorio donde se guardará "reaches.shp"       

    """
    
    fshape = gpd.GeoDataFrame(reachesvec)
    
    
    ids   = fshape['ID'].values
    celsi = fshape['ID_i'].values
    celsf = fshape['ID_des'].values
    
    geom = []
    for nn in range(len(ids)):
        celi = int(celsi[nn])
        celf = int(celsf[nn])
        
        celda = celi
        
        geo = []
        while celda!=celf:
            xstr = round(flowvec.loc[celda,'X'],4)
            ystr = round(flowvec.loc[celda,'Y'],4)
            
            geo.append(Point(xstr,ystr))
            
            celda = int(flowvec.loc[celda,'IDD'])
        
        xstr = round(flowvec.loc[celda,'X'],4)
        ystr = round(flowvec.loc[celda,'Y'],4)
        geo.append(Point(xstr,ystr))
        
        # if len(geo) < 2:
        #     geo.append(Point(geo[-1].x-1,geo[-1].y-1))
            
        line = LineString(geo).wkt
        geom.append(line)

    fshape['geometry'] = geom
    fshape['geometry'] = fshape['geometry'].apply(wkt.loads)
    fshape.set_geometry(col='geometry', inplace=True)
        
    fshape.to_file(os.path.join(output_pth,'reaches.shp'))
#==============================================================================
#-- Parámetros de confinamiento
#==============================================================================
def DicCellTramo(FlowVectorA,TramosA):
    u"""
    Crea un diccionario con los diferentes tramos y sus respectivas celdas

    Entradas
    ----------
    FlowVector : DataFrame
        flow vector with the next information on each column:
                  [0] Id cell
                  [1] Id destination cell
                  [2] Cumulate Area
                  [3] Non cumulate Length
                  [4] Elevation (DEM)
                  [5] X: coordinate
                  [6] Y: coordinate
                  [7] i
                  [8] j
                  [9] Slope
                  [10] Cell type [0: Hillslope; 1:Drainage; 2:Ephemeral drainage;
                       3:reservoir]
                  [11] Nodes
                  [12] New corrected elevation in the drainage network
                  [13] corresponding Link number of the drainge-type cell


    Tramos:     DataFrame
                [0] Link ID
                [1] Link length
                [2] Link Slope
                [3] Link mean Basin Area
                [4] Id of destination link
                [5] initial cell Id
                [6] end cell Id
                [7] destination link initial cell id
                [8] corresponding Hillslope Area

    Salidas
    -------
    diccelltramo: diccionario
                    Key: Número de tramo
                    Value: Array con las celdas del tramo

    """
    diccelltramo={}
    for t in range(TramosA.shape[0]):
        diccelltramo[t] = FlowVectorA.loc[np.bitwise_and(FlowVectorA['Tipo']>0.0,FlowVectorA['ID_tramo']==t),'ID'].astype(int)
    return diccelltramo
#==============================================================================
def EjeSec(XYl,Rtrenzado,dx):
    u"""
    Genera eje secundario

    Entradas
    ----------
    XYl: Array
        Coordinates of each point drainage network


    Rtrenzado: Float
               Width of braided river with same area

    dx: float
        DEM resolution

    Salidas
    -------
    XY2: array
         Coordinates of new axis

    """
    XYl=RePOLY(XYl,min(0.2*Rtrenzado,dx))
    n=len(XYl)
    XY2=np.empty([n,2])-9999
    for i in range (n):
        anillo=((XYl[:,0]-XYl[i,0])**2+(XYl[:,1]-XYl[i,1])**2)**0.5
        k=np.where(anillo<=Rtrenzado)[0]
        XY2[i,0]=np.mean(XYl[k,0])
        XY2[i,1]=np.mean(XYl[k,1])
    return XY2
#==============================================================================
def RePOLY(XYl,malla):
    u"""
    Axis resample

    Entradas
    ----------
    XY1: Array
        Coordinates of each point in the axis

    malla: Float
           Resample size

    Salidas
    -------
    XY2: array
         Resample coordinates

    """
    S,L=sinuosidad(XYl)
    XYl=np.insert(XYl,2,0,axis=1)  #Agrega columna de long acumulada
    LongNoAcum=((XYl[:-1,0]-XYl[1:,0])**2+(XYl[:-1,1]-XYl[1:,1])**2)**0.5
    XYl[1:,2]=LongNoAcum.cumsum()
#    for i in range(XYl.shape[0]-1):
#        delta=((XYl[i,0]-XYl[i+1,0])**2+(XYl[i,1]-XYl[i+1,1])**2)**0.5
#        XYl[i+1,2]=XYl[i,2]+delta
    dist=malla
    XY2=[]
    XY2.append([XYl[0,0], XYl[0,1], XYl[0,2]])
    fila=1
    while dist<L*S:
        k=np.where(XYl[:,2]<=dist)[0]
        i=len(k)-1
        j=i + 1
        Xinicio=XYl[i,0];Yinicio=XYl[i,1];
        Xfinal=XYl[j,0];Yfinal=XYl[j,1];

        d = dist-XYl[i,2]
        teta = np.arctan(abs(Yfinal - (Yinicio + 0.01)) / abs(Xfinal - (Xinicio + 0.01)))

        if (Xfinal >= Xinicio) & (Yfinal >= Yinicio):
            x = Xfinal - d * np.cos(teta)
            y = Yfinal - d * np.sin(teta)

        if (Xfinal <= Xinicio) & (Yfinal >= Yinicio):
            x = Xfinal + d * np.cos(teta)
            y = Yfinal - d * np.sin(teta)

        if (Xfinal <= Xinicio) & (Yfinal <= Yinicio):
            x = Xfinal + d * np.cos(teta)
            y = Yfinal + d * np.sin(teta)

        if (Xfinal >= Xinicio) & (Yfinal <= Yinicio):
            x = Xfinal - d * np.cos(teta);
            y = Yfinal + d * np.sin(teta);

        XY2.append([x,y,dist])
        dist=dist+malla
        fila+=1

    XY2.append([XYl[-1,0], XYl[-1,1], XYl[-1,2]])
    XY2=np.array(XY2)

    return XY2
#==============================================================================
def sinuosidad(XYl):
    """
    Calcula la sinuosidad

    Entradas
    ----------
    XYl : array_like
        Coordenadas de cada punto de la red de drenaje.

    Salidas
    -------
    S : float
        Sinuosidad calculada como la longitud del tramo dividida entre la longitud recta
    L : float
        Distancia plana entre el punto inicial y final del tramo

    """
    if XYl.shape[0]<=1:
        S=0.0
        L=0.0
    else:
        LongNoAcum=((XYl[:-1,0]-XYl[1:,0])**2+(XYl[:-1,1]-XYl[1:,1])**2.0)**0.5
        dist=LongNoAcum.sum()
        L=(((XYl[0,0]-XYl[-1,0])**2+(XYl[0,1]-XYl[-1,1])**2.0)**0.5)
        S=dist/L
    return S,L

#==============================================================================
#-- Definir Laderas aferentes de cada tramo
#==============================================================================
def DefinirLaderas(ArregloFlujo):
    u"""
    Enumera las laderas con el mismo numero del tramo al que tributan sus aguas

    Entradas
    --------
    ArregloFlujo: DataFrame
       con los atributos de la cuenca
       [0] Identificador de la celda de partida
       [1] Identificador de la celda de destino
       [2] Área de drenaje acumulada
       [3] Longitud no acumulada en la dirección de drenaje
       [4] Elevación sobre el nivel del mar
       [5] Coordenada Este del centro de pixel
       [6] Coordenada Norte del centro de pixel
       [7] Índice i (fila) de la celda en el esquema ráster
       [8] Índice j (columna) de la celda en el esquema ráster
       [9] Pendiente en la dirección de drenaje
       [10] Tipo de celda [0: Ladera; 1: Cauce]
       [11] Nodos hidrológicos y topográficos
       [12] Elevación corregida sobre la red de drenajes
       [13] Identificador de tramo al que pertenece cada celda de cauce

    Salidas
    -------
    ArregloFlujo:  DataFrame
       con los atributos de la cuenca. Modifica las siguientes columnas:
       ['ID_tramo'] Asigna a cada celda (cauce y ladera) el identificador de tramo al que drenan
    """
    divisorias = np.setdiff1d(ArregloFlujo['ID'],ArregloFlujo['IDD']).astype(int)
    for cabecera in divisorias:
        celda    = cabecera
        trayecto = []
        while ArregloFlujo.loc[celda,'ID_tramo'] == -9999.0 and ArregloFlujo.loc[celda,'ID'] != ArregloFlujo.loc[celda,'IDD']:
            trayecto.append(celda)
            celda = int(ArregloFlujo.loc[celda,'IDD'])
        ArregloFlujo.loc[trayecto,'ID_tramo'] = ArregloFlujo.loc[celda,'ID_tramo']
    return ArregloFlujo


#==============================================================================
def FloodPlane(FlowVectorA,TramosA, dem_rst,red_rst):
    u"""
    Define Llanura de inundación

    Entradas
    ----------
    FlowVector : DataFrame
                con los atributos de la cuenca
                  [0] Id cell
                  [1] Id destination cell
                  [2] Cumulate Area
                  [3] Non cumulate Length
                  [4] Elevation (DEM)
                  [5] X: coordinate
                  [6] Y: coordinate
                  [7] i
                  [8] j
                  [9] Slope
                  [10] Cell type [0: Hillslope; 1:Drainage; 2:Ephemeral drainage;
                       3:reservoir]
                  [11] Nodes
                  [12] New corrected elevation in the drainage network
                  [13] corresponding Link number of the cell


    Tramos:     DataFrame
                [0] Link ID
                [1] Link length
                [2] Link Slope
                [3] Link mean Basin Area
                [4] Id of destination link
                [5] initial cell Id
                [6] end cell Id
                [7] destination link initial cell id
                [8] corresponding Hillslope Area

    ddem_rst: Raster 
              con el modelo de elevación digital

    red_rst: Raster 
             con la red de drenaje de la cuenca

    Salidas
    -------
    Tramos:     DataFrame
                Se crean dos nuevas columnas
                'Confinamiento': Grado de confinamiento
                'Qb'           : Caudal banca llena

    MapGeoLLANURA: llanura potenial

    MapGeoTRENZADO: Ancho cinturon de meandros o rio trenzado

    MapGeoVALLE: llanura topografica

    MapRedVirtual: representa red virtual

    """

    TramosA['Confinamiento'] = -9999  #Grado de confinamiento
    TramosA['Qb']            = -9999  #Caudal Banca llena

    Qb  = 1.4796 * np.power(TramosA['Area_cuenca_prom'],0.7484)         #caudal de banca llena para Tr=2.33 anos (UPME-UNALMED): Region Colombia
    WB1 = 17.748 * np.power(TramosA['Area_cuenca_prom'],0.3508)         #ancho medio de banca llena en rio Trenzado, equivalete al ancho de cinturon de meandros
    WB2 = 1.2399 * np.power(TramosA['Area_cuenca_prom'],0.5451)         #Ancho sistema migratorio


    radioTrenz = WB1/2                                 #Radio para rio trenzado
    NcTrenzado = np.ceil(radioTrenz/dem_rst['clsz'])
    factor1    = 6                                       #factor para radio sistema trenzado
    factor2    = 5.1132                                  #factor para radio sistema migratorio
    radio      = np.maximum(np.ceil(factor1*radioTrenz/dem_rst['clsz']),np.ceil(factor2*WB2/dem_rst['clsz']))
    #Ncllanura=np.ceil((radio)/dem_rst['clsz'])+1 #Llanura potencial (Dodov y Foufoula)
    Ncllanura  = np.ceil(radio)
    dic        = DicCellTramo(FlowVectorA,TramosA)

    # Crear Variables para almacenar mapas
    GeoLLANURA    = np.ones([dem_rst['nrows'],dem_rst['ncols']])*dem_rst['nodt']
    GeoTRENZADO   = np.ones([dem_rst['nrows'],dem_rst['ncols']])*dem_rst['nodt']
    GeoVALLE      = np.ones([dem_rst['nrows'],dem_rst['ncols']])*dem_rst['nodt']
    redSECUNDARIA = np.ones([dem_rst['nrows'],dem_rst['ncols']])*dem_rst['nodt']

    # Para cada Tramo calcular ancho Trenzado, ancho valle y ancho llanura
    for t in range(TramosA.shape[0]):

        CeldasValle    = np.ones([dem_rst['nrows'],dem_rst['ncols']])*-9999.   #Matriz temporal para almacenar celdas de la llanura topografica
        CeldasTrenzado = np.ones([dem_rst['nrows'],dem_rst['ncols']])*-9999.  #Matriz temporal para almacenar celdas del cinturon de meandros (rio Trenzado)
        CeldasLLanura  = np.ones([dem_rst['nrows'],dem_rst['ncols']])*-9999.   #Matriz temporal para almacenar celdas de la llanura potencial
        CeldasLadera   = np.ones([dem_rst['nrows'],dem_rst['ncols']])*-9999.
        iLadera        = FlowVectorA.loc[FlowVectorA['ID_tramo']==t,'i'].astype(int)
        jLadera        = FlowVectorA.loc[FlowVectorA['ID_tramo']==t,'j'].astype(int)
        CeldasLadera[iLadera,jLadera] = 1.0


        XY = np.array(FlowVectorA.loc[dic[t],['X','Y']])

        if radioTrenz[t]<2*dem_rst['clsz']:
            newXY=XY
        elif XY.shape[0]<=1:
            newXY=XY
        else:
            newXY=EjeSec(XY,radioTrenz[t],dem_rst['clsz']) # se debe formular eje auxiliar esta pendiente !!!!!

        ii,jj=XYtoij(newXY[:,0],newXY[:,1],dem_rst)
        redSECUNDARIA[ii,jj]=1.0

        # calcular limites de busqueda ancho trenzado para cada celda de cauce en el tramo
        a = np.array(ii-NcTrenzado[t]).astype(int)
        b = np.array(ii+NcTrenzado[t]).astype(int)
        c = (jj-NcTrenzado[t]).astype(int)
        d = (jj+NcTrenzado[t]).astype(int)
#        a[np.where(a<0)]=0
#        b[np.where(b > dem_rst['nrows']-1)]=dem_rst['nrows']-1
#        c[np.where(c<0)]=0
#        d[np.where(d > dem_rst['ncols']-1)]=dem_rst['ncols']-1

        # calcular indices dentro limites ancho trenzado
        aux_itren = np.array([np.arange(a[m],b[m]+1) for m in range(len(a))])
        aux_itren = np.repeat(aux_itren,int(2*NcTrenzado[t]+1),axis=0)
        aux_itren = aux_itren.reshape(aux_itren.size)
        aux_itren[np.where(aux_itren < 0)[0]] = 0
        aux_itren[np.where(aux_itren > dem_rst['nrows']-1)[0]] = dem_rst['nrows']-1

        aux_jtren = np.array([np.arange(c[m],d[m]+1) for m in range(len(c))])
        aux_jtren = np.repeat(aux_jtren,int(2*NcTrenzado[t]+1),axis=1)
        aux_jtren = aux_jtren.reshape(aux_jtren.size)
        aux_jtren[np.where(aux_jtren < 0)[0]] = 0
        aux_jtren[np.where(aux_jtren > dem_rst['ncols']-1)[0]] = dem_rst['ncols']-1

        # Marcar celdas de ancho trenzado en tramo
        CeldasTrenzado[aux_itren,aux_jtren]=1
        CeldasTrenzado[CeldasLadera != 1.0]=dem_rst['nodt']
        GeoTRENZADO[np.where(CeldasTrenzado==1.0)]=t

        # calcular limites de busqueda de ancho llanura
        a1=ii-Ncllanura[t].astype(int)
        b1=ii+Ncllanura[t].astype(int)
        c1=jj-Ncllanura[t].astype(int)
        d1=jj+Ncllanura[t].astype(int)
#        a1[np.where(a1<0)]=0
#        b1[np.where(b1 > dem_rst['nrows']-1)]=dem_rst['nrows']-1
#        c1[np.where(c1<0)]=0
#        d1[np.where(d1 > dem_rst['ncols']-1)]=dem_rst['ncols']-1

        # calcular indices dentro limites de ancho llanura
        aux_illan=np.array([np.arange(a1[m],b1[m]+1) for m in range(len(a1))])
        aux_illan=np.repeat(aux_illan,int(2*Ncllanura[t]+1),axis=0)
        aux_illan=aux_illan.reshape(aux_illan.size)
        aux_illan[np.where(aux_illan < 0)[0]]=0
        aux_illan[np.where(aux_illan > dem_rst['nrows']-1)[0]]=dem_rst['nrows']-1

        aux_jllan=np.array([np.arange(c1[m],d1[m]+1) for m in range(len(c1))])
        aux_jllan=np.repeat(aux_jllan,int(2*Ncllanura[t]+1),axis=1)
        aux_jllan=aux_jllan.reshape(aux_jllan.size)
        aux_jllan[np.where(aux_jllan < 0)[0]]
        aux_jllan[np.where(aux_jllan > dem_rst['ncols']-1)[0]]=dem_rst['ncols']-1

        # Identificar celdas de la potencial llanura
        CeldasLLanura[aux_illan,aux_jllan]=1
        CeldasLLanura[CeldasLadera != 1.0]=dem_rst['nodt']
        GeoLLANURA[np.where(CeldasLLanura==1.0)]=t

        # muestrear elevaciones en los limites de ancho de llanura
        #discos=np.reshape(dem_rst['mtrx'][aux_illan,aux_jllan],(int(2*Ncllanura[t]+1)**2,len(c1)),'F')

        # Muestrear elevaciones en los limites del ancho del rio trenzado
        discos=np.reshape(dem_rst['mtrx'][aux_itren,aux_jtren],(int(2*NcTrenzado[t]+1)**2,len(c)),'F')

        # Crear arreglo para comparar altura de referencia dentro del limite en valle
        Zperc=np.percentile(discos,10,axis=0)
        frac=0.0
        Zref=Zperc + frac
        Zref=np.repeat(Zref,int(2*Ncllanura[t]+1)**2,axis=0)
        Zref=np.reshape(Zref,Zref.size)

        # Identificar celdas del cinturon geomorfologico
        is_valle=np.bitwise_or((dem_rst['mtrx'][aux_illan,aux_jllan]<=Zref),(red_rst['mtrx'][aux_illan,aux_jllan]>=0),(redSECUNDARIA[aux_illan,aux_jllan]==1))
        CeldasValle[aux_illan[is_valle],aux_jllan[is_valle]]=1
        CeldasValle[CeldasLadera != 1.0]=dem_rst['nodt']
        GeoVALLE[np.where(CeldasValle==1.0)]=t


        areaTrenzado=np.sum(CeldasTrenzado==1)*(dem_rst['clsz']**2)      #area de la zona potencialmente trenzada
        areaValle=np.sum(CeldasValle==1)*(dem_rst['clsz']**2)             #area del las celdas del corredor fluvial
        #areaLLanura=len(np.where(CeldasLLanura==1)[0])*(dem_rst['clsz']**2)         #area del las celdas del corredor fluvial

        TramosA.loc[t,'Confinamiento'] =areaValle/areaTrenzado;     #Grado de confinamiento [9]
        TramosA.loc[t,'Qb']=Qb[t];                                  #Caudal de banca llena  [10]

    MapGeoLLANURA={'ncols':dem_rst['ncols'],'nrows':dem_rst['nrows'],'nclls':dem_rst['nclls'],'xll':dem_rst['xll'],'yll':dem_rst['yll'],'xur':dem_rst['xur'],'yur':dem_rst['yur'],'clsz':dem_rst['clsz'],'nodt':-9999,'mtrx':GeoLLANURA,'prj':''}
    MapGeoTRENZADO={'ncols':dem_rst['ncols'],'nrows':dem_rst['nrows'],'nclls':dem_rst['nclls'],'xll':dem_rst['xll'],'yll':dem_rst['yll'],'xur':dem_rst['xur'],'yur':dem_rst['yur'],'clsz':dem_rst['clsz'],'nodt':-9999,'mtrx':GeoTRENZADO,'prj':''}
    MapGeoVALLE={'ncols':dem_rst['ncols'],'nrows':dem_rst['nrows'],'nclls':dem_rst['nclls'],'xll':dem_rst['xll'],'yll':dem_rst['yll'],'xur':dem_rst['xur'],'yur':dem_rst['yur'],'clsz':dem_rst['clsz'],'nodt':-9999,'mtrx':GeoVALLE,'prj':''}
    MapRedVirtual={'ncols':dem_rst['ncols'],'nrows':dem_rst['nrows'],'nclls':dem_rst['nclls'],'xll':dem_rst['xll'],'yll':dem_rst['yll'],'xur':dem_rst['xur'],'yur':dem_rst['yur'],'clsz':dem_rst['clsz'],'nodt':-9999,'mtrx':redSECUNDARIA,'prj':''}

    return TramosA, MapGeoLLANURA, MapGeoTRENZADO, MapGeoVALLE, MapRedVirtual


#==============================================================================
#--Clasificación de tramos de acuerdo a Flores y Beechie
#==============================================================================
def GetMorpho(FlowVectorA,TramosA):

    u"""
    Clasifica los tramos de acuerdo a su morfología (Flores, 2006) (Beachie, 2006)

    Entradas
    ----------
    FlowVector :  DataFrame
                DataFrame con los atributos de la cuenca
    Tramos:     DataFrame
                [0] Link ID
                [1] Link length
                [2] Link Slope
                [3] Link mean Basin Area
                [4] Id of destination link
                [5] initial cell Id
                [6] end cell Id
                [7] destination link initial cell id
                [8] corresponding Hillslope Area
                [9] Grado de confinamiento
                [10] Bank full streamflow
                [11] Nodos
                [12] Elevación corregida
                [13] Número de tramo

    Salidas
    -------
    FlowVector :  DataFrame
                 Se crean columnas:
                 ['Beechie']: Patrón de alineamiento (Beechie, 2006)
                 ['Flores']:  Forma del lecho (Flores et al., 2006)
    Tramos:     DataFrame
                 Se crean columnas:
                 ['Potencia']: Potencia del tramo
                 ['Beechie']: Patrón de alineamiento (Beechie, 2006)
                 ['Flores']:  Forma del lecho (Flores et al., 2006)

    """
    TramosA['Potencia'] = -9999  #Potencia [11]
    TramosA['Beechie']  = -9999  #Beachie [12]
    TramosA['Flores']   = -9999  #Flores [13]
    
    FlowVectorA['Beechie']  = -9999 
    FlowVectorA['Flores']   = -9999 
    
    #W=coA^b; Q=c1A^d
    d=0.8
    b=0.4
    coef=d-b
    TramosA['Potencia'] = np.abs(TramosA['Pendiente']*np.power(TramosA['Area_cuenca_prom'],coef))

    #Patrones de alineamiento
    # Trenzado          : Patron =5 ;Braided S>0.1*Q**-0.42
    TramosA.loc[TramosA['Pendiente']>(0.1*(TramosA['Qb']**(-0.42))),'Beechie'] = 5.0
    # Island-Braided    : Patron =4 ;Island-Braided Q>15 , 0.05*Q**-0.61<S<0.1*0.1*Q**-0.42
    TramosA.loc[(TramosA['Qb']>15) * (TramosA['Pendiente']>0.05*TramosA['Qb']**-0.61)*(TramosA['Pendiente']<0.1*TramosA['Qb']**-.42),'Beechie'] = 4.0
    # Meandrico         : Patron =3 ;Meandering Q>15 , S<0.05*Q**-0.61
    TramosA.loc[(TramosA['Qb']>15) * (TramosA['Pendiente']<0.05*TramosA['Qb']**-0.61),'Beechie'] = 3.0
    # Recto             : Patron =2 ;Straight Q<15 S<0.1*Q**-0.42
    TramosA.loc[(TramosA['Qb']<15) * (TramosA['Pendiente']<0.1*TramosA['Qb']**-0.42),'Beechie']  = 2.0
    # Confinado         : Patron =1 ;Confined GC<=1
    #UmbralConf=(4*1.2399*TramosA[:,3]**(0.5451-0.3823))/(15.336)
    TramosA.loc[TramosA['Confinamiento']<=1.0,'Beechie'] = 1.0
    #TramosA[(TramosA[:,9]<=UmbralConf),len(TramosA[0])-2]=1


    #Formas del del lecho
    #Limitados por capacidad
    pool_riffle = np.logical_and(np.abs(TramosA['Pendiente'])<=0.025,TramosA['Potencia']<=0.055)
    plane_bed   = np.logical_and(np.abs(TramosA['Pendiente'])<=0.025,TramosA['Potencia']>0.055)
    step_pool   = np.logical_and(np.abs(TramosA['Pendiente'])>0.025, TramosA['Potencia']<=0.055)
    cascade     = np.logical_and(np.abs(TramosA['Pendiente'])>0.025, TramosA['Potencia']>0.055)
    TramosA.loc[pool_riffle,'Flores'] = 4
    TramosA.loc[plane_bed,'Flores']   = 3
    TramosA.loc[step_pool,'Flores']   = 2
    TramosA.loc[cascade,'Flores']     = 1

    #Para colocar en flow Vector patron de alineamiento
    ReachtoFvec(TramosA,FlowVectorA,'Beechie')
#
#    aux=zip(TramosA[:,len(TramosA[0])-2],TramosA[:,0])
#    dic={}
#    dic[-9999]=-9999
#    map(lambda x:fillDic(x,dic),aux)
#    FlowVectorA[:,14]=np.array(map(lambda x: dic[x],FlowVectorA[:,13]))

    #Para colocar en flow Vector Patron formas del lecho
    ReachtoFvec(TramosA,FlowVectorA,'Flores')
#
#    aux=zip(TramosA[:,len(TramosA[0])-1],TramosA[:,0])
#    dic={}
#    dic[-9999]=-9999
#    map(lambda x:fillDic(x,dic),aux)
#    FlowVectorA[:,15]=np.array(map(lambda x: dic[x],FlowVectorA[:,13]))

    return FlowVectorA,TramosA

#==============================================================================
# Geometría hidráulica
#==============================================================================
def GeoHidraulica (FlowVectorA,TramosA):
    u"""
    Geometría hidráulica

    Entradas
    ----------
    FlowVector :  DataFrame
                Atributos de la cuenca
    Tramos:     DataFrame
                [0] Link ID
                [1] Link length
                [2] Link Slope
                [3] Link mean Basin Area
                [4] Id of destination link
                [5] initial cell Id
                [6] end cell Id
                [7] destination link initial cell id
                [8] corresponding Hillslope Area
                [9] Grado de confinamiento
                [10] Bank full streamflow
                [11] link power
                [12] Alignment Patterns Beachie
                [13] Bed forms Flores

    Salidas
    -------
    Tramos:     DataFrame
                ['Wb'] Ancho banca llena

    FlowVector :  DataFrame
                  ['Qb'] Caudal banca llena
                  ['Wb'] Ancho banca llena

   """

    TramosA['Wb']     = -9999 #Ancho Banca Llena [14]
    FlowVectorA['Qb'] = -9999
    FlowVectorA['Wb'] = -9999
    
    #Limitados por suministro
    TramosA.loc[np.logical_and(TramosA['Beechie']==1,np.logical_or(TramosA['Flores']==1,TramosA['Flores']==2)),'Wb'] = \
            7.5789*TramosA.loc[np.logical_and(TramosA['Beechie'] == 1,np.logical_or(TramosA['Flores']==1,TramosA['Flores']==2)),'Area_cuenca_prom']**0.2025
    #Limitados por capacidad
    TramosA.loc[np.logical_and(TramosA['Beechie']==1,np.logical_or(TramosA['Flores']==3,TramosA['Flores']==4)),'Wb'] = \
            1.7712*TramosA.loc[np.logical_and(TramosA['Beechie']==1,np.logical_or(TramosA['Flores']==3,TramosA['Flores']==4)),'Area_cuenca_prom']**0.4217
    #Limitados por capacidad Libres
    TramosA.loc[TramosA['Beechie']==2,'Wb'] = 1.2399*TramosA.loc[TramosA['Beechie']==2,'Area_cuenca_prom']**0.5451
    #Trenzados y cinturon de meandros
    TramosA.loc[np.logical_or(TramosA['Beechie']==3,np.logical_or(TramosA['Beechie']==4,TramosA['Beechie']==5)),'Wb'] = \
            15.336*TramosA.loc[np.logical_or(TramosA['Beechie']==3,np.logical_or(TramosA['Beechie']==4,TramosA['Beechie']==5)),'Area_cuenca_prom']**0.3823

    #Para colocar en flow Vector caudales banca llena
#    aux=zip(TramosA[:,10],TramosA[:,0])
#    dic={}
#    dic[-9999]=-9999
#    map(lambda x:fillDic(x,dic),aux)
#    FlowVectorA[:,16]=np.array(map(lambda x: dic[x],FlowVectorA[:,13]))
    mascara = FlowVectorA['ID_tramo'] != -9999.0
    filas   = FlowVectorA.loc[mascara,'ID_tramo'].astype(int)
    FlowVectorA.loc[mascara,'Qb'] = TramosA.loc[filas,'Qb'].values

    #Para colocar en flow Vector ancho banca llena
#    aux=zip(TramosA[:,14],TramosA[:,0])
#    dic={}
#    dic[-9999]=-9999
#    map(lambda x:fillDic(x,dic),aux)
#    FlowVectorA[:,17]=np.array(map(lambda x: dic[x],FlowVectorA[:,13]))
    FlowVectorA.loc[mascara,'Wb'] = TramosA.loc[filas,'Wb'].values

    return TramosA, FlowVectorA

#==============================================================================
# ETP cenicafe
#==============================================================================
def ETPcenicafe(elev):
    """
    Evapotranspiración potencial de acuerdo a Jaramillo (2006)

    Entradas
    ----------
    elev : array
        Elevación (metros sobre nivel del mar)

    Salidas
    -------
    ETP: array
        Evapotranspiración potencial en m/dia
    """
    ETP=4.37*np.exp((-0.0002*elev))/1000.
    return ETP

#==============================================================================
# Manning on the drainage network as a function of the slope
#==============================================================================
def ManningCauce(So):
    """
    Se calcula coeficiente de rugosidad de Manning dada la pendiente del tramo, de acuerdo con Yochum, Comiti, Whol, David & Mao (2014)

    Entradas
    ----------
    So : array
        Pendiente

    Salidas
    -------
    n : array
        Manning's roughness coefficient or array of coefficients
    """
    # n de Manning para la condición de aguas bajas
    nl = 1.85*So**0.79
    
    # n de Manning para la condición de banca llena 
    nb = 0.5005052393*nl**0.8793182984
    
    return nb

#==============================================================================
# Propiedades de tramo a flowvector
#==============================================================================

def ReachtoFvec(Tramos,FlowVectorA,key):
    draincells     = FlowVectorA.index[FlowVectorA['Tipo'] > 0.0]
    draincells_IDT = list(FlowVectorA.loc[draincells,'ID_tramo'].astype(int))
    FlowVectorA.loc[draincells,key] = Tramos.loc[draincells_IDT,key].values
    return FlowVectorA

#==============================================================================
# Calcula profundidad de banca llena
#==============================================================================
def BankfullH(Qb, So, Wb, n):
    """
    Radio hidráulico del canal para una sección transversal parabólica, dada el área de flujo, el y la geometría en condiciones de 
    banca llena (ancho y profundidad).

    Entradas
    ----------
    Qb : array
        Caudal de banca llena

    So :
        Pendiente del tramo

    Wb: array
        Ancho en condiciones de banca llena

    n: array
        Coeficiente de rugosidad de Manning en condiciones de banca llena

    Salidas
    -------
    Hb: array
        Profundidad en condiciones de banca llena
    """
    
    num_tra=np.size(Wb)
    a4=(-1.)*((3.**3.)*2.*(Qb**3.)*(n**3.)*(So**(-3./2.))*(Wb**(-7.)))
    a2=(-1.)*((3.**4.)*(2.**(-1.))*(Qb**3.)*(n**3.)*(So**(-3./2.))*(Wb**(-5.)))
    a0=(-1.)*((3.**5.)*(2.**(-5.))*(Qb**3.)*(n**3.)*(So**(-3./2.))*(Wb**(-3.)))
    p=[[1, a4[i], 0, a2[i], 0, a0[i]] for i in range(np.size(Qb))]
    Hb=np.zeros(num_tra)
    for i in range(0,num_tra):
        raices=np.roots(p[i][:])
        mascara=np.isreal(raices)
        Hb[i]=np.real(raices[mascara])
    
#    Hb=np.copy(Wb)*0.0
#    for i in range(Qb.size):
#        numiter=1
#        Hbsem=((Qb[i]*n[i])/(Wb[i]*(So[i]**0.5)))**(3.0/5.0)
#        err=1000
#        while err > 0.001 and numiter < 20:
#            LadoIzq=((Qb[i]**3.0)*(n[i]**3.0)/(So[i]**1.5))*((Wb[i]**2.0)+(16.0*(Hbsem**2.0)/3.0)+((8.0*(Hbsem**2)/(3.0*Wb[i]))**2.0))
#            HbDer=(LadoIzq**0.2)/((2.0/3.0)*Wb[i])
#            err=np.abs(HbDer-Hbsem)
#            Hbsem=HbDer
#            numiter=+1
#        Hb[i]=Hbsem
    return Hb
#==============================================================================
# Calculo de alfa y beta
#==============================================================================
def RadioHidraulicoArea(Hb,Wb):

    u"""
    Encuentra los coeficientes y exponentes de una relación potencial entre
    el radio hidráulico y el área mojada de una sección parabólica

    Entradas
    ----------
    Hb : Arreglo de reales,
        Profundidades de banca llena
    Wb : Arreglo de reales,
        Anchos de banca llena

    Salidas
    -------
    alfa : Arreglo de reales,
        Coeficientes de la relación potencial
    beta : Arreglo de reales,
        Exponentes de la relación potencial
    """
    num_pun=100
    num_sec=np.size(Hb)
    lin_int=np.zeros(num_sec)
    beta=np.zeros(num_sec)
    for i in range(0,num_sec):
        k=4.0*Hb[i]/(Wb[i]**2.0)
        y=np.linspace(Hb[i]/num_pun,Hb[i],num_pun)
        W=2.0*np.sqrt(y/k)
        x=4.0*y/W
        Pm=(W/2.0)*(np.sqrt(1.0+np.power(x,2.0))+np.log(x+np.sqrt(1.0+np.power(x,2.0)))/x)
        Am=2.0*W*y/3.0
        Rh=Am/Pm
        Am=np.log10(Am)
        Rh=np.log10(Rh)
        mtrx=np.vstack([Am,np.ones(len(Am))]).T
        beta[i],lin_int[i]=np.linalg.lstsq(mtrx,Rh)[0]
    alfa=np.power(10.0,lin_int)
    return alfa,beta

#==============================================================================
# ESTIMACIÓN DE PARÁMETROS DEL MODELO MDLC-ADZ
#==============================================================================

def ParametrosMDLCADZ(Hb,Wb,S0_in,nb,L_in,Qb,beta=0.0,num_pun=200):
    
    u"""
    Calcula los coeficientes y exponentes de relaciones Us=alfa1*Q^beta1 y
    maxUs=alfa2*Q^beta2, donde Us=L/tm y maxUs=L/tau utilizando el modelo
    MDLC-ADZ (Camacho y Lees, 2000)
    
    Entradas
    ----------
    Hb : Arreglo de reales,
        Profundidades de banca llena
    Wb : Arreglo de reales,
        Anchos de banca llena
    S0_in : Arreglo de reales
        Pendiente de tramos
    nb : Arreglo de reales
        Coeficiente de rugosidad de Manning para la condición de banca llena
    L_in : Arreglo de reales
        Longitud de los tramos
    beta : Real,
        Coeficiente de retraso del soluto
    [num_pun] : Entero,
        Número de puntos a considerar en los diagramas de dispersión    
    
    Salidas
    -------
    alfa1: Arreglo de reales,
        Coeficientes de la relación potencial U=alfa1*Am**beta1
    beta1: Arreglo de reales,
        Exponentes de la relación potencial U=alfa1*Q**beta1
    alfa2: Arreglo de reales,
        Coeficientes de la relación potencial Us=alfa2*Q**beta2
    beta2: Arreglo de reales,
        Exponentes de la relación potencial Us=alfa2*Q**beta2
    alfa3: Arreglo de reales,
        Coeficientes de la relación potencial maxUs=alfa3*Q**beta3
    beta3: Arreglo de reales,
        Exponentes de la relación potencial maxUs=alfa3*Q**beta3
    """
    
    # Calados para los que se construirá el diagrama de dispersión
    H=(np.linspace(1.0,0.0,num_pun,False)[:,None]*2.0*Hb).T
    
    # Cálculo de la geometría hidráulica transversal para una sección parabólica
    k=4.0*Hb/np.square(Wb)              # Coeficiente de geometría hidráulica de la parábola
    k=np.array([k,]*num_pun).T
    W=2.0*np.sqrt(H/k)                  # Arreglo de anchos correspondientes a los calados guardados en H
    Am=2.0/3.0*W*H                      # Arreglo de areas mojadas correspondientes a los calados guardados en H
    x=4.0*H/W
    Pm=(W/2.0)*(np.sqrt(1.0+np.power(x,2.0))+np.log(x+np.sqrt(1.0+np.power(x,2.0)))/x)
    Rh=Am/Pm
    
    # Cálculo de los n_Manning correspondientes a los calados guardados en H
    # de acuerdo con lo encontrado por Yochum. El cálculo se acota por la
    # condición de aguas bajas
    nl=2.197080245069931*nb**1.1372446153111921
    Sz=Hb*(nb/0.41)**(1.0/0.69)
    Sz=Sz[:,None]
    n=0.41*(H/Sz)**-0.69
    for i in range(0,np.size(n,0)):
        n[i,:][n[i,:]>nl[i]]=nl[i]
        n[i,:][n[i,:]<nb[i]]=nb[i]
    
    S0=np.array([S0_in,]*num_pun).T
    Q=1.0/n*Am*np.power(Rh,2.0/3.0)*np.sqrt(S0) # m3/s
    U=Q/Am # m/s
    Dh=Am/W
    
    # Cálculo de la celeridad de la onda cinemática (ver Lança, 2000)
    a=np.power(n*np.power(Pm,2.0/3.0)/np.sqrt(S0),3.0/5.0)
    b=3.0/5.0
    c=1.0/(a*b*np.power(Q,b-1.0))
    
    # Estimación de parámetros del modelo MDLC-ADZ (ver Rogeliz, 2011)
    m=c/U
    theta=np.arctan(S0)
    Fr2=np.square(U/np.sqrt(9.81*Dh*np.cos(theta))) # Número de Froude al cuadrado
    L=np.array([L_in,]*num_pun).T
    t1=1.0+(m-1.0)*Fr2
    t2=H/(S0*L)
    t3=L/(m*U)
    t4=1.0-np.square(m-1.0)*Fr2
    K=3.0*t1*t2*t3/(2.0*m)
    n=4.0*m*t4/(9.0*np.square(t1)*t2)
    tau_fl=t3*(1.0-2.0*t4/(3.0*t1))
    tm=m*(n*K+tau_fl)*(1.0+beta)
    tau=m*tau_fl*(1.0+beta)
    
    # Conversión de unidades
    U=U*86400.0         # De m/s a m/día
    Us=L/tm*86400.0     # De m/s a m/día
    maxUs=L/tau*86400.0 # De m/s a m/día
    Q=Q*86400.0         # De m3/s a m3/día
    
    # Ajuste de las curvas potenciales U=alfa1*Am**beta1,
    # Us=alfa2*Q^beta2 y maxUs=alfa3*Q^beta3
    x1=np.log(Am)
    x2=np.log(Q)
    y1=np.log(U)
    y2=np.log(Us)
    y3=np.log(maxUs)
    num_sec=np.size(Hb)
    alfa1=np.zeros(num_sec)
    beta1=np.zeros(num_sec)
    alfa2=np.zeros(num_sec)
    beta2=np.zeros(num_sec)
    alfa3=np.zeros(num_sec)
    beta3=np.zeros(num_sec)
    for i in range(0,num_sec):
        matriz1=np.vstack([x1[i,:],np.ones(num_pun)]).T
        matriz2=np.vstack([x2[i,:],np.ones(num_pun)]).T
        beta1[i],alfa1[i]=np.linalg.lstsq(matriz1,y1[i,:])[0]
        beta2[i],alfa2[i]=np.linalg.lstsq(matriz2,y2[i,:])[0]
        beta3[i],alfa3[i]=np.linalg.lstsq(matriz2,y3[i,:])[0]
    alfa1=np.exp(alfa1)
    alfa2=np.exp(alfa2)
    alfa3=np.exp(alfa3)

    return alfa1,beta1,alfa2,beta2,alfa3,beta3
    

#==============================================================================
# Create Basin
#==============================================================================
def GetBasin(flowvecA,TramosA):
    u"""
    Builds first Topolgy (without any scenario charged) as a list of dictionaries

    Entradas
    ----------
    flowvec : array_like
        Each row corresponds to a cell of the basin. Each column has a property of each cell according with the following indices:
        [0] Id cell.
        [1] Id destination cell.
        [2] Cumulate Area.
        [3] Non cumulate Length.
        [4] Elevation (DEM).
        [5] X: coordinate.
        [6] Y: coordinate.
        [7] i.
        [8] j.
        [9] Slope.
        [10] Cell type [0: Hillslope; 1:Drainage; 2:Ephemeral drainage; 3:reservoir].
        [11] Nodes.
        [12] New corrected elevation in the drainage network.
        [13] corresponding Link number of the drainge-type cell.
        [14] Alignment Patterns Beachie.
        [15] Bed forms Flores.
        [16] Bank full streamflow.
        [17] Width full bank.

    Tramos : array_like
        Each row corresponds to a link of the basin. Each column has a property of each link according with the following indices:
        [0] Link ID
        [1] Link length
        [2] Link Slope
        [3] Link mean Basin Area
        [4] Id of destination link
        [5] initial cell Id
        [6] end cell Id
        [7] destination link initial cell id
        [8] corresponding Hillslope Area
        [9] Grado de confinamiento
        [10] Bank full streamflow
        [11] link power
        [12] Alignment Patterns Beachie
        [13] Bed forms Flores
        [14] Width full bank

    Salidas
    -------
    Basin :  list of dictionaries
        List of 4 dictionaries with each Basin considered in the model according the indices:
        [0] Hillslope cells.
        [1] Drainage network cells.
        [2] Ephemeral drainage network cells.
        [3] Reservoir cells.
    """
    flowvecLadera= flowvecA[np.where(flowvecA[:,10]==0)][:,0]
    aux1=zip(range(len(flowvecLadera)),(flowvecLadera))
    flowvecCauce= flowvecA[np.where(flowvecA[:,10]==1)][:,0]
    aux2=zip(range(len(flowvecCauce)),(flowvecCauce))
    flowvecCarcava=flowvecA[np.where(flowvecA[:,10]==2)][:,0]
    aux3=zip(range(len(flowvecCarcava)),(flowvecCarcava))
    flowvecEmbalse=flowvecA[np.where(flowvecA[:,10]>=3.0)][:,0]
    aux4=zip(range(len(flowvecEmbalse)),(flowvecEmbalse))
    dic1={}
    map(lambda x:fillDic(x,dic1),aux1)
    map(lambda x:fillDic(x,dic1),aux2)
    map(lambda x:fillDic(x,dic1),aux3)
    map(lambda x:fillDic(x,dic1),aux4)
    ind_dic=np.shape(flowvecA)[1]
    aux0=flowvecA[:,3] * 0. - 9999.
    aux1=np.array(aux0, ndmin=2).T
    flowvecA=np.append(flowvecA,aux1,axis=1) #Dic donde drena
    ind_idd=np.shape(flowvecA)[1]
    flowvecA=np.append(flowvecA,aux1,axis=1) #id donde drena
    flowvecA[:,ind_dic]=np.floor(flowvecA[flowvecA[:,1].astype(int),10])
    flowvecA[:,ind_idd]=np.array(map(lambda x:dic1[x],flowvecA[:,1]))

    Ladera={'Id_dic':flowvecA[np.where(flowvecA[:,10]==0)][:,ind_dic],\
            'Id_dren': flowvecA[np.where(flowvecA[:,10]==0)][:,ind_idd],\
            'X-Y':flowvecA[np.where(flowvecA[:,10]==0)][:,5:7],\
            'i,j':flowvecA[np.where(flowvecA[:,10]==0)][:,7:9],\
            'dx':flowvecA[np.where(flowvecA[:,10]==0)][:,3],\
            'A':flowvecA[np.where(flowvecA[:,10]==0)][:,2],\
            'So':flowvecA[np.where(flowvecA[:,10]==0)][:,9],\
            'ETP':flowvecA[np.where(flowvecA[:,10]==0)][:,18],\
            'n_Mann':np.nan,\
            'Hi':np.nan,\
            'Hu1':np.nan,\
            'Hu3':np.nan,\
            'Ks':np.nan,\
            'Kp':np.nan,\
            'P':np.nan,\
            'C':np.nan,\
            'K':np.nan,\
            'Arenas':np.nan,\
            'Limos':np.nan,\
            'Arcillas':np.nan}

    Cauce={'Id_dic':flowvecA[np.where(flowvecA[:,10]==1)][:,ind_dic],\
            'Id_dren': flowvecA[np.where(flowvecA[:,10]==1)][:,ind_idd],\
            'Id_Tramo':flowvecA[np.where(flowvecA[:,10]==1)][:,13],\
            'Id_Control':np.where(flowvecA[np.where(flowvecA[:,10]==1)][:,11]==4)[0],\
            'X-Y':flowvecA[np.where(flowvecA[:,10]==1)][:,5:7],\
            'i,j':flowvecA[np.where(flowvecA[:,10]==1)][:,7:9],\
            'dx':flowvecA[np.where(flowvecA[:,10]==1)][:,3],\
            'A':flowvecA[np.where(flowvecA[:,10]==1)][:,2],\
            'So':flowvecA[np.where(flowvecA[:,10]==1)][:,9],\
            'Qb':flowvecA[np.where(flowvecA[:,10]==1)][:,16],\
            'Wb':flowvecA[np.where(flowvecA[:,10]==1)][:,17],\
            'ETP':flowvecA[np.where(flowvecA[:,10]==1)][:,18],\
            'n_Mann':flowvecA[np.where(flowvecA[:,10]==1)][:,19],\
            'Hb':flowvecA[np.where(flowvecA[:,10]==1)][:,20],\
            'Hi':np.nan,\
            'Hu1':np.nan,\
            'Hu3':np.nan,\
            'Ks':np.nan,\
            'Kp':np.nan,\
            'P':np.nan,\
            'C':np.nan,\
            'K':np.nan,\
            'Arenas':np.nan,\
            'Limos':np.nan,\
            'Arcillas':np.nan,\
            'E_banca':np.nan}

    Carcava={'Id_dic':flowvecA[np.where(flowvecA[:,10]==2)][:,ind_dic],\
            'Id_dren': flowvecA[np.where(flowvecA[:,10]==2)][:,ind_idd],\
            'Id_Tramo':flowvecA[np.where(flowvecA[:,10]==2)][:,13],\
            'Id_Control':np.where(flowvecA[np.where(flowvecA[:,10]==2)][:,11]==4)[0],\
            'X-Y':flowvecA[np.where(flowvecA[:,10]==2)][:,5:7],\
            'i,j':flowvecA[np.where(flowvecA[:,10]==2)][:,7:9],\
            'dx':flowvecA[np.where(flowvecA[:,10]==2)][:,3],\
            'A':flowvecA[np.where(flowvecA[:,10]==2)][:,2],\
            'So':flowvecA[np.where(flowvecA[:,10]==2)][:,9],\
            'Qb':flowvecA[np.where(flowvecA[:,10]==2)][:,16],\
            'Wb':flowvecA[np.where(flowvecA[:,10]==2)][:,17],\
            'ETP':flowvecA[np.where(flowvecA[:,10]==2)][:,18],\
            'n_Mann':flowvecA[np.where(flowvecA[:,10]==2)][:,19],\
            'Hb':flowvecA[np.where(flowvecA[:,10]==2)][:,20],\
            'Hi':np.nan,\
            'Hu1':np.nan,\
            'Hu3':np.nan,\
            'Ks':np.nan,\
            'Kp':np.nan,\
            'P':np.nan,\
            'C':np.nan,\
            'K':np.nan,\
            'Arenas':np.nan,\
            'Limos':np.nan,\
            'Arcillas':np.nan}

    Embalse={'Id_dic':flowvecA[np.where(flowvecA[:,10]>=3)][:,ind_dic],\
            'Id_dren': flowvecA[np.where(flowvecA[:,10]>=3)][:,ind_idd],\
            'Id_Tramo':flowvecA[np.where(flowvecA[:,10]>=3)][:,13],\
            'Id_Control':np.where(flowvecA[np.where(flowvecA[:,10]>=3)][:,11]==4)[0],\
            'X-Y':flowvecA[np.where(flowvecA[:,10]>=3)][:,5:7],\
            'i,j':flowvecA[np.where(flowvecA[:,10]>=3)][:,7:9],\
            'dx':flowvecA[np.where(flowvecA[:,10]>=3)][:,3],\
            'A':flowvecA[np.where(flowvecA[:,10]>=3)][:,2],\
            'So':flowvecA[np.where(flowvecA[:,10]>=3)][:,9],\
            'Qb':flowvecA[np.where(flowvecA[:,10]>=3)][:,16],\
            'Wb':flowvecA[np.where(flowvecA[:,10]>=3)][:,17],\
            'type':(flowvecA[np.where(flowvecA[:,10]>=3)][:,10]-3)*10,\
            'ETP':flowvecA[np.where(flowvecA[:,10]>=3)][:,18],\
            'n_Mann':flowvecA[np.where(flowvecA[:,10]>=3)][:,19],\
            'Hb':flowvecA[np.where(flowvecA[:,10]==3)][:,20],\
            'Hi':np.nan,\
            'Hu1':np.nan,\
            'Hu3':np.nan,\
            'Ks':np.nan,\
            'Kp':np.nan}
    Embalse['n_Mann'][np.where(Embalse['type']==0)]=0.01
    Basin=[Ladera,Cauce,Carcava,Embalse]

    return Basin

#==============================================================================
# Create Basin Table
#==============================================================================
def GetBasinTable(flowvec, basin_scen, tipo_cuenca):
    u"""
    Construye la primera topología (sin ningún escenario) como DataFrame

    Entradas
    ----------
    flowvec : DataFrame
        Cada fila corresponde a una celda de la cuenca. Cada columna describe una propiedad de la cuenca de la siguiente manera:
        ID         : Id de la celda de la cuenca
        IDD        : ID de la celda de destino
        Area_Acum  : Área acumulada
        Lon_NoAcum : Longitud no acumulada. Igual 2^0.5*dx si el flujo es en diagonal (dir = 1,3,6 u 8) o dx si el flujo es recto (dir=)
        Elevacion  : Elevación de la celda en metros sobre el nivel del mar
        X 		   : Coordenada Este en algún sistema de referencia de coordenadas planas
        Y 		   : Coordenada Norte en algún sistema de referencia de coordenadas planas
        i 		   : indice de fila de la matriz del raster del modelo e elevación digital
        j          : indice de columna de la matriz del raster del modelo e elevación digital
        Pendiente  :  valor entre 0 y 1 correspondiente a la pendiente del terreno
        Tipo       : Tipo de celda [0: Ladera; 1:Drenaje; 2:drenaje efímero; 3:embalse, 4: celda de acumulación(presa)]
        Nodos      : 0: cauce , 1: cabeceras, 2: confluencias o nodos hidrológicos , 3: Nodos topográficos
        Elev_corr  : Elevación corregida en la red de drenaje
        ID_tramo   : ID de tramo
        Beechie    : Clasificación de las corrientes de acuerdo al patrón de alineamiento
        Flores     : Clasificación de las corrientes de acuerdo a la forma del lecho
        Qb         : Caudal Banca Llena
        Wb         : Ancho Banca Llena
        ETP        : Evapotranspiración Potencial
        manning    : Coeficiente de rugosidad
        Hb         : Altura Banca Llena
        alfa1      : Coeficiente de la relación Área Mojada - Velocidad media de flujo U=alfa1*Am**beta1
        beta1      : Exponente de la relación Área Mojada - Velocidad media de flujo U=alfa1*Am**beta1
        alfa2      : Coeficiente de la relación Caudal - Velocidad media del soluto Us=alfa2*Q**beta2
        beta2      : Exponente de la relación Caudal - Velocidad media del soluto Us=alfa2*Q**beta2
        alfa3      : Coeficiente de la relación Caudal - Velocidad máxima del soluto maxUs=alfa3*Q**beta3
        beta3      : Exponente de la relación Caudal - Velocidad máxima del soluto maxUs=alfa3*Q**beta3
        Qmedio     : Caudal medio de largo plazo 



    Tramos : DataFrame
        Cada fila corresponde a un tramo de la cuenca. Cada columna describe una propiedad del tramo de la siguiente manera:
        ID               : ID tramo
        Longitud         : Longitud del tramo
        Pendiente        : Pendiente del tramo
        Area_cuenca_prom : Area promedio de cuenca en el tramo
        IDD              : ID del tramo de destino
        ID_i             : ID celda de inicio del tramo
        ID_f             : ID celda de fin del tramo
        ID_des           : ID de la celda inicial del tramo de destino
        Area_ladera      : Área correpondiente de ladera
        Confinamiento    : Grado de confinamiento
        Qb               : Caudal Banca Llena
        Potencia 		 : Potencia del tramo
        Beechie          : Clasificación de las corrientes de acuerdo al patrón de alineamiento
        Flores           : Clasificación de las corrientes de acuerdo a la forma del lecho
        Wb         		 : Ancho Banca Llena
        manning   		 : Coeficiente de rugosidad
        Hb       	     : Altura Banca Llena
        alfa1    		 : Coeficiente de la relación Área Mojada - Velocidad media de flujo U=alfa1*Am**beta1
        beta1   	     : Exponente de la relación Área Mojada - Velocidad media de flujo U=alfa1*Am**beta1
        alfa2 		     : Coeficiente de la relación Caudal - Velocidad media del soluto Us=alfa2*Q**beta2
        beta2     		 : Exponente de la relación Caudal - Velocidad media del soluto Us=alfa2*Q**beta2
        alfa3   	     : Coeficiente de la relación Caudal - Velocidad máxima del soluto maxUs=alfa3*Q**beta3
        beta3   	     : Exponente de la relación Caudal - Velocidad máxima del soluto maxUs=alfa3*Q**beta3
        Qmedio 		     : Caudal medio de largo plazo 

    Salidas
    -------
    Basin :  DataFrame
    """
    #Clasifica las celdas entre 0: ladera, 1: cauce con drenado y 2: cauce con acumulación
    cauce = np.logical_or(flowvec['Tipo']==1,flowvec['Tipo']==2)
    cauce = np.logical_or(cauce,flowvec['Tipo']==3.1)
    cauce = np.logical_or(cauce,flowvec['Tipo']==3.2)
    cauce = np.logical_or(cauce,flowvec['Tipo']==4)

    basin_scen.loc[cauce,'tipo']              = 1.0
    basin_scen.loc[flowvec['Tipo']==4,'tipo'] = 2.0   # Cauce con acumulación
    basin_scen['embalse']                     = (flowvec['Tipo']>=3.0).astype(int)

    # Llena el objeto cuenca con el escenario de hidrología / sedimentos
    basin_scen.loc[:,['destino','tramo','X','Y','Z','L','S','alfa1','beta1','H5b','W5b']] = flowvec.loc[:,['IDD','ID_tramo','X','Y','Elevacion','Lon_NoAcum','Pendiente','alfa1','beta1','Hb','Wb']].values
    basin_scen.loc[:,'Q5b'] = flowvec['Qb'].values*86400.

    # Reproyecta las coordenadas de los centros de celda hacia el sistema WGS84
    # para obtener latitud y longitud

    # Importa el SRC Origen Nacional desde WKT
    from pyproj.crs import CRS
    proj_crs = CRS.from_wkt('PROJCS["MAGNA-SIRGAS / Origen-Nacional",GEOGCS["MAGNA-SIRGAS",DATUM["Marco_Geocentrico_Nacional_de_Referencia",SPHEROID["GRS 1980",6378137,298.257222101],TOWGS84[0,0,0,0,0,0,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4686"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",4],PARAMETER["central_meridian",-73],PARAMETER["scale_factor",0.9992],PARAMETER["false_easting",5000000],PARAMETER["false_northing",2000000],UNIT["metre",1],AUTHORITY["EPSG","9377"]]')
    magna=proj_crs
    wgs84=Proj(init='epsg:4326')
    long,lat=transform(magna,wgs84,flowvec['X'].values,flowvec['Y'].values)
    basin_scen.loc[:,'lat'] = lat; basin_scen.loc[:,'lon'] = long

    # Cambia el ID para que comience en 1, de acuerdo con la convención de FORTRAN
    basin_scen['destino'] = basin_scen['destino'] + 1

    if (tipo_cuenca == 'calidad'):
        # Llena el objeto cuenca con el escenario de calidad
        basin_scen.loc[:,['alfa2','beta2','alfa3','beta3']] = flowvec.loc[:,['alfa2','beta2','alfa3','beta3']].values
#        #Pone la elevación corregida en los cauces
#        basin_scen[cauce,50]=flowvec[cauce,12]
    return basin_scen

#==============================================================================
# Find closest channel network cell
#==============================================================================
def CeldaCauceCercana(flowvec, rad_bus = 300, coordenadas = False):
    u""" Funcion que busca puntos más cercanos a red de drenaje 
    teniendo en cuenta 2 criterios: distancia y direccion de drenaje 
    """
    flowvec['IDCCC'] = -9999
    for index, row in flowvec.iterrows():
        x    = row['X']
        y    = row['Y']
        id   = row['ID']
        tipo = row['Tipo']
        if tipo == 1:
            flowvec.loc[index,'IDCCC'] = id
        else:
            #Encontrar punto en la red de drenaje al que drenan el punto con coordenadas definidas 
            id_d = id
            while tipo==0:
                id_d = flowvec.loc[id_d,'IDD']
                tipo = flowvec.loc[id_d,'Tipo']
            dis_dre = np.sqrt((flowvec.X[id_d]-x)**2 + (flowvec.Y[id_d]-y)**2)
            #Encontrar punto en la red de drenaje con menor distancia a punto seleccionado en un radio determinado
            id = np.sqrt((flowvec.X-x)**2 + (flowvec.Y-y)**2) < rad_bus
            id = id[id==True].index
            id = flowvec.loc[id][flowvec.loc[id]['Tipo']==1].index 
            # definicion inicial celda donde drena 
            id_cauce = id_d
            try:
                dis_rad = np.sqrt((flowvec.loc[id].X-x)**2 + (flowvec.loc[id].Y-y)**2)
                id_r = dis_rad.loc[dis_rad == dis_rad.min()].index[0]
                dis_rad  = dis_rad.min()
                # si tienen misma distancia se toma punto con mayor area acumulada 
                if dis_dre == dis_rad:
                # se toma punto con mayor area si los dos estan a la misma distancia 
                    if flowvec.loc[id_d,'Area_Acum'] >= flowvec.loc[id_r,'Area_Acum'] :
                        id_cauce = id_d 
                    else:
                        id_cauce = id_r
                    distancia = dis_dre
                # se toma punto con menor distancia 
                elif dis_rad < dis_dre:
                    id_cauce = id_r
                    distancia = dis_rad
                else: 
                    id_cauce = id_d
                    distancia = dis_dre
            except:
                id_r = -9999
                distancia = dis_dre
            flowvec.loc[index,'IDCCC'] = id_cauce
    print('Columna IDCCC calculada exitosamente.')
    return flowvec

#==============================================================================
# Raster a objeto cuenca
#==============================================================================
def SampleRSTtoBSN(basin, raster_fp, name, types=[0,1,2,3]):
    """
    Cargar ráster en un objeto cuenca

    Parameters
    ----------
    basin:
        Objeto cuenca con columnas 'X', 'Y'

    raster_fp:
        file path of the raster to sample

    name:
        name of the new Basin field

    types:
        list of codes of the elements to sample. the codes are [0]hillslope, [1]drainage, [2]ephimeral drainage and [3]reservoir. default [0,1,2,3].
    """
    raster=ReadRaster(raster_fp)
    for n in types:
        X=basin[n]['X']
        Y=basin[n]['Y']
        i,j=XYtoij(X,Y,raster)
        if type(raster['mtrx'][0,0]) is int:
            basin[n][name]=raster['mtrx'][i,j].astype(np.int32)
        if type(raster['mtrx'][0,0]) is float:
            basin[n][name]=raster['mtrx'][i,j].astype(np.float32)
        else:
            basin[n][name]=raster['mtrx'][i,j]
    return basin

#==============================================================================
# Objeto cuenca a ráster
#==============================================================================
def SampleBSNtoRST(basin, baseRaster_fp, name, types=[0,1,2,3]):
    """
    Load raster data into a basin object

    Parameters
    ----------
    basin:
        Basin object with the basin 'x,y' field

    raster_fp:
        file path of the raster to sample

    name:
        name of the new Basin field

    types:
        list of codes of the elements to sample. the codes are [0]hillslope, [1]drainage, [2]ephimeral drainage and [3]reservoir. default [0,1,2,3].
    """
    baseRaster=ReadRaster(baseRaster_fp)
    if type(baseRaster['mtrx'][0,0]) is int:
        mtrx=(baseRaster['mtrx']*0 - 9999).astype(np.int32)
    if type(baseRaster['mtrx'][0,0]) is float:
        mtrx=(baseRaster['mtrx']*0. - 9999.).astype(np.float32)
    for n in types:
        X=basin[n]['X-Y'][:,0]
        Y=basin[n]['X-Y'][:,1]
        i,j=XYtoij(X,Y,baseRaster)
        mtrx[i,j]=basin[n][name]
    raster=BuildRaster(xll=baseRaster['xll'],yll=baseRaster['yll'],clsz=baseRaster['clsz'],nodt=baseRaster['nodt'],mtrx=mtrx)
    return raster

#==============================================================================
# Objeto cuenca a ráster
#==============================================================================
def SampleTBLtoRST(table, baseRaster, valcol,col_tipo= 'Tipo',types=[0,1,2,3.0,3.1,3.2], return_raster=False,raster_pth='',prcsn=gdal.GDT_Float32):
    """
    Guarda tabla a ráster

    Entradas
    ----------
    DataFrame:
        DataFrame con campos 'X' y 'Y'

    baseRaster:
        Ráster base

    valcol:
        columna con valores a muestrear
    col_tipo: str
        columna con la información de tipos de celda
    types:
        Lista de códigos, los códigos son 0: Ladera; 1:Drenaje; 2:drenaje efímero; 3:embalse, 4: celda de acumulación(presa)

    return_raster: Boolean
                Por defecto False. Si True el ráster se entrega

    raster_pth: string
                Nombre de archivo y directorio donde se va a guardar el ráster

    Salidas
    -------

    raster:
        Raster con los valores de la columna muestreada
    """

    mtrx = baseRaster['mtrx']*0 + baseRaster['nodt']
    for n in types:
        indices   = table.index[table[col_tipo]==n]
        X         = table.loc[indices,'X'].values
        Y         = table.loc[indices,'Y'].values
        i,j       = XYtoij(X,Y,baseRaster)
        mtrx[i,j] = table.loc[indices,valcol]
    raster = BuildRaster(xll=baseRaster['xll'],yll=baseRaster['yll'],clsz=baseRaster['clsz'],nodt=baseRaster['nodt'],mtrx=mtrx)
    if raster_pth != '':
        WriteRaster(raster_pth,raster,prcsn=prcsn)
    if return_raster==True:
        return raster
    else:
        return

#==============================================================================
# Ráster a DataFrame
#==============================================================================
def SampleRSTtoTBL(basin, raster_fp, xycols=[2,3], types=[0.,1.], coltypes=0,colind='',return_column=False):
    """
    Carga raster en columna de objeto cuenca

    Entradas
    ----------
    basin:
        DataFrame

    raster_fp:
        Dirección del ráster a ser usado

    xycols:
        Coordenadas 'X', 'Y' (no debería necesitarse)

    types:
        Lista de códigos, los códigos son 0: Ladera; 1:Drenaje; 2:drenaje efímero; 3:embalse, 4: celda de acumulación(presa)

    coltypes:
        columna donde están guardados lo índices con el tipo de celda.
    
    colind:
        columna del objeto cuenca a ser llenada con valores del ráster

    return_column:
        Si True, devuelve una columna con los datos del ráster
    """
    raster = ReadRaster(raster_fp)
    ChangeNoData(raster,-9999.0)
#     if colind == '':
#         if return_column==False:
#             colind  = basin.shape[1]
#             basin[] = np.insert(basin,colind,-9999.,axis=None)
    #Sample on table
    if return_column == False:
        for n in types:
            # coordinates to sample
            X = basin.loc[basin['Tipo'] == n,'X'].values
            Y = basin.loc[basin['Tipo'] == n,'Y'].values
            # calculate row and column indexes to sample
            i,j = XYtoij(X,Y,raster)
            # Sample0
            basin.loc[basin['Tipo'] == n,colind] = raster['mtrx'][i,j]
        return basin
    else:
        # coordinates to sample
        print(basin.shape)
        X   = basin['X'].values
        Y   = basin['Y'].values
        # calculate row and column indexes to sample
        i,j = XYtoij(X,Y,raster)
        return raster['mtrx'][i,j]

#==============================================================================
# Index where cummulate drainage Area changes in Basin Table
#==============================================================================

def IndUniqArea(basin):
    """
    Create an array with indices where accumulated area values changes
    """
    aux1=np.unique(basin[:,-1],return_counts=True)[1]
    indexes=aux1.cumsum()
    return indexes

#==============================================================================
# Crea el Objeto cuenca
#==============================================================================
def BuildBasin(areas,dir,dem,xy_c,raster = False,Uarea=1.,Ulong=300.,Uareap=.25,
               NoCellTramo=5,ventana=11,bsnmask_pth='',resmask_pth='',
               red_pth='',xyA_ini='',bsnobj_pth='',return_table=False,
               output_pth='',rivhead_pth='',sbn_pth='',percent=0.0,
               output_callback=None):
    u"""
    Construye un objeto cuenca como un DataFrame con las características
    mofológicas básicas para el estudio de la cuenca. Cada fila del DataFrame
    representa una celda de la cuenca.

    Entradas
    ----------
    areas: String
        Directorio o ráster donde se encuentra el raster de Areas de flujo acumulado

    dir: String
        Directorio o ráster donde se encuentra el raster de Direcciones de flujo

    dem: String
        Directorio o ráster donde se encuentra el raster del Modelo de Elevación Digital

    raster: Boolean
        True si se están ingresando los ráster ya leídos, False si se ingresan directorios para leer los ráster
        
    xy_c: array
        Coordenadas de cierre de la cuenca

    Uarea: float
        Umbral de área que define el drenaje

    Ulong: float
        Umbral de longitud para tramos de primer orden

    Uareap: float
        Umbral de área que distingue entre drenaje permanente y efímero

    NoCellTramo: integer
        Distancia mínima entre dos nodos topográficos.

    ventana: integer
        Amplitud de la media y mediana movil para suavizar los perfiles del canal

    bsnmask_pth: string
        Directorio donde se encuentra el raster con la máscara de la cuenca.
        Si '', la cuenca se traza usando xy_c.

    resmask_pth: string
        Directorio donde se encuentra el raster con la máscara de embalses.
        Si '', se asume que no existen embalses en la cuenca.

    red_pth: string
        Directorio donde se encuentra el raster la red de drenaje. Si '', la red se traza usando Uarea.
    
    xyA_ini: array
        Límites superiores de drenaje para cuencas incompletas. Cada fila representa una salida de subcuencas aguas arriba,
        con la primera columna siendo la coordenada X, la segunda la coordenada Y y la tercera es el area de la cuenca aferente,
        si la corrección de area no se necesita, la tercera columna será -9999.

    output_pth: string
        Directorio donde se guardan los archivos de salida
    
    rivhead_pth : str
        Directorio donde se encuentra el shapefile con las cabeceras
    
    sbn_pth : str
        Directorio donde se encuentra el raster de Areas de Interes, definidas como 1
    
    percent : Float
        Umbral del percentil de los tramos de longitud de primer orden

    Salidas
    -------
    Basin :  DataFrame con atributos de la cuenca

    outputs:
        Archivos ASCII con los mapas calculados de la morfología de la cuenca, guardados en output_pth.
    """
    #--------------------------------------------------------------------------
    #--DEFINE PRINT FUNCTION
    #--------------------------------------------------------------------------
    if output_callback is None:
        show_message = print_message
    else:
        show_message = output_callback

    #--------------------------------------------------------------------------
    #--Read disk data and check data sufficiency
    #--------------------------------------------------------------------------
    if raster == False:
        areas_rst = ReadRaster(areas)
        dir_rst   = ReadRaster(dir)
        dem_rst   = ReadRaster(dem)
    else:
        areas_rst = areas
        dir_rst   = dir
        dem_rst   = dem

    show_message ('Se cargaron los mapas')
    if red_pth != '':
        red_rst=ReadRaster(red_pth)
    elif red_pth == '' and Uarea == '':
        show_message('El Umbral de área es requerido si no se tiene ráster de red de drenaje')
        return
    if bsnmask_pth == '' and len(xy_c)==0:
        show_message('Las coordenadas de salida de la cuenca se requieren, si no se tiene un ráster con una máscara de cuenca')
        return

    #--------------------------------------------------------------------------
    #--Create Basin mask
    #--------------------------------------------------------------------------
    # Load Basin mask
    if bsnmask_pth != '':
        bsnmask_rst = ReadRaster(bsnmask_pth)
        bsnmask_rst = SampleRastToRast(rasterA=bsnmask_rst,rasterB=dir_rst)
        show_message('Se cargó la máscara de cuenca')

    # Create Basin mask
    elif bsnmask_pth == '':
        show_message('Se empieza la creación de la máscara de cuenca')
        bsnmask_rst = TraceBasin(X=xy_c[0,0],Y=xy_c[0,1],dir_rst=dir_rst,return_raster=False,return_mask=True,XY_ini=xyA_ini)
        print('Coordenadas cierre de la cuenca')
        print(xy_c[0,0],xy_c[0,1])
        show_message('Se trazó la máscara de cuenca')

    # Calculate number of cells inside the basin
    numbsncells = (bsnmask_rst['mtrx'] != bsnmask_rst['nodt']).sum()
    if numbsncells > 20:
        show_message('Una máscara de cuenca con ' + str(numbsncells) + ' celdas es usada')
    else:
        show_message('ADVERTENCIA: Traced Basin tiene menos de 20 celdas. Por favor revisa que las coordenadas de salida de la cuenca están sobre la red de drenaje. De igual manera, la cuenca fue trazada.')
    import matplotlib.pyplot as plt
    plt.imshow(bsnmask_rst['mtrx'])

    # Save basin mask
    if output_pth != '':
        if bsnmask_pth == '':
            WriteRaster(path=os.path.join(output_pth,'bsnmask.tif'),rst=bsnmask_rst)
            show_message('Se guardó la máscara de cuenca')
            
    #--------------------------------------------------------------------------
    #--Se crea la máscara de embalse
    #--------------------------------------------------------------------------
    if resmask_pth != '':
        resmask_rst = ReadRaster(resmask_pth)
        resmask_rst = SampleRastToRast(rasterA=resmask_rst,rasterB=dem_rst)
        WriteRaster(path=os.path.join(output_pth,'resmask.tif'),rst=resmask_rst) # Guarda la mascara remuestreada
        show_message('Se cargó máscara de embalse')
    else:
        show_message('ADVERTENCIA: Máscara de embalse no fue encontrada. Se asume que no existen embalses en la cuenca')

    #--------------------------------------------------------------------------
    #--Create flow vector
    #--------------------------------------------------------------------------
    flowvec, AuxReajuste = FlowVector(dem_rst=dem_rst,dir_rst=dir_rst,areas_rst=areas_rst,bsnmask_rst=bsnmask_rst,Uarea=Uarea,xyA_ini=xyA_ini)
    show_message('Se creó DataFrame con los atributos de la cuenca')

    #--------------------------------------------------------------------------
    #--Corrects negative slopes in the whole raster
    #--------------------------------------------------------------------------
    flowvec = SlopeFit(FlowVectorA=flowvec)
    show_message('Se corrige la pendiente en todas las celdas de la cuenca')
    
    #--------------------------------------------------------------------------
    #--Trazar red de drenaje desde las cabeceras cargadas desde el shapefile 
    #--------------------------------------------------------------------------
    if rivhead_pth != '':
        flowvec = TraceRiverheads(flowvec,rivhead_pth,bsnmask_rst)
        show_message('Se trazó red de drenaje desde las cabeceras')
    
    #--------------------------------------------------------------------------
    #--Calcular nodos hidrológicos
    #--------------------------------------------------------------------------
    #Identificar Cabeceras
    flowvec = IdCab(FlowVectorA=flowvec)
    show_message('Se marcaron las cabeceras de la cuenca como 1')
    
    # Marca como confluencias los puntos donde inicia la red de la cuenca , si es
    # una cuenca incompleta
    if len(xyA_ini) >0:
        for n in range(xyA_ini.shape[0]):
            ij_des = IJ_des(xyA_ini[n,:],dir_rst)
            inicio = flowvec.index[(flowvec['i'] == ij_des[0]) & (flowvec['j'] == ij_des[1])]
            flowvec.loc[inicio,'Nodos'] = 2.0
    
    #Identificar confluencias
    flowvec                   = HydrologicNodes(FlowVectorA=flowvec)
    flowvec['Nodos'].iloc[-1] = 2.0 # Marca como confluencia la última celda
    show_message('Se marcaron confluencias')
    
    #--------------------------------------------------------------------------
    #--Corrección de la red de drenaje de acuerdo al umbral de longitud
    #--------------------------------------------------------------------------
    if percent == 0.0:
        flowvec = UmbralLong(flowvec,Ulong)
    else:
        flowvec = PercUmbralLong(flowvec, bsnmask_rst,Ulong,percent,sbn_pth)
    
    # - Después de la corrección por umbral de longitud se deben volver a
    # identificar las cabeceras y las confluencias en caso de que haya
    # confluencias de más de 2 corrientes y nuevas cabeceras donde antes había
    # nodos hidrológicos
    # Borrar confluencias (Nodos hidrológicos) viejas
    flowvec.loc[flowvec['Nodos'] == 2,'Nodos'] = 0.0
    # Recalcular cabeceras y confluencias
    flowvec = IdCab(FlowVectorA=flowvec)
    flowvec = HydrologicNodes(FlowVectorA=flowvec)
    show_message('Se corrigió la red de drenaje con el umbral de longitud')
    
    #--------------------------------------------------------------------------
    #--Calculo de nodos topográficos
    #--------------------------------------------------------------------------
    # cabeceras donde se corta la cuenca
    if len(xyA_ini) > 0:
        for n in range(xyA_ini.shape[0]):
            # calcular indices de ráster de celda de destino del inicio de cuenca
            ij_des=IJ_des(xyA_ini[n,:],dir_rst)
            # calcular id de celda de inicio
            celda = flowvec.index[(flowvec['i'] == ij_des[0]) & (flowvec['j'] == ij_des[1])]
            inicio = flowvec.loc[celda,'ID']
            # Marcar como cabecera el inicio de la cuenca
            flowvec.loc[inicio,'Nodos']=1.0
            
    # Nodos Topográficos (Knick nodes)
    flowvec = KnickNodes(FlowVectorA=flowvec,NoCllTram=NoCellTramo,ventana=ventana,peso_media=0.8,tol=0.01)
    show_message('Se calcularon los nodos topográficos y se corrigió la elevación de las celdas de la cuenca')
    #--------------------------------------------------------------------------
    #--Mark ephimeral drainage cells
    #--------------------------------------------------------------------------
    #flowvec=Carcavas(flowvec,Uarea,Uareap)
    #--------------------------------------------------------------------------
    #--Se calculan las propiedades de Tramos
    #--------------------------------------------------------------------------
    # Calcular tramos y propiedades de Tramos
    flowvec,reachesvec = GetReaches(flowvec)
    show_message('Función GetReaches finalizó')
    # Cambiar pendientes en cauce por pendiente tramo
    flowvec = ReachtoFvec(reachesvec,flowvec,'Pendiente')
    show_message('Se cambiaron las pendientes de las celdas de de drenaje por pendientes de los tramos')
    # calcular laderas aferentes
    flowvec = DefinirLaderas(ArregloFlujo=flowvec)
    show_message('Función DefinirLaderas finalizó')
    # confinamiento
    red_rst = DrainageMap(dir_rst,flowvec,AuxReajuste,'Nodos',os.path.join(output_pth,'red.asc'))
    show_message('Se construyó ráster de red')
    reachesvec, MapGeoLLANURA, MapGeoTRENZADO, MapGeoVALLE, MapRedVirtual = FloodPlane(FlowVectorA=flowvec,TramosA=reachesvec,dem_rst=dem_rst,red_rst=red_rst)
    show_message('Función FloodPlane finalizó')

    # Clasificación de corrientes de acuerdo a Beechie y Flores
    flowvec,reachesvec=GetMorpho(flowvec,reachesvec)
    show_message('Función GetMorpho finalizó')
    #--------------------------------------------------------------------------
    #-- Geometría hidráulica
    #--------------------------------------------------------------------------
    reachesvec,flowvec=GeoHidraulica(flowvec,reachesvec)
    show_message('Función GeoHidraulica finalizó')
    #--------------------------------------------------------------------------
    #--Calculo de ETP
    #--------------------------------------------------------------------------
    #Se agrega ETP a flowvec
    flowvec['ETP'] = ETPcenicafe(elev=flowvec['Elevacion'])
    show_message('Función ETPcenicafe finalizó')
    #--------------------------------------------------------------------------
    #--Se calcula Manning en la red de drenaje
    #--------------------------------------------------------------------------
    #Se agrega el coeficiente de manning en el DataFrame de tramos
    reachesvec['manning'] = ManningCauce(reachesvec['Pendiente'])
    show_message('Función ManningCauce finalizó')
    #Se agrega el coeficiente de manning en el DataFrame flowvec
    flowvec['manning'] = -9999.0
    flowvec=ReachtoFvec(reachesvec,flowvec,'manning')
    show_message('Se agrega el coeficiente de Manning al DataFrame de atributos de la cuenca')
    
    #--------------------------------------------------------------------------
    #--Calcular Hb
    #--------------------------------------------------------------------------
	#Se agrega profundidad de banca llena en el DataFrame de tramos reachesvec
    reachesvec['Hb'] = BankfullH(Qb=reachesvec['Qb'].values, So=reachesvec['Pendiente'].values, \
                                 Wb=reachesvec['Wb'].values, n=reachesvec['manning'].values)
    show_message('Función BankfullH finalizó')
	#Se agrega profundidad de banca llena en el DataFrame flowvec
    flowvec = ReachtoFvec(reachesvec,flowvec,'Hb')
    show_message('Se completan parámetros de banca llena en el DataFrame de atributos de la cuenca')
    
    #--------------------------------------------------------------------------
    #--Cálculos de geometría hidráulica transversal y parámetros MDLC-ADZ
    #--------------------------------------------------------------------------
    reachesvec['alfa1'],reachesvec['beta1'],reachesvec['alfa2'],reachesvec['beta2'],reachesvec['alfa3'],reachesvec['beta3'] = \
    ParametrosMDLCADZ(Hb=reachesvec['Hb'].values,Wb=reachesvec['Wb'].values,S0_in=reachesvec['Pendiente'].values,nb=reachesvec['manning'].values,L_in=reachesvec['Longitud'].values,Qb=reachesvec['Qb'].values)#,num_pun=len(reachesvec))
    
    # alfa1 y beta1 : Coeficiente y exponente de la relación Área Mojada - Velocidad media de flujo
    flowvec = ReachtoFvec(reachesvec,flowvec,'alfa1')
    flowvec = ReachtoFvec(reachesvec,flowvec,'beta1')
    
    # alfa2 y beta 2 : Coeficiente y exponente de la relación Caudal - Velocidad media del soluto
    flowvec = ReachtoFvec(reachesvec,flowvec,'alfa2')
    flowvec = ReachtoFvec(reachesvec,flowvec,'beta2')
    
    # alfa3 y beta 3 : Coeficiente y exponente de la relación Caudal - Velocidad máxima del soluto
    flowvec = ReachtoFvec(reachesvec,flowvec,'alfa3')
    flowvec = ReachtoFvec(reachesvec,flowvec,'beta3')  
    
    #--------------------------------------------------------------------------
    #--Cálculos del caudal medio de largo plazo
    #--------------------------------------------------------------------------
    reachesvec['Qmedio'] = Qmedio(reachesvec['Area_cuenca_prom'].values)
    flowvec = ReachtoFvec(reachesvec,flowvec,'Qmedio')
    
    flowvec = flowvec.fillna(-9999.0)

    #--------------------------------------------------------------------------
    #--Se guardan los tramos como shapefile
    #--------------------------------------------------------------------------    
    if output_pth != '':
        # Guardar tramos como shape
        Reaches2Shp(flowvec,reachesvec,output_pth=output_pth)
        show_message('reachesvec se guardó como shapefile')
        
    #--------------------------------------------------------------------------
    #--Se marcan las celdas tipo embalse como tipo '3.'
    #--------------------------------------------------------------------------
    if resmask_pth != '':
        aux1=SampleRSTtoTBL(basin=flowvec,raster_fp=os.path.join(output_pth,'resmask.tif'), types=[0.,1.],colind='raster',return_column=True)
        flowvec.loc[aux1>0.0,'Tipo'] = (flowvec.loc[aux1>0.0,'Tipo']/10.0)+3.0
        print(np.unique(flowvec['Tipo']))
       
        if len(flowvec[flowvec['Tipo']>=3.0])>0 and resmask_pth != '':
            # Hace que la celda de pie de presa drene a sí misma
            flowvec = CellReservoir(flowvec,resmask_rst)
            # Hace que las celdas de cauce entre celdas de embalse se marquen como celdas de embalse
            flowvec = CorrReservoir(flowvec)
            show_message('Se marcan celdas de embalse en el DataFrame de atributos de la cuenca')
            print(np.unique(flowvec['Tipo']))

    #--------------------------------------------------------------------------
    #--Se adiciona un campo que almacene la celda de cauce más cercana
    #--------------------------------------------------------------------------
    flowvec = CeldaCauceCercana(flowvec, rad_bus = 300, coordenadas = False)

    #--------------------------------------------------------------------------
    #--Build Basin
    #--------------------------------------------------------------------------
    #if return_table==True:
    #    basin=GetBasinTable(flowvec,reachesvec)
    #    show_message('basin table created')
    #else:
    #    basin=GetBasin(flowvec,reachesvec)
    #    show_message('basin object created')
    
    #--------------------------------------------------------------------------
    #--Return and Save outputs
    #--------------------------------------------------------------------------
    # Se guardan mapas y archivos de salida
    if output_pth != '':
        print(np.unique(flowvec['Tipo']))
        # Save drainage network, nodes, reaches, flores class, beechie class
        DrainageMap(dir_rst,FlowVectorA=flowvec,AuxReajuste=AuxReajuste,topic='Tipo',path=os.path.join(output_pth,'net.tif'))
        show_message('Se guardó mapa de red de drenaje')
        DrainageMap(dir_rst,FlowVectorA=flowvec,AuxReajuste=AuxReajuste,topic='Nodos',path=os.path.join(output_pth,'nodes.tif'))
        show_message('Se guardó mapa con información de los nodos')
        DrainageMap(dir_rst,FlowVectorA=flowvec,AuxReajuste=AuxReajuste,topic='ID_tramo',path=os.path.join(output_pth,'reaches.tif'))
        show_message('Se guardó mapa con información de los tramos')
        DrainageMap(dir_rst,FlowVectorA=flowvec,AuxReajuste=AuxReajuste,topic='Beechie',path=os.path.join(output_pth,'beechie.tif'))
        show_message('Se guardó mapa con clasificación Beechie')
        DrainageMap(dir_rst,FlowVectorA=flowvec,AuxReajuste=AuxReajuste,topic='Flores',path=os.path.join(output_pth,'flores.tif'))
        show_message('Se guardó mapa con clasificación Flores')
        # Save MapGeoLLANURA, MapGeoTRENZADO, MapGeoVALLE, MapRedVirtual
        WriteRaster(os.path.join(output_pth,'MapGeoLLANURA.tif'),MapGeoLLANURA)
        show_message('Se guardó mapa MapGeoLLANURA')
        WriteRaster(os.path.join(output_pth,'MapGeoTRENZADO.tif'),MapGeoTRENZADO)
        show_message('Se guardó mapa MapGeoTRENZADO')
        WriteRaster(os.path.join(output_pth,'MapGeoVALLE.tif'),MapGeoVALLE)
        show_message('Se guardó mapa MapGeoVALLE')
        WriteRaster(os.path.join(output_pth,'MapRedVirtual.tif'),MapRedVirtual)
        show_message('Se guardó mapa MapRedVirtual')
        SampleTBLtoRST(table=flowvec,baseRaster=dem_rst,valcol='ID_tramo',raster_pth=os.path.join(output_pth,'watersheds.tif'),prcsn=gdal.GDT_Int32)
        show_message('Se guardó mapa watersheds')
        
        # Se adiciona el campo correspondiente al ID de la celda de cauce más cercana en flowvec
        flowvec = CeldaCauceCercana(flowvec, rad_bus = 300, coordenadas = False)
         
        # Se guarda flowvec
        flowvec.to_csv(os.path.join(output_pth,'basin.txt'),index=False)
        show_message('Se guardó archivo con los atributos de la cuenca como:')
        show_message(os.path.join(output_pth,'basin.txt'))
        print(np.unique(flowvec['Tipo']))
       

        # Se guarda reachesvec
        reachesvec.to_csv(os.path.join(output_pth,'reachesvec.txt'),index=False)
        show_message('Se guardó archivo con los atributos de los tramos como:')
        show_message(os.path.join(output_pth,'reachesvec.txt'))

        
        return flowvec

#==============================================================================
# Load scenary
#==============================================================================
def Load_scen(output_pth,lista_rutas_mapas,tipo_cuenca,basin='',version="",bsninput_pth='',output_callback = None):
    """
    Carga un escenario al objeto de cuenca de la morfología base, y calcula la
    Erodabilidad de Banca del escenario, si se crea un objeto de cuenca del tipo
    "sedimentos" o "calidad".
    
    Entradas
    ----------
    
    basin : DataFrame
        Objeto de cuenca que contiene la morfología base, si está cargado en memoria.
        Contiene los siguientes campos (columnas) para cada pixel (fila):
            https://docs.google.com/spreadsheets/d/1Uaqax-U_87OueKhfPQC0EAzdjo9NKz_J36HCur_JgZs/edit?usp=sharing
    
    output_pth : string
        Ruta donde se guarda el archivo de texto con el objeto cuenca y el mapa de erodabilidad de la banca Eb5
    
    bsninput_pth : string
        Ruta del archivo de cuenca de la morfología Base si se quiere cargar uno preexistente
    
    lista_rutas_mapas : Diccionario
        Lista de rutas de los mapas de variación de los parámetros del modelo para un
        escenario. Define si se creará un objeto de cuenca hidroligia, sedimentos o calidad,
        de acuerdo con el número de mapas de entrada.
        Recibe una ruta para el ráster de MapGeoVALLE.tif y de direcciones

    tipo_cuenca : str    
        hidrologia = [SU0, SU1, SU3, fks, ks, kp, perd, n_manning]
        sedimentos = hidrologia + [Arenas, Limos, Arcillas, K_usle, C_usle, P_usle, MapGeovalle]
        calidad = hidrologia + sedimentos + [ccuINCA,cargaCT,cargaNO3,cargaNH4,cargaNO]
    
    output_callback : 
        Define la forma en que se imprimen las advertencias de la simulación en la interfaz.
        Por defecto es igual a "None".
    
    Return
    ------
    
    basin_scen : DataFrame
        Array con las propiedades del objeto de cuenca, cuyas columnas representan las
        variables descritas en siguiente enlace:
            https://docs.google.com/spreadsheets/d/1NrePOh36NKMUHjfhpYR-Pab8urBJjdeUep_jhi6eeJY/edit?usp=sharing
    
    """
    #==========================================================================
    # Define la forma en que se muestran los mensajes en la interfaz
    #==========================================================================
    if output_callback is None:
        show_message = print_message
    else:
        show_message = output_callback
    
    #==========================================================================
	# Cargar objeto cuenca
    #==========================================================================
    if np.all(basin == ''):
        if bsninput_pth == '':
            # No se cargó información suficiente
            show_message("No se encuentra archivo de cuenca de Morfología Base")
        else:
            # Cargar archivo de cuenca de Morfología Base
            basin = pd.read_csv(bsninput_pth)
            show_message("Lectura de archivo de cuenca completa")
    else:
        show_message("Archivo de cuenca en memoria")
    
    #==========================================================================
    # Revisar validez del objeto cuenca base recibido
    #==========================================================================
    try:
        num_col_base = basin.shape[1]
    except:
        # El archivo de cuenca cargado no es válido
        show_message("El archivo de cuenca cargado no es válido")
    
    if num_col_base != 28:
        # El archivo de cuenca cargado no es válido
        show_message("El archivo de cuenca cargado no es válido")
    else:
        # El archivo de cuenca es válido
        show_message("El archivo de cuenca cargado es válido")
    
    #==========================================================================
    # Crear arreglo de numpy del archivo cuenca con escenario y llenar datos base
    #==========================================================================
    #Calcular número de celdas
    num_celdas = basin.shape[0]
    
    dtypes={'tipo':'int32','destino':'int32','tramo':'int32','llanura':'int32','embalse':'int32','X':'float64',
                    'Y':'float64','Z':'float64','lat':'float64','lon':'float64','L':'float64','S':'float64','D':'float64',
                    'alfa1':'float64','beta1':'float64','S0':'float64','S1':'float64','S2':'float64','S3':'float64',
                    'S4':'float64','S5':'float64','H5b':'float64','W5b':'float64','Q5b':'float64','HU':'float64','LAI':'float64'}
    sediDtypes = {'arcS2S':'float64','limS2S':'float64','areS2S':'float64','arcS2D':'float64','limS2D':'float64',
                'areS2D':'float64','arcS5S':'float64','limS5S':'float64','areS5S':'float64','arcS5D':'float64',
                'limS5D':'float64','areS5D':'float64'}
    calDtypes={'alfa2':'float64','beta2':'float64','alfa3':'float64','beta3':'float64','S0EC':'float64','S1ECp':'float64',
                'S1ECs':'float64','S2EC':'float64','S3EC':'float64','S4EC':'float64','S0NO':'float64','S1NOp':'float64',
                'S1NOs':'float64','S2NO':'float64','S3NO':'float64','S4NO':'float64','S0NH4':'float64','S1NH4p':'float64',
                'S1NH4s':'float64','S2NH4':'float64','S3NH4':'float64','S4NH4':'float64','S0NO3':'float64','S1NO3p':'float64',
                'S1NO3s':'float64','S2NO3':'float64','S3NO3':'float64','S4NO3':'float64','S0PO':'float64','S1POp':'float64',
                'S1POs':'float64','S2PO':'float64','S3PO':'float64','S4PO':'float64','S0PI':'float64','S1PIp':'float64',
                'S1PIs':'float64','S2PI':'float64','S3PI':'float64','S4PI':'float64','S0PO_fb':'float64','S1POp_fb':'float64',
                'S1POs_fb':'float64','S2PO_fb':'float64','S3PO_fb':'float64','S0PI_fb':'float64','S1PIp_fb':'float64',
                'S1PIs_fb':'float64','S2PI_fb':'float64','S3PI_fb':'float64','OD':'float64','CDBO':'float64','CE':'float64',
                'EC':'float64','NO3':'float64','NH4':'float64','NO':'float64','PO':'float64','PI':'float64','PT':'float64',
                'CA':'float64','alk':'float64','pH':'float64'}
    #Crear cuenca base
    if tipo_cuenca == 'hidrologia':
        # crear DataFrame del tamaño necesario
        basin_scen = pd.DataFrame(0,columns = ['tipo','destino','tramo','llanura','embalse','X','Y','Z','lat','lon','L','S','D',\
                                             'alfa1','beta1','S0','S1','S2','S3','S4','S5','H5b','W5b','Q5b','HU','LAI'],
                                    index = np.arange(num_celdas))
        show_message("Inicia parametrización hidrología")
    elif tipo_cuenca == 'sedimentos':
        # se modifica la operación para python 3.8, verificar en otros entornos
        #dtypes = dtypes|sediDtypes
        dtypes.update(sediDtypes)
        # crear DataFrame del tamaño necesario
        basin_scen = pd.DataFrame(0,columns = ['tipo','destino','tramo','llanura','embalse','X','Y','Z','lat','lon','L','S',\
                                             'D','alfa1','beta1','S0','S1','S2','S3','S4','S5','H5b','W5b','Q5b','HU','LAI',\
                                             'arcS2S','limS2S','areS2S','arcS2D','limS2D','areS2D','arcS5S','limS5S','areS5S',\
                                             'arcS5D','limS5D','areS5D'],
                                    index = np.arange(num_celdas))
        show_message("Inicia parametrización sedimentos")
    elif tipo_cuenca == 'calidad':
        # se modifica la operación para python 3.8, verificar en otros entornos
        #dtypes= dtypes|sediDtypes
        #dtypes= dtypes|calDtypes
        dtypes.update(sediDtypes)
        dtypes.update(calDtypes)
        # crear DataFrame del tamaño necesario
        basin_scen = pd.DataFrame(0, columns= ['tipo','destino','tramo','llanura','embalse','X','Y','Z','lat','lon','L','S','D',\
                                           'alfa1','beta1','S0','S1','S2','S3','S4','S5','H5b','W5b','Q5b','HU','LAI','arcS2S',\
                                           'limS2S','areS2S','arcS2D','limS2D','areS2D','arcS5S','limS5S','areS5S','arcS5D',\
                                           'limS5D','areS5D','alfa2','beta2','alfa3','beta3','S0EC','S1ECp','S1ECs','S2EC',\
                                           'S3EC','S4EC','S0NO','S1NOp','S1NOs','S2NO','S3NO','S4NO','S0NH4','S1NH4p','S1NH4s',\
                                           'S2NH4','S3NH4','S4NH4','S0NO3','S1NO3p','S1NO3s','S2NO3','S3NO3','S4NO3','S0PO',\
                                           'S1POp','S1POs','S2PO','S3PO','S4PO','S0PI','S1PIp','S1PIs','S2PI','S3PI','S4PI',\
                                           'S0PO_fb','S1POp_fb','S1POs_fb','S2PO_fb','S3PO_fb','S0PI_fb','S1PIp_fb','S1PIs_fb',\
                                           'S2PI_fb','S3PI_fb','OD','CDBO','CE','EC','NO3','NH4','NO','PO','PI','PT','CA','alk','pH'], 
                                   index = np.arange(num_celdas))
        show_message("Inicia parametrización calidad")
    else:
        #- El tipo de cuenca no es válido
        show_message("Tipo de cuenca no es válido")
    basin_scen = GetBasinTable(basin, basin_scen, tipo_cuenca)
    #==========================================================================
    # Agregación de datos de escenario que se leen como mapas desde la interfaz
    #==========================================================================
    #Cargar mapa de llanura
    basin_scen['llanura'] = ReadFromRSTIJ(basin,lista_rutas_mapas['llanura'])

    #Gradiente direccion topográfico DGT
    ddir=ReadFromRSTIJ(basin,lista_rutas_mapas['dir'])
    vgdt = {
            1: 225.,
            2: 270.,
            3: 315.,
            4: 180.,
            6: 0.,
            7: 135.,
            8: 90.,
            9: 45.,
            }
    for zz in range(len(ddir)):
        ddir[zz]=vgdt[ddir[zz]]
    basin_scen['D'] = ddir

    #Ajustar para que las variables exclusivas de cauce queden con valor -9999 en ladera
    basin_scen.loc[basin_scen['tipo'] == 0, ['tramo','alfa1','beta1','H5b','W5b','Q5b']] = -9999.0
    print("Se configuran con valor de -9999 en ladera las variables exclusivas de cauce")
    
    #==========================================================================
    # Escritura del archivo cuenca con escenario
    #==========================================================================
    basin_pth = output_pth+'.txt'
    # Write table
    dA = np.min(basin_scen['L'])**2.
    # Encabezado con metadatos
    if tipo_cuenca == 'hidrologia':
        header='[NÚMERO DE CELDAS]\n'+str(num_celdas)+'\n\n[ÁREA DE LAS CELDAS]\n'+str(dA)+'\n\n[TIPO DE TOPOLOGÍA]\n'+version+'\n\n[MATRIZ DE VARIABLES]\n'+'tipo destino tramo llanura embalse X Y Z lat lon L S D alfa1 beta1 S0 S1 S2 S3 S4 S5 H5b W5b Q5b HU LAI'
    elif tipo_cuenca == 'sedimentos':
        header='[NÚMERO DE CELDAS]\n'+str(num_celdas)+'\n\n[ÁREA DE LAS CELDAS]\n'+str(dA)+'\n\n[VERSIÓN DE TOPOLOGÍA]\n'+version+'\n\n[MATRIZ DE VARIABLES]\n' +'tipo destino tramo llanura embalse X Y Z lat lon L S D alfa1 beta1 S0 S1 S2 S3 S4 S5 H5b W5b Q5b HU LAI arcS2S limS2S areS2S arcS2D limS2D areS2D arcS5S limS5S areS5S arcS5D limS5D areS5D'
    elif tipo_cuenca == 'calidad':
        header='[NÚMERO DE CELDAS]\n'+str(num_celdas)+'\n\n[ÁREA DE LAS CELDAS]\n'+str(dA)+'\n\n[TIPO DE TOPOLOGÍA]\n'+version+'\n\n[MATRIZ DE VARIABLES]\n'+'tipo destino tramo llanura embalse X Y Z lat lon L S D alfa1 beta1 S0 S1 S2 S3 S4 S5 H5b W5b Q5b HU LAI arcS2S limS2S areS2S arcS2D limS2D areS2D arcS5S limS5S areS5S arcS5D limS5D areS5D alfa2 beta2 alfa3 beta3 S0EC S1ECp S1ECs S2EC S3EC S4EC S0NO S1NOp S1NOs S2NO S3NO S4NO S0NH4 S1NH4p S1NH4s S2NH4 S3NH4 S4NH4 S0NO3 S1NO3p S1NO3s S2NO3 S3NO3 S4NO3 S0PO S1POp S1POs S2PO S3PO S4PO S0PI S1PIp S1PIs S2PI S3PI S4PI S0PO_fb S1POp_fb S1POs_fb S2PO_fb S3PO_fb S0PI_fb S1PIp_fb S1PIs_fb S2PI_fb S3PI_fb OD CDBO CE EC NO3 NH4 NO PO PI PT CA alk pH'
    # Cuerpo con los valores
    file = open(basin_pth,'w')
    file.write(header)

    #Transformación para que los datos se escriban con el menor número posible de cifras significativas
    #(así el archivo de la cuenca pesa menos y fortran se demora menos leyendola)
    basin_scen = basin_scen.astype(dtypes)
    basin_scen = basin_scen.astype(str)
    for index in basin_scen.index:
        row = basin_scen.loc[index,:].tolist()
        line = '\n'+' '.join(row)
        file.write(line)
    file.close()

    return basin_scen

#==============================================================================
# Corrección de áreas
#==============================================================================

def IJ_des(xyA_ini,dir_rst):
    """
    Cálcula la fila y columna de la celda destino

    Parámetros
    ----------
    xyA_ini : array
        Array de coordenadas x,y de la celda actual

    dir_rst : ráster
        ráster de direcciones

    Retorna
    -------
    ij_des : array
        Array con la fila y columna de la celda destino
    """
    aux_ij=np.array([[1,1,1,0,0,0,-1,-1,-1],[-1,0,1,-1,0,1,-1,0,1]])
    i,j=XYtoij(X=xyA_ini[0],Y=xyA_ini[1],raster=dir_rst)
    i_new=i+aux_ij[0,dir_rst['mtrx'][i,j].astype(int)-1]
    j_new=j+aux_ij[1,dir_rst['mtrx'][i,j].astype(int)-1]
    ij_des=np.array([i_new,j_new])
    return ij_des
#==============================================================================
def AreasFit(FlowVectorA,xyA_ini,dir_rst,areas_rst):
    """
    Corrige las áreas acumuladas de una cuenca incompleta, tomando el offset del valor que tiene y el que realmente 
    debería tener de acuerdo a lo ingresado por el usuario 
    """
    for n in range(xyA_ini.shape[0]):
        i_fin,j_fin = XYtoij(X=xyA_ini[n,0],Y=xyA_ini[n,1],raster=areas_rst) #Transforma x,y a i,j
        ij_des      = IJ_des(xyA_ini[n,:],dir_rst)
        inicio      = FlowVectorA.index[(FlowVectorA['i']==ij_des[0]) & (FlowVectorA['j']==ij_des[1])]
        FlowVectorA.loc[inicio,'Nodos'] = 2.0   # nodo hidrológico
        
        if xyA_ini[n,2] != -9999:
            corr     = xyA_ini[n,2] - areas_rst['mtrx'][i_fin,j_fin]
            cell     = inicio
            nextcell = int(FlowVectorA.loc[cell,'IDD'])
            while cell != -9999: # En este punto la celda de cierre de cuenca no se ha marcado como que drene a sí misma, si no que tiene 'IDD' = -9999
                FlowVectorA.loc[cell,'Area_Acum'] = FlowVectorA.loc[cell,'Area_Acum']+corr
                cell = int(FlowVectorA.loc[cell,'IDD'])
                nextcell = int(FlowVectorA.loc[cell,'IDD'])
            FlowVectorA.loc[cell,'Area_Acum'] = FlowVectorA.loc[cell,'Area_Acum']+corr
            
    return FlowVectorA
#==============================================================================
# Ebanca
#==============================================================================
def EbancaRST(MapGeoVALLE_pth,C_pth,K_pth,Erst_pth,return_Erst=False):
    MapGeoVALLE=ReadRaster(MapGeoVALLE_pth)
    C=ReadRaster(C_pth)
    K=ReadRaster(K_pth)
    valles=np.unique(MapGeoVALLE['mtrx'])[1:]
    Erst=copy.copy(MapGeoVALLE)
    Erst['mtrx']=Erst['mtrx']*0.0 + Erst['nodt']
    # Remuestrear raster al mismo tamaño y resolución
    C=SampleRastToRast(C,Erst)
    K=SampleRastToRast(K,Erst)
    # Calcular erodabilidad
    for k in valles:
        i,j=np.where(MapGeoVALLE['mtrx']==k)
        Erst['mtrx'][i,j]=np.mean(C['mtrx'][i,j])*np.mean(K['mtrx'][i,j])
    WriteRaster(Erst_pth, Erst)
    if return_Erst==False:
        return
    else:
        return Erst
#==============================================================================
# Ancho seccion parabolica en funcion del area y parametros de banca (Wb - Hb)
#==============================================================================

def wofA(A, Wb, Hb):
    """
    Channel width for a parabolic cross-section given the flow area and the geometry of bankfull conditions (width and depth).

    Entradas
    ----------
    A : array
        Area or array of areas of the speciffic condition to calculate channel width

    Wb: array
        Bankfull conditions width or array of widths

    Hb: array
        Bankfull conditions depth or array of depths

    Salidas
    -------
    w :
        Width or array of widths of channel surface for the given area and Entradas
    """
    w=(1.5*(Wb**2.)*(1./Hb)*A)**(1./3.)
    return w

#==============================================================================
# Profundidad seccion parabolica en funcion del area y parametros de banca (Wb - Hb)
#==============================================================================

def hofA(A, Wb, Hb):
    """
    Channel depth for a parabolic cross-section given the flow area and the geometry of bankfull conditions (width and depth).

    Entradas
    ----------
    A : array
        Area or array of areas of the speciffic conditions to calculate channel depth

    Wb: array
        Bankfull conditions width or array of widths

    Hb: array
        Bankfull conditions depth or array of depths

    Salidas
    -------
    h : array
        Depth or array of depths of the channel flow section for the given area and Entradas
    """
    h=((1.5*(Hb**0.5)*(1./Wb)*A)**(1./2.))**(1./3.)
    return h

#==============================================================================
# Radio hidraulico seccion parabolica en funcion del area y parametros de banca (Wb - Hb)
#==============================================================================

def RhofA(A, w, Wb, Hb):
    """
    Channel hydraulic radius for a parabolic cross-section given the flow area and width and the geometry of bankfull conditions (width and depth).

    Entradas
    ----------
    A : array
        Area or array of areas of the speciffic conditions to calculate channel width

    w :
        Width or array of widths of channel surface of the speciffic conditions

    Wb: array
        Bankfull conditions width or array of widths

    Hb: array
        Bankfull conditions depth or array of depths

    Salidas
    -------
    Rh : array
        hydraulic radius or array of hydraulic radius of channel surface for the speciffic conditions
    """
    Rh=A/((w/2.)*(((1+(4*Hb*(Wb**(-2))*w)**2)**0.5)+(1/(4*Hb*(Wb**(-2))*w))*(np.log((4*Hb*(Wb**(-2))*w)+((1+(4*Hb*(Wb**(-2))*w)**2)**0.5)))))
    return Rh

#==============================================================================
# Corrección celdas tipo embalse
#==============================================================================
def CorrReservoir(flowvec):
    """
    Corrección en las celdas de la red de drenaje que deben ser embalse

    Entradas
    ----------
    flowvec : DataFrame
        DataFrame con atributos de la cuenca

    Salidas
    -------
    flowvec : DataFrame
        Cambios en columna 'Tipo'
    """
    
    #Retiene las celdas de pie de presa
    pie_presa = flowvec['Tipo']==4
    
    Vini = flowvec['Tipo'].astype(float) #tipo de celdas -> embalse
    #las q deberian ser embalse por que les drena embalse
    Vend = copy.copy(Vini)
    uniq = 2
    while uniq==2:
        ind_embalse  = np.where(Vend==3.1)[0]
        dren_embalse = flowvec.loc[ind_embalse,'IDD'].astype(int)
        uniq         = len(np.unique(Vend[dren_embalse]))
        Vend[np.unique(dren_embalse)] = 3.1
    flowvec['Tipo']    = Vend
    
    #Remarca las celdas de pie de presa con el valor 4
    flowvec.loc[pie_presa,'Tipo'] = 4
    
    return flowvec

#==============================================================================
# Ultima celda de embalse se drena a ella misma
#==============================================================================
#import matplotlib.pyplot as plt
def CellReservoir(flowvec,resmask_rst):
    """
    Correccón última celda del embalse, se drena a sí misma

    Entradas
    ----------
    flowvec : DataFrame
        DataFrame con atributos de la cuenca
    resmask_rst: RASTER
        Raster con máscara de cada embalse

    Salidas
    -------
    flowvec : DataFrame
        Cambios en columna 'Tipo'
    """
    id_res = flowvec.index[flowvec['Tipo']>=3.0] # identificar celdas tipo embalse
    coord  = np.array(list(zip(flowvec.loc[id_res,'X'],flowvec.loc[id_res,'Y']))) # coordenadas XY
    # min_dis=5*np.max(flowvec[:,3]) #definición mínima distancia
    i,j    = XYtoij(coord[:,0],coord[:,1],resmask_rst)
    fcluster = resmask_rst['mtrx'][i,j] #cluster.hierarchy.fclusterdata(coord,min_dis,criterion='distance')
    print(fcluster )
    for k in set(fcluster): # encontrar ultima celda del embalse
        pos_clu    = np.where(fcluster==k)[0]
        pos_clu_fv = id_res[pos_clu]
        maxres     = np.max(flowvec.loc[pos_clu_fv,'Area_Acum'])
        ind        = flowvec.index[(flowvec['Area_Acum']==maxres)&(flowvec['Tipo']==3.1)]
        flowvec.loc[ind,'IDD']  = flowvec.loc[ind,'ID'] #Hace la que la celda de pie de presa se drene a sí misma
        flowvec.loc[ind,'Tipo'] = 4              #Marca la celda de pie de presa como tipo de celda 4: celda de acumulación
    return flowvec


def ReadCuenca (path_basin):
    with open(path_basin,'r') as f:
        datos = f.readlines()

    NCells=int(datos[1])
    AreaCells=float(datos[4])
    mtrx=[]
    
    
    #for i in datos[7:]:
    
    #print i; mtrx.append([float(K) for K in i.split()])
    for l in datos[7:]:
        #print 'oe'
        mtrx.append(np.array([float(i) for i in l.split()]))
        
    #mtrx = pd.read_csv(path_basin, sep=' ', skiprows=7).values
    mtrx = np.array(mtrx)
    
    #print mtrx.shape
    #print 'fin'
    f.close
    return NCells,AreaCells,mtrx

def WriteCuenca (new_path_basin,NCells,AreaCells,mtrx):
    header='[NÚMERO DE CELDAS]\n'+str(NCells)+'\n\n[ÁREA DE LAS CELDAS]\n'+str(AreaCells)+'\n\n[MATRIZ DE VARIABLES] %tip_cel %id_des %X %Y %etp %Hui %Hu1 %Hu3 %ks %kp %perd %lon_dre %pen_dre %n_mann %alfa %beta %S0 %S1 %S2 %S3 %S4 %S5 %are_acu %por_are %por_lim %por_arc %K_usle %C_usle %P_usle %anc_ban %pro_ban %Eb %Embalse %Q_med'
    file=open(new_path_basin,'w')
    file.write(header)
    for row in mtrx:
      line='\n'+' '.join(row.astype('str'))
      file.write(line)
    file.close()
    return ('Write new basin')

def escorrentia_to_basin(path_basin, path_escorrentia, path_basin_sed):
    # leer cuenca base
    NCells,AreaCells,mtrx=ReadCuenca (path_basin)
    print(mtrx.shape, NCells)
    mtrx_new=copy.copy(mtrx)
    mtrx_new=np.c_[mtrx_new,np.zeros(mtrx_new.shape[0])]
    Q = SampleRSTtoTBL(mtrx,path_escorrentia,xycols=[2, 3], types=[0.0, 1.0], coltypes=0, colind='', return_column=True)
    mtrx_new[:,33]=Q
    WriteCuenca (path_basin_sed,NCells,AreaCells,mtrx_new)
    return ('Write Basin Sed')

#==============================================================================
# Escribir ráster con una variable de cuenca
#==============================================================================
def Basin2Raster(FlowVectorA,columna,ruta_escritura = '', return_raster = True, prcsn = gdal.GDT_Float32, solo_red = False):
    """
    Crea un ráster con una variable del objeto cuenca
    
    Entradas
    ----------
    FlowVectorA: Objeto cuenca
        tabla de cuenca con el campo a escribir en la columna "column"
        FlowVectorA[:,5] : X
        FlowVectorA[:,6] : Y
        FlowVectorA[:,column] : Valores
    
    columna: integer
        índice de la columna con el valor que se quiere volver ráster
        
    ruta_escritura: string
        Opcional. Ruta donde se quiere escribir el 
        
    prcsn: gdal-type
        Opcional. Por defecto gdal.GDT_Float32. Precisión con la que se quieren
        guardar los valores del ráster
        
    solo_red: Bool
        True si se quiere que el ráster sólo tenda valores en la red de drenaje
    
    Salidas
    -------
    raster: raster
        Diccionario raster con la matriz georreferenciada de los valores de la
        tabla de cuenca en la columna column
    """
    # coordenadas x y y en el objeto de cuenca
    coord_x = np.unique(FlowVectorA['X'])
    coord_y = np.unique(FlowVectorA['Y'])
    # propiedades raster
    clsz = np.abs(np.mean(coord_x[1:] - coord_x[:-1]))
    xll = np.min(coord_x) - (2.5* clsz)
    yll = np.min(coord_y) - (2.5* clsz)
    nodt = -9999.
    rows = coord_y.size + 4
    columns = coord_x.size + 4
    mtrx = np.zeros((rows,columns)) - 9999.
    # raster vacío con 2 filas más que el mínimo de la cobertura
    raster = BuildRaster(xll,yll,clsz,nodt,mtrx)
    # indices de matriz raster para el objeto cuenca
    ies,jotas = XYtoij(FlowVectorA['X'].values,FlowVectorA['Y'].values,raster)
    # llenar raster
    raster['mtrx'][ies,jotas] = FlowVectorA[columna].values
    if solo_red:
        celdas_ladera = FlowVectorA.index[FlowVectorA['Tipo'] == 0.0]
        ies,jotas = XYtoij(FlowVectorA.loc[celdas_ladera,'X'].values,FlowVectorA.loc[celdas_ladera,'Y'].values,raster)
        raster['mtrx'][ies,jotas] = FlowVectorA[columna] = nodt
    if ruta_escritura != '':
        WriteRaster(ruta_escritura,raster,prcsn)
    if return_raster:
        return raster
    else:
        return


#==============================================================================
# Calcular caudal medio
#==============================================================================
def Qmedio(area):
    """
    Calcula el caudal medio de largo plazo en Base a la relación potencial
    encontrada por Rendón-Álvarez (2019)
    
    Entradas
    ----------
    
    area : array_like
        Arreglo unidimensional que contiene las áreas acumuladas [km2] para las
        que se va a calcular el caudal.
    
    Salidas
    -------
    
    Qm : array_like
        Arreglo unidimensional con los caudales [m3/s] correspondientes a cada
        area acumulada.
        
    """
    Qm = 1.533*area**0.7418 # primer momento estadístico de caudal (UPME-UNALMED): Region Colombia
    return Qm

def ReadFromRSTIJ(basin,val_rst):
    """Rutina para leer un raster y arrojar un vector con las celdas pertenecientes a la cuenca
    Entradas
    ---------
    
    basin: Objeto cuenca
    
    ruta_rst: Ruta del raster a cargar
    
      
    """
    if type(val_rst) == str:
        raster=ReadRaster(val_rst)
    else:
        raster = val_rst

    num_celdas = basin.shape[0]
    variable=np.zeros(num_celdas, dtype=float)
    for s in range(0,num_celdas):
        variable[s]=raster['mtrx'][int(basin.loc[s,'i']),int(basin.loc[s,'j'])]
    
    return variable
