
import os, sys
PathProject = os.path.split(os.getcwd())[0]
sys.path.append(os.path.join(PathProject,'functions'))
import spotpy
import numpy as np
import pandas as pd
import json
from .ficheros import *
from .Raster import *
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
import scipy
import rasterio
import numpy as np
from pysheds.grid import Grid


def get_latest_folder_Piori(directory, nuevo_nombre=None):
    """
    Esta función encuentra el último directorio modificado y opcionalmente cambia su nombre.

    Parameters:
        directory (str): Ruta donde se almacenan las simulaciones.
        nuevo_nombre (str, opcional): Nuevo nombre para la carpeta más reciente. Si no se proporciona, no se renombra.

    Returns:
        str: Nombre de la carpeta más reciente (renombrada si se especificó).
    """

    # Obtiene una lista de todas las subcarpetas en el directorio especificado
    folders = [f for f in Path(directory).iterdir() if f.is_dir()]

    if not folders:
        print("No se encontraron carpetas en el directorio.")
        return None

    # Encuentra la carpeta más reciente basada en la fecha de modificación
    latest_folder = max(folders, key=lambda f: f.stat().st_mtime)

    # Si se proporciona un nuevo nombre, renombra la carpeta
    if nuevo_nombre:
        nuevo_path = latest_folder.parent / nuevo_nombre

        # Verificar si ya existe una carpeta con el nuevo nombre
        if nuevo_path.exists():
            print(f"Error: Ya existe una carpeta con el nombre '{nuevo_nombre}'.")
            return latest_folder.name

        # Renombrar la carpeta
        latest_folder.rename(nuevo_path)
        print(f"La carpeta '{latest_folder.name}' ha sido renombrada a '{nuevo_nombre}'.")
        return nuevo_nombre

    # Devuelve solo el nombre de la carpeta si no se renombra
    return latest_folder.name

def ExeSIGAModel_Priori(ProjectPath, NameSIGAFolder, NameExe, PathExeSIGA,escenario):
    from subprocess import Popen
    """
    Esta función ejecuta el modelo SIGA utilizando Cygwin.

    Parameters:
        ProjectPath:
                Ruta del proyecto definido por el usuario.
        NameSIGAFolder:
                Nombre de la carpeta que corresponde a una configuración de SIGA-CAL.
        NameExe:
                Nombre del escenario de simulación.
        PathExeSIGA:
                Ruta del ejecutable de SIGA.
    """
    salida = os.path.join(ProjectPath, NameSIGAFolder, 'salidas', escenario)
    os.makedirs(salida, exist_ok=True)
    # Construir las rutas necesarias
    siga_exe = rf'{PathExeSIGA}'.replace('C:', 'C').replace(os.sep, '/')

    # Obtener el directorio actual y el ejecutable SIGA
    current_dir = os.path.dirname(__file__)

    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    dir_siga_exe = os.path.join(parent_dir, siga_exe)

    # Convertir la ruta de Windows a la ruta de Cygwin
    dir_siga_exe_cygwin = dir_siga_exe.replace('C:', '/cygdrive/c').replace('\\', '/')

    # Construir el comando de ejecución
    path_eje = rf'{ProjectPath}/{NameSIGAFolder}'.replace(os.sep, '/')
    comando = f'. /etc/profile; {dir_siga_exe_cygwin} -t {path_eje} -n {NameExe} -a true'

    # Ejecutar el modelo usando Popen
    print('IMPORTANTE:', comando)
    p = Popen(['C:/cygwin64/bin/bash.exe', '-c', comando])
    # Esperar a que el proceso termine
    stdout, stderr = p.communicate()

    if p.returncode == 0:
        print("El modelo SIGA se ejecutó correctamente.")
    else:
        print(f"Error en la ejecución de SIGA. Código de salida: {p.returncode}")
        print(f"Salida del proceso: {stdout}")
        print(f"Errores del proceso: {stderr}")

    get_latest_folder_Piori(salida, nuevo_nombre='Simulacion_'+escenario)

    return p.returncode

def calcular_promedio_rasters(rasters):
    """Calcula el promedio de una lista de archivos raster."""
    if not rasters:
        return None

    with rasterio.open(rasters[0]) as src:
        perfil = src.profile
        datos = src.read(1).astype(float)
        nodata = perfil.get('nodata', None)  # Obtener el valor de nodata si existe

        # Reemplazar valores nodata por NaN
        if nodata is not None:
            datos[datos == nodata] = np.nan

    for raster in rasters[1:]:
        with rasterio.open(raster) as src:
            nodata = perfil.get('nodata', None)  # Obtener el valor de nodata si existe

            # Reemplazar valores nodata por NaN
            if nodata is not None:
                datos[datos == nodata] = np.nan

            datos += src.read(1).astype(float)

    promedio = datos / len(rasters)
    return promedio, perfil

def muestrearMapas(rutaMapas, listaMapas, cuenca, escenario):
    """ Muestrea los mapas de interés en formato sgabr

    Entradas:
    -----------
    - rutaMapas: string
        string con la ruta a la carpeta mapas

    - listaMapas: list
        Lista de strings cuyos elementos son los nombres de los mapas que se desean muestrear

    Retorna:
    --------
    - rasterDict: diccionario
        Diccionario cuyas elemento son los mapas muestreados


    """
    for mapa in listaMapas:
        rutaCompleta = rutaMapas + mapa + '/' + mapa + '_media.sgabr'
        rasterTemp = Raster(rutaCompleta)
        cuenca[mapa + f'_{escenario}'] = rasterTemp.getValueByCoordinates(cuenca.X.values, cuenca.Y.values)

    return cuenca


def celdasCauceAguasArriba(cuenca):
    '''Devuelve una lista de las celdas de cauce a las que solo les drenan celdas de cauce'''

    idDrenan, conteo = np.unique(cuenca['destino'], return_counts=True)
    celdasDrenaUno = set(idDrenan[conteo == 1])
    idCauce = set(cuenca[cuenca['tipo'] == 1]['idCelda'])
    celdasCauceDrenaUno = list(celdasDrenaUno.intersection(idCauce))

    return (celdasCauceDrenaUno)


def rastrearEntradas(variable, cuenca, tipoEntradas=''):
    """Rastrea las entradas a una celda de acuerdo con la información de la celda de destino
    (destino) en el objeto cuenca

    Entradas:
    -----------
    - variable: str
        string con el nombre de la variable de interés, a la cual se quiere calcular las entradas
        a las celdas de destino.

    - cuenca: pandas.DataFrame
        pandas.DataFrame cuyos índices son el id de la celda correspondiente
        del objeto cuenca y contiene las siguientes columnas:
            - 'destino': id de la celda de destino

    Retorna:
    --------
    - cuenca: pandas.DataFrame
        pandas.Dataframe de entrada con la columna adicional de las entradas a cada celda de la
        variable de interés.

    """
    # en qué tipo de celdas se quieren acumular entradas
    if tipoEntradas == 'cauce':
        cuencasuma = cuenca[cuenca['tipo'] == 1].groupby(by=['destino']).sum()
    elif tipoEntradas == 'ladera':
        cuencasuma = cuenca[cuenca['tipo'] == 0].groupby(by=['destino']).sum()
    else:
        cuencasuma = cuenca.groupby(by=['destino']).sum()

    # guardar entradas en la cuenca
    nombreColumna = f'{variable}_ent{tipoEntradas}'
    cuenca.loc[:, nombreColumna] = 0.0
    cuenca.loc[cuencasuma.index, nombreColumna] = cuencasuma[variable]
    return cuenca

def sumarFracciones(cuenca, listaVariables, listaFracciones, nombreEscenario, total='', primeroFraccion = False, borrarOriginal=True):
    """
    Suma fracciones de una misma variable cuyo nombre está compuesto por
    {variable}{fraccion}_{nombreEscenario} si primeroFraccion es False o por
    {fraccion}{variable}_{nombreEscenario} si primeroFraccion es True. Guarda el resultado en la columna
    {nombreTotal}_{nombreEscenario} del objeto cuenca, y borra las variables iniciales si borrarOriginal=True
    """
    if primeroFraccion:
        for variable in listaVariables:
            lista = [f'{fraccion}{variable}_{nombreEscenario}' for fraccion in listaFracciones]
            cuenca[f'{variable}{total}_{nombreEscenario}'] = cuenca[lista].sum(axis=1)
            if borrarOriginal:
                cuenca = cuenca.drop(lista,axis=1)
    else:
        for variable in listaVariables:
            lista = [f'{variable}{fraccion}_{nombreEscenario}' for fraccion in listaFracciones]
            cuenca[f'{variable}{total}_{nombreEscenario}'] = cuenca[lista].sum(axis=1)
            if borrarOriginal:
                cuenca = cuenca.drop(lista,axis=1)
    return cuenca


def calcularRC(entradas, salidas):
    """
    Función para calcular el coeficiente de retención de la variable cuyas entradas y salidas en
    cada pixel se proveen.

    entradas:
    ----------

    - entradas: pandas.Series
        pandas.Series cuyos índices son el id de la celda correspondiente
        del objeto cuenca y sus valores son las entradas de la variable.

    - salidas: pandas.Series
        pandas.Series cuyos índices son el id de la celda correspondiente
        del objeto cuenca y sus valores son las salidas de la variable.

    Retorna:
    --------

    - RC: pandas.Series
        pandas.Series cuyos índices son el id de la celda correspondiente
        del objeto cuenca y sus valores son el coeficiente de retención de la variable.
    """
    RC = (entradas - salidas) / entradas  # operaciones vectoriales

    return RC


def calcularDk(cuenca, cabeceras, variable, escenarioBase, escenarioRestauracion):
    """
    Calcula la entrega (Dk) de la variable

    Entradas:
    -----------
    - cuenca: pandas.DataFrame
        pandas.DataFrame cuyos índices son el id de la celda correspondiente
        debe tener las siguientes columnas:
            - I4_lineaBase : Variable I4 para el escenario lineaBase
            - I4_restauración : Variable I4 para el escenario restauración
            - RC_flujoBase : columna del coeficiente de retención

    - cabeceras: lista
        lista con las celdas de cabecera

    - tipo: integer
        1: cálculo para flujo base
        2: cálculo para sedimentos
        3: cálculo para calidad

    - variable: str
        nombre de la variable para la cual se está calculando Dk

    Retorna:
    --------

    - cuenca: pandas.DataFrame
        pandas.DataFrame de entrada con la variable Dk_{variable} adicionada como columna

    """
    # Calcular la trayectoria de flujo
    for celda in cabeceras:

        celdaEval = celda
        celdaDest = cuenca.loc[celdaEval, 'destino']
        listaCamino = [celdaEval]
        cauces = 0

        # Se itera hasta la primera celda de cauce, es decir, la segunda celda cuyo destino sea un cauce
        while cauces < 2:
            listaCamino.append(celdaDest)
            celdaEval = celdaDest
            celdaDest = cuenca.loc[celdaEval, 'destino']
            tipoDest = cuenca.loc[celdaEval, 'tipo']
            embalse = cuenca.loc[celdaEval, 'embalse']
            if (tipoDest == 1.0) or (embalse == 1.0):
                cauces += 1

        # Filtrar cuenca en el camino de la celda de cabecera evaluada invertir el orden
        cuencaFilt = cuenca[cuenca['idCelda'].isin(listaCamino)][::-1]

        # calcular la productoria
        cuencaFilt['prod'] = cuencaFilt['1-RC'].cumprod()

        # Cálculo de dk para red de drenaje de celda de cabecera evaluada
        if variable == 'flujoBase':
            cuencaFilt[f'Dk_{variable}'] = -1 * cuencaFilt['prod'] * (
                        cuencaFilt[f'I4_{escenarioBase}'] - cuencaFilt[f'I4_{escenarioRestauracion}'])

        if variable == 'sedimentos':
            cuencaFilt['prod'].iloc[:-1] = cuencaFilt['prod'].iloc[1:]
            cuencaFilt[f'Dk_{variable}'] = cuencaFilt['prod'] * (
                        cuencaFilt[f'E2E_{escenarioBase}'] - cuencaFilt[f'E2E_{escenarioRestauracion}'])

        if (variable == 'N' or variable == 'P' or variable == 'EC'):
            cuencaFilt[f'Dk_{variable}'] = cuencaFilt['prod'] * (
                        cuencaFilt[f'Q{variable}_{escenarioBase}'] - cuencaFilt[f'Q{variable}_{escenarioRestauracion}'])

        cuencaFilt = cuencaFilt[['idCelda', 'prod', f'Dk_{variable}']]

        # Merge del archivo cuenca original con el de la red filtrada
        cuenca.loc[cuencaFilt.index, 'prod'] = cuencaFilt['prod']
        cuenca.loc[cuencaFilt.index, f'Dk_{variable}'] = cuencaFilt[f'Dk_{variable}']

    # cuenca[f'Dk_{variable}'].fillna(value=0.0, inplace=True)
    cuenca = cuenca.rename(columns={"prod": f"{variable}Prod"})
    cuenca = cuenca.drop(['1-RC'], axis=1)
    return cuenca

def calcularRCLaderaCauce(entradasLadera, entradasCauce, salidasLadera, salidasCauce, variable, cuenca):
    """

    """

    # Separar la cuenca en ladera y cauce
    cuencaCauce = cuenca[cuenca['tipo'] == 1.0]
    cuencaLadera = cuenca[cuenca['tipo'] == 0.0]

    # Calcular RC para ladera y cauce
    cuencaCauce[f'RC_{variable}'] = calcularRC(entradasCauce, salidasCauce)
    cuencaLadera[f'RC_{variable}'] = calcularRC(entradasLadera, salidasLadera)

    # Donde las entradas son cero, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuencaCauce.loc[(entradasCauce <= umbral).values, [f'RC_{variable}']] = 0.0
    cuencaLadera.loc[(entradasLadera <= umbral).values, [f'RC_{variable}']] = 0.0

    # filtrar donde las entradas son iguales a 0.0 que genera RC inf
    cuencaCauce.loc[(entradasCauce == 0.0).values, [f'RC_{variable}']] = 0.0
    cuencaLadera.loc[(entradasLadera == 0.0).values, [f'RC_{variable}']] = 0.0

    # Ensamblar ladera y cauce
    cuencaSedim = pd.concat([cuencaCauce, cuencaLadera])
    cuencaSedim = cuencaSedim[['idCelda', f'RC_{variable}']]

    # Ensamblar resultados a cuenca original
    cuenca = pd.merge(left=cuenca, right=cuencaSedim, left_on='idCelda', right_on='idCelda', how='left')

    # cálculo 1-RC
    cuenca['1-RC'] = 1 - cuenca[f'RC_{variable}']

    return cuenca


def muestreo_escenarios(rutaCuenca,rutaMapasBase,rutaMapasRestauracion,name_scenario,new_escenario):
    # ===================================================== Leer cuenca =====================================================
    cuenca =  pd.read_csv(rutaCuenca, sep=" ", skiprows=10)
    cuenca = cuenca[['tipo', 'destino', 'embalse', 'X', 'Y']]
    cuenca[['tipo', 'destino', 'embalse']] = cuenca[['tipo', 'destino', 'embalse']].astype(int)
    cuenca.index =  cuenca.index+1
    cuenca['idCelda'] = np.array(cuenca.index)


    # ===================================================== Muestreo escenario lineaBase======================================================================================
    # Muestrear mapas del escenario base
    # listaMapas=['I4','E4']

    listaMapas=['I4','E4','O4','arcE2S','limE2S','areE2S','arcE2D','limE2D','areE2D','arcE2E','limE2E','areE2E','arcE5S','limE5S','areE5S','arcE5D','limE5D','areE5D',
                'arcE5E','limE5E','areE5E','EC','QEC','E2EC','E3EC','I4EC','NO','NH4','NO3','QNO','QNH4','QNO3','E2NO','E2NH4','E2NO3','E3NO','E3NH4','E3NO3','I4NO3','I4NH4',
                'I4NO','PI','PO','QPI','QPO','E2PI','E2PO','E3PI','E3PO','I4PI','I4PO']
    cuenca = muestrearMapas(rutaMapasBase,listaMapas,cuenca,name_scenario)



    # Sumar fracciones de sedimentos
    listaVariables = ['E2S','E2D','E2E','E5S','E5D','E5E']
    listaFracciones = ['arc','lim','are']
    cuenca = sumarFracciones(cuenca, listaVariables, listaFracciones, name_scenario, primeroFraccion = True)

    # Sumar variables calidad de agua
    listaVariables = ['Q','E2','E3','I4']
    listaTotales = ['N','P']
    listaFracciones = {'N':['NO','NO3','NH4'], 'P':['PI','PO']}
    for total in listaTotales:
        cuenca = sumarFracciones(cuenca, listaVariables, listaFracciones[total], name_scenario,total=total)


    cuenca[f'N_{name_scenario}'] = cuenca[[f'NH4_{name_scenario}',f'NO_{name_scenario}',f'NO3_{name_scenario}']].sum(axis=1)
    cuenca[f'P_{name_scenario}'] = cuenca[[f'NH4_{name_scenario}',f'NO_{name_scenario}',f'NO3_{name_scenario}']].sum(axis=1)



    #========================================================= Muestreo escenario restauración ==========================================================================
    listaMapas=['I4','arcE2E','limE2E','areE2E','QEC','QNO','QNH4','QNO3','QPI','QPO']
    cuenca = muestrearMapas(rutaMapasRestauracion,listaMapas,cuenca,new_escenario)

    # sumar fracciones de sedimentos
    listaVariables = ['E2E']
    listaFracciones = ['arc','lim','are']
    cuenca = sumarFracciones(cuenca, listaVariables, listaFracciones, new_escenario, primeroFraccion = True)

    # Sumar variables calidad de agua
    listaVariables = ['Q']
    listaTotales = ['N','P']
    listaFracciones = {'N':['NO','NO3','NH4'],'P':['PI','PO']}
    for total in listaTotales:
        cuenca = sumarFracciones(cuenca, listaVariables, listaFracciones[total], new_escenario,total=total)

    return cuenca

def proceso_flujo_base(cuenca,name_scenario,new_escenario):
    #=================================================================== flujo base =========================================================================================
    sustancia = 'flujoBase'
    # rastrear entradas necesarias
    cuenca = rastrearEntradas(f'E4_{name_scenario}', cuenca)

    # calcular RC para flujo base
    entradasFlujoBase = cuenca[f'I4_{name_scenario}'] + cuenca[f'E4_{name_scenario}_ent']
    salidasFlujoBase = cuenca[f'E4_{name_scenario}'] + cuenca[f'O4_{name_scenario}']
    cuenca[f'RC_{sustancia}'] = calcularRC(entradasFlujoBase, salidasFlujoBase)

    # Donde las entradas son cero, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasFlujoBase<=umbral).values,[f'RC_{sustancia}']] = 0.0

    #Llenar pixeles urbanos con 1.0
    entradasSalidasNulas = entradasFlujoBase[salidasFlujoBase==0.0]
    cuenca.loc[entradasSalidasNulas.index,f'RC_{sustancia}'] = 1.0

    # Calculo Dk flujo base
    # calcular 1-RC
    cuenca['1-RC'] = 1 - cuenca[f'RC_{sustancia}']
    cuenca[f'Dk_{sustancia}'] = np.zeros(len(cuenca))*np.nan

    # calcular celdas de cabeceras
    celdasAll = set(cuenca['idCelda'])
    celdasDest = set(cuenca['destino'])
    cabeceras = list(celdasAll.difference(celdasDest))
    celdasCauce = celdasCauceAguasArriba(cuenca)


    cabeceras = cabeceras + celdasCauce
    cuenca = calcularDk(cuenca,cabeceras,sustancia,name_scenario,new_escenario)


    # enmascarar Dk para celdas de embalse
    cuenca[f'Dk_{sustancia}'][cuenca['embalse']==1]=0.0

    return cuenca

def proceso_Sedimentos(cuenca,name_scenario,new_escenario):

    # ================================================ SEDIMENTOS =============================================================================
    #Crear entradas para sedimentos
    cuenca = rastrearEntradas(f'E2S_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E2D_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E2E_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E2S_{name_scenario}', cuenca, tipoEntradas= 'ladera')
    cuenca = rastrearEntradas(f'E2D_{name_scenario}', cuenca, tipoEntradas= 'ladera')
    cuenca = rastrearEntradas(f'E2E_{name_scenario}', cuenca, tipoEntradas= 'ladera')
    cuenca = rastrearEntradas(f'E5S_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E5D_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E5S_{name_scenario}', cuenca, tipoEntradas= 'cauce')
    cuenca = rastrearEntradas(f'E5D_{name_scenario}', cuenca, tipoEntradas= 'cauce')
    cuenca = rastrearEntradas(f'E5E_{name_scenario}', cuenca, tipoEntradas= 'cauce')


    # calcular entradas para sedimentos
    sustancia = 'sedimentos'
    # Lo que reciben las celdas de la entrega de celdas tipo ladera
    listaEntradasSedimentos = [f'{variable}_{name_scenario}_entladera' for variable in ['E2S','E2D','E2E']]
    entradasSedimentos = cuenca[listaEntradasSedimentos].sum(axis=1)

    # se añaden las entradas de cauce a cauce
    listaEntradasSedimentos = [f'{variable}_{name_scenario}_entcauce' for variable in ['E5S','E5D','E5E']]
    entradasSedimentos[cuenca['tipo']==1.0] = entradasSedimentos[cuenca['tipo']==1.0] + cuenca[cuenca['tipo']==1.0][listaEntradasSedimentos].sum(axis=1)

    # se añaden las entradas de la erosión en la misma celda sólo en cauce
    entradasSedimentos[cuenca['tipo']==1.0] = entradasSedimentos[cuenca['tipo']==1.0] + cuenca[cuenca['tipo']==1.0][f'E2E_{name_scenario}']

    # calcular salidas sedimentos
    # calcular salidas para sedimentos en celdas tipo ladera
    listaSalidas = [f'{variable}_{name_scenario}' for variable in ['E2S','E2D']]
    salidasSedimentos = cuenca[listaSalidas].sum(axis=1)

    # calcular salidas para sedimentos en celdas tipo cauce
    listaSalidas = [f'{variable}_{name_scenario}' for variable in ['E5S','E5D']]
    salidasSedimentos[cuenca['tipo']==1.0] = cuenca[cuenca['tipo']==1.0][listaSalidas].sum(axis=1)

    # Calcular RC sedimentos
    cuenca[f'RC_{sustancia}'] = calcularRC(entradasSedimentos, salidasSedimentos)

    # Donde las entradas son cero, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasSedimentos<=umbral).values,[f'RC_{sustancia}']] = 0.0

    # Donde las entradas y las salidas son casi iguales, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasSedimentos-salidasSedimentos<=umbral).values,[f'RC_{sustancia}']] = 0.0

    # calcular celdas de cabeceras
    celdasAll = set(cuenca['idCelda'])
    celdasDest = set(cuenca['destino'])
    cabeceras = list(celdasAll.difference(celdasDest))
    celdasCauce = celdasCauceAguasArriba(cuenca)

    cabeceras = cabeceras + celdasCauce

    # ================================================= Calculo Dk sedimentos ============================================================================
    # Calculo Dk sedimentos
    # calcular 1-RC
    cuenca['1-RC'] = 1 - cuenca[f'RC_{sustancia}']
    cuenca = calcularDk(cuenca,cabeceras,sustancia,name_scenario,new_escenario)

    # enmascarar Dk para celdas de embalse
    cuenca[f'Dk_{sustancia}'][cuenca['embalse']==1]=0.0

    return cuenca

def proceso_Ecoli(cuenca,name_scenario,new_escenario):

    # ================================================ ECOLI =============================================================================
    #Crear entradas para EC
    cuenca = rastrearEntradas(f'E2EC_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E3EC_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'QEC_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'EC_{name_scenario}', cuenca,'cauce')

    # calcular entradas para EC
    sustancia = 'EC'
    # Lo que reciben las celdas de la entrega de celdas tipo ladera
    listaEntradasEC = [f'{variable}_{name_scenario}_ent' for variable in ['E2EC','E3EC','QEC']]
    entradasEC = cuenca[listaEntradasEC].sum(axis=1)

    # se añaden las entradas de cauce a cauce
    listaEntradasEC = [f'{variable}_{name_scenario}_entcauce' for variable in ['EC']]
    entradasEC[cuenca['tipo']==1.0] = entradasEC[cuenca['tipo']==1.0] + cuenca[cuenca['tipo']==1.0][listaEntradasEC].sum(axis=1)


    # calcular salidas sedimentos
    # calcular salidas para EC en celdas tipo ladera
    listaSalidasEC = [f'{variable}_{name_scenario}' for variable in ['I4EC','E2EC','E3EC']]
    salidasEC = cuenca[listaSalidasEC].sum(axis=1)

    # calcular salidas para sedimentos en celdas tipo cauce
    listaSalidasEC = [f'{variable}_{name_scenario}' for variable in ['EC','I4EC']]
    salidasEC[cuenca['tipo']==1.0] = cuenca[cuenca['tipo']==1.0][listaSalidasEC].sum(axis=1)

    # Calcular RC EC
    cuenca[f'RC_{sustancia}'] = calcularRC(entradasEC, salidasEC)

    # Donde las entradas son cero, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasEC<=umbral).values,[f'RC_{sustancia}']] = 0.0

    # # Donde las entradas y las salidas son casi iguales, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasEC-salidasEC<=umbral).values,[f'RC_{sustancia}']] = 0.0

    # calcular celdas de cabeceras
    celdasAll = set(cuenca['idCelda'])
    celdasDest = set(cuenca['destino'])
    cabeceras = list(celdasAll.difference(celdasDest))
    celdasCauce = celdasCauceAguasArriba(cuenca)

    cabeceras = cabeceras + celdasCauce
    # ================================================= Calculo EC  ============================================================================
    # Calculo Dk EC
    # calcular 1-RC
    cuenca['1-RC'] = 1 - cuenca[f'RC_{sustancia}']
    cuenca = calcularDk(cuenca,cabeceras,sustancia,name_scenario,new_escenario)

    # enmascarar Dk para celdas de embalse
    cuenca[f'Dk_{sustancia}'][cuenca['embalse']==1]=0.0

    return cuenca

def proceso_Nitrogeno(cuenca,name_scenario,new_escenario):

    # ============================================ NITRÓGENO =============================================================================
    #Crear entradas para N
    cuenca = rastrearEntradas(f'E2N_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E3N_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'QN_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'N_{name_scenario}', cuenca,'cauce')

    sustancia = 'N'
    # calcular entradas para N
    # Lo que reciben las celdas de la entrega de celdas tipo ladera
    listaEntradasN = [f'{variable}_{name_scenario}_ent' for variable in ['E2N','E3N','QN']]
    entradasN = cuenca[listaEntradasN].sum(axis=1)

    # se añaden las entradas de cauce a cauce
    listaEntradasN = [f'{variable}_{name_scenario}_entcauce' for variable in ['N']]
    entradasN[cuenca['tipo']==1.0] = entradasN[cuenca['tipo']==1.0] + cuenca[cuenca['tipo']==1.0][listaEntradasN].sum(axis=1)


    # calcular salidas sedimentos
    # calcular salidas para EC en celdas tipo ladera
    listaSalidasN = [f'{variable}_{name_scenario}' for variable in ['I4N','E2N','E3N']]
    salidasN = cuenca[listaSalidasN].sum(axis=1)

    # calcular salidas para sedimentos en celdas tipo cauce
    listaSalidasN = [f'{variable}_{name_scenario}' for variable in ['N','I4N']]
    salidasN[cuenca['tipo']==1.0] = cuenca[cuenca['tipo']==1.0][listaSalidasN].sum(axis=1)

    # Calcular RC_N
    cuenca[f'RC_{sustancia}'] = calcularRC(entradasN, salidasN)

    # Donde las entradas son cero, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasN<=umbral).values,[f'RC_{sustancia}']] = 0.0

    # # Donde las entradas y las salidas son casi iguales, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasN-salidasN<=umbral).values,[f'RC_{sustancia}']] = 0.0


    # calcular celdas de cabeceras
    celdasAll = set(cuenca['idCelda'])
    celdasDest = set(cuenca['destino'])
    cabeceras = list(celdasAll.difference(celdasDest))
    celdasCauce = celdasCauceAguasArriba(cuenca)

    cabeceras = cabeceras + celdasCauce
    # ================================================= Calculo EC  ============================================================================
    # Calculo Dk_EC
    # calcular 1-RC
    cuenca['1-RC'] = 1 - cuenca[f'RC_{sustancia}']
    cuenca = calcularDk(cuenca,cabeceras,sustancia,name_scenario,new_escenario)

    # enmascarar Dk para celdas de embalse
    cuenca[f'Dk_{sustancia}'][cuenca['embalse']==1]=0.0

    return cuenca

def proceso_fosforo(cuenca,name_scenario,new_escenario):


    # ============================================ FÓSFORO =============================================================================
    #Crear entradas para P
    cuenca = rastrearEntradas(f'E2P_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'E3P_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'QP_{name_scenario}', cuenca)
    cuenca = rastrearEntradas(f'P_{name_scenario}', cuenca,'cauce')

    sustancia = 'P'
    # calcular entradas para P
    # Lo que reciben las celdas de la entrega de celdas tipo ladera
    listaEntradasP = [f'{variable}_{name_scenario}_ent' for variable in ['E2P','E3P','QP']]
    entradasP = cuenca[listaEntradasP].sum(axis=1)

    # se añaden las entradas de cauce a cauce
    listaEntradasP = [f'{variable}_{name_scenario}_entcauce' for variable in ['P']]
    entradasP[cuenca['tipo']==1.0] = entradasP[cuenca['tipo']==1.0] + cuenca[cuenca['tipo']==1.0][listaEntradasP].sum(axis=1)


    # calcular salidas fósforo
    # calcular salidas para fósforo en celdas tipo ladera
    listaSalidasP = [f'{variable}_{name_scenario}' for variable in ['I4P','E2P','E3P']]
    salidasP = cuenca[listaSalidasP].sum(axis=1)

    # calcular salidas para fósforo en celdas tipo cauce
    listaSalidasP = [f'{variable}_{name_scenario}' for variable in ['P','I4P']]
    salidasP[cuenca['tipo']==1.0] = cuenca[cuenca['tipo']==1.0][listaSalidasP].sum(axis=1)

    # Calcular RC_P
    cuenca[f'RC_{sustancia}'] = calcularRC(entradasP, salidasP)

    # Donde las entradas son cero, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasP<=umbral).values,[f'RC_{sustancia}']] = 0.0

    # # Donde las entradas y las salidas son casi iguales, llenar RC con 0.0 (para estabilidad numérica en RC)
    umbral = 1e-05
    cuenca.loc[(entradasP-salidasP<=umbral).values,[f'RC_{sustancia}']] = 0.0

    # calcular celdas de cabeceras
    celdasAll = set(cuenca['idCelda'])
    celdasDest = set(cuenca['destino'])
    cabeceras = list(celdasAll.difference(celdasDest))
    celdasCauce = celdasCauceAguasArriba(cuenca)

    cabeceras = cabeceras + celdasCauce

    # ================================================= Calculo P ============================================================================
    # Calculo Dk_P
    # calcular 1-RC
    cuenca['1-RC'] = 1 - cuenca[f'RC_{sustancia}']
    cuenca = calcularDk(cuenca,cabeceras,sustancia,name_scenario,new_escenario)

    # enmascarar Dk para celdas de embalse
    cuenca[f'Dk_{sustancia}'][cuenca['embalse'] == 1] = 0.0

    return cuenca


def export_result(cuenca, path_result, Name_Basin='Basin'):
    """
    Exporta resultados de cuenca a Excel y rasters individuales.

    Parámetros:
    ----------
    cuenca : pandas.DataFrame
        DataFrame con todos los indicadores calculados
    path_result : str
        Ruta base donde se guardarán los resultados (temp/05-Prioritization/)
    Name_Basin : str
        Nombre de la cuenca (default: 'Basin')
    """
    cuencaExportar = cuenca[
        ['idCelda', 'tipo', 'X', 'Y', 'RC_flujoBase', 'RC_sedimentos', 'RC_EC', 'RC_N', 'RC_P', 'flujoBaseProd',
         'sedimentosProd', 'ECProd',
         'NProd', 'PProd', 'dI4', 'dE2E', 'dQEC', 'dQN', 'dQP', 'Dk_flujoBase', 'Dk_sedimentos', 'Dk_EC', 'Dk_N',
         'Dk_P']]

    cuencaExportar = cuencaExportar.loc[:, ~cuencaExportar.columns.duplicated()].copy()

    # Crear carpeta 01_Result si no existe
    result_folder = os.path.join(path_result, '01_Result')
    os.makedirs(result_folder, exist_ok=True)

    # Exportar Excel
    excel_path = os.path.join(result_folder, f'cuencaDk_{Name_Basin}.xlsx')
    cuencaExportar.to_excel(excel_path, index=None)
    print(f"  ✓ Excel exportado: {excel_path}")

    # Exportar rasters individuales
    for nombreRaster in cuencaExportar.columns.values.tolist():
        dictRaster = {'x': cuenca['X'].values, 'y': cuenca['Y'].values,
                      'values': cuencaExportar[nombreRaster].values}

        raster = Raster(dictRaster)
        raster_path = os.path.join(result_folder, f'{nombreRaster}.tif')
        raster.exportRaster(raster_path)

    print(f"  ✓ {len(cuencaExportar.columns)} rasters exportados en 01_Result/")


# ========================================================================================================
# PARTE 2: AGREGACIÓN ESPACIAL POR POLÍGONOS (02_Priorizacion_Esquema_part_2.py)
# ========================================================================================================

import geopandas as gpd
from shapely.geometry import Polygon
import copy


def getGeometry(obj, dx, dy):
    """
    Genera una geometría poligonal con forma de celda ráster
    a partir de las coordenadas X e Y del centro de la celda

    Parámetros:
    ----------
    obj : objeto con atributos X e Y
        Objeto que tenga entre sus atributos un X y un Y
        que representen las coordenadas horizontales y verticales del centro de una celda raster
    dx : flotante
        Tamaño horizontal de celda dividido por 2
    dy : flotante
        Tamaño vertical de celda dividido por 2

    Retorna:
    -------
    geom : geometría poligonal de Shapely
        Polígono generado
    """
    xmin = obj.X - dx
    xmax = obj.X + dx
    ymin = obj.Y - dy
    ymax = obj.Y + dy
    coords = [
        [xmin, ymax],
        [xmax, ymax],
        [xmax, ymin],
        [xmin, ymin]
    ]
    geom = Polygon(coords)
    return geom


def transformacionRangoUnitario(x, minimo=None):
    """
    Calcula una transformación de rango unitario para un arreglo

    Parámetros:
    ----------
    x : arreglo numpy
        Arreglo de valores a transformar
    minimo : flotante (opcional)
        Valor sobre el que se desea que la transformación pivote (ejemplo: 0)

    Retorna:
    -------
    y : arreglo numpy
        Arreglo de valores transformados
    """
    maximo = np.nanmax(x)
    if minimo is None:
        minimo = np.nanmin(x)
    rango = maximo - minimo

    y = (x - minimo) / rango

    return y


def procesar_priorizacion_completa(
    rutaCuenca,
    rutaMapasBase,
    rutaMapasRestauracion,
    name_scenario,
    new_escenario,
    Name_Basin,
    path_result,
    rutaPoligonos,
    nombreID,
    rutaCoberturasLineaBase,
    rutaCoberturasRestauracion,
    fcPonderacion,
    percIntervencion,
    umbral=0.5,
    procesar_variables=None,
    lookup_idQ=None,
    rutaWatershed=None
):
    """
    Función principal que ejecuta todo el proceso de priorización completo.
    Integra extracción de resultados y agregación espacial.

    Parámetros:
    ----------
    rutaCuenca : str
        Ruta al archivo de topología (Topology_*.txt)
    rutaMapasBase : str
        Ruta a mapas del escenario base (carpeta mapas/)
    rutaMapasRestauracion : str
        Ruta a mapas del escenario NBS (carpeta mapas/)
    name_scenario : str
        Nombre del escenario base (Step 2)
    new_escenario : str
        Nombre del escenario NBS (Step 3)
    Name_Basin : str
        Nombre de la cuenca
    path_result : str
        Ruta donde se guardarán resultados (temp/05-Prioritization/)
    rutaPoligonos : str
        Ruta al shapefile de polígonos de agregación
    nombreID : str
        Nombre del campo ID en shapefile
    rutaCoberturasLineaBase : str
        Ruta a CLC.tif del escenario base
    rutaCoberturasRestauracion : str
        Ruta a CLC.tif del escenario NBS
    fcPonderacion : dict
        Factores de ponderación {'flujoBase': 1.0, 'sedimentos': 1.0, 'EC': 0.2, 'N': 0.6, 'P': 0.7}
    percIntervencion : list
        Lista de porcentajes de intervención [20, 40, 80]
    umbral : float
        Umbral de área intervenida (default 0.5)
    procesar_variables : dict
        Variables a procesar {'flujo_base': True, 'sedimentos': True, 'ecoli': True, 'nitrogeno': True, 'fosforo': True}
    lookup_idQ : dict (opcional)
        Diccionario de conversión LULC → idQ
    rutaWatershed : str (opcional)
        Ruta al archivo Watershed.shp para obtener el CRS del proyecto.
        Si no se proporciona, se asumirá EPSG:32619

    Retorna:
    -------
    cuenca : pandas.DataFrame
        DataFrame con todos los indicadores calculados
    poligonos : geopandas.GeoDataFrame
        GeoDataFrame con polígonos y prioridades
    agrupacionPoligonalDk : pandas.DataFrame
        DataFrame con agregación por polígonos
    """

    print("="*80)
    print("INICIANDO PROCESO DE PRIORIZACIÓN COMPLETO")
    print("="*80)

    # ========== PASO 1: Muestreo de escenarios ==========
    print("\n[PASO 1/8] Muestreando escenarios...")
    cuenca = muestreo_escenarios(rutaCuenca, rutaMapasBase, rutaMapasRestauracion, name_scenario, new_escenario)
    print(f"  ✓ Cuenca cargada: {len(cuenca)} celdas")

    # ========== PASO 2: Procesamiento de indicadores ==========
    print("\n[PASO 2/8] Procesando indicadores...")

    if procesar_variables is None:
        procesar_variables = {
            'flujo_base': True,
            'sedimentos': True,
            'ecoli': True,
            'nitrogeno': True,
            'fosforo': True
        }

    if procesar_variables.get('flujo_base', False):
        print("  - Procesando flujo base...")
        cuenca = proceso_flujo_base(cuenca, name_scenario, new_escenario)
        print("    ✓ Dk_flujoBase calculado")

    if procesar_variables.get('sedimentos', False):
        print("  - Procesando sedimentos...")
        cuenca = proceso_Sedimentos(cuenca, name_scenario, new_escenario)
        print("    ✓ Dk_sedimentos calculado")

    if procesar_variables.get('ecoli', False):
        print("  - Procesando E.coli...")
        cuenca = proceso_Ecoli(cuenca, name_scenario, new_escenario)
        print("    ✓ Dk_EC calculado")

    if procesar_variables.get('nitrogeno', False):
        print("  - Procesando Nitrógeno...")
        cuenca = proceso_Nitrogeno(cuenca, name_scenario, new_escenario)
        print("    ✓ Dk_N calculado")

    if procesar_variables.get('fosforo', False):
        print("  - Procesando Fósforo...")
        cuenca = proceso_fosforo(cuenca, name_scenario, new_escenario)
        print("    ✓ Dk_P calculado")

    # ========== PASO 3: Calcular deltas ==========
    print("\n[PASO 3/8] Calculando deltas...")
    cuenca['dI4'] = -1.0 * (cuenca[f'I4_{name_scenario}'] - cuenca[f'I4_{new_escenario}'])
    cuenca['dE2E'] = cuenca[f'E2E_{name_scenario}'] - cuenca[f'E2E_{new_escenario}']
    cuenca['dQEC'] = cuenca[f'QEC_{name_scenario}'] - cuenca[f'QEC_{new_escenario}']
    cuenca['dQN'] = cuenca[f'QN_{name_scenario}'] - cuenca[f'QN_{new_escenario}']
    cuenca['dQP'] = cuenca[f'QP_{name_scenario}'] - cuenca[f'QP_{new_escenario}']
    print("  ✓ Deltas calculados: dI4, dE2E, dQEC, dQN, dQP")

    # ========== PASO 4: Exportar resultados de cuenca ==========
    print("\n[PASO 4/8] Exportando resultados de cuenca...")
    os.makedirs(os.path.join(path_result, '01_Inputs'), exist_ok=True)

    cuencaExportar = cuenca[
        ['idCelda', 'tipo', 'X', 'Y', 'RC_flujoBase', 'RC_sedimentos', 'RC_EC', 'RC_N', 'RC_P',
         'flujoBaseProd', 'sedimentosProd', 'ECProd', 'NProd', 'PProd',
         'dI4', 'dE2E', 'dQEC', 'dQN', 'dQP',
         'Dk_flujoBase', 'Dk_sedimentos', 'Dk_EC', 'Dk_N', 'Dk_P']
    ]

    cuencaExportar = cuencaExportar.loc[:, ~cuencaExportar.columns.duplicated()].copy()
    excel_path = os.path.join(path_result, '01_Inputs', f'cuencaDk_{Name_Basin}.xlsx')
    cuencaExportar.to_excel(excel_path, index=None)
    print(f"  ✓ Excel exportado: {excel_path}")

    # Exportar rasters individuales
    os.makedirs(os.path.join(path_result, '02_Rasters'), exist_ok=True)
    for nombreRaster in cuencaExportar.columns.values.tolist():
        dictRaster = {'x': cuenca['X'].values, 'y': cuenca['Y'].values,
                      'values': cuencaExportar[nombreRaster].values}
        raster = Raster(dictRaster)
        raster.exportRaster(os.path.join(path_result, '02_Rasters', f'{nombreRaster}.tif'))
    print(f"  ✓ {len(cuencaExportar.columns)} rasters exportados")

    # ========== PASO 5: Agregación espacial por polígonos ==========
    print("\n[PASO 5/8] Procesando agregación espacial...")

    # Determinar CRS del proyecto
    if rutaWatershed and os.path.exists(rutaWatershed):
        watershed = gpd.read_file(rutaWatershed)
        project_crs = watershed.crs
        if project_crs is None:
            print("  ⚠ WARNING: Watershed.shp no tiene CRS. Se usará EPSG:32619")
            project_crs = 'EPSG:32619'
    else:
        print("  ⚠ WARNING: No se proporcionó rutaWatershed. Se usará EPSG:32619")
        project_crs = 'EPSG:32619'

    print(f"  ✓ CRS del proyecto: {project_crs}")

    # Leer shapefile de polígonos
    poligonos = gpd.read_file(rutaPoligonos)

    # Verificar y convertir CRS
    if poligonos.crs is None:
        raise ValueError(
            f"El shapefile de polígonos no tiene sistema de coordenadas definido.\n"
            f"Por favor, asigne un sistema de coordenadas al shapefile.\n"
            f"El proyecto está configurado con: {project_crs}")
    else:
        # Convertir al CRS del proyecto
        poligonos = poligonos.to_crs(project_crs)

    poligonos['AreaTotal'] = poligonos.area
    poligonos = poligonos.dropna(subset=nombreID, axis=0)
    poligonos = poligonos.dissolve(by=nombreID, aggfunc='first')
    poligonos = poligonos.reset_index()
    print(f"  ✓ Polígonos cargados: {len(poligonos)}")

    # Rastreo de celdas con cambio
    coberturasLineaBase = Raster(rutaCoberturasLineaBase)
    coberturasRestauracion = Raster(rutaCoberturasRestauracion)

    cuenca['cobLineaBase'] = coberturasLineaBase.getValueByCoordinates(cuenca.X, cuenca.Y)
    cuenca['cobRestauracion'] = coberturasRestauracion.getValueByCoordinates(cuenca.X, cuenca.Y)

    cuenca['cambio'] = 0
    cuenca.loc[cuenca['cobLineaBase'] != cuenca['cobRestauracion'], 'cambio'] = 1
    celdas_con_cambio = cuenca['cambio'].sum()
    print(f"  ✓ Celdas con cambio identificadas: {celdas_con_cambio}")

    # Inferencia de tamaños de celda
    clsz = []
    clsz.append(np.nanmean(np.sort(cuenca.X.unique())[1:] - np.sort(cuenca.X.unique())[:-1]))
    clsz.append(np.nanmean(np.sort(cuenca.Y.unique())[1:] - np.sort(cuenca.Y.unique())[:-1]))
    clsz = tuple(clsz)
    dx = clsz[0] / 2.
    dy = clsz[1] / 2.
    print(f"  ✓ Tamaño de celda: {clsz[0]:.2f} x {clsz[1]:.2f} m")

    # Crear GeoDataFrame con celdas con cambio
    cuencaCambio = cuenca.loc[cuenca['cambio'] == 1, :].copy()
    geometry = cuencaCambio.apply(lambda x: getGeometry(x, dx, dy), axis=1)
    cuencaCambio = gpd.GeoDataFrame(cuencaCambio, geometry=geometry, crs=project_crs)

    # Intersección
    poligonosInterCuenca = poligonos.overlay(cuencaCambio, how='intersection')
    print(f"  ✓ Intersección calculada: {len(poligonosInterCuenca)} polígonos")

    # Cálculo de áreas y ponderación
    areaCelda = cuencaCambio.area.mean()
    poligonosInterCuenca['Area_m2'] = poligonosInterCuenca.area
    poligonosInterCuenca['PondArea'] = poligonosInterCuenca['Area_m2'] / areaCelda

    # Agregación de área con cambio
    areaConCambio = poligonosInterCuenca.groupby(nombreID)['Area_m2'].agg('sum')
    poligonos = poligonos.join(areaConCambio, on=nombreID, how='left')
    poligonos = poligonos.rename(columns={'Area_m2': 'AreaCambio'})
    poligonos['PercCambio'] = poligonos['AreaCambio'] / poligonos['AreaTotal']

    # ========== PASO 6: Cálculo de índice de priorización ==========
    print("\n[PASO 6/8] Calculando índice de priorización...")

    # Ponderación y normalización de Dk
    columnasDk = [f'Dk_{key}' for key in fcPonderacion.keys()]

    for Dk in columnasDk:
        poligonosInterCuenca[Dk] = poligonosInterCuenca[Dk] * poligonosInterCuenca['PondArea']
        poligonosInterCuenca[f'{Dk}_tru'] = transformacionRangoUnitario(poligonosInterCuenca[Dk], minimo=0)

    # Dk total
    columnasDkTRU = [f'Dk_{key}_tru' for key in fcPonderacion.keys()]
    poligonosInterCuenca['Dk_total'] = 0.

    for DkTRU in columnasDkTRU:
        poligonosInterCuenca['Dk_total'] = poligonosInterCuenca['Dk_total'] + \
            poligonosInterCuenca[DkTRU] * fcPonderacion[DkTRU.split('_')[1]]

    # Agrupación por polígonos
    columnasAnalisis = columnasDk + ['Dk_total', 'Area_m2']
    agrupacionPoligonalDk = poligonosInterCuenca.groupby(nombreID)[columnasAnalisis].agg('sum')
    agrupacionPoligonalDk = agrupacionPoligonalDk.reset_index()

    # Ordenación y área acumulada
    agrupacionPoligonalDk = agrupacionPoligonalDk.sort_values(by='Dk_total', ascending=False)
    areaTotalIntervenible = agrupacionPoligonalDk['Area_m2'].sum()
    agrupacionPoligonalDk['AreaAcu_m2'] = agrupacionPoligonalDk['Area_m2'].cumsum()
    agrupacionPoligonalDk['AreaAcuPer'] = (agrupacionPoligonalDk['AreaAcu_m2'] / areaTotalIntervenible) * 100.

    # Transformación de rango unitario y acumulados
    for Dk in columnasDk + ['Dk_total']:
        agrupacionPoligonalDk[f'{Dk}_tru'] = transformacionRangoUnitario(agrupacionPoligonalDk[Dk], minimo=0)
        agrupacionPoligonalDk[f'{Dk}_tru_acum'] = \
            transformacionRangoUnitario(agrupacionPoligonalDk[f'{Dk}_tru'].cumsum())

    print(f"  ✓ Índice de priorización calculado para {len(agrupacionPoligonalDk)} polígonos")

    # ========== PASO 7: Generar visualizaciones ==========
    print("\n[PASO 7/8] Generando visualizaciones...")
    os.makedirs(os.path.join(path_result, '03_Prioritization'), exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_title(Name_Basin.capitalize(), fontdict={'fontname': 'Calibri', 'fontweight': 'bold',
                                                     'fontsize': '20', 'va': 'bottom', 'color': 'darkslategrey'})
    ax.set_xlabel('Porcentaje acumulado de área intervenida',
                   fontdict={'fontname': 'Calibri', 'fontsize': '18', 'color': 'darkslategrey'})
    ax.set_ylabel('Entrega (Dk) en rango unitario',
                   fontdict={'fontname': 'Calibri', 'fontsize': '18', 'color': 'darkslategrey'})

    plt.xticks(fontname='Calibri', color='darkslategrey', size=14)
    plt.yticks(fontname='Calibri', color='darkslategrey', size=14)

    for col in columnasDkTRU + ['Dk_total_tru']:
        agrupacionPoligonalDk.plot.line(x='AreaAcuPer', y=f'{col}_acum', ax=ax, label=col.replace('_tru', ''))

    ax.grid(color='whitesmoke')
    legend = ax.legend(prop={'family': 'Calibri', 'size': '10'}, loc=4, framealpha=0.5)
    for text in legend.get_texts():
        text.set_color("darkslategrey")

    plt.savefig(os.path.join(path_result, '03_Prioritization', f'Dk_vs_AreaAcumPer_{Name_Basin}.png'), dpi=300)
    plt.close()
    print(f"  ✓ Gráfico guardado: Dk_vs_AreaAcumPer_{Name_Basin}.png")

    # Exportar entregas
    entregas_path = os.path.join(path_result, '03_Prioritization', f'entregas_{Name_Basin}.xlsx')
    agrupacionPoligonalDk[['AreaAcu_m2', 'AreaAcuPer'] + [f'{col}_acum' for col in columnasDkTRU] + ['Dk_total_tru_acum']].to_excel(
        entregas_path, index=False)
    print(f"  ✓ Excel entregas guardado: entregas_{Name_Basin}.xlsx")

    # ========== PASO 8: Generar escenarios de intervención ==========
    print("\n[PASO 8/8] Generando escenarios de intervención...")
    os.makedirs(os.path.join(path_result, '04_Scenarios'), exist_ok=True)

    # Asignar prioridad
    agrupacionPoligonalDk['prioriz'] = np.arange(len(agrupacionPoligonalDk['Dk_total_tru_acum'])) + 1

    # Asignar escenarios de intervención
    for perc in percIntervencion:
        nameCol = f'escen{perc:03}'
        agrupacionPoligonalDk[nameCol] = 0
        agrupacionPoligonalDk.loc[agrupacionPoligonalDk['AreaAcuPer'] <= perc, nameCol] = 1

    # Adecuar índices
    poligonosInterCuenca.index = poligonosInterCuenca[nombreID]
    poligonosInterCuenca.index.name = None
    agrupacionPoligonalDk.index = agrupacionPoligonalDk[nombreID]
    agrupacionPoligonalDk.index.name = None

    # Join del DataFrame de agrupación sobre polígonos intersectados
    poligonosInterCuenca = poligonosInterCuenca.join(agrupacionPoligonalDk, how='left', lsuffix='_polig', rsuffix='_agrup')
    poligonosInterCuenca = poligonosInterCuenca.reset_index()

    # Asignar escenarios a celdas
    for perc in percIntervencion:
        nameCol = f'escen{perc:03}'
        cuenca[nameCol] = 0
        poligonosInterCuencaFilt = poligonosInterCuenca.loc[poligonosInterCuenca[nameCol] == 1, ['idCelda', 'PondArea']]
        poligonosInterCuencaFilt = poligonosInterCuencaFilt.groupby('idCelda').agg('sum')
        indicesCeldas = poligonosInterCuencaFilt[poligonosInterCuencaFilt['PondArea'] >= umbral].index
        cuenca.loc[indicesCeldas, nameCol] = 1

    # Generar coberturas de intervención
    cuencaCoberturas = copy.deepcopy(cuenca[['idCelda', 'tipo', 'X', 'Y', 'cobLineaBase', 'cobRestauracion'] + \
                                            [f'escen{perc:03}' for perc in percIntervencion]])

    for perc in percIntervencion:
        # Asignar cobertura inicial por defecto
        cuencaCoberturas[f'cob{perc:03}'] = cuencaCoberturas['cobLineaBase']
        # Asignar cobertura de restauración a celdas a intervenir
        cuencaCoberturas.loc[cuencaCoberturas[f'escen{perc:03}'] == 1, f'cob{perc:03}'] = \
            cuencaCoberturas['cobRestauracion']

        # Almacenamiento de cobertura en objeto Raster
        sourceDict = {
            'x': cuencaCoberturas['X'],
            'y': cuencaCoberturas['Y'],
            'values': cuencaCoberturas[f'cob{perc:03}']
        }
        cobertura = Raster(source=sourceDict)
        cobertura.exportRaster(os.path.join(path_result, '04_Scenarios', f'CLC_{Name_Basin}_escen{perc:03}.tif'))
        cobertura.exportRaster(os.path.join(path_result, '04_Scenarios', f'CLC_{Name_Basin}_escen{perc:03}.sgabr'))

        # Crear mapa idQ si se proporciona lookup
        if lookup_idQ is not None:
            idQ = cobertura.lookupRaster(lookup=lookup_idQ)
            idQ.exportRaster(os.path.join(path_result, '04_Scenarios', f'idQ_{Name_Basin}_escen{perc:03}.tif'))
            idQ.exportRaster(os.path.join(path_result, '04_Scenarios', f'idQ_{Name_Basin}_escen{perc:03}.sgabr'))

    print(f"  ✓ {len(percIntervencion)} escenarios generados: {percIntervencion}")

    # Join con polígonos
    poligonos.index = poligonos[nombreID]
    poligonos.index.name = None
    poligonos = poligonos.join(agrupacionPoligonalDk, how='left', lsuffix='', rsuffix='_agrup')

    # Exportar shapefile
    shp_path = os.path.join(path_result, '04_Scenarios', f'poligonos_{Name_Basin}.shp')
    poligonos.to_file(filename=shp_path, driver="ESRI Shapefile")
    print(f"  ✓ Shapefile guardado: poligonos_{Name_Basin}.shp")

    # Exportar Excel
    excel_path = os.path.join(path_result, '04_Scenarios', f'poligonos_{Name_Basin}.xlsx')
    poligonos.to_excel(excel_path)
    print(f"  ✓ Excel polígonos guardado: poligonos_{Name_Basin}.xlsx")

    # Exportar cuenca con cambios
    cuenca_shp_path = os.path.join(path_result, '04_Scenarios', f'cuencaCambio_{Name_Basin}.shp')
    cuencaCambio.to_file(filename=cuenca_shp_path, driver="ESRI Shapefile")
    print(f"  ✓ Shapefile cuenca cambio guardado: cuencaCambio_{Name_Basin}.shp")

    print("\n" + "="*80)
    print("PROCESO DE PRIORIZACIÓN COMPLETADO EXITOSAMENTE")
    print("="*80)
    print(f"\nResultados guardados en: {path_result}")
    print(f"  - 01_Inputs/: DataFrame cuenca y rasters individuales")
    print(f"  - 02_Rasters/: Rasters de indicadores")
    print(f"  - 03_Prioritization/: Gráficos y entregas")
    print(f"  - 04_Scenarios/: Shapefiles y escenarios de intervención")

    return cuenca, poligonos, agrupacionPoligonalDk