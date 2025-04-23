import pandas as pd
import numpy as np
# from datetime import datetime
import os

def create_outputs_folders(path_out, name_project, scenario=None):
    main_folders = ['entradas', 'salidas', 'temp']
    inside_folders = ['calibracion', 'ciclos','condicionesfrontera', 'ejecucion', 'mapas', 'movmasa', 'muestras', 'series', 'topologia',
                      'vertimientos']

    project_path = os.path.join(path_out, name_project)

    for main_folder in main_folders:
        main_folder_path = os.path.join(project_path, main_folder)
        os.makedirs(main_folder_path, exist_ok=True)

        if main_folder == 'entradas':
            for inside_folder in inside_folders:
                inside_folder_path = os.path.join(main_folder_path, inside_folder)
                os.makedirs(inside_folder_path, exist_ok=True)

                if inside_folder == 'mapas':
                    scenario_folder = 'BaseLine' if scenario is None else scenario
                    scenario_path = os.path.join(inside_folder_path, scenario_folder)
                    os.makedirs(scenario_path, exist_ok=True)


def escribe_archivo(archivo,encabezados,dic):
    """
    Función para escribir un archivo tipo fichero SIGA
    
    Entradas
    ----------
        archivo : str
                nombre con el que el fichero se va a guardar
    encabezados : lista de strings
                lista con nombre de cada encabezado de fichero, entiéndase por encabezado las líneas que están definidas por []. 
                Nota: No debería variar a no ser que se modifique estructura de archivos
            dic : diccionario
                diccionario que contiene una llave por encabezado, con lo que se llenan los ficheros. 
                Nota: Variable de acuerdo a configuración de la corrida
    """
    archivo=open(archivo,'w')

    for i in range(len(dic)):
        archivo.write(encabezados[i] + '\n' + dic[list(dic.keys())[i]] + '\n\n')
    archivo.close()
    
def diccionario_factores_calibracion(dic_archivos,ruta_csv):
    """
    Función que crea diccionario de factores de calibración con información con la que se va a llenar el fichero. Lee como dataframe los 
    diferentes archivos csv de calibración que están relacionados a cada encabezado del fichero y complementa campos que no se llenan con 
    un dataframe si no con un solo valor, como el caso de [NÚMEROS DE FILAS DE LOS FACTORES MATRICIALES].
    
    Entradas
    ----------
    dic_archivos : diccionario
                diccionario que contiene nombre de archivo csv relacionado con cada encabezado del fichero.
        ruta_csv : str
                    ruta donde se encuentran los archivos csv para construir fichero de calibración.
    Salidas
    ----------
    dic_factores_calibracion : diccionario
                diccionario que contiene inofrmación para llenar fichero de acuerdo a cada encabezado del fichero de calibración.
        
    """
    dic_factores_calibracion = dic_archivos.copy()

    for i in dic_factores_calibracion.keys():
        if dic_factores_calibracion[i].endswith('.csv'):
            df  = pd.read_csv(ruta_csv + dic_factores_calibracion[i], encoding='utf8')
            print('Se leyó archivo: ' + ruta_csv + dic_factores_calibracion[i])
            if i == 'Restauracion':
                n_rest = len(df)
            if i == 'ModeloFenologico':
                n_fenol = len(df)
            if i == 'CargasDifusas':
                n_carga = len(df)
            df_string = df.to_string(index=False)
            df_string = '\n'.join(' '.join(line.split()) for line in df_string.split('\n'))
            dic_factores_calibracion[i] = df_string 

    dic_factores_calibracion['n_filas_matrices'] = 'n_rest n_fenol n_carga\n' + \
                                                   ' '.join(np.array([n_rest,n_fenol,n_carga]).astype(str))
    return dic_factores_calibracion

def diccionario_muestreador(dic_archivos,ruta_csv):
    """
    Función que crea diccionario del muestreador con información con la que se va a llenar el fichero. Lee como dataframe los 
    diferentes archivos csv del muestreador que están relacionados a cada encabezado del fichero y complementa campos que no se llenan con 
    un dataframe si no con un solo valor, como el caso de [NÚMERO DE PUNTOS DE CONTROL].
    
    Entradas
    ----------
    dic_archivos : diccionario
                diccionario que contiene nombre de archivo csv relacionado con cada encabezado del fichero.
        ruta_csv : str
                    ruta donde se encuentran los archivos csv para construir fichero muestreador.
    Salidas
    ----------
    dic_muestreador : diccionario
                diccionario que contiene inoformación para llenar fichero de acuerdo a cada encabezado del fichero muestreador
        
    """
    dic_muestreador = dic_archivos.copy()

    for i in dic_muestreador.keys():
        if dic_muestreador[i].endswith('.csv'):
            df  = pd.read_csv(ruta_csv + dic_muestreador[i], encoding='latin-1')
            print('Se leyó archivo: ' + ruta_csv + dic_muestreador[i])
            if i == 'PuntosControl':
                n_ptos_control = str(len(df))
            if i == 'TramosControl':
                n_tramos_control = str(len(df))

            if df.shape[0]==0 :
                df_string = df.columns.values.tolist()
                df_string = ' '.join(df_string)
            else: 
                df_string = df.to_string(index=False)
                df_string = '\n'.join(' '.join(line.split()) for line in df_string.split('\n'))
                
            dic_muestreador[i] = df_string 

    dic_muestreador['n_ptos_control']   = n_ptos_control
    dic_muestreador['n_tramos_control'] = n_tramos_control
    return dic_muestreador

def diccionario_solo_df(dic_df):
    """
    Esta función se usa para la creación de los diccionarios asociados a ficheros que solo se llenan con información de un DataFrame,
    este es el caso de vertimientos, condiciones de frontera y ciclo anual de cargas difusas.
    
    Entradas
    ----------
    dic_df : diccionario
                diccionario que contiene información asociada a los diferentes encabezados de cada fichero, donde una llave se asocia a
                'datos' que contiene una variable tipo DataFrame
    Salidas
    ----------
    dic_df : diccionario
                diccionario que contiene información asociada a los diferentes encabezados de cada fichero lista para ser escrita.   
    """
    dic_df['datos'] = dic_df['datos'].to_string(index=False)
    dic_df['datos'] = '\n'.join(' '.join(line.split()) for line in dic_df['datos'].split('\n'))

    return dic_df

import time

def diccionario_series(df,df_metadatos,datos_faltantes_series=-9999):
    """
    Esta función se usa para la creación de los diccionarios asociados a ficheros de series de tiempo.
    
    Entradas
    ----------
    df : DataFrame
        DataFrame con series de tiempo de la variable de interés, donde las filas corresponden a fecha y las columnas al punto de monitoreo. Las
        columnas tienen como nombre el código del punto de monitoreo.
    df_metadatos : DataFrame
                DataFrame con los metadatos de df, debe tener por lo menos las columnas 'codigo', asociado a los nombres de las columnas de 
                df, 'x','y', 'z' con las coordenadas de los puntos de monitoreo.
    datos_faltantes_series: float
                            valor de faltantes en la serie de tiempo, se asume por defecto que es -999
    Salidas
    ----------
    dic_df : diccionario
                diccionario que contiene información asociada a los diferentes encabezados de cada fichero lista para ser escrita.   
    """
    df           = df.sort_index(axis=1)
    df_metadatos = df_metadatos.sort_values(by= 'codigo').reset_index(drop=True)

    dic_series = {'num_ser' : str(np.shape(df)[1]),
                   'num_per' : str(np.shape(df)[0]),
                   'datos_faltantes' : str(datos_faltantes_series),
                   'x'       : ' '.join(df_metadatos['x'].values.astype(str)),
                   'y'       : ' '.join(df_metadatos['y'].values.astype(str)),
                   'z'       : ' '.join(df_metadatos['z'].values.astype(str)),
                   'datos'   : df.round(4)
                  }

    dic_series['datos'].insert(loc=0, column='Año', value= dic_series['datos'].index.year.astype(int))
    dic_series['datos'].insert(loc=1, column='Mes', value= dic_series['datos'].index.month.astype(int))
    dic_series['datos'].insert(loc=2, column='Día', value= dic_series['datos'].index.day.astype(int))

    dic_series['datos'] = dic_series['datos'].fillna(datos_faltantes_series)
 
    texto = [' '.join(dic_series['datos'].columns.astype(str))]
    datos = dic_series['datos'].values
    # for i in range(0,len(datos)):
    #     linea=' '.join(datos[i,:].astype(str))
    #     texto.append(linea)

    for i in range(0, len(datos)):
        # Convertir las primeras tres columnas a enteros y luego unirlos con el resto de las columnas
        linea = ' '.join(
            [str(int(datos[i, 0])),  # Año
             str(int(datos[i, 1])),  # Mes
             str(int(datos[i, 2]))]  # Día
            + [f'{x:g}' for x in datos[i, 3:]]  # Convertir los demás valores sin notación científica
        )
        texto.append(linea)
    dic_series['datos'] = '\n'.join(texto)

    return dic_series


class Ficheros(dict):
    """
    Clase para escribir ficheros de entrada del SIGA
    """
    def __init__(self, *args, **kwargs):
        super(Ficheros, self).__init__(*args, **kwargs)
        self.__dict__ = self
        print(self)
        


    def crear_archivo_muestreador(self,dic_archivos,ruta_csv):

        dic_muestreador = diccionario_muestreador(dic_archivos,ruta_csv)


        encabezados =  ['[NÚMERO DE PUNTOS DE CONTROL]', '[LISTADO DE VARIABLES PARA PUNTOS DE CONTROL]', \
                        '[LISTADO DE PUNTOS DE CONTROL]', '[NÚMERO DE TRAMOS DE CONTROL]', \
                        '[LISTADO DE VARIABLES PARA TRAMOS DE CONTROL]', '[LISTADO DE TRAMOS DE CONTROL]', \
                        '[LISTADO DE VARIABLES PARA MAPAS]', '[ETIQUETAS DE VARIABLES]', '[UNIDADES DE VARIABLES]', \
                        '[FACTORES DE CONVERSION DE UNIDADES]', '[LISTADO DE VARIABLES PARA MAPAS SOLO EN CAUCE]', \
                        '[LISTADO DE VARIABLES PARA ENMASCARAR MAPAS EN CAUCE]']

        if os.path.exists(self.ruta_cuenca_SIGA + 'muestras/'):
            pass
        else:
            os.makedirs(self.ruta_cuenca_SIGA + 'muestras/')

        escribe_archivo(self.ruta_cuenca_SIGA + 'muestras/'+ self.archivo_muestreador,encabezados,dic_muestreador)

        print('Aviso de crear_archivo_muestreador: El archivo '+self.archivo_muestreador+' se creó correctamente')



    def crear_archivo_movmasa(self, dic_movmasa):
        if type(dic_movmasa['datos']) != str:
            dic_movmasa = diccionario_solo_df(dic_movmasa)
        encabezados = ['[NÚMERO DE DESLIZAMIENTOS]', '[VALOR DE DATOS FALTANTES]', '[MATRIZ DE DATOS]']
        os.makedirs(self.ruta_cuenca_SIGA +'movmasa', exist_ok=True)
        escribe_archivo(self.ruta_cuenca_SIGA +'movmasa/MovMasa.txt', encabezados, dic_movmasa)
        print('El archivo MovMasa.txt se creó correctamente')
    def crear_archivo_vertimientos(self, dic_df):
        dic_vertimientos = diccionario_solo_df(dic_df)
        encabezados = ['[NÚMERO DE PUNTOS]', '[VALOR DE DATOS FALTANTES]', '[MATRIZ DE DATOS]']
        os.makedirs(self.ruta_cuenca_SIGA +'vertimientos', exist_ok=True)
        escribe_archivo(self.ruta_cuenca_SIGA +'vertimientos/vertimientos.txt', encabezados, dic_vertimientos)
        print('El archivo vertimientos_base_med.txt se creó correctamente')
    def crear_archivo_ciclo_cargas(self, dic_ciclos):
        if type(dic_ciclos['datos']) != str:
            dic_ciclos = diccionario_solo_df(dic_ciclos)
        encabezados = ['[VARIABLE]', '[PERIODO EN DÍAS DEL CICLO]', '[NÚMERO DE CICLOS REPRESENTADOS](2 para un ciclo lineal, cualquier número entero para cualquier otro tipo de ciclo)',
                       '[RESOLUCIÓN TEMPORAL EN DÍAS]', '[MODELO A PARTIR DE VARIBLE SECUNDARIA]%ninguno %lineal', '[CICLO DE LA VARIABLE]']
        os.makedirs(self.ruta_cuenca_SIGA +'ciclos', exist_ok=True)
        escribe_archivo(self.ruta_cuenca_SIGA +'ciclos/ciclos.txt', encabezados, dic_ciclos)
        print('El archivo ciclo_anual_cargasdifusas.txt se creó correctamente')
    def ExportarSeriesSIGA(self,df,df_metadatos,nombre_archivo_serie):

        """
        Exporta datos a un archivo de series de tiempo del modelo SIGA

        Entradas
        --------
        ruta: Cadena de caracteres,
            Ruta para escribir el archivo en disco
        series: Diccionario,
            Contiene, al menos, las siguientes claves
            'num_ser' Entero, Número de series/puntos/estaciones
            'num_per' Entero, Número de periodos
            'X'       Matriz de reales, Coordenadas de localización de los puntos
            'codigos' Códigos de las series a exportar
            'datos'   pandas Dataframe con los datos

        Autor
        -----
        Rendón-Álvarez, J. P., jprendona@gmail.com

        """
        start = time.time()
        dic_series = diccionario_series(df,df_metadatos)
        start2 = time.time()
        print('Tiempo diccionario: ',start2-start)

        encabezados =  ['[NÚMERO DE SERIES]', '[NÚMERO DE REGISTROS]','[VALOR DE DATOS FALTANTES]', '[COORDENADAS X]', \
                        '[COORDENADAS Y]', '[COORDENADAS Z]','[MATRIZ DE DATOS]']

        # Abre el archivo y escribe los metadatos
        if os.path.exists(self.ruta_cuenca_SIGA + 'series/'):
            pass
        else:
            os.makedirs(self.ruta_cuenca_SIGA + 'series/')

        escribe_archivo(self.ruta_cuenca_SIGA + 'series/' + nombre_archivo_serie,encabezados,dic_series)

        start3 = time.time()
        print('Tiempo escritura: ',start3-start2)

        print('{Aviso de ExportarSeriesSIGA: Las series ' + nombre_archivo_serie + ' se exportaron correctamente')


    def crear_archivo_factores_calibracion(self,dic_archivos,ruta_csv):

        dic_factores_calibracion = diccionario_factores_calibracion(dic_archivos,ruta_csv)

        encabezados = ['[NÚMEROS DE FILAS DE LOS FACTORES MATRICIALES]', '[FACTORES DEL MODELO METEOROLOGICO]', \
                       '[FACTORES DEL MODELO DE RESTAURACION DE LA VEGETACION]', '[FACTORES DEL MODELO FENOLOGICO]', \
                       '[FACTORES DEL MODELO HIDROLOGICO]', '[FACTORES DEL MODELO SEDIMENTOLOGICO]', \
                       '[FACTORES DEL MODELO GEOTÉCNICO]', '[FACTORES DE APLICACION DE CARGAS DIFUSAS]', \
                       '[NUMERO DE REZAGOS TEMPORALES A ALMACENAR]',  '[FACTORES DEL MODELO DE DBO CARBONACEA]', \
                       '[FACTORES DEL MODELO DE REAIREACION]', '[FACTORES DEL MODELO DE ESCHERICHIA COLI]', \
                       '[FACTORES DEL MODELO DE NITROGENO]', '[FACTORES DEL MODELO DE FOSFORO]', \
                       '[FACTORES DEL MODELO DE COLIFORMES TOTALES]', '[FACTORES DEL MODELO DE TURBIEDAD]',\
                       '[FACTORES DEL MODELO DE CONDUCTIVIDAD ELECTRICA]', \
                       '[FACTORES PARA EL MODELO DE SOLIDOS SUSPENDIDOS TOTALES]', \
                       '[FACTORES THETA PARA LA CORRECCION DE TASAS POR TEMPERATURA]', \
                       '[FACTORES PARA TEMPERATURA DEL SUELO Y DEL AGUA]','[FACTORES DEL MODELO ADZ-QUASAR]', \
                       '[ESTADO DE LA CUENCA PARA LA EVALUACION DEL ICA-EPM]','[FACTORES PARA MODELO DE pH]']


        if os.path.exists(self.ruta_cuenca_SIGA + 'calibracion/'):
            pass
        else:
            os.makedirs(self.ruta_cuenca_SIGA + 'calibracion/')

        escribe_archivo(self.ruta_cuenca_SIGA + 'calibracion/'+ self.archivo_fact_cal,encabezados,dic_factores_calibracion)

        print('Aviso de crear_archivo_factores_calibracion: El archivo '+self.archivo_fact_cal + ' se creó correctamente')
    def crear_archivo_condiciones_frontera(self,dic_df):

        dic_vertimientos = diccionario_solo_df(dic_df)

        encabezados =  ['[NÚMERO DE PUNTOS]','[VALOR DE DATOS FALTANTES]','[TIPO CONDICIONES DE FRONTERA]','[MATRIZ DE DATOS]']
        os.makedirs(self.ruta_cuenca_SIGA + 'condicionesfrontera', exist_ok=True)
        escribe_archivo(self.ruta_cuenca_SIGA + 'condicionesfrontera/'+ 'BoundaryCondition.txt',encabezados,dic_vertimientos)
        print('Aviso de crear_archivo_vertimientos: El archivo '+'BoundaryCondition.txt'+' se creó correctamente')


    def crear_archivo_ejecucion(self):

        encabezados =  ['[NOMBRE DEL ESCENARIO]', '[FECHA DE INICIO DE LA SIMULACION] Año Mes Día', '[FECHA DE FINALIZACION DE LA SIMULACION] Año Mes Día',\
                        '[TAMAÑO DE PASO TEMPORAL DE SIMULACION EN DÍAS]', '[NOMBRE DEL ARCHIVO DE TOPOLOGÍA]',\
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE TEMPERATURA MÍNIMA]', \
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE TEMPERATURA MÁXIMA]', \
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE HUMEDAD RELATIVA]', \
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE RADIACION SOLAR]', \
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE PRECIPITACION]',\
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE VIENTO ZONAL]', \
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE VIENTO MERIDIONAL]',\
                        '[NOMBRE DEL ARCHIVO DE CAUDAL LÍQUIDO DESCARGADO]', '[NOMBRE DEL ARCHIVO DE CAUDAL LÍQUIDO CAPTADO]', \
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE ARCILLAS MEZCLADAS EN DESCARGAS]',\
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE LIMOS MEZCLADAS EN DESCARGAS]', \
                        '[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE ARENAS MEZCLADAS EN DESCARGAS]',\
                        '[NOMBRE DEL ARCHIVO DE MOVIMIENTOS EN MASA]', '[NOMBRE DEL ARCHIVO DE FACTORES DE CALIBRACION]',\
                        '[NOMBRE DEL ARCHIVO MUESTREADOR DE RESULTADOS]','[NOMBRE DEL ARCHIVO DE SERIES DE TIEMPO DE TASAS DE RESTAURACION]',\
                        '[NOMBRE DEL ARCHIVO DE PUNTOS DE VERTIMIENTO]', '[NOMBRE DEL ARCHIVO DE CICLOS ANUALES DE CARGAS DIFUSAS]',\
                        '[NOMBRE DEL ARCHIVO DE CONDICIONES DE FRONTERA]','[TIPO DE EJECUCION DEL MODELO]', '[EVALUAR MODELO GEOTÉCNICO]',\
                        '[EVALUAR pH CON CUENCA CALIENTE]']

        

        variables = [self.escenario,' '.join(self.fecha_inicio.split('/')),' '.join(self.fecha_fin.split('/')),  \
                     self.paso_temporal_dias, self.archivo_topologia,self.archivo_Tmin, self.archivo_Tmax, self.archivo_HR, \
                     self.archivo_RS, self.archivo_PV, self.archivo_Ux, self.archivo_Uy, self.archivo_Qdesc, \
                     self.archivo_Qcapt, self.archivo_Qarc,self.archivo_Qlim, self.archivo_Qare, self.archivo_movmasa, \
                     self.archivo_fact_cal, self.archivo_muestreador, self.archivo_tasrest, self.archivo_vertimientos,\
                     self.archivo_ciclo_cargas,self.archivo_condiciones_frontera,self.tipo_ejecucion,self.modelo_geotecnico,self.pH_cuencacaliente]

        tit = np.arange(0,len(variables))
        d1  = dict(zip(tit,variables))
        
        if os.path.exists(self.ruta_cuenca_SIGA + 'ejecucion/'):
            pass
        else:
            os.makedirs(self.ruta_cuenca_SIGA + 'ejecucion/')


        escribe_archivo(self.ruta_cuenca_SIGA + 'ejecucion/' + self.archivo_ejecucion,encabezados,d1)


        print('Aviso de crear_archivo_ejecucion: El archivo '+self.archivo_ejecucion + ' se creó correctamente')



