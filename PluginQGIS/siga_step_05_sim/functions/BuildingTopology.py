# -*- coding: utf-8 -*-
"""
 Python 3.8
--------------------------------------------------------------------------
                           Información Basica
--------------------------------------------------------------------------
 Autores       : Jonathan Nogales Pimentel
                 Miguel Angel Caños Ramos
 Email         : jonathan.nogales@tnc.org
                 miguel.canon@tnc.org
 Fecha         : Abril-2024

--------------------------------------------------------------------------
 Este programa es de uso libre: Usted puede redistribuirlo y/o modificarlo
 bajo los términos de la licencia pública general GNU. El autor no se hace
 responsable de los usos que pueda tener. Para mayor información revisar
 http://www.gnu.org/licenses/

 -------------------------------------------------------------------------
 Descripción del Código
 -------------------------------------------------------------------------
 Este código permite la construcción del archivo topológico que utiliza
 SIGA-CAL para interpretar una cuenca.
 Se parte del código original desarrollado por Gotta Ingeniería S.A.S
    Created on Thursday Oct 13 09:43:24 2022
    @author: Manantiales

 -------------------------------------------------------------------------
 Entradas
 -------------------------------------------------------------------------
 Raster [tif]

 -------------------------------------------------------------------------
 Salidas
 -------------------------------------------------------------------------
"""

# ------------------------------------------------------------------------
# Importar librerías
# ------------------------------------------------------------------------

import os
import sys
import time
import numpy as np
from .Morphology import BuildBasin, Load_scen
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)


def build_topology(dem_path, flow_dir_path, flow_accum_path, shp_embalse, shp_cabeceras, raster_nbs, x_coor, y_coor, escenario,name_project, output_folder, criterio):
    # Ruta del directorio de trabajo
    ProjectPath = os.path.split(os.getcwd())[0]

    # Entradas Obligatorias
    ruta_mde = dem_path
    ruta_dir = flow_dir_path
    ruta_areas = flow_accum_path
    ruta_embalse = shp_embalse
    ruta_cab = shp_cabeceras
    ruta_sbn = raster_nbs
    umbral_longitud = 50.0
    celdas_tramo = 3
    ventana = 3
    percentil = 75.0
    Uarea = 1.0

    tip_top = 'calidad'
    version = "SIGA_CAL_V1.0"

    cierres = {name_project: np.array([[x_coor, y_coor]])}

    for zona in cierres.keys():
        # Inicio cronómetro de ejecución
        ts = time.time()

        print('--------------------------------------------------------------------------------')
        print(f'El procesamiento de la zona "{zona}" se ha iniciado')
        print('--------------------------------------------------------------------------------')

        # Coordenadas de cierre de la cuenca
        XY_cierre = cierres[zona]

        # Coordenadas de las cabeceras - Aplica para cuenca incompleta
        XY_ini = ''

        # Check para creación de carpeta de resultados
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # ----------------------------------------------------------------------------------------------
        # Construcción del objeto cuenca
        #   Criterio 1 - Sencillo (solo umbrales de área y longitud)
        # ----------------------------------------------------------------------------------------------
        if criterio == 1:
            # Construcción de objeto
            cuenca = BuildBasin(areas=ruta_areas,
                                dir=ruta_dir,
                                dem=ruta_mde,
                                xy_c=XY_cierre,
                                xyA_ini=XY_ini,
                                Uarea=Uarea,
                                Ulong=umbral_longitud,
                                Uareap=.25,
                                NoCellTramo=celdas_tramo,
                                ventana=ventana,
                                return_table=True,
                                output_pth=output_folder,
                                bsnobj_pth=os.path.join(output_folder, 'basin.txt'))

        # ----------------------------------------------------------------------------------------------
        # Construcción del objeto cuenca
        #   Criterio 2 - Intermedio (solo umbrales de área y longitud, trazando desde cabeceras)
        # ----------------------------------------------------------------------------------------------
        elif criterio == 2:
            cuenca = BuildBasin(areas=ruta_areas,
                                dir=ruta_dir,
                                dem=ruta_mde,
                                xy_c=XY_cierre,
                                xyA_ini=XY_ini,
                                Uarea=Uarea,
                                Ulong=umbral_longitud,
                                NoCellTramo=celdas_tramo,
                                ventana=ventana,
                                return_table=True,
                                output_pth=output_folder,
                                bsnobj_pth=os.path.join(output_folder, 'basin.txt'),
                                resmask_pth=ruta_embalse)

        # ----------------------------------------------------------------------------------------------
        # Construcción del objeto cuenca
        #   Criterio 3 - Completo (trazado desde cabeceras, teniendo en cuenta percentil y SBN)
        # ----------------------------------------------------------------------------------------------
        elif criterio == 3:
            cuenca = BuildBasin(areas=ruta_areas,
                                dir=ruta_dir,
                                dem=ruta_mde,
                                xy_c=XY_cierre,
                                xyA_ini=XY_ini,
                                Uarea=Uarea,
                                Ulong=umbral_longitud,
                                NoCellTramo=celdas_tramo,
                                ventana=ventana,
                                return_table=True,
                                output_pth=output_folder,
                                bsnobj_pth=os.path.join(output_folder, 'basin.txt'),
                                resmask_pth=ruta_embalse,
                                rivhead_pth=ruta_cab, percent=1.0)

        # ----------------------------------------------------------------------------------------------
        # Construcción del objeto cuenca
        #   Criterio 4 - Completo (trazado desde cabeceras, teniendo en cuenta percentil y SBN)
        # ----------------------------------------------------------------------------------------------
        elif criterio == 4:
            cuenca = BuildBasin(areas_pth=ruta_areas,
                                dir_pth=ruta_dir,
                                dem_pth=ruta_mde,
                                xy_c=XY_cierre,
                                xyA_ini=XY_ini,
                                Ulong=umbral_longitud,
                                NoCellTramo=celdas_tramo,
                                ventana=ventana,
                                return_table=True,
                                output_pth=output_folder,
                                bsnobj_pth=os.path.join(output_folder, 'basin.txt'),
                                resmask_pth=ruta_embalse,
                                rivhead_pth=ruta_cab,
                                percent=percentil,
                                sbn_pth=ruta_sbn)

        # ----------------------------------------------------------------------------------------------
        # Creación de archivo de topología
        # ----------------------------------------------------------------------------------------------
        # Direcciones de los mapas necesarios para crear escenarios
        rutas_mapas = {'dir': ruta_dir,
                       'llanura': os.path.join(output_folder, 'MapGeoVALLE.tif')}

        # Creación del archivo cuenca (topología)
        cuenca2 = Load_scen(output_pth=os.path.join(output_folder, f'Topology_{zona}'),
                            lista_rutas_mapas=rutas_mapas,
                            tipo_cuenca=tip_top,
                            basin=cuenca,
                            version=version)

        # ----------------------------------------------------------------------------------------------
        # Impresión del resumen de procesamiento
        # ----------------------------------------------------------------------------------------------
        print('--------------------------------------------------------------------------------')
        print(f'La zona "{zona}" ha sido procesada exitosamente en {(time.time() - ts) / 60.:.2f} minutos')
        print('--------------------------------------------------------------------------------')
        print('::::: RESUMEN :::::')
        print('--------------------------------------------------------------------------------')
        print(f'Nombre del escenario: {escenario}')
        print(f'Tipo de topología: {tip_top}')
        print(f'Versión de procesamiento: {version}')
        print(f'Coordenadas de cierre: {XY_cierre[0][0]} E,{XY_cierre[0][1]} N')
        print(f'Carpeta de Salida: ...{output_folder[len(ProjectPath):]}')
        print('--------------------------------------------------------------------------------')
        print('')