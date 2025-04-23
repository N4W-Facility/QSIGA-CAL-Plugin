# - * - coding: utf - 8 -*-
"""
 Python 3.8
--------------------------------------------------------------------------
                           Información Basica
--------------------------------------------------------------------------
 Autores       : Jonathan Nogales Pimentel
                 Miguel Angel Caños Ramos
 Email         : jonathan.nogales@tnc.org
                 miguel.canon@tnc.org
 Fecha         : Julio-2024

--------------------------------------------------------------------------
 Este programa es de uso libre: Usted puede redistribuirlo y/o modificarlo
 bajo los términos de la licencia publica general GNU. El autor no se hace
 responsable de los usos que pueda tener. Para mayor información revisar
 http://www.gnu.org/licenses/

 -------------------------------------------------------------------------
 Descripción del Código
 -------------------------------------------------------------------------
 Este código contiene la configuración para calibración de todos los
 modelos que integra SIGA-CAL, de acuerdo con la estructura definida por
 la librería SPOTPY. Con esta configuración se realiza la calibración del
 modelo.

"""

# ----------------------------------------------------------------------------------------------------------------------
# Importar librerías
# ----------------------------------------------------------------------------------------------------------------------
import os, sys
PathProject = os.path.split(os.getcwd())[0]
sys.path.append(os.path.join(PathProject,'functions'))
import spotpy
import numpy as np
import pandas as pd
import json
from .ficheros import *

import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
import scipy



# ----------------------------------------------------------------------------------------------------------------------
# Objeto para calibración del módulo hidrológico de SIGA-CAL
# ----------------------------------------------------------------------------------------------------------------------
class Cal_SIGA(object):
    def __init__(self, Param_Min, Param_Max, ProjectPath, JSONPath, NameSIGAFolder,
                 NameCalibration, NameExe, NameSceExe, TS_Obs, OmitObs, FObjFun,
                 TypeFobjCal, PathExeSIGA, NameSIGAModel='Hy', NameFunObj='NASH'):
        """
        Inicializa la clase para la calibración del módulo hidrológico de SIGA-CAL
        Parameters
            Param_Min:
                    Valores minimos de los parámetros de calivrbación
            Param_Max:
            ProjectPath
            JSONPath
            NameSIGAFolder,
            NameCalibration
            NameExe
            NameSceExe
            TS_Obs
            OmitObs
            FObjFun,
            TypeFobjCal
            NameSIGAModel

        Parameters - Modelo hidrológico
            fc_PV:
                    Factor de calibración de la lámina de precipitación vertical
            fc_PH:
                    Factor de calibración de la lámina de precipitación horizontal
            fc_SU0:
                    Factor de calibración del almacenamiento máximo foliar (SU0)
            fc_SU1:
                    Factor de calibración del almacenamiento máximo capilar (SU1)
            fc_SU3:
                    Factor de calibración del almacenamiento máximo gravitacional (SU3)
            fc_ks:
                    Factor de calibración de la permeabilidad subsuperficial (ks)
            fc_kp:
                    Factor de calibración de la permeabilidad subterránea (kp)
            fc_ki:
                    Factor de calibración de la permeabilidad de intercambio subterráneo (ki)
            fc_U2:
                    Factor de calibración de la velocidad de flujo superficial en ladera (U2)
            fc_U3:
                    Factor de calibración de la velocidad de flujo gravitacional (U3)
            fc_U4:
                    Factor de calibración de la velocidad de flujo subterránea (U4)
            fc_U5:
                    Factor de calibración de la velocidad de flujo superficial en cauce (U5)
            us:
                    Umbral [0,1] de retención en páramo - Valor por defecto 0.75
            ui:
                    Umbral [0,1] de liberación en páramo - Valor por defecto 0.45

        Parameters - Modelo Sedimentológico
            fc_E2Smax:
                    Factor de calibración de la capacidad de transporte en ladera (Ad) [0, 10]
            fc_E5Smax:
                    Factor de calibración de la capacidad de transporte en cauce/banca (Ad) [0, 10]
            fc_Bx:
                    Factor de calibración de la erosión en banca (Ad). Si es cero, el proceso de erosión en banca se desactiva [0, 10]
            b:
                    Exponente de desagregación del flujo en cauces (Ad) [1,2]
            Us_arc:
                    Velocidad de sedimentación de las arcillas (m.d-1) [0, 10000]
            Us_lim:
                    Velocidad de sedimentación de los limos (m.d-1) [0, 10000]
            Us_are:
                    Velocidad de sedimentación de las arenas (m.d-1) [0, 10000]
            Ds_arc:
                    Diámetro medio de las arcillas (m) [0, 0.01]
            Ds_lim:
                    Diámetro medio de los limos (m) [0, 0.01]
            Ds_are:
                    Diámetro medio de las arenas (m) [0, 0.01]
            Gs:
                    Gravedad específica del sedimento (Ad) [0, 3]
            rhoS:
                    Densidad de las partículas de suelo (kg.m-3) [0, 3000]

        Parameters - Modelo de Nitrógeno
            ts_desni:
                    Tasa de desnitrificación a 20 °C
            ts_fija:
                    Tasa de fijación a 20 °C
            ts_pltNO3:
                    Tasa de absorción de nitratos por parte de las plantas a 20 °C
            ts_nitri:
                    Tasa de nitrificación a 20 °C
            ts_miner:
                    Tasa de mineralización a 20 °C
            ts_inmov:
                    Tasa de inmovilización a 20 °C
            ts_pltNH4:
                    Tasa de absorción de nitrógeno amoniacal por parte de las plantas a 20 °C
            ts_hidNO:
                    Tasa de hidrólisis del nitrógeno orgánico a 20 °C
            Us_NO:
                    Velocidad de sedimentación del nitrógeno orgánico
            K_odn:
                    Coeficiente exponencial para la inhibición de la desnitrificación por la disminución de la
                    concentración de oxígeno disuelto
            fc_QNO:
                    Factor de calibración que multiplica la carga de nitrógeno orgánico
            fc_QNO3:
                     Factor de calibración que multiplica la carga de nitratos
            fc_QNH4:
                    Factor de calibración que multiplica la carga de nitrógeno amoniacal
            th_N:
                    Base del modelo potencial para la corrección de tasas del modelo de nitrógeno diferentes a la
                    tasa de hidrólisis
            th_hidN:
                    Base del modelo potencial para la corrección de la tasa de hidrólisis

        Parameters - Modelo de Fósforo
            ts_minerP:
                    Tasa de mineralización del fósforo a 20 °C
            ts_inmovP:
                    Tasa de inmovilización del fósforo a 20 °C
            ts_pltPO:
                    Tasa de absorción de fósforo orgánico por parte de las plantas a 20 °C
            ts_pltPI:
                    Tasa de absorción de fósforo inorgánico por parte de las plantas a 20 °C
            ts_fbPOin:
                    Tasa de adherencia de fósforo orgánico
            ts_fbPOout:
                    Tasa de liberación de fósforo orgánico
            ts_fbPIin:
                    Tasa de adherencia de fósforo inorgánico
            ts_fbPIout:
                    Tasa de liberación de fósforo inorgánico
            ts_hidPO:
                    Tasa de hidrólisis del fósforo orgánico
            Us_PO:
                    Velocidad de sedimentación del fósforo orgánico
            Us_PI:
                    Velocidad de sedimentación del fósforo inorgánico
            fc_QPO:
                    Factor de calibración que multiplica la carga de fósforo orgánico
            fc_QPI:
                    Factor de calibración que multiplica la carga de fósforo inorgánico
            th_hidP:
                    Base del modelo potencial para la corrección de la tasa de hidrólisis de fósforo orgánico
            th_P:
                    Base del modelo potencial para la corrección de las tasas de fósforo en ladera

        Parameters - DBO
            ts_d:
                    Tasa de descomposición de la DBO carbonácea a 20°C
            Us_cdbo:
                    Velocidad de sedimentación de la DBO carbonácea
            th_cdbo:
                    Base del modelo potencia para la corrección de tasa de descomposición de la CDBO

        Parameters - Oxígeno Disuelto
            mod_rair:
                    Modelo para el cálculo de la tasa de reaireación
                    (1: Tsivoglou y Neal, 2: Oconnor y Dobbins, 3: Churchill, 4: Owen y Gibss, 5: Covar)
            fc_ts_ra:
                    Factor de calibración que multiplica a la tasa de reaireación
            th_rair:
                    Base del modelo potencial para la corrección de la tasa de reaireación

        Parameters - EColi
            ts_regen:
                    Tasa neta de muerte y recrecimiento de las bacterias
            mrt_min:
                    Muerte mínima de bacterias
            Us_sst:
                    Velocidad de sedimentación de los SST
            fr_bs:
                    Fracción de bacterias adheridas a los sedimentos
            fc_QEC:
                    Factor de calibración que multiplica la carga de Escherichia coli
            th_bact:
                    Base del modelo potencial para la corrección de la tasa de muerte y recrecimiento de bacterias

        Parameters - Turbiedad
            K_turb:
                    Coeficiente de la relación potencial entre turbiedad (UNT) y sólidos suspendidos totales (kg/m3)
            exp_turb:
                    Exponente de la relación potencial entre turbiedad (UNT) y sólidos suspendidos totales (kg/m3)

        Parameters - Solidos Suspendidos Totales
            ts_d:
                    Tasa de descomposición de la DBO carbonácea a 20°C
            Us_cdbo:
                    Velocidad de sedimentación de la DBO carbonácea
            mod_rair:
                    Modelo para el cálculo de la tasa de reaireación

        Parameters - Conductividad eléctrica
            b_CE:
                    Conductividad eléctrica base de las corrientes a 25 °C
            alfa:
                    Factor de corrección de la conductividad eléctrica por temperatura. Su valor varía entre
                    0.019 y 0.023 °C-1

        Parameters - Temperatura
            t_ss:
                    Coeficiente de la temperatura del suelo para el cálculo de la temperatura del suelo
            t_sa:
                    Coeficiente de la temperatura del aire para el cálculo de la temperatura del suelo
            t_ww:
                    Coeficiente de la temperatura del agua para el cálculo de la temperatura del agua
            t_wa:
                    Coeficiente de la temperatura del aire para el cálculo de la temperatura del agua

        Parameters - pH
            pCO2:

            roc:

            fc_alk:
        """

        if NameSIGAModel == 'Hy':
            self.params = [spotpy.parameter.Uniform('fc_PV',        Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('fc_PH',        Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('fc_SU0',       Param_Min[2],   Param_Max[2]),
                           spotpy.parameter.Uniform('fc_SU1',       Param_Min[3],   Param_Max[3]),
                           spotpy.parameter.Uniform('fc_SU3',       Param_Min[4],   Param_Max[4]),
                           spotpy.parameter.Uniform('fc_ks',        Param_Min[5],   Param_Max[5]),
                           spotpy.parameter.Uniform('fc_kp',        Param_Min[6],   Param_Max[6]),
                           spotpy.parameter.Uniform('fc_ki',        Param_Min[7],   Param_Max[7]),
                           spotpy.parameter.Uniform('fc_U2',        Param_Min[8],   Param_Max[8]),
                           spotpy.parameter.Uniform('fc_U3',        Param_Min[9],   Param_Max[9]),
                           spotpy.parameter.Uniform('fc_U4',        Param_Min[10],  Param_Max[10]),
                           spotpy.parameter.Uniform('fc_U5',        Param_Min[11],  Param_Max[11]),
                           spotpy.parameter.Uniform('us',           Param_Min[12],  Param_Max[12]),
                           spotpy.parameter.Uniform('ui',           Param_Min[13],  Param_Max[13]),
                           ]

        elif NameSIGAModel == 'Sed':
            self.params = [spotpy.parameter.Uniform('fc_E2Smax',    Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('fc_E5Smax',    Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('fc_Bx',        Param_Min[2],   Param_Max[2]),
                           spotpy.parameter.Uniform('b',            Param_Min[3],   Param_Max[3]),
                           spotpy.parameter.Uniform('Us_arc',       Param_Min[4],   Param_Max[4]),
                           spotpy.parameter.Uniform('Us_lim',       Param_Min[5],   Param_Max[5]),
                           spotpy.parameter.Uniform('Us_are',       Param_Min[6],   Param_Max[6]),
                           spotpy.parameter.Uniform('Ds_arc',       Param_Min[7],   Param_Max[7]),
                           spotpy.parameter.Uniform('Ds_lim',       Param_Min[8],   Param_Max[8]),
                           spotpy.parameter.Uniform('Ds_are',       Param_Min[9],   Param_Max[9]),
                           spotpy.parameter.Uniform('Gs',           Param_Min[10],  Param_Max[10]),
                           spotpy.parameter.Uniform('rhoS',         Param_Min[11],  Param_Max[11]),
                           ]

        elif (NameSIGAModel == 'WQ_NT')|(NameSIGAModel == 'WQ_NO')|(NameSIGAModel == 'WQ_NH4')|(NameSIGAModel == 'WQ_NO3'):
            self.params = [spotpy.parameter.Uniform('ts_desni',     Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('ts_fija',      Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('ts_pltNO3',    Param_Min[2],   Param_Max[2]),
                           spotpy.parameter.Uniform('ts_nitri',     Param_Min[3],   Param_Max[3]),
                           spotpy.parameter.Uniform('ts_miner',     Param_Min[4],   Param_Max[4]),
                           spotpy.parameter.Uniform('ts_inmov',     Param_Min[5],   Param_Max[5]),
                           spotpy.parameter.Uniform('ts_pltNH4',    Param_Min[6],   Param_Max[6]),
                           spotpy.parameter.Uniform('ts_hidNO',     Param_Min[7],   Param_Max[7]),
                           spotpy.parameter.Uniform('Us_NO',        Param_Min[8],   Param_Max[8]),
                           spotpy.parameter.Uniform('K_odn',        Param_Min[9],   Param_Max[9]),
                           spotpy.parameter.Uniform('fc_QNO',       Param_Min[10],  Param_Max[10]),
                           spotpy.parameter.Uniform('fc_QNO3',      Param_Min[11],  Param_Max[11]),
                           spotpy.parameter.Uniform('fc_QNH4',      Param_Min[12],  Param_Max[12]),
                           spotpy.parameter.Uniform('th_N',         Param_Min[13],  Param_Max[13]),
                           spotpy.parameter.Uniform('th_hidN',      Param_Min[14],  Param_Max[14]),
                           ]

        elif (NameSIGAModel == 'WQ_PT')|(NameSIGAModel == 'WQ_PO')|(NameSIGAModel == 'WQ_PI'):
            self.params = [spotpy.parameter.Uniform('ts_minerP',    Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('ts_inmovP',    Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('ts_pltPO',     Param_Min[2],   Param_Max[2]),
                           spotpy.parameter.Uniform('ts_pltPI',     Param_Min[3],   Param_Max[3]),
                           spotpy.parameter.Uniform('ts_fbPOin',    Param_Min[4],   Param_Max[4]),
                           spotpy.parameter.Uniform('ts_fbPOout',   Param_Min[5],   Param_Max[5]),
                           spotpy.parameter.Uniform('ts_fbPIin',    Param_Min[6],   Param_Max[6]),
                           spotpy.parameter.Uniform('ts_fbPIout',   Param_Min[7],   Param_Max[7]),
                           spotpy.parameter.Uniform('ts_hidPO',     Param_Min[8],   Param_Max[8]),
                           spotpy.parameter.Uniform('Us_PO',        Param_Min[9],   Param_Max[9]),
                           spotpy.parameter.Uniform('Us_PI',        Param_Min[10],  Param_Max[10]),
                           spotpy.parameter.Uniform('fc_QPO',       Param_Min[11],  Param_Max[11]),
                           spotpy.parameter.Uniform('fc_QPI',       Param_Min[12],  Param_Max[12]),
                           spotpy.parameter.Uniform('th_hidP',      Param_Min[13],  Param_Max[13]),
                           spotpy.parameter.Uniform('th_P',         Param_Min[14],  Param_Max[14]),
                           ]

        elif NameSIGAModel == 'WQ_DBO':
            self.params = [spotpy.parameter.Uniform('ts_d',         Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('Us_cdbo',      Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('th_cdbo',      Param_Min[2],   Param_Max[2]),
                           ]

        elif NameSIGAModel == 'WQ_OD':
            self.params = [spotpy.parameter.Uniform('mod_rair',     Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('fc_ts_ra',     Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('th_rair',      Param_Min[2],   Param_Max[2]),
                           ]

        elif NameSIGAModel == 'WQ_EColi':
            self.params = [spotpy.parameter.Uniform('ts_regen',     Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('mrt_min',      Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('Us_sst',       Param_Min[2],   Param_Max[2]),
                           spotpy.parameter.Uniform('fr_bs',        Param_Min[3],   Param_Max[3]),
                           spotpy.parameter.Uniform('fc_QEC',       Param_Min[4],   Param_Max[4]),
                           spotpy.parameter.Uniform('th_bact',      Param_Min[5],   Param_Max[5]),
                           ]

        elif NameSIGAModel == 'WQ_CT':
            self.params = [spotpy.parameter.Uniform('K_ct',         Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('exp_ct',       Param_Min[1],   Param_Max[1]),
                           ]

        elif NameSIGAModel == 'WQ_Tur':
            self.params = [ spotpy.parameter.Uniform('K_turb',      Param_Min[0],   Param_Max[0]),
                            spotpy.parameter.Uniform('exp_turb',    Param_Min[1],   Param_Max[1]),
                           ]

        elif NameSIGAModel == 'WQ_SST':
            self.params = [spotpy.parameter.Uniform('fr_arc',       Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('fr_lim',       Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('fr_are',       Param_Min[2],   Param_Max[2]),
                           ]

        elif NameSIGAModel == 'WQ_Co':
            self.params = [spotpy.parameter.Uniform('b_CE',         Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('alfa',         Param_Min[1],   Param_Max[1]),
                           ]

        elif NameSIGAModel == 'WQ_T':
            self.params = [spotpy.parameter.Uniform('t_ss',         Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('t_sa',         Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('t_ww',         Param_Min[2],   Param_Max[2]),
                           spotpy.parameter.Uniform('t_wa',         Param_Min[3],   Param_Max[3]),
                           ]

        elif NameSIGAModel == 'WQ_pH':
            self.params = [spotpy.parameter.Uniform('pCO2',         Param_Min[0],   Param_Max[0]),
                           spotpy.parameter.Uniform('roc',          Param_Min[1],   Param_Max[1]),
                           spotpy.parameter.Uniform('fc_alk',       Param_Min[2],   Param_Max[2]),
                           ]

        # Nombre de la variable para función objetivo
        self.NameVarSIGA    = Select_NameVarSIGA(NameSIGAModel)
        # Ruta de la ruta de trabajo
        self.ProjectPath    = ProjectPath
        # Ruta del archivo JSON con la info de archivos del modelo
        self.JSONPath       = JSONPath
        # Nombre de la carpeta que corresponde a una configuración de SIGA-CAL
        self.NameSIGAFolder = NameSIGAFolder
        # Serie de tiempo de caudales para calibración
        self.TS_Obs         =  TS_Obs
        # Nombre del archivo de ejecución para calibración
        self.NameExe        = NameExe
        # Nombre del escenario de ejecución
        self.NameSceExe     = NameSceExe
        # Nombre del archivo de parámetros
        self.NameCal        = NameCalibration
        # Cantidad de datos iniciales que se omiten para la evaluación de la función objetivo
        self.OmitObs        = OmitObs
        # Factor multiplicador de la función objetivo
        self.FacObjF        = FObjFun
        # Tipo de calibración [Time Series = 0 | t-Student = 1]
        self.TypeFobjCal    = TypeFobjCal
        # Modelo a calibrar [Hy, Sed, WQ_N, WQ_P, WQ_DBO, WQ_OD, WQ_CT, WQ_EColi, WQ_T, WQ_Co, WQ_pH]
        self.NameSIGAModel  = NameSIGAModel
        # Path del ejecutable de SIGA
        self.PathExeSIGA    = PathExeSIGA
        # Nombre dela función obejtivo
        self.NameFunObj     = NameFunObj
        # Pendiente del tanque 4
        self.S_S4           = 1

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # Configurar vector de parámetros
        Params = vector[:]

        self.TS_Sim, self.S_S4, TS_S4 =  simulationSIGA(Params,
                                                        self.ProjectPath,
                                                        self.JSONPath,
                                                        self.NameSIGAFolder,
                                                        self.NameSIGAModel,
                                                        self.NameVarSIGA,
                                                        self.NameSceExe,
                                                        self.NameExe,
                                                        self.PathExeSIGA,
                                                        self.NameCal,
                                                        self.TS_Obs)
        return self.TS_Sim.iloc[:, 0].values

    def evaluation(self):
        return self.TS_Obs.iloc[:,0].values

    def objectivefunction(self, simulation, evaluation):
        return OF(self.OmitObs,self.FacObjF,self.TS_Sim,self.TS_Obs,self.TypeFobjCal, self.NameFunObj, self.S_S4)


def Select_NameVarSIGA(NameSIGAModel):
    if NameSIGAModel == 'Hy':
        NameVarSIGA = 'Q5(m3.d-1)'

    elif NameSIGAModel == 'Sed':
        NameVarSIGA = 'ESDET(m3)'

    elif (NameSIGAModel == 'WQ_NT') | (NameSIGAModel == 'WQ_NO') | (NameSIGAModel == 'WQ_NH4') | (NameSIGAModel == 'WQ_NO3'):
        if (NameSIGAModel == 'WQ_NO'):
            NameVarSIGA = 'NO(kg.m-3)'
        elif (NameSIGAModel == 'WQ_NH4'):
            NameVarSIGA = 'NH4(kg.m-3)'
        elif (NameSIGAModel == 'WQ_NO3'):
            NameVarSIGA = 'NO3(kg.m-3)'

    elif (NameSIGAModel == 'WQ_PT') | (NameSIGAModel == 'WQ_PO') | (NameSIGAModel == 'WQ_PI'):
        if (NameSIGAModel == 'WQ_PO'):
            NameVarSIGA = 'PO(kg.m-3)'
        elif (NameSIGAModel == 'WQ_PI'):
            NameVarSIGA = 'PI(kg.m-3)'

    elif NameSIGAModel == 'WQ_DBO':
        NameVarSIGA = 'CDBO5(kg.m-3)'

    elif NameSIGAModel == 'WQ_OD':
        NameVarSIGA = 'OD(kg.m-3)'

    elif NameSIGAModel == 'WQ_EColi':
        NameVarSIGA = 'EC(ufc.m-3)'

    elif NameSIGAModel == 'WQ_CT':
        NameVarSIGA = 'CT(ufc.m-3)'

    elif NameSIGAModel == 'WQ_Tur':
        NameVarSIGA = 'TUR(NTU)'

    elif NameSIGAModel == 'WQ_SST':
        NameVarSIGA = 'SST(kg.m-3)'

    elif NameSIGAModel == 'WQ_Co':
        NameVarSIGA = 'CE(uS.cm-1)'

    elif NameSIGAModel == 'WQ_T':
        NameVarSIGA = ''

    elif NameSIGAModel == 'WQ_pH':
        NameVarSIGA = 'pH(Unidades_pH)'

    return NameVarSIGA

def simulationSIGA(Params, ProjectPath, JSONPath, NameSIGAFolder, NameSIGAModel,
                   NameVarSIGA,NameSceExe, NameExe, PathExeSIGA, NameCal, TS_Obs):
    # Factor de pendiente inicial
    TS_S4   = 0
    S_S4    = 1

    # Modificar y escribir archivo de parámetros
    if NameSIGAModel == 'Hy':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='Hy')
    elif NameSIGAModel == 'Sed':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='Sed')
    elif (NameSIGAModel == 'WQ_NT') | (NameSIGAModel == 'WQ_NO') | (NameSIGAModel == 'WQ_NH4') | (NameSIGAModel == 'WQ_NO3'):
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[:-2], TypeModel='N')
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[-2:], TypeModel='Theta', Posi=[2, 3])
    elif (NameSIGAModel == 'WQ_PT') | (NameSIGAModel == 'WQ_PO') | (NameSIGAModel == 'WQ_PI'):
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[:-2], TypeModel='P')
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[-2:], TypeModel='Theta', Posi=[5, 6])
    elif NameSIGAModel == 'WQ_DBO':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[:-1], TypeModel='DBO')
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[-1], TypeModel='Theta', Posi=[0])
    elif NameSIGAModel == 'WQ_OD':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[:-1], TypeModel='OD')
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[-1], TypeModel='Theta', Posi=[4])
    elif NameSIGAModel == 'WQ_EColi':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[:-1], TypeModel='EColi')
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params[-1], TypeModel='Theta', Posi=[1])
    elif NameSIGAModel == 'WQ_CT':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='CT')
    elif NameSIGAModel == 'WQ_Tur':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='Tur')
    elif NameSIGAModel == 'WQ_SST':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='SST')
    elif NameSIGAModel == 'WQ_T':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='T')
    elif NameSIGAModel == 'WQ_Co':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='Co')
    elif NameSIGAModel == 'WQ_pH':
        Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='pH')

    # Crear archivo de calibración en formato SIGA
    Create_CalibrationFile(ProjectPath, JSONPath, NameSIGAFolder, NameCal)

    # Ejecutar modelo
    ExeSIGAModel(ProjectPath, NameSIGAFolder, NameExe, PathExeSIGA)

    # Leer serie de tiempo
    FolderExe = get_latest_folder(os.path.join(ProjectPath, NameSIGAFolder, 'salidas', NameSceExe))

    TS_Sim = pd.DataFrame(index=TS_Obs.index)

    for NF in TS_Obs.columns:
        # Nombre del archivo de salida del modelo
        FilePathTS = os.path.join(ProjectPath, NameSIGAFolder, 'salidas', NameSceExe, FolderExe,
                                  'series', NF + '_valor.txt')

        # Leer serie de tiempo
        TS_Sim_i = ReadTimeSeries(FilePathTS)

        # m3/día -> m3/s
        Factor1 = 1 / (3600 * 24)

        # Pasar de kg/m3 a mg/l
        Factor2 = 1000

        # pasar 1/m3 a 1/l
        Factor3 = 1 / 1000

        if NameSIGAModel == 'Hy':
            # Ajustar periodo temporal
            TS_Sim_i = pd.merge(pd.DataFrame(index=TS_Obs.index), TS_Sim_i[NameVarSIGA], left_index=True,right_index=True, how='outer')
            TS_Sim[NF] = TS_Sim_i.values * Factor1

        elif NameSIGAModel == 'Sed':
            # Ajustar series
            TS_Sim_i = pd.merge(pd.DataFrame(index=TS_Obs.index), TS_Sim_i[[NameVarSIGA, 'Q5(m3.d-1)']],left_index=True, right_index=True, how='outer')
            TS_Sim[NF] = ((TS_Sim_i[NameVarSIGA].values * Params[-1]) / TS_Sim_i['Q5(m3.d-1)'].values) * Factor2

        elif NameSIGAModel == 'WQ_NT':
            TS_Sim_i = pd.merge(pd.DataFrame(index=TS_Obs.index),TS_Sim_i[['NO3(kg.m-3)', 'NH4(kg.m-3)', 'NO(kg.m-3)']], left_index=True,right_index=True, how='outer')
            TS_Sim[NF] = TS_Sim_i.sum(axis=1).values * Factor2

        elif NameSIGAModel == 'WQ_PT':
            TS_Sim_i = pd.merge(pd.DataFrame(index=TS_Obs.index),TS_Sim_i[['NO3(kg.m-3)', 'NH4(kg.m-3)', 'NO(kg.m-3)']], left_index=True,right_index=True, how='outer')
            TS_Sim[NF] = TS_Sim_i.sum(axis=1).values * Factor2

        elif ((NameSIGAModel == 'WQ_NO') | (NameSIGAModel == 'WQ_NH4') | (NameSIGAModel == 'WQ_NO3') | (NameSIGAModel == 'WQ_PO') | \
              (NameSIGAModel == 'WQ_PO') | (NameSIGAModel == 'WQ_DBO') | (NameSIGAModel == 'WQ_OD') | (NameSIGAModel == 'WQ_SST')):
            TS_Sim_i = pd.merge(pd.DataFrame(index=TS_Obs.index), TS_Sim_i[NameVarSIGA], left_index=True,right_index=True, how='outer')
            TS_Sim[NF] = TS_Sim_i.values * Factor2

        elif (NameSIGAModel == 'WQ_Co') | (NameSIGAModel == 'WQ_Tur') | (NameSIGAModel == 'WQ_pH') | (NameSIGAModel == 'WQ_T'):
            TS_Sim_i = pd.merge(pd.DataFrame(index=TS_Obs.index), TS_Sim_i[NameVarSIGA], left_index=True,right_index=True, how='outer')
            TS_Sim[NF] = TS_Sim_i.values

        elif (NameSIGAModel == 'WQ_EColi') | (NameSIGAModel == 'WQ_CT'):
            TS_Sim_i = pd.merge(pd.DataFrame(index=TS_Obs.index), TS_Sim_i[NameVarSIGA], left_index=True,right_index=True, how='outer')
            TS_Sim[NF] = TS_Sim_i.values * Factor3

    if NameSIGAModel == 'Hy':
        # Nombre del archivo de salida del modelo
        FilePathTS = os.path.join(ProjectPath, NameSIGAFolder, 'salidas', NameSceExe, FolderExe, 'series', 'RiverMounthBasin_media.txt')

        # Leer serie de tiempo
        Tmp = ReadTimeSeries(FilePathTS)

        # Ajustar periodo temporal
        TS_S4 = pd.merge(pd.DataFrame(index=TS_Obs.index), Tmp['S4(m)'], left_index=True, right_index=True,how='outer') #TANQUE 2, TANQUE 3

        # Ajuste lineal
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.arange(0, len(TS_S4['S4(m)']), 1),TS_S4['S4(m)'] * 1000)
        S_S4 = 1 + abs(slope)

    return TS_Sim, S_S4, TS_S4

def OF(OmitObs,FacObjF,TS_Sim,TS_Obs,TypeFobjCal,NameFunObj, OtherFactor):
    """
    Esta función cálcula la metríca de desempeño para la calibración del modelo

    Parameters
        OmitObs:
                Número de datos iniciales que se omiten para la estimación de la metríca
        FacObjF:
                Factor multiplicador para cambiar el sentido de la metríca de calibración.
                -1 Maximización | 1 Minimización
        TS_Sim:
                Valores simulados por el modelo
        TS_Obs:
                Valores observador
        TypeFobjCal:
                Tipo de metríca a estimar. 1: t-Student, 0: Serie de tiempo
    """

    # Cantidad de datos iniciales que se omiten para la evaluación de la función objetivo
    n = OmitObs

    ObjFun = 0
    for i in range(len(TS_Sim.columns)):
        # Serie de tiempo simulada
        Sim = TS_Sim.iloc[n:, i].values

        # Serie de tiempo observada
        Obs = TS_Obs.iloc[n:, i].values

        # Tomar solo los valores que no son NaN
        id = np.where(~np.isnan(Obs) & ~np.isnan(Sim))
        Obs = Obs[id]
        Sim = Sim[id]

        if TypeFobjCal == 0:
            if NameFunObj == 'NASH':
                # Metrica 1 - Coefficient of determination
                Part_1 = np.nansum((Obs - Sim) ** 2)
                Part_2 = np.nansum((Obs - np.nanmean(Obs)) ** 2)
                ObjFun += FacObjF * (Part_1 / Part_2)
            elif NameFunObj == 'IoAd':
                Part_1 = np.nansum((Obs - Sim) ** 2)
                Part_2 = np.nansum((np.abs((Obs - np.nanmean(Obs))) + np.abs((Sim - np.nanmean(Obs)))) ** 2)
                ObjFun += FacObjF * (Part_1 / Part_2)
            elif NameFunObj == 'EAN':
                Part_1 = np.nansum(np.abs(Obs - Sim))
                Part_2 = np.nansum(np.abs(Obs - np.nanmedian(Obs)))
                ObjFun += FacObjF*(Part_1 / Part_2)
            elif NameFunObj == 'r':
                ObjFun += FacObjF*spotpy.objectivefunctions.correlationcoefficient(Obs, Sim)
            elif NameFunObj == 'R2':
                ObjFun += FacObjF*spotpy.objectivefunctions.rsquared(Obs, Sim)
            elif NameFunObj == 'MSE':
                ObjFun += FacObjF*spotpy.objectivefunctions.mse(Obs, Sim)
            elif NameFunObj == 'RMSE':
                ObjFun += FacObjF*spotpy.objectivefunctions.rmse(Obs, Sim)
            elif NameFunObj == 'RRMSE':
                ObjFun += FacObjF*spotpy.objectivefunctions.rrmse(Obs, Sim)
            elif NameFunObj == 'COMPOS':
                # Metrica 1 - Coefficient of determination
                Part_1 = np.nansum((Obs - np.nanmean(Obs)) * (Sim - np.nanmean(Sim)))
                Part_2 = np.nansum((Obs - np.nanmean(Obs)) ** 2) * np.nansum((Sim - np.nanmean(Sim)) ** 2)
                RSqr = (Part_1 / np.sqrt(Part_2)) ** 2

                # Metrica 2 - Nash–Sutcliffe coefficient
                Part_1 = np.nansum((Obs - Sim) ** 2)
                Part_2 = np.nansum((Obs - np.nanmean(Obs)) ** 2)
                NSE = (Part_1 / Part_2)

                # Metrica 3 - Index of agreement
                Part_2 = np.nansum((np.abs((Obs - np.nanmean(Obs))) + np.abs((Sim - np.nanmean(Obs)))) ** 2)
                IoAd = (Part_1 / Part_2)

                # Metrica 4 - Error absoluto normalizado
                Part_1 = np.nansum(np.abs(Obs - Sim))
                Part_2 = np.nansum(np.abs(Obs - np.nanmedian(Obs)))
                NAE = (Part_1 / Part_2)

                # La función objetivo es una combinación lineal de las métricas anteriores donde los pesos que se asignan son
                # los siguientes:
                # RSqr:  0.2
                # NSE:   0.4
                # IoAd:  0.2
                # NAE:   0.2
                ObjFun += FacObjF * ((0.2 * (1 - RSqr)) + (0.4 * NSE) + (0.2 * IoAd) + (0.2 * NAE))
        else:
            t_stat, p_valor = scipy.stats.ttest_ind(Obs, Sim)
            ObjFun += FacObjF * t_stat

    return (ObjFun / len(TS_Sim.columns))*OtherFactor

def get_latest_folder(directory):
    """
    Esta función estima el último directorio

    Parameters
        directory:
                Ruta donde se almacenan las simulaciones
    """

    # Obtiene una lista de todas las subcarpetas en el directorio especificado
    folders = [f for f in Path(directory).iterdir() if f.is_dir()]

    if not folders:
        return None

    # Ordena las carpetas por la fecha de modificación (del más reciente al más antiguo)
    latest_folder = max(folders, key=lambda f: f.stat().st_mtime)

    # Devuelve solo el nombre de la carpeta
    return latest_folder.name

def Modify_csvHyParms(ProjectPath, JSONPath, NameSIGAFolder, Params, TypeModel='Hy',Posi=[]):
    """
    Esta función modifica el archivo csv de parámetros hidrológicos con el cual se configura
    el archivo txt de calibración de SIGA

    Parameters
        ProjectPath:
                Ruta del proyecto definido por el usuario
        JSONPath:
                Ruta donde se localiza el archivo JSON con la configuración de archivos
        NameSIGAFolder:
                Nombre de carpeta del proyecto
        Params:
                Vector de parámetros
        TypeModel:
                Sigla del modelo al cual pertenecen los parámetros
        Posi:
                Posición específica para modificar los parámetros
    """

    # Leer estructura de datos en json
    with open(JSONPath) as json_file:
        DicFicCal = json.load(json_file)

    if TypeModel == 'Hy':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["ModeloHidrologico"]
    elif TypeModel == 'Sed':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["ModeloSedimentologico"]
    elif TypeModel == 'Geo':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["ModeloGeotecnico"]
    elif TypeModel == 'CaDi':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["CargasDifusas"]
    elif TypeModel == 'DBO':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["DBOCarbonicea"]
    elif TypeModel == 'OD':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["Reaireacion"]
    elif TypeModel == 'EColi':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["EColi"]
    elif TypeModel == 'N':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["Nitrogeno"]
    elif TypeModel == 'P':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["Fosforo"]
    elif TypeModel == 'CT':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["ColiformesTotales"]
    elif TypeModel == 'Tur':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["Turbiedad"]
    elif TypeModel == 'Co':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["ConductividadElectrica"]
    elif TypeModel == 'SST':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["SolidosSuspendidosTotales"]
    elif TypeModel == 'Theta':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["CorreccionTemperatura"]
    elif TypeModel == 'T':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["Temperatura"]
    elif TypeModel == 'ADZ':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["FactoresADZQUASAR"]
    elif TypeModel == 'pH':
        # Nombre del archivo de parámetros
        FileName = DicFicCal["dic_factores_calibracion"]["pH"]

    # Ruta del archivo CSV de parámetros del Modelo SIGA-CAL

    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    PathFile = os.path.join(parent_dir,'Data','archivos_csv_Cal', 'calibracion', FileName)

    # Leer parámetros
    HyParams = pd.read_csv(PathFile)

    # Remplazar parámetros
    if not Posi:
        HyParams.iloc[:]    =  Params
    else:
        HyParams.iloc[0,Posi] = Params

    # Guardar archivo de parámetros
    HyParams.to_csv(PathFile,index=False)


def Modify_csvMuestreador_CP( JSONPath, NameSIGAFolder, Catalog, NameSIGAModel, X, Y):
    """
    Esta función modifica el archivo csv de puntos de control con el cual se configura
    el archivo txt de muestreador de SIGA

    Parameters
        ProjectPath:
                Ruta del proyecto definido por el usuario
        JSONPath:
                Ruta donde se localiza el archivo JSON con la configuración de archivos
        NameSIGAFolder:
                Nombre de carpeta del proyecto
        Catalog:
                Catálogo de estaciones con las coordenadas X y Y de los puntos de control y su respectivo nombre
        NameSIGAModel:
                Nombre del modelo del cual se quieren extraer los resultados
    """

    # Leer estructura de datos en json
    with open(JSONPath) as json_file:
        DicFicCal = json.load(json_file)

    # Nombre del archivo de parámetros hidrológicos
    FileName = DicFicCal["dic_muestreador"]["PuntosControl"]

    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    # Ruta del archivo CSV de parámetros del Modelo SIGA-CAL
    PathFile = os.path.join(parent_dir,'Data','archivos_csv_Cal', 'muestreador', FileName)

    # Identificar puntos de control para la calibración del modelo hidrológico
    id = Catalog[NameSIGAModel].values == 1

    # Crear archivo
    ControlPoints = pd.DataFrame(data=np.vstack([Catalog.iloc[id, [1, 2, 0]].values, [X, Y, 'RiverMounthBasin']]),
                                 columns=['coord_x_or_aju', 'coord_y_or_aju', 'Nombre'])
    ControlPoints['Estadistico'] = 'valor'
    ControlPoints.iloc[-1, -1] = 'media'

    # Guardar archivo
    ControlPoints.to_csv(PathFile,index=False)


def Create_CalibrationFile(ProjectPath, JSONPath, NameSIGAFolder, NameCalibration):
    """
    Esta función crea el archivo txt de calibración de SIGA

    Parameters
        ProjectPath:
                Ruta del proyecto definido por el usuario
        JSONPath:
                Ruta donde se localiza el archivo JSON con la configuración de archivos
        NameSIGAFolder:
                Nombre de carpeta del proyecto
        NameCalibration:
                Nombre del archivo de calibración
    """

    # Ruta de la carpeta de entradas en el Modelo SIGA-CAL
    SIGAPath = os.path.join(ProjectPath,NameSIGAFolder,'entradas' + os.path.sep)

    # Ruta del archivo CSV de parámetros del Modelo SIGA-CAL

    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    CsvParamPath = os.path.join(parent_dir,'Data','archivos_csv_Cal','calibracion' + os.path.sep)

    # Leer estructura de datos en json
    with open(JSONPath) as json_file:
        DicFicCal = json.load(json_file)

    # Configurar clase de ficheros
    DicFic = Ficheros(ruta_cuenca_SIGA=SIGAPath,archivo_fact_cal=NameCalibration)

    # Crear archivo txt en el formato del Modelo SIGA-CAL
    DicFic.crear_archivo_factores_calibracion(DicFicCal['dic_factores_calibracion'],CsvParamPath)


def Create_SamplerFile(ProjectPath, JSONPath, NameSIGAFolder, NameSampler):
    """
    Esta función crea el archivo txt de muestreador de SIGA

    Parameters
        ProjectPath:
                Ruta del proyecto definido por el usuario
        JSONPath:
                Ruta donde se localiza el archivo JSON con la configuración de archivos
        NameSIGAFolder:
                Nombre de carpeta del proyecto
        NameSampler:
                Nombre del archivo de muestreo
    """

    # Ruta de la carpeta de entradas en el Modelo SIGA-CAL
    SIGAPath = os.path.join(ProjectPath,NameSIGAFolder,'entradas' + os.path.sep)


    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    # Ruta del archivo CSV de parámetros del Modelo SIGA-CAL
    CsvParamPath = os.path.join(parent_dir,'Data','archivos_csv_Cal','muestreador' + os.path.sep)

    # Leer estructura de datos en json
    with open(JSONPath) as json_file:
        DicFicCal = json.load(json_file)

    # Configurar clase de ficheros
    DicFic = Ficheros(ruta_cuenca_SIGA=SIGAPath,
                      archivo_muestreador=NameSampler)

    # Crear archivo txt en el formato del Modelo SIGA-CAL
    DicFic.crear_archivo_muestreador(DicFicCal['dic_muestreador'],CsvParamPath)


def ReadTimeSeries(FilePathTS):
    """
    Esta función leer una serie de tiempo en el formato de SIGA

    Parameters
        ProjectPath:
                Ruta del archivo de serie de tiempo
    """

    # Leer archivo
    with open(FilePathTS, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # Extraer número de registros y valor de datos faltantes
    num_records = int(lines[4].strip().split()[-1])
    missing_value = float(lines[7].strip().split()[-1])

    # Encontrar el índice donde comienzan los datos
    start_idx = 19  # Los datos comienzan en la línea 20, índice 19 (0 basado)

    # Leer los nombres de las columnas
    columns = lines[start_idx].strip().split()

    # Leer los datos
    data = []
    for line in lines[start_idx + 1:]:
        if line.strip():
            data.append(list(map(float, line.strip().split())))

    # Crear el DataFrame
    TS = pd.DataFrame(data, columns=columns)

    # Reemplazar los valores de datos faltantes
    TS.replace(missing_value, pd.NA, inplace=True)

    # Crear la columna Date
    TS['Date'] = pd.to_datetime(TS[['Año', 'Mes', 'Día']].astype(int).astype(str).agg('-'.join, axis=1),
                                format='%Y-%m-%d').dt.strftime('%d/%m/%Y')

    # Reordenar las columnas para que Date esté al inicio
    TS = TS[['Date'] + [col for col in TS.columns if col not in ['Date']]]

    #
    TS['Date'] = pd.to_datetime(TS['Date'].values, format='%d/%m/%Y')
    TS.set_index('Date', inplace=True)

    return TS


def ExeSIGAModel4(ProjectPath, NameSIGAFolder, NameExe, PathExeSIGA):
    """
    Esta función ejecuta el modelo SIGA

    Parameters
        ProjectPath:
                Ruta del proyecto definido por el usuario
        NameSIGAFolder:
                Nombre de la carpeta que corresponde a una configuración de SIGA-CAL
        NameExe:
                Nombre del escenario de simulación
        PathExeSIGA:
                Ruta del ejecutable de SIGA
    """

    # Nombre
    ErrorFile   = os.path.join(ProjectPath,NameSIGAFolder,'Log_Calibration.txt')

    # Linea de código para ejecutar SIGA-CAL
    ProcessSIGA = PathExeSIGA + ' -t "' + os.path.join(ProjectPath, NameSIGAFolder) + '" -n ' + NameExe + ' -a true 2>&1 | tee ' + ErrorFile + ' &'
    print('IMPORNTATE:',ProcessSIGA)
    # Ejecutar modelo
    result = subprocess.run(ProcessSIGA, shell=True)

    return result

def ExeSIGAModel(ProjectPath, NameSIGAFolder, NameExe, PathExeSIGA):
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

    # Definir la ruta del archivo de log
    ErrorFile = os.path.join(ProjectPath, NameSIGAFolder, 'Log_Calibration.txt').replace(os.sep, '/')

    # Comando con redirección de salida a archivo de log usando tee
    comando = f'. /etc/profile; {dir_siga_exe_cygwin} -t {path_eje} -n {NameExe} -a true 2>&1 | tee {ErrorFile}'

    # comando = f'. /etc/profile; {dir_siga_exe_cygwin} -t {path_eje} -n {NameExe} -a true'

    # Ejecutar el modelo usando Popen
    print('IMPORTANTE:', comando)
    p = Popen(['C:/cygwin64/bin/bash.exe', '-c', comando])
    # Esperar a que el proceso termine,shell= True
    stdout, stderr = p.communicate()

    timeStarted = time.time()                                 # Se guarda tiempo de inicio



    if p.returncode == 0:
        print("El modelo SIGA se ejecutó correctamente.")
    else:
        print(f"Error en la ejecución de SIGA. Código de salida: {p.returncode}")
        print(f"Salida del proceso: {stdout}")
        print(f"Errores del proceso: {stderr}")

    p.communicate()
    end_t = time.time()
    print('Simulación se demoró:')
    print(end_t- timeStarted)

    return p.returncode




