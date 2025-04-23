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
 Este código contiene la configuración del modelo de crecimiento fenológico
 De acuerdo con la estructura definida por la librería SPOTPY. Con esta
 configuración se realiza la calibración del modelo.

 -------------------------------------------------------------------------
 Entradas
 -------------------------------------------------------------------------
 Para inicializar la configuración se deben ingresar los siguientes
 Param_Min:
    Array [1,12] con los valores mínimos que pueden tomar los parámetros
    del modelo
 Param_Max:
    Array [1,12] con los valores máximos que pueden tomar los parámetros
    del modelo
 M_SOS1:
    Número del mes de inicio de transición de lluvias
 M_SOS2
    Número del mes de finalización de transición de lluvias
 Month_Data
    Array [n,1] de los meses de cada día en el periodo de simulación
 temperatures_min:
    Array [n,1] de las temperaturas díarias mínimas del aire
 temperatures_max
    Array [n,1] de las temperaturas díarias máximas del aire
 evaluation_data
    Array [n,1] de los LAI díarios descargados de MODIS
 SMI_data
    Array [n,1] de las humedades del suelo díarias
 OmitObs
    Número de días que se omiten al inicio de las series para la evaluación
    de la función objetivo

 -------------------------------------------------------------------------
 Salidas
 -------------------------------------------------------------------------
 Best_Param:
    Mejor conjunto de parámetros para el modelo
"""

# ----------------------------------------------------------------------------------------------------------------------
# Importar librerías
# ----------------------------------------------------------------------------------------------------------------------
import spotpy
import numpy as np

# Definir el modelo de crecimiento de la vegetación
class VegetationGrowthModel:
    def __init__(self, Param_Min, Param_Max, M_SOS1, M_SOS2, Month_Data, temperatures_min,
                 temperatures_max, evaluation_data, SMI_data, OmitObs, FObjFun):
        """
        Inicializa la clase del modelo de crecimiento fenológico

        Parameters
            t_inf:
                    Temperatura mínima para la cual ocurre desarrollo vegetativo [°C]
            tf_sup:
                    Factor para determinar la temperatura máxima para la cual ocurre desarrollo vegetativo [°C]
                    Nota: Para que no se traslapen las temperaturas en la generación aleatoria de parámetros,
                    la Tmax se determina como una proporción de la Tmin
            lai_min:
                    Índice de área foliar mínimo que alcanza la cobertura [m2/m2]
            f_lai:
                    Índice de área foliar máximo que alcanza la cobertura [m2/m2]
                    Nota: Para que no se traslapen los LAI en la generación aleatoria de parámetros, la LAI_max se
                    determina como una proporción de la LAI_min
            l1:
                    Factor de forma 1 de la relación sigmoide de fracción máxima del LAI [Adimensional]
            l2:
                    Factor de forma 2 de la relación sigmoide de fracción máxima del LAI [Adimensional]
            l3:
                    Factor de forma 1 de la relación sigmoide de LAI [Adimensional]
            phu:
                    Potencial de unidades de calor que puede acumular la cobertura durante su fase de crecimiento [°C]
            fr_phu_sen:
                    Fracción del potencial de unidades de calor con el cual inicia la fase de senescencia [Adimensional]
            SMIc:
                    Umbral crítico del índice de humedad del suelo con el cual se detona el crecimiento [Ad]
            nSMI:
                    Número de días consecutivos en los que se debe superar el umbral crítico de humedad del suelo para
                    desatar una nueva fase de crecimiento
            gdd_c:
                    Condición inicial de unidades de calor acumuladas
            SOS1:
                    Mes de inicio de la transición a la temporada de lluvias
            SOS2:
                    Mes de finalización de la transición a la temporada de lluvias
            SMI:
                    Humedad del suelo
            OmitObs:
                    Cantidad de datos iniciales que se omiten para la evaluación de la función objetivo
            FObjFun:
                    Factor multiplicador de la función objetivo
        """

        # Generar parámetros utilizando una distribución uniforme
        self.params = [
            spotpy.parameter.Uniform('t_inf',       Param_Min[0],  Param_Max[0]),
            spotpy.parameter.Uniform('tf_sup',      Param_Min[1],  Param_Max[1]),
            spotpy.parameter.Uniform('lai_min',     Param_Min[2],  Param_Max[2]),
            spotpy.parameter.Uniform('f_lai',       Param_Min[3],  Param_Max[3]),
            spotpy.parameter.Uniform('l1',          Param_Min[4],  Param_Max[4]),
            spotpy.parameter.Uniform('l2',          Param_Min[5],  Param_Max[5]),
            spotpy.parameter.Uniform('l3',          Param_Min[6],  Param_Max[6]),
            spotpy.parameter.Uniform('phu',         Param_Min[7],  Param_Max[7]),
            spotpy.parameter.Uniform('fr_phu_sen',  Param_Min[8],  Param_Max[8]),
            spotpy.parameter.Uniform('SMIc',        Param_Min[9],  Param_Max[9]),
            spotpy.parameter.Uniform('nSMI',        Param_Min[10], Param_Max[10]),
            spotpy.parameter.Uniform('gdd_c',       Param_Min[11], Param_Max[11]),
        ]

        # Mes de inicio de transición de lluvias
        self.SOS1       = M_SOS1
        # Mes de finalización de la transición de lluvias
        self.SOS2       = M_SOS2
        # Meses de la serie de tiempo
        self.Month      = Month_Data
        # Temperatura mínima en el día (°C)
        self.Tmin       = temperatures_min
        # Temperatura máxima en el día (°C)
        self.Tmax       = temperatures_max
        # Serie de tiempo de datos LAI - Datos MODIS
        self.ObsData    = evaluation_data
        # Serie de humedad del suelo
        self.SMI        = SMI_data
        # Cantidad de datos iniciales que se omiten para la evaluación de la función objetivo
        self.OmitObs    = OmitObs
        # Factor multiplicador de la función objetivo
        self.FacObjF    = FObjFun

    def parameters(self):
        # Generación aleatoria de parámetros
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # Guarda los parámetros de la iteración i
        self.Params_i = vector

        # Ejecuta el modelo de crecimiento con los parámetros de la iteración i
        return PhenologicalModel(vector, self.SOS1, self.SOS2, self.Month, self.Tmin, self.Tmax, self.SMI)

    def evaluation(self):
        # Retorna los datos observados
        return self.ObsData

    def objectivefunction(self, simulation, evaluation):
        # Cantidad de datos iniciales que se omiten para la evaluación de la función objetivo
        n       = self.OmitObs

        # Serie de tiempo simulada
        Sim     = np.array(simulation[n+1:])

        # Serie de tiempo observada
        Obs     = np.array(evaluation[n:])

        # Metrica 1 - Coefficient of determination
        Part_1  = np.sum((Obs - np.nanmean(Obs))*(Sim - np.nanmean(Sim)))
        Part_2  = np.sum((Obs - np.nanmean(Obs))**2)*np.sum((Sim - np.nanmean(Sim))**2)
        RSqr    = (Part_1/np.sqrt(Part_2))**2

        # Metrica 2 - Nash–Sutcliffe coefficient
        Part_1  = np.sum((Obs - Sim)**2)
        Part_2  = np.sum((Obs - np.nanmean(Obs)) ** 2)
        NSE     = (Part_1 / Part_2)

        # Metrica 3 - Index of agreement
        Part_2  = np.sum((np.abs((Obs - np.nanmean(Obs))) + np.abs((Sim - np.nanmean(Obs))))**2)
        IoAd    = (Part_1 / Part_2)

        # Metrica 4 - Error absoluto normalizado
        Part_1  = np.sum(np.abs(Obs - Sim))
        Part_2  = np.sum(np.abs(Obs - np.nanmedian(Obs)))
        NAE     = (Part_1 / Part_2)

        # La función objetivo es una combinación lineal de las métricas anteriores donde los pesos que se asignan son
        # los siguientes:
        # RSqr:  0.2
        # NSE:   0.4
        # IoAd:  0.2
        # NAE:   0.2
        ObjFun  = self.FacObjF*((0.2*(1-RSqr)) + (0.4*NSE) + (0.2*IoAd) + (0.2*NAE))

        return ObjFun

def calculate_gdd(temperature_min, temperature_max, t_inf, t_sup):
    """
       Calcula las unidades térmicas acumuladas (GDD).

       Parameters:
       temperature_min (array-like): Temperaturas mínimas diarias.
       temperature_max (array-like): Temperaturas máximas diarias.
       t_inf (float): Temperatura mínima para el desarrollo vegetativo.
       t_sup (float): Temperatura máxima para el desarrollo vegetativo.

       Returns:
       array-like: Unidades térmicas acumuladas (GDD) diarias.
   """
    t_avg = (temperature_min + temperature_max) / 2
    td = np.minimum(np.maximum(t_avg - t_inf, 0), t_sup - t_inf)
    return td

def lai_growth_phase(fr_lai_max_j_prev,fr_phu_j, lai_max, lai_prev, l1, l2):
    """
    Calcula el incremento del LAI durante la fase de crecimiento.

    Parameters:
    fr_phu_j (float): Fracción de unidades térmicas acumuladas (PHU) en el día j.
    lai_max (float): Máximo índice de área foliar (LAI).
    lai_prev (float): Índice de área foliar (LAI) en el día anterior (j-1).
    l1 (float): Parámetro de forma de la relación sigmoide.
    l2 (float): Parámetro de forma de la relación sigmoide.

    Returns:
    float: Incremento del índice de área foliar (LAI) en el día j.
    """

    fr_lai_max_j = fr_phu_j / (fr_phu_j + np.exp(l1 - (l2 * fr_phu_j)))
    delta_lai = (fr_lai_max_j - fr_lai_max_j_prev) * lai_max * (1 - np.exp(5 * (lai_prev - lai_max)))
    return delta_lai, fr_lai_max_j

def lai_senescence_phase(fr_phu_j, lai_min, lai_opt, l3, fr_phu_sen):
    """
    Calcula el LAI durante la fase de senescencia.

    Parameters:
    fr_phu_j (float): Fracción de unidades térmicas acumuladas (PHU) en el día j.
    lai_min (float): Mínimo índice de área foliar (LAI).
    lai_opt (float): Índice de área foliar (LAI) óptimo.
    l3 (float): Parámetro de forma para la fase de senescencia.
    fr_phu_sen (float): Fracción de PHU para inicio de senescencia.

    Returns:
    float: Índice de área foliar (LAI) en el día j durante la fase de senescencia.
    """
    t_j = 12 * (((1 - fr_phu_j) / (1 - fr_phu_sen)) - 0.5)
    lai = lai_min + ((lai_opt - lai_min) / (1 + np.exp(-l3 * t_j)))
    return lai

def PhenologicalModel(vector, SOS1, SOS2, Month, temperatures_min, temperatures_max, SMI):
    '''
    En esta función se evalua el modelo de crecimiento fenológico con el juego de parámetros definidos por el algoritmo
    de optimización
    '''

    # Definir parámetros
    t_inf, ft_sup, lai_min, f_lai, l1, l2, l3, phu, fr_phu_sen, SMIc, nSMI, gdd_cumulative = vector

    # Temperatura superior
    t_sup   = t_inf*ft_sup

    # LAI Máximo
    lai_max = lai_min*f_lai

    # Cálculo de unidades térmicas acumuladas (GDD)
    gdd_daily = calculate_gdd(temperatures_min, temperatures_max, t_inf, t_sup)

    # Se define como condición inicial de LAi el LAI_min
    lai_values      = [lai_min]
    # Se define como condición de la fracción de LAI máximo como cero
    fr_lai_max_j    = 0
    # Se define como condición inicial de LAI el LAI_min
    lai             = lai_min
    # Se define como condición inicial de LAI_Opt el LAI_max
    lai_opt         = lai_max
    # Condición inicial de número de días acumulados para detonar el crecimineto
    nSMI_c          = 0
    # Estatus para activar el crecimiento al finalizar el mes de lluvias y no ha iniciado el crecimiento
    StatusEnd       = False

    for gdd, m, SMI_i in zip(gdd_daily, Month, SMI):
        # Unidades de calor acumuladas
        gdd_cumulative += gdd

        # Fracción del potencial de unidades de calor
        fr_phu_j = gdd_cumulative / phu

        # Modelo de crecimiento
        if fr_phu_j < fr_phu_sen:
            # Fase de crecimiento
            if lai < lai_max:
                delta_lai, fr_lai_max_j = lai_growth_phase(fr_lai_max_j, fr_phu_j, lai_max, lai_values[-1], l1, l2)
                lai = lai_values[-1] + delta_lai
            else:
                lai = lai_max
            # Asignar el LAI del paso de tiempo como LAI_opt
            lai_opt = lai
        else:
            # Fase de senescence
            lai = lai_senescence_phase(fr_phu_j, lai_min, lai_opt, l3, fr_phu_sen)

            # Si el día de simulación 𝑗 se encuentra entre el primer día del mes de inicio de la transición a la
            # temporada de lluvias (𝑆𝑂𝑆1) y el último día del mes de finalización de la transición a la temporada
            # de lluvias (𝑆𝑂𝑆2), y no se ha iniciado un nuevo ciclo de crecimiento, se evalúa el índice de humedad
            # del suelo 𝑆𝑀𝐼𝑗
            if (m >= SOS1)&(m <= SOS2):
                StatusEnd = True

                # Contador del número de días en el que se excede el umbral crítico de la humedad del suelo
                if SMI_i >= SMIc:
                    nSMI_c += 1

            # Si se cumple que 𝑆𝑀𝐼_𝑖 ≥ 𝑆𝑀𝐼𝑐 por un número 𝑛𝑆𝑀𝐼 consecutivo de días, donde 𝑆𝑀𝐼𝑐 es un umbral crítico,
            # se inicia una nueva fase de crecimiento para la vegetación, detonada por la humedad en el suelo.
            # Entonces, el valor de 𝐻𝑈𝑗 se iguala a cero y el 𝐼𝐴𝐹 se iguala a 𝐼𝐴𝐹𝑚í𝑛
            if nSMI_c >= nSMI:
                # Reiniciar contador de número de días acumulados
                nSMI_c          = 0
                # Reiniciar unidades de calor acumuladas
                gdd_cumulative  = 0
                # Reiniciar fracción máxima
                fr_lai_max_j    = 0
                # Definir LAI como LAI_min
                lai             = lai_min

            # Si el último día del mes 𝑆𝑂𝑆2 no se ha cumplido la condición anterior, se fuerza el inicio de la fase
            # de crecimiento. Entonces, el valor de 𝐻𝑈𝑗 se iguala a cero y el 𝐼𝐴𝐹 se iguala a 𝐼𝐴𝐹𝑚í𝑛.
            if (StatusEnd)&(m > SOS2):
                # Cambiar estatus para activar el crecimiento al finalizar el mes de lluvias y no ha iniciado el
                # crecimiento
                StatusEnd       = False
                # Reiniciar unidades de calor acumuladas
                gdd_cumulative  = 0
                # Reiniciar fracción máxima
                fr_lai_max_j    = 0
                # Definir LAI como LAI_min
                lai             = lai_min

            # Si el LAI es menor al LAI_min se iguala al LAI_min y se
            if lai <= lai_min:
                # Reiniciar unidades de calor acumuladas
                gdd_cumulative  = 0
                # Reiniciar fracción máxima
                fr_lai_max_j    = 0
                # Definir LAI como LAI_min
                lai             = lai_min

        lai_values.append(lai)

    return lai_values