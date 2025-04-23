
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
import rasterio
import numpy as np
from pysheds.grid import Grid




def ExeSIGAModel_General(ProjectPath, NameSIGAFolder, NameExe, PathExeSIGA,escenario):
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


    return p.returncode

