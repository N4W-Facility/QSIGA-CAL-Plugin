# Common Functions - QSIGA_CAL

Esta carpeta contiene todas las funciones compartidas entre los diferentes subplugins de QSIGA_CAL.

## Archivos Compartidos

### Funciones Principales (Duplicadas previamente en múltiples Steps)

- **BuildingTopology.py** - Construcción de topología de cuencas
- **Functions_Raster.py** - Operaciones con rasters (versión más reciente de Step_05)
- **Morphology.py** - Análisis morfológico de cuencas
- **ShpToRaster.py** - Conversión de shapefile a raster
- **ficheros.py** - Manejo de archivos de configuración SIGA (versión más reciente de Step_03)
- **spotpy_setup_SIGA.py** - Configuración para calibración SIGA con SPOTPY (versión más reciente de Step_03)

### Funciones Específicas por Módulo

- **spotpy_setup_LAI.py** - Calibración de LAI (Step_02)
- **sensitivity_analysis_complete.py** - Análisis de sensibilidad (Step_03)
- **Fuctions_Priority.py** - Funciones de priorización (Step_04)
- **Raster.py** - Operaciones adicionales con rasters (Step_04)
- **Fuctions_Simulations.py** - Funciones para simulaciones (Step_05)

### Archivos de Configuración

- **Diccionario_ficheros.json** - Diccionario de archivos de configuración

## Uso

Todos los subplugins ahora importan desde esta carpeta central usando rutas relativas:

```python
# Desde los subplugins (Step_01, Step_02, etc.)
from ..common_functions.ficheros import Ficheros
from ..common_functions import spotpy_setup_SIGA as spotpy_setup
from ..common_functions.Functions_Raster import resample_to_match
```

## Versiones Prioritarias

Cuando había múltiples versiones de un archivo, se priorizaron las siguientes:

1. **spotpy_setup_SIGA.py** - Step_03 (Nov 20, 2024) - Versión más reciente y completa
2. **Functions_Raster.py** - Step_05 (Abril 2, 2025) - Versión más actualizada
3. **ficheros.py** - Step_03 (Nov 13, 2024) - Versión con últimas correcciones
4. **BuildingTopology.py** - Step_01 (Dic 11, 2024) - Versión más reciente

## Notas

- Los archivos duplicados en las carpetas `functions` de cada Step deben ser eliminados por el usuario una vez verificado el correcto funcionamiento.
- Esta estructura centralizada facilita el mantenimiento y evita inconsistencias entre versiones.
