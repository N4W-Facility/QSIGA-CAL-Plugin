import os, sys
import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_bounds
import geopandas as gpd


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

    if PathRaster != '':
        # Convertir shp a raster
        rasterize_rivernet = rasterize(
            [(shape, 1) for shape in df['geometry']],
            out_shape=shape,
            transform=transform,
            fill=0,
            all_touched=True,
            dtype=rasterio.uint8)
    else:
        # Convertir shp a raster
        rasterize_rivernet = rasterize(
            [(shape, 1) for shape in df['geometry', Value]],
            out_shape=shape,
            transform=transform,
            fill=0,
            all_touched=True,
            dtype=rasterio.uint8)

    # Guardar raster
    with rasterio.open(
        PathOutput, 'w',
        driver='GTiff',
        dtype=rasterio.uint8,
        count=1,
        width=shape[1],
        height=shape[0],
        transform=transform
    ) as dst:
        dst.write(rasterize_rivernet, indexes=1)