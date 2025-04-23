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
This module contains classes to manipulate raster files

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
import matplotlib.pyplot as plt, matplotlib.cm as cm
import struct
from osgeo import gdal
from netCDF4 import Dataset
import os
import copy
import time
import rasterio
import warnings

class Raster():

    u"""
    Create an object to handle raster single-band data

    """
    _undefinedProjID = ('unknown', 'undefined', '')

    def __init__(self, source=None, nodt=None):
        
        u"""
        Constructor of a Raster object from a file

        Parameters
        ----------
        source : str or dict, optional
            Raster data source. It can be a raster file direction in the hard disk (string)
            or a dictionary of one dimensional arrays containing X coordinates, Y coordinates
            and values for the construction of a Raster object. In the last case, the value
            corresponding to each (x, y) pair will be stored in the corresponding raster cell
                'x': X coordinates array
                'y': Y coordinates array
                'values': values array
            If no data source is specified an empty Raster object is created
        nodt : int, optional
            Custom NoData value for the object
        """
        self.path = None
        self.proj = 'Unknown'

        if source is None:
            self.createRaster()
            print(f'An empty Raster object has been created successfully',
                end='\r', flush=True)
            
        elif isinstance(source, str):
            self.path = source
            if os.path.exists(self.path):
                self.loadRaster()
                if isinstance(nodt, int):  # Change the NoData value if specified
                    self.changeNoDataValue(nodt)
                    print(f'The NoData value of the Raster object has ben changed to {self.nodt}',
                        end='\r', flush=True)
                elif nodt is not None:
                    warnings.warn(f'The value {nodt} is not valid as a NoData marker, NoData value will not be changed')
                print(f'The file {self.path.split("/")[-1]} has been loaded successfully',
                    end='\r', flush=True)

            else:
                warnings.warn(f'The file {self.path.split("/")[-1]} could not be loaded, the path does not exist')
                return None

        elif isinstance(source, dict):
            if isinstance(nodt, int):  # Change the NoData value if specified
                self.setNoDataValue(nodt)
            else:
                self.setNoDataValue()
                if nodt is not None:
                    warnings.warn(f'The value {nodt} is not valid as a NoData marker, NoData value will be set to -9999')
                else:
                    print(f'NoData value will be set automatically to -9999',
                        end='\r', flush=True)
            self.createFromDict(source)
            print(f'A Raster object has been created successfully from a dictionary',
                end='\r', flush=True)
    
    def __str__(self):
        if self.path is not None:
            text =  f'\n------------------------------------------------------------------------' \
                    f'\nRaster object loaded from the file:' \
                    f'\n-----------------------------------' \
                    f'\n{self.path.split("/")[-1]}' \
                    f'\n------------------------------------------------------------------------' \
                    f'\n' \
                    f'\n------------------------------------------------------------------------' \
                    f'\nMetadata:' \
                    f'\n---------' \
                    f'\nNumber of rows: {self.nrows}' \
                    f'\nNumber of columns: {self.ncols}' \
                    f'\nLower left corner coordinates: ({self.xll}, {self.yll})' \
                    f'\nUpper right corner coordinates: ({self.xur}, {self.yur})' \
                    f'\nCell size (x, y): ({self.clsz[0]}, {self.clsz[1]})' \
                    f'\nNoData value: {self.nodt}' \
                    f'\nCoordinate reference system: '
        else:
            text =  f'\n------------------------------------------------------------------------' \
                    f'\nRaster object not loaded from any source file' \
                    f'\n------------------------------------------------------------------------' \
                    f'\nMetadata:' \
                    f'\n---------' \
                    f'\nNumber of rows: {self.nrows}' \
                    f'\nNumber of columns: {self.ncols}' \
                    f'\nLower left corner coordinates: ({self.xll}, {self.yll})' \
                    f'\nUpper right corner coordinates: ({self.xur}, {self.yur})' \
                    f'\nCell size (x, y): ({self.clsz[0]}, {self.clsz[1]})' \
                    f'\nNoData value: {self.nodt}' \
                    f'\nCoordinate reference system: '
        
        if self.proj.lower() in self._undefinedProjID:
            text = f'{text}Unknown'
        else:
            try:
                text = f"""{text}{self.proj.split(',')[-2].split('"')[1]}:""" \
                    f"""{self.proj.split(',')[-1].split('"')[1]}"""
            except:
                text = f"""{text}Coordinate system not recognized"""
        
        return f'{text}' \
               f'\n------------------------------------------------------------------------' \
               f'\n'
    
    def createRaster(self):
        u"""
        Create an empty Raster object
        """
        self.ncols = None
        self.nrows = None
        self.nclls = None
        self.xll = None
        self.yll = None
        self.xur = None
        self.yur = None
        self.clsz = None
        self.nodt = None
        self.mtrx = None
    
    def createFromDict(self, rasterDict):
        u"""
        Create a Raster object from a dictionary of one dimensional arrays with the
        same length

        Parameters
        ----------
        rasterDict : dict
            Dictionary of one dimensional arrays containing X coordinates, Y coordinates and values
            for the construction of a Raster object. The value corresponding to each (x, y) pair
            will be stored in the corresponding raster cell
                'x': X coordinates array
                'y': Y coordinates array
                'values': values array
        """
        try:
            X = rasterDict['x']
        except:
            raise Exception(f'Key "x" not found in the source dictionary')
        
        try:
            Y = rasterDict['y']
        except:
            raise Exception(f'Key "y" not found in the source dictionary')
        
        try:
            values = rasterDict['values']
        except:
            raise Exception(f'Key "values" not found in the source dictionary')

        xCoords = np.sort(np.unique(X))
        yCoords = np.sort(np.unique(Y))

        self.ncols = len(xCoords)
        self.nrows = len(yCoords)
        self.nclls = self.nrows * self.ncols
        clszx = np.mean(xCoords[1:] - xCoords[:-1])
        clszy = np.mean(yCoords[1:] - yCoords[:-1])
        self.clsz = (clszx, clszy)
        self.xll = np.min(xCoords) - self.clsz[0]/2
        self.yll = np.min(yCoords) - self.clsz[1]/2
        self.xur = np.max(xCoords) + self.clsz[0]/2
        self.yur = np.max(yCoords) + self.clsz[1]/2
        self.mtrx = np.ones((self.nrows, self.ncols)) * self.nodt

        for x, y, val in zip(X, Y, values):
            i, j = self.xy2ij(x, y)
            self.mtrx[i, j] = val
    
    def loadRaster(self):

        u"""
        Load a raster after identifiying the file extension
        """
        loadMethod = {'asc'  : self.loadAsciiRaster,
                      'tif'  : self.loadTifRaster,
                      'tiff' : self.loadTifRaster,
                      'sgabr': self.loadSgabrRaster
                      }
        
        # Identify the file extension and load the raster
        self._originFormat = self.path.lower().split('.')[-1]
        try:
            loadMethod[self._originFormat]()
        except KeyError:
            print(f'{self._originFormat} is not a recognized file extension. File not loaded.')
            return None
    
    def exportRaster(self, exportPath, prcsn=gdal.GDT_Float32):

        u"""
        Export a Raster object to a file

        Parameters
        ----------
        path : text
            Raster file direction in the hard disk
        """
        exportMethod = {'asc'  : self.exportAsciiRaster,
                        'tif'  : self.exportTifRaster,
                        'tiff' : self.exportTifRaster,
                        'sgabr': self.exportSgabrRaster
                        }

        # Identify the desired file extension and write the raster
        ext = exportPath.lower().split('.')[-1]
        try:
            exportMethod[ext](exportPath, prcsn)
        except TypeError:
            exportMethod[ext](exportPath)
        except KeyError:
            print(f'{ext} is not a recognized file extension. File not exported.')
            return None
        
        print(f'A Raster object has been exported successfully into the {exportPath} file.',
              end='\r', flush=True)

    def loadAsciiRaster(self):

        u"""
        Load an Arc/ASCII raster file from the hard disk
        """
        # Read the header information as text
        asc_file = open(self.path,'r')
        asc_line = asc_file.readline()
        temp,ncols = str.split(asc_line)
        asc_line = asc_file.readline()
        temp,nrows = str.split(asc_line)
        asc_line = asc_file.readline()
        temp,xll = str.split(asc_line)
        asc_line = asc_file.readline()
        temp,yll = str.split(asc_line)
        asc_line = asc_file.readline()
        temp,clszx = str.split(asc_line)
        asc_line = asc_file.readline()
        temp,clszy = str.split(asc_line)
        asc_line = asc_file.readline()
        temp,nodt = str.split(asc_line)
        asc_file.close()

        # Convert header information from text to numbers
        self.ncols = int(ncols)
        self.nrows = int(nrows)
        self.xll = float(xll)
        self.yll = float(yll)
        self.clsz = (float(clszx), float(clszy))
        self.nodt = float(nodt)

        # Read the matrix information
        self.mtrx = np.loadtxt(self.path, skiprows = 7)

        # Verification assertions
        assert self.nrows == np.size(self.mtrx,0)
        assert self.ncols == np.size(self.mtrx,1)

        # Calculate additional values
        self.nclls = self.nrows*self.ncols
        self.xur = self.xll+self.ncols*self.clsz[0]
        self.yur = self.yll+self.nrows*self.clsz[1]
    
    def loadTifRaster(self):

        u"""
        Load a GeoTiff raster file from the hard disk
        """
        # Read the georeferenciation properties and the data matrix
        tif = gdal.Open(self.path)
        self.mtrx = tif.GetRasterBand(1).ReadAsArray().astype(float)
        self.proj = tif.GetProjection()
        georef = tif.GetGeoTransform()
        if georef[1]<0:
            self.xll = georef[0]+georef[1]*self.mtrx.shape[1]#np.size(mtrx,1)
        else:
            self.xll = georef[0]
        if georef[5]<0:
            self.yll = georef[3]+georef[5]*self.mtrx.shape[0]#np.size(mtrx,0)
        else:
            self.yll = georef[3]
        clszx = np.abs(georef[1])
        clszy = np.abs(georef[5])
        self.clsz = (clszx, clszy)
        self.nodt = tif.GetRasterBand(1).GetNoDataValue()
        self.nrows = np.size(self.mtrx,0)
        self.ncols = np.size(self.mtrx,1)

        # Calculate additional values
        self.nclls = self.nrows*self.ncols
        self.xur = self.xll+self.ncols*self.clsz[0]
        self.yur = self.yll+self.nrows*self.clsz[1]
        tif = None
    
    def loadSgabrRaster(self):
    
        u"""
        Load a SGABR raster file from the hard disk
        """
        # Read the .sgabr file as a binary
        bin_file = open(self.path,'rb')
        bin_content = bin_file.read()
        bin_file.close()

        # Rastreo de la información entre los datos binarios
        (ncols, nrows, self.xll, self.yll, clsz, self.nodt) = \
            struct.unpack('ffffff',bin_content[0:24])
        self.ncols = int(ncols)
        self.nrows = int(nrows)
        self.nclls = self.ncols*self.nrows
        self.clsz = (clsz, clsz)
        self.xur = self.xll+self.ncols*self.clsz[0]
        self.yur = self.yll+self.nrows*self.clsz[1]
        self.mtrx = np.reshape(struct.unpack('f'*int(self.nclls),bin_content[24:]), \
                                            (self.nrows,self.ncols))
        del bin_content
    
    def exportAsciiRaster(self, exportPath):
    
        u"""
        Save a Raster object in the hard disk using the Arc/ASCII format

        Parameters
        ----------
        exportPath: text
            Direction in the hard disk in which the Arc/ASCII file will be written
        """
        # Build the header information string
        hdr='ncols        '+str(self.ncols)\
         +'\nnrows        '+str(self.nrows)\
         +'\nxllcorner    '+str(self.xll)\
         +'\nyllcorner    '+str(self.yll)\
         +'\ndx           '+str(self.clsz[0])\
         +'\ndy           '+str(self.clsz[1])\
         +'\nnodata_value '+str(self.nodt)

        # Save the data in the hard disk and finishes
        file = open(exportPath,'w')
        file.write(hdr)
        for row in self.mtrx:
            line = '\n'+' '.join(row.astype('str'))
            file.write(line)
        file.close()

    def exportTifRaster(self, exportPath, prcsn=gdal.GDT_Float32):

        u"""
        Save a Raster object in the hard disk using the GeoTiff format

        Parameters
        ----------
        exportPath: text
            Direction in the hard disk in which the GeoTiff file will be written
        prcsn: optional, integer
            Data precision. By default gdal.GDT_Float32, but any other data type
            can be specified
        """
        # Create the driver and the file
        driver = gdal.GetDriverByName("GTiff")
        tif = driver.Create(exportPath, self.ncols, self.nrows, 1, prcsn)

        # Set the metadata, band data and finishes
        tif.SetGeoTransform((self.xll, self.clsz[0], 0.0, self.yur, 0.0, -self.clsz[1]))
        tif.SetProjection(self.proj)
        band = tif.GetRasterBand(1)
        band.SetNoDataValue(self.nodt)
        band.WriteArray(self.mtrx, 0, 0)
        band.FlushCache()
        tif = None; band = None
    
    def exportSgabrRaster(self, exportPath):
    
        u"""
        Save a Raster object in the hard disk using the SIGA SGABR format

        Parameters
        ----------
        exportPath: text
            Direction in the hard disk in which the SGABR file will be written
        """
        # Set matrix values as float
        self.mtrx = self.mtrx.astype('f')

        # Pack the metadata into a binary variable
        bin_mtdt = struct.pack('f'*6, self.ncols, self.nrows, self.xll, self.yll,
                               0.5*(self.clsz[0]+self.clsz[1]), self.nodt)
        
        # Pack the data matrix into a binary variable
        bin_mtrx = self.mtrx.tobytes()

        # Write the metadata into the new SGABR file
        new_rst = open(exportPath,'wb')
        new_rst.write(bin_mtdt)
        new_rst.close()

        # Write the data matrix into the new SGABR file
        new_rst = open(exportPath,'ab')
        new_rst.write(bin_mtrx)
        new_rst.close()
        
        del bin_mtdt, bin_mtrx
    
    def changeNoDataValue(self, newNodt = -9999.0):
        
        u"""
        Change the no data value of a raster matrix to a custom value,
        both in the nodt and mtrx attributes

        Parameters
        ----------
        newNodt : number
            New no data value (default = -9999.0)
        """
        # Replace values equivalent to current NoData to the new NoData value
        self.replaceNoDataValue(newNodt)
        # Set the nodt attribute to the new NoData value
        self.setNoDataValue(newNodt)
    
    def replaceValue(self, oldValue, newValue):
        u"""
        Set the cells that have certain value to a new different
        value
        
        Parameters
        ----------
        oldValue: number
            Value to be replaced
        newValue: number
            New value to assign to the oldValue cells
        """
        self.mtrx[np.where(self.mtrx == oldValue)] = newValue

    def replaceNoDataValue(self, newValue):
        
        u"""
        Set the values equivalent to the current nodt value to a
        different value. When this method is used without using
        setNoDataValue, you can turn current NoData values into
        data values

        Parameters
        ----------
        newValue : number
            New value to replace current NoData values
        """
        self.replaceValue(self.nodt, newValue)
    
    def setNoDataValue(self, newNodt = -9999.0):
        
        u"""
        Set the nodt attribute of a Raster object to a new value.
        When this method is used without using replaceNoDataValue,
        all current NoData values in the raster matrix will become
        valid values automatically

        Parameters
        ----------
        newNodt : number
            New no data value (default = -9999.0)
        """
        self.nodt = newNodt
    
    def xy2ij(self, x, y):
        
        u"""
        Transform X Y coordinates to i j indices

        Parameters
        ----------
        x: float, numpy array or list of floats
            X coordinate
        y: float, numpy array or list of floats
            Y coordinate

        Returns
        -------
        i: int or numpy array of ints
            Row index
        j: int or numpy array of ints
            Colum index
        """
        x = np.array(x)
        y = np.array(y)
        i = (self.nrows - np.ceil((y - self.yll) / self.clsz[1])).astype(int)
        j = (np.floor((x - self.xll) / self.clsz[0])).astype(int)
        return i, j
    
    def getValueByIndex(self, i, j):
        
        u"""
        Get a raster cell value given a row (i) column (j) index

        Parameters
        ----------
        i: int, numpy array or list of ints
            Row index
        j: int, numpy array or list of ints
            Colum index

        Returns
        -------
        value: float or numpy array of floats
            Raster cell value
        """
        try:
            type(i) == type(j)
        except:
            raise Exception(f"Type mismatch between the specified indexes " \
                            f"({type(i)} versus {type(j)})")
        
        try:
            iter(i)
            iter(j)
            if len(i) != len(j):
                shortest = np.min((len(i), len(j)))
                warnings.warn(f"Length mismatch between the specified indexes " \
                              f"({len(i)} versus {len(j)}), the output will have " \
                              f"the length of the shortest one ({shortest})")
            if np.min(i) >= 0 and np.min(j) >= 0 and np.max(i) < self.nrows and np.max(j) < self.ncols:
                value = self.mtrx[i, j]
            else:
                warnings.warn(f"Some of the specified indexes are out of bounds, " \
                              f"the value for these points will be set to np.nan")
                value = np.zeros(len(i)) * np.nan
                for n, (iInd, jInd) in enumerate(zip(i, j)):
                    if iInd >= 0 and jInd >= 0 and iInd < self.nrows and jInd < self.ncols:
                        value[n] = self.mtrx[iInd, jInd]
            value[value == self.nodt] = np.nan
        except:
            if i >= 0 and j >= 0 and i < self.nrows and j < self.ncols:
                value = self.mtrx[i, j]
                if value == self.nodt: value = np.nan
            else:
                warnings.warn(f"The specified index pair ({i}, {j}) is out of bounds, " \
                              f"the value for this point will be set to np.nan")
                value = np.nan
        
        return value
    
    def getValueByCoordinates(self, x, y):
        
        u"""
        Get a raster cell value given X Y coordinates

        Parameters
        ----------
        x: float, numpy array or list of floats
            X coordinate
        y: float, numpy array or list of floats
            Y coordinate

        Returns
        -------
        value: float or numpy array of floats
            Raster cell value
        """
        i, j = self.xy2ij(x, y)
        value = self.getValueByIndex(i, j)
        return value
    
    def stats(self):
        u"""
        Get a summary of the main stats of a Raster matrix

        Returns
        -------
        statsSummary: dict
            Raster stats stored in a dictionary
        """
        statsSummary = {}

        statsSummary['max'] = self.max()
        statsSummary['min'] = self.min()
        statsSummary['mean'] = self.mean()
        statsSummary['sum'] = self.sum()

        if self.path is not None:
            text =  f'\n------------------------------------------------------------------------' \
                    f'\nRaster object stats summary:' \
                    f'\n-----------------------------------' \
                    f'\nFile: {self.path.split("/")[-1]}' \
                    f'\n------------------------------------------------------------------------' \
                    f'\n' \
                    f'\n------------------------------------------------------------------------' \
                    f'\nMaximum: {statsSummary["max"]}' \
                    f'\nMinimum: {statsSummary["min"]}' \
                    f'\nMean: {statsSummary["mean"]}' \
                    f'\nSum: {statsSummary["sum"]}'
        else:
            text =  f'\n------------------------------------------------------------------------' \
                    f'\nRaster object stats summary:' \
                    f'\n-----------------------------------' \
                    f'\nFile: None' \
                    f'\n------------------------------------------------------------------------' \
                    f'\n' \
                    f'\n------------------------------------------------------------------------' \
                    f'\nMaximum: {statsSummary["max"]}' \
                    f'\nMinimum: {statsSummary["min"]}' \
                    f'\nMean: {statsSummary["mean"]}' \
                    f'\nSum: {statsSummary["sum"]}'
        
        print(text)

        return statsSummary
    
    def max(self):
        
        u"""
        Get the maximum value of the raster matrix

        Returns
        -------
        max: float
            Raster max value
        """
        values = self.mtrx[self.mtrx != self.nodt]
        max = np.max(values)
        return(max)
    
    def min(self):
        
        u"""
        Get the minimum value of the raster matrix

        Returns
        -------
        min: float
            Raster min value
        """
        values = self.mtrx[self.mtrx != self.nodt]
        min = np.min(values)
        return(min)
    
    def mean(self):
        
        u"""
        Get the average value of the raster matrix

        Returns
        -------
        mean: float
            Raster mean value
        """
        values = self.mtrx[self.mtrx != self.nodt]
        mean = np.mean(values)
        return(mean)
    
    def sum(self):
        
        u"""
        Get the sum of the raster matrix values

        Returns
        -------
        sum: float
            Raster sum value
        """
        values = self.mtrx[self.mtrx != self.nodt]
        sum = np.sum(values)
        return(sum)

    def plot(self):
        
        u"""
        Generate a simple plot representing the raster matrix
        """
        fig, ax = plt.subplots(figsize = (7,7))
        plotMtrx = self.mtrx.copy().astype(float)
        plotMtrx[plotMtrx == self.nodt] = np.nan
        plt.imshow(plotMtrx)
        plt.colorbar()
    
    def copy(self):

        u"""
        Generate a separate copy of the raster object
        """
        return copy.deepcopy(self)
    
    def rasterGrid(self, missing=False):

        X = np.linspace(self.xll+0.5*self.clsz[0], self.xur-0.5*self.clsz[0], self.ncols)
        Y = np.linspace(self.yur-0.5*self.clsz[1], self.yll+0.5*self.clsz[1], self.nrows)
        X, Y = np.meshgrid(X, Y)
        X = np.resize(X, self.nclls)
        Y = np.resize(Y, self.nclls)
        if not missing:
            V = np.resize(self.mtrx, self.nclls)
            available = V!=self.nodt
            X = X[available]
            Y = Y[available]
        return np.vstack((X,Y)).T
    
    def lookupRaster(self, lookup):
        """
        Replaces the values of a categorical raster into a new
        categorical raster, based on a user defined lookup table.

        Parameters:
        -----------
        lookup: dictionary
            Lookup table as a dictionary. Keys are the raster
            existent categories, values are the desired
            categories
        """
        
        rst2 = self.copy()
        rows = rst2.nrows
        cols = rst2.ncols
        
        ts = time.time()
        
        for key in lookup.keys():
            rst2.replaceValue(key, lookup[key])
        
        print(f"Tiempo de ejecución: {str(time.time()-ts)} segundos")
        return rst2

#------------------------------------------------------------------------

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
    georef=tif.GetGeoTransform()
    with rasterio.open(path, 'r') as src:
        mtrx = src.read()[0]
    src =None
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
    tif = None
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
    rst['mtrx'][0,:]=rst['nodt']
    rst['mtrx'][rst['nrows']-1,:]=rst['nodt']
    rst['mtrx'][:,0]=rst['nodt']
    rst['mtrx'][:,rst['ncols']-1]=rst['nodt']

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

    Parameters
    ----------
    dir_rst: Flow direction map
    FlowVector :  numpy array
                  [0] Id celda
                  [1] Id celda destino
                  [2] Area acumulada
                  [3] Long No acum
                  [4] Cota del DEM
                  [5] X: coordenada
                  [6] Y: coordenada
                  [7] i
                  [8] j
                  [9] pendiente
                  [10] tipo celda [0: ladera; 1:cauce; 2:embalse]
                  [11] Nodos [0:es cauce; 1:cabecera; 2:hidrologico; 3:topografico;
                            4:antropico, pie de presa]
                  [12] Cota corregida cauce
                  [13] Numero de tramo al que corresponde cada celda que
                  sea cauce
                  [14] Clasificacion de Beechie para celda tipo cauce
                  [15] Clasificacion de Flores para celda tipo cauce

    topic: int
            Contenido del mapa que se va a hacer asi:
            10:Dibuja red de drenaje
            11: Grafica los nodos en la red de drenaje y la red
            13: Grafica los tramos de la red
            14: Grafica la clasificacion segun flores
            15: Grafica la clasificacion segun Beechie

    path: string
            Donde se va a almacenar el mapa creado

    Returns
    -------
    Map : ASC
            Devuelve un mapa en formato ASC

    """
    #Para guardar en raster
    mapT1=np.zeros((dir_rst['nrows'],dir_rst['ncols']))-9999.0
    isnetwork=np.where(np.logical_or(FlowVectorA.loc[:,'Tipo']==1.0, FlowVectorA.loc[:,'Tipo']>3.0))[0]
    i=FlowVectorA.loc[isnetwork,'i'].astype(int)
    j=FlowVectorA.loc[isnetwork,'j'].astype(int)
    mapT1[i,j]=FlowVectorA.loc[isnetwork,topic]
    map_rst_1={'ncols':dir_rst['ncols'],'nrows':dir_rst['nrows'],'nclls':dir_rst['nclls'],'xll':dir_rst['xll'],'yll':dir_rst['yll'],'xur':dir_rst['xur'],'yur':dir_rst['yur'],'clsz':dir_rst['clsz'],'nodt':-9999,'mtrx':mapT1,'prj':''}
    if path != '':
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
    if path_rst1 != '':
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


def RasterGrid(rst,missing=False):
    X=np.linspace(rst['xll']+0.5*rst['clsz'],rst['xur']-0.5*rst['clsz'],rst['ncols'])
    Y=np.linspace(rst['yur']-0.5*rst['clsz'],rst['yll']+0.5*rst['clsz'],rst['nrows'])
    X,Y=np.meshgrid(X,Y)
    X=np.resize(X,rst['nclls'])
    Y=np.resize(Y,rst['nclls'])
    if not missing:
        V=np.resize(rst['mtrx'],rst['nclls'])
        available=V!=rst['nodt']
        X=X[available]
        Y=Y[available]
    return np.vstack((X,Y)).T

def WarpRaster(path_in,path_out,crs_in,crs_out):
    a=gdal.Warp(path_out,path_in,srcSRS=crs_in,dstSRS=crs_out,dstNodata=-9999.0)
    print(path_in+' successfully warped')
    return


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