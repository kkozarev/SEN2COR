#!/usr/bin/env python
'''
Created on Feb 24, 2012
@author: umuellerwilm
'''
import fnmatch
import subprocess
import sys, os
import glob
import Image

from tables import *
from numpy import *
from tables.description import *
from distutils.dir_util import copy_tree, mkpath
from distutils.file_util import copy_file
from L2A_Config import L2A_Config
from L2A_Library import rectBivariateSpline, stdoutWrite, showImage
from lxml import etree, objectify
from L2A_XmlParser import L2A_XmlParser
from L2A_Borg import Borg

try:
    from osgeo import gdal,osr
    from osgeo.gdalconst import *
    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal,osr
    from gdalconst import *
# SIITBX-47: to suppress user warning due to the fact that 
# http://trac.osgeo.org/gdal/ticket/5480 is not implemented
# in the current openJPEG driver for windows used by ANACONDA:
gdal.PushErrorHandler('CPLQuietErrorHandler')

class Particle(IsDescription):
    bandName = StringCol(8)
    projectionRef = StringCol(512)
    geoTransformation = Int32Col(shape=6)
    rasterXSize = UInt16Col()
    rasterYSize = UInt16Col()
    rasterCount = UInt8Col()


class L2A_Tables(Borg):
    def __init__(self, config, L1C_TILE_ID):
        self.config = config
        L1C_UP_MASK = '*1C_*'
        if(fnmatch.fnmatch(config.workDir, L1C_UP_MASK) == False):
            self.config.tracer.fatal('identifier "*1C_*" is missing for the Level-1C_User_Product')
            self.config.exitError()
            return False

        AUX_DATA = '/AUX_DATA'
        IMG_DATA = '/IMG_DATA'
        QI_DATA = '/QI_DATA'
        GRANULE = '/GRANULE/'

        if os.name == 'posix':
            self._DEV0 = ' > /dev/null'
        else:
            self._DEV0 = ' > NUL'

        # Resolution:
        self._resolution = int(self.config.resolution)
        if(self._resolution == 10):
            self._bandIndex = [1,2,3,7]
            self._ownChannelIndex = [1,2,3,7]
            self._resChannelIndex = None
            self._nBands = 4
            bandDir = '/R10m'
        elif(self._resolution == 20):
            self._bandIndex = [1,2,3,4,5,6,8,11,12]
            self._ownChannelIndex = [4,5,6,8,11,12]
            self._resChannelIndex = [0,1,2,3,9,10]
            self._nBands = 9
            bandDir = '/R20m'
        elif(self._resolution == 60):
            self._bandIndex = [0,1,2,3,4,5,6,8,9,11,12]
            self._ownChannelIndex = [0,9,10]
            self._resChannelIndex = [1,2,3,4,5,6,8,11,12]
            self._nBands = 11
            bandDir = '/R60m'

        BANDS = bandDir
        #Creation_Date = self.config.creationDate
        # generate new Tile ID:
        L2A_TILE_ID = L1C_TILE_ID[:4] + 'USER' + L1C_TILE_ID[8:]
        #L2A_TILE_ID = L2A_TILE_ID[:25] + Creation_Date + L2A_TILE_ID[40:]
        L2A_TILE_ID = L2A_TILE_ID.replace('L1C_', 'L2A_')
        self.config.L2A_TILE_ID = L2A_TILE_ID
        L2A_TILE_ID_SHORT = '/' + L2A_TILE_ID[:55]
        L1C_TILE_ID = config.workDir + GRANULE + L1C_TILE_ID
        L2A_TILE_ID = config.L2A_UP_DIR + GRANULE + L2A_TILE_ID

        if(os.path.exists(L2A_TILE_ID) == False):
            os.mkdir(L2A_TILE_ID)
            copy_tree(L1C_TILE_ID + QI_DATA, L2A_TILE_ID + QI_DATA)

        self.config.logger.info('new working directory is: ' + L2A_TILE_ID)

        filelist = sorted(os.listdir(L1C_TILE_ID))
        found = False
        for filename in filelist:
            if(fnmatch.fnmatch(filename, L1C_UP_MASK) == True):
                found = True
                break
        if found == False:
            self.config.tracer.fatal('No metadata in tile')
            self.config.exitError()

        L1C_TILE_MTD_XML = L1C_TILE_ID + '/' + filename
        L2A_TILE_MTD_XML = filename
        L2A_TILE_MTD_XML = L2A_TILE_MTD_XML[:4] + 'USER' + L2A_TILE_MTD_XML[8:]
        L2A_TILE_MTD_XML = L2A_TILE_MTD_XML.replace('L1C_', 'L2A_')
        L2A_TILE_MTD_XML = L2A_TILE_ID + '/' + L2A_TILE_MTD_XML
        self.config.L1C_TILE_MTD_XML = L1C_TILE_MTD_XML
        self.config.L2A_TILE_MTD_XML = L2A_TILE_MTD_XML        

        if(self._resolution == 60):
            copy_file(L1C_TILE_MTD_XML, L2A_TILE_MTD_XML)
            xp = L2A_XmlParser(self.config, 'T2A')
            if(xp.convert() == False):
                self.tracer.fatal('error in converting tile metadata to level 2A')
                self.exitError()
            
            #update tile id in ds metadata file.
            xp = L2A_XmlParser(self.config, 'DS2A')
            ti = xp.getTree('Image_Data_Info', 'Tiles_Information')
            Tile = objectify.Element('Tile', tileId = self.config.L2A_TILE_ID)
            ti.Tile_List.append(Tile)
            xp.export()

        L1C_ImgDataDir = L1C_TILE_ID + IMG_DATA
        self._L2A_ImgDataDir = L2A_TILE_ID + IMG_DATA

        self._L1C_bandDir = L1C_ImgDataDir
        self._L2A_bandDir = self._L2A_ImgDataDir + BANDS

        if(os.path.exists(self._L2A_bandDir) == False):
            mkpath(self._L2A_bandDir)

        self._L1C_QualityMasksDir = L1C_TILE_ID + QI_DATA
        self._L2A_QualityDataDir = L2A_TILE_ID + QI_DATA
        self._L2A_AuxDataDir = L2A_TILE_ID + AUX_DATA

        if(os.path.exists(self._L2A_AuxDataDir) == False):
            mkpath(self._L2A_AuxDataDir)
            # copy configuration to AUX dir:
            dummy, basename = os.path.split(config._L2A_TILE_MTD_XML)
            fnAux = basename.replace('_MTD_', '_GIP_')
            target = self._L2A_AuxDataDir + '/' + fnAux
            copy_file(config.configFn, target)

        if(os.path.exists(self._L2A_QualityDataDir) == False):
            mkpath(self._L2A_QualityDataDir)
        
        '''
        if(self.config.traceLevel == 'DEBUG'):
            self._testdir = L2A_TILE_ID + '/TESTS_' + str(self._resolution) + '/'
            if(os.path.exists(self._testdir) == False):
                mkpath(self._testdir)
        '''
        #
        # the File structure:
        #-------------------------------------------------------
        pre = L2A_TILE_ID_SHORT[:9]
        post = L2A_TILE_ID_SHORT[13:]
        self._L2A_Tile_BND_File = self._L2A_bandDir + L2A_TILE_ID_SHORT + '_BXX_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_VIS_File = self._L2A_bandDir        + pre + '_VIS' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_AOT_File = self._L2A_bandDir        + pre + '_AOT' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_WVP_File = self._L2A_bandDir        + pre + '_WVP' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_SDW_File = self._L2A_bandDir        + pre + '_SDW' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_SLP_File = self._L2A_bandDir        + pre + '_SLP' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_ASP_File = self._L2A_bandDir        + pre + '_ASP' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_DEM_File = self._L2A_AuxDataDir     + pre + '_DEM' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_CLD_File = self._L2A_QualityDataDir + pre + '_CLD' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_SNW_File = self._L2A_QualityDataDir + pre + '_SNW' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_SCL_File = self._L2A_ImgDataDir     + pre + '_SCL' + post + '_' + str(self._resolution) + 'm.jp2'
        self._L2A_Tile_PVI_File = self._L2A_QualityDataDir + pre + '_PVI' + post + '.png'

        self._ImageDataBase = self._L2A_bandDir + '/.database.h5'
        self._TmpFile = self._L2A_bandDir + '/.tmpfile.tif'
        self._TmpDemFile = self._L2A_bandDir + '/.tmpdem.tif'
        self._acMode = False # default setting for scene classification

        # Geodata from image metadata:
        self._cornerCoordinates = None
        self._geoTransformation = None
        self._geoExtent = None
        self._projectionRef = None

        # Band Names:
        self._bandNames = ['B01','B02','B03','B04','B05','B06','B07','B08','B8A',\
                        'B09','B10','B11','B12','DEM','SCL','SNW','CLD','AOT',\
                        'WVP','VIS','SCM','PRV','ILU','SLP','ASP','HAZ','SDW',\
                        'DDV','HCW','ELE']

        # the mapping of the channels and bands
        self._B01 = 0
        self._B02 = 1
        self._B03 = 2
        self._B04 = 3
        self._B05 = 4
        self._B06 = 5
        self._B07 = 6
        self._B08 = 7
        self._B8A = 8
        self._B09 = 9
        self._B10 = 10
        self._B11 = 11
        self._B12 = 12
        self._DEM = 13
        self._SCL = 14
        self._SNW = 15
        self._CLD = 16
        self._AOT = 17
        self._WVP = 18
        self._VIS = 19
        self._SCM = 20
        self._PRV = 21
        self._ILU = 22
        self._SLP = 23
        self._ASP = 24
        self._HAZ = 25
        self._SDW = 26
        self._DDV = 27
        self._HCW = 28
        self._ELE = 29

        self.config.tracer.debug('Module L2A_Tables initialized with resolution %d' % self._resolution)

        return


    def get_ac_mode(self):
        return self._acMode


    def set_ac_mode(self, value):
        self._acMode = value


    def del_ac_mode(self):
        del self._acMode


    def get_corner_coordinates(self):
        return self._cornerCoordinates


    def get_geo_extent(self):
        return self._geoExtent


    def get_projection(self):
        return self._projection


    def set_corner_coordinates(self, value):
        self._cornerCoordinates = value


    def set_geo_extent(self, value):
        self._geoExtent = value


    def set_projection(self, value):
        self._projection = value


    def del_corner_coordinates(self):
        del self._cornerCoordinates


    def del_geo_extent(self):
        del self._geoExtent


    def del_projection(self):
        del self._projection


    def getBandNameFromIndex(self, index):
        return self._bandNames[index]


    def get_band_index(self):
        return self._bandIndex


    def get_n_bands(self):
        return self._nBands


    def get_db_name(self):
        return self._dbName


    def set_band_index(self, value):
        self._bandIndex = value


    def set_n_bands(self, value):
        self._nBands = value


    def set_db_name(self, value):
        self._dbName = value


    def del_band_index(self):
        del self._bandIndex


    def del_n_bands(self):
        del self._nBands


    def del_db_name(self):
        del self._dbName


        # end mapping of channels and bands

    def __exit__(self):
        sys.exit(-1)


    def __del__(self):
        self.config.tracer.debug('Module L2A_Tables deleted')


    def get_config(self):
        return self._config


    def set_config(self, value):
        self._config = value


    def del_config(self):
        del self._config


    def get_b01(self):
        return self._B01


    def get_b02(self):
        return self._B02


    def get_b03(self):
        return self._B03


    def get_b04(self):
        return self._B04


    def get_b05(self):
        return self._B05


    def get_b06(self):
        return self._B06


    def get_b07(self):
        return self._B07


    def get_b08(self):
        return self._B08


    def get_b8a(self):
        return self._B8A


    def get_b09(self):
        return self._B09


    def get_b10(self):
        return self._B10


    def get_b11(self):
        return self._B11


    def get_b12(self):
        return self._B12


    def get_dem(self):
        return self._DEM


    def get_scl(self):
        return self._SCL


    def get_qsn(self):
        return self._SNW


    def get_qcl(self):
        return self._CLD


    def get_aot(self):
        return self._AOT


    def get_wvp(self):
        return self._WVP


    def get_vis(self):
        return self._VIS


    def get_scm(self):
        return self._SCM


    def get_prv(self):
        return self._PRV


    def set_b01(self, value):
        self._B01 = value


    def set_b02(self, value):
        self._B02 = value


    def set_b03(self, value):
        self._B03 = value


    def set_b04(self, value):
        self._B04 = value


    def set_b05(self, value):
        self._B05 = value


    def set_b06(self, value):
        self._B06 = value


    def set_b07(self, value):
        self._B07 = value


    def set_b08(self, value):
        self._B08 = value


    def set_b8a(self, value):
        self._B8A = value


    def set_b09(self, value):
        self._B09 = value


    def set_b10(self, value):
        self._B10 = value


    def set_b11(self, value):
        self._B11 = value


    def set_b12(self, value):
        self._B12 = value


    def set_dem(self, value):
        self._DEM = value


    def set_scl(self, value):
        self._SCL = value


    def set_qsn(self, value):
        self._SNW = value


    def set_qcl(self, value):
        self._CLD = value


    def set_aot(self, value):
        self._AOT = value


    def set_wvp(self, value):
        self._WVP = value


    def set_vis(self, value):
        self._VIS = value


    def set_scm(self, value):
        self._SCM = value


    def set_prv(self, value):
        self._PRV = value


    def del_b01(self):
        del self._B01


    def del_b02(self):
        del self._B02


    def del_b03(self):
        del self._B03


    def del_b04(self):
        del self._B04


    def del_b05(self):
        del self._B05


    def del_b06(self):
        del self._B06


    def del_b07(self):
        del self._B07


    def del_b08(self):
        del self._B08


    def del_b8a(self):
        del self._B8A


    def del_b09(self):
        del self._B09


    def del_b10(self):
        del self._B10


    def del_b11(self):
        del self._B11


    def del_b12(self):
        del self._B12


    def del_dem(self):
        del self._DEM


    def del_scl(self):
        del self._SCL


    def del_qsn(self):
        del self._SNW


    def del_qcl(self):
        del self._CLD


    def del_aot(self):
        del self._AOT


    def del_wvp(self):
        del self._WV


    def del_vis(self):
        del self._VIS


    def del_scm(self):
        del self._SCM


    def del_prv(self):
        del self._PRV


    def get_ilu(self):
        return self._ILU


    def get_slp(self):
        return self._SLP


    def get_asp(self):
        return self._ASP


    def set_ilu(self, value):
        self._ILU = value


    def set_slp(self, value):
        self._SLP = value


    def set_asp(self, value):
        self._ASP = value


    def del_ilu(self):
        del self._ILU


    def del_slp(self):
        del self._SLP


    def del_asp(self):
        del self._ASP


    def get_sdw(self):
        return self._SDW


    def set_sdw(self, value):
        self._SDW = value


    def del_sdw(self):
        del self._SDW


    def get_ddv(self):
        return self._DDV


    def set_ddv(self, value):
        self._DDV = value


    def del_ddv(self):
        del self._DDV

    def get_hcw(self):
        return self._HCW


    def get_ele(self):
        return self._ELE


    def set_hcw(self, value):
        self._HCW = value


    def set_ele(self, value):
        self._ELE = value


    def del_hcw(self):
        del self._HCW


    def del_ele(self):
        del self._ELE


    B01 = property(get_b01, set_b01, del_b01, "B01's docstring")
    B02 = property(get_b02, set_b02, del_b02, "B02's docstring")
    B03 = property(get_b03, set_b03, del_b03, "B03's docstring")
    B04 = property(get_b04, set_b04, del_b04, "B04's docstring")
    B05 = property(get_b05, set_b05, del_b05, "B05's docstring")
    B06 = property(get_b06, set_b06, del_b06, "B06's docstring")
    B07 = property(get_b07, set_b07, del_b07, "B07's docstring")
    B08 = property(get_b08, set_b08, del_b08, "B08's docstring")
    B8A = property(get_b8a, set_b8a, del_b8a, "B8A's docstring")
    B09 = property(get_b09, set_b09, del_b09, "B09's docstring")
    B10 = property(get_b10, set_b10, del_b10, "B10's docstring")
    B11 = property(get_b11, set_b11, del_b11, "B11's docstring")
    B12 = property(get_b12, set_b12, del_b12, "B12's docstring")
    DEM = property(get_dem, set_dem, del_dem, "DEM's docstring")
    SCL = property(get_scl, set_scl, del_scl, "SCL's docstring")
    SNW = property(get_qsn, set_qsn, del_qsn, "SNW's docstring")
    CLD = property(get_qcl, set_qcl, del_qcl, "CLD's docstring")
    AOT = property(get_aot, set_aot, del_aot, "AOT's docstring")
    WVP = property(get_wvp, set_wvp, del_wvp, "WVP's docstring")
    VIS = property(get_vis, set_vis, del_vis, "VIS's docstring")
    SCM = property(get_scm, set_scm, del_scm, "SCM's docstring")
    PRV = property(get_prv, set_prv, del_prv, "PRV's docstring")
    ILU = property(get_ilu, set_ilu, del_ilu, "ILU's docstring")
    SLP = property(get_slp, set_slp, del_slp, "SLP's docstring")
    SDW = property(get_sdw, set_sdw, del_sdw, "SDW's docstring")
    ASP = property(get_asp, set_asp, del_asp, "ASP's docstring")
    DDV = property(get_ddv, set_ddv, del_ddv, "DDV's docstring")
    HCW = property(get_hcw, set_hcw, del_hcw, "HCW's docstring")
    ELE = property(get_ele, set_ele, del_ele, "ELE's docstring")
    config = property(get_config, set_config, del_config, "config's docstring")
    bandIndex = property(get_band_index, set_band_index, del_band_index, "bandIndex's docstring")
    nBands = property(get_n_bands, set_n_bands, del_n_bands, "nBands's docstring")
    dbName = property(get_db_name, set_db_name, del_db_name, "dbName's docstring")
    cornerCoordinates = property(get_corner_coordinates, set_corner_coordinates, del_corner_coordinates, "cornerCoordinates's docstring")
    geoExtent = property(get_geo_extent, set_geo_extent, del_geo_extent, "geoExtent's docstring")
    projection = property(get_projection, set_projection, del_projection, "projection's docstring")
    acMode = property(get_ac_mode, set_ac_mode, del_ac_mode, "acMode's docstring")

    def ReprojectCoords(self,coords,src_srs,tgt_srs):
        ''' Reproject a list of x,y coordinates.
            @type geom:     C{tuple/list}
            @param geom:    List of [[x,y],...[x,y]] coordinates
            @type src_srs:  C{osr.SpatialReference}
            @param src_srs: OSR SpatialReference object
            @type tgt_srs:  C{osr.SpatialReference}
            @param tgt_srs: OSR SpatialReference object
            @rtype:         C{tuple/list}
            @return:        List of transformed [[x,y],...[x,y]] coordinates
        '''
        trans_coords=[]
        transform = osr.CoordinateTransformation( src_srs, tgt_srs)
        for x,y in coords:
            x,y,z = transform.TransformPoint(x,y)
            trans_coords.append([x,y])
        return trans_coords


    def GetExtent(self, gt,cols,rows):
        ''' Return list of corner coordinates from a geotransform

            @type gt:   C{tuple/list}
            @param gt: geotransform
            @type cols:   C{int}
            @param cols: number of columns in the dataset
            @type rows:   C{int}
            @param rows: number of rows in the dataset
            @rtype:    C{[float,...,float]}
            @return:   coordinates of each corner
        '''
        ext=[]
        xarr=[0,cols]
        yarr=[0,rows]

        for px in xarr:
            for py in yarr:
                x=gt[0]+(px*gt[1])+(py*gt[2])
                y=gt[3]+(px*gt[4])+(py*gt[5])
                ext.append([x,y])
            yarr.reverse()
        return ext


    def J2kToH5(self):
        # convert JPEG-2000 input files to H5 file format
        # initialize H5 database for usage:
        workDir = self._L1C_bandDir
        os.chdir(workDir)
        tmpfile = self._TmpFile
        database = self._ImageDataBase
        command = 'gdal_translate '
        option1 = ' -ot UInt16 -outsize '
        rasterX = None

        if(os.path.isfile(database)):
            os.remove(database)
            self.config.tracer.info('Old database removed')
        self.config.timestamp('L2A_Tables: start import')
        dirs = sorted(os.listdir(workDir))
        # own channels first
        channels = self._ownChannelIndex
        for i in channels:
            for filename in dirs:
                bandName = self.getBandNameFromIndex(i)
                filemask = '*_L1C_*_%3s.jp2' % bandName
                if fnmatch.fnmatch(filename, filemask):
                    if(rasterX == None):
                        # get the target resolution and metadata for the resampled bands below:
                        indataset = gdal.Open(filename, GA_ReadOnly)
                        rasterX = indataset.RasterXSize
                        rasterY = indataset.RasterYSize
                        xp = L2A_XmlParser(self.config, 'T2A')           
                        tg = xp.getTree('Geometric_Info', 'Tile_Geocoding')
                        ulx = tg.Geoposition[0].ULX
                        uly = tg.Geoposition[0].ULY
                        res = float32(self.config.resolution)
                        self._geoTransformation = [ulx,res,0.0,uly,0.0,-res]
                        projectionRef = indataset.GetProjectionRef()
                        extent = self.GetExtent(self._geoTransformation, rasterX, rasterY)
                        self._cornerCoordinates = asarray(extent)
                        src_srs = osr.SpatialReference()
                        src_srs.ImportFromWkt(indataset.GetProjection())
                        tgt_srs = src_srs.CloneGeogCS()
                        if(tgt_srs != None):
                            geo_extent = self.ReprojectCoords(extent,src_srs,tgt_srs)
                            self._geoExtent = asarray(geo_extent)
                            self._projectionRef = src_srs
                        indataset = None
                        
                    callstr = command + filename + ' ' + tmpfile + self._DEV0
                    if(subprocess.call(callstr, shell=True) != 0):
                        self.config.tracer.fatal('shell execution error using gdal_translate')
                        self.config.exitError()
                        return False

                    self.ToH5_step2(i)
                    break

        # now the channels with different resolutions:
        channels = self._resChannelIndex
        # no resampling for 10m bands:
        if(channels != None):
            for i in channels:
                for filename in dirs:
                    bandName = self.getBandNameFromIndex(i)
                    filemask = '*_L1C_*_%3s.jp2' % bandName
                    if fnmatch.fnmatch(filename, filemask):
                        callstr = command + option1 + str(rasterX) + ' ' + str(rasterY) + ' ' + filename + ' ' + tmpfile + self._DEV0
                        if(subprocess.call(callstr, shell=True) != 0):
                            self.config.tracer.fatal('shell execution error using gdal_translate')
                            self.config.exitError()
                            return False
                        self.ToH5_step2(i)
                        break

        upsampling = False
        # 20m bands only: perform an up sampling of VIS from 60 m channels to 20
        if(self._resolution == 20):
            workDir = self._L2A_bandDir.replace('R20m', 'R60m')
            srcResolution = '_60m'
            channels = [17,19]
            upsampling = True
            self.config.tracer.info('perform up sampling of AOT and VIS from 60m channels to 20m')

        # 10m bands only: perform an up sampling of SCL, AOT, WVP, and VIS from 20 m channels to 10
        elif(self._resolution == 10):
            workDir = self._L2A_bandDir.replace('R10m', 'R20m')
            srcResolution = '_20m'
            channels = [14,17,18,19]
            upsampling = True
            self.config.tracer.info('perform up sampling of SCL, WVP, AOT and VIS from 20m channels to 10m')

        if(upsampling == True):
            os.chdir(workDir)
            dirs = sorted(os.listdir(workDir))
            for i in channels:
                for filename in dirs:
                    bandName = self.getBandNameFromIndex(i)
                    if (bandName == 'VIS') or (bandName == 'SCL') or bandName == 'AOT':
                        filemask = '*_' + bandName + '_L2A_*' + srcResolution + '.jp2'                    
                    else:
                        filemask = '*_L2A_*_' + bandName + srcResolution + '.jp2'
                    if fnmatch.fnmatch(filename, filemask):
                        if (bandName == 'VIS') or (bandName == 'SCL'):
                            # SCL and VIS are Byte Type:
                            option = ' -ot Byte -outsize '
                        else:
                            option = ' -ot UInt16 -outsize '
                        callstr = command + option + str(rasterX) + ' ' + str(rasterY) + ' ' + filename + ' ' + tmpfile + self._DEV0
                        if(subprocess.call(callstr, shell=True) != 0):
                            self.config.tracer.fatal('shell execution error using gdal_translate')
                            self.config.exitError()
                            return False
                        self.ToH5_step2(i)
                        break

        if(self.gdalDEM() == True):
            # generate hill shadow, slope and aspect using DEM:
            if(self.gdalDEM_Shade() == False):
                self.config.tracer.fatal('shell execution error generating DEM shadow')
                self.config.exitError()
                return False

            if(self.gdalDEM_Slope() == False):
                self.config.tracer.fatal('shell execution error generating DEM slope')
                self.config.exitError()
                return False

            if(self.gdalDEM_Aspect() == False):
                self.config.tracer.fatal('shell execution error generating DEM aspect')
                self.config.exitError()
                return False

            os.remove(self._TmpDemFile)

        self.config.timestamp('L2A_Tables: stop import')
        return True


    def get_utm_zone(self, longitude):
        return (int(1+(longitude+180.0)/6.0))


    def is_northern(self, latitude): # Determines if given latitude is a northern for UTM
            if (latitude < 0.0):
                return 0
            else:
                return 1


    def transform_utm_to_wgs84(self, easting, northing, zone1, zone2):
        utm_coordinate_system = osr.SpatialReference()
        utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set geographic coordinate system to handle lat/lon
        zone = zone1
        hemi = zone2
        # SIITBX-48
        if(hemi == 'N'): # N is Northern Hemisphere
            utm_coordinate_system.SetUTM(zone, 1)
        else:
            utm_coordinate_system.SetUTM(zone, 0)
        wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system
        # create transform component
        utm_to_wgs84_geo_transform = osr.CoordinateTransformation(utm_coordinate_system, wgs84_coordinate_system)
        return utm_to_wgs84_geo_transform.TransformPoint(easting, northing, 0) # returns lon, lat, altitude


    def transform_wgs84_to_utm(self, lon, lat):
        utm_coordinate_system = osr.SpatialReference()
        utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set geographic coordinate system to handle lat/lon
        utm_coordinate_system.SetUTM(self.get_utm_zone(lon), self.is_northern(lat))
        wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system
        # create transform component
        wgs84_to_utm_geo_transform = osr.CoordinateTransformation(wgs84_coordinate_system, utm_coordinate_system)
        return wgs84_to_utm_geo_transform.TransformPoint(lon, lat, 0) # returns easting, northing, altitude


    def gdalDEM(self):
        import urllib
        import zipfile

        demDir = self.config.getStr('Common_Section', 'DEM_Directory')
        if demDir == 'false':
            self.config.logger.info('DEM dir does not exist, flat surface is used')
            return False

        self.config.logger.info('Start DEM alignment for tile')
        workDir = self.config.home + '/' + demDir
        if(os.path.exists(workDir) == False):
            mkpath(workDir)
        os.chdir(workDir)

        xy = self.cornerCoordinates
        xp = L2A_XmlParser(self.config, 'T2A')
        tg = xp.getTree('Geometric_Info', 'Tile_Geocoding')
        hcsName = tg.HORIZONTAL_CS_NAME.text
        zone = hcsName.split()[4]
        zone1 = int(zone[:-1])
        zone2 = zone[-1:].upper()
        lonMin, latMin, dummy = self.transform_utm_to_wgs84(xy[1,0], xy[1,1], zone1, zone2)
        lonMax, latMax, dummy = self.transform_utm_to_wgs84(xy[3,0], xy[3,1], zone1, zone2)

        lonMinId = int((-180-lonMin)/(-360)*72)+1
        lonMaxId = int((-180-lonMax)/(-360)*72)+1
        latMinId = int((60-latMax)/(120)*24)+1 # this is inverted by intention
        latMaxId = int((60-latMin)/(120)*24)+1 # this is inverted by intention

        if(lonMinId < 1) or (lonMaxId > 72) or (latMinId < 1) or (latMaxId > 24):
            self.config.logger.info('no SRTM dataset available for this tile, flat surface is used')
            return False

        for i in range(lonMinId, lonMaxId+1):
            for j in range(latMinId, latMaxId+1):
                zipFn = 'srtm_{:0>2d}_{:0>2d}.zip'.format(i,j)
                if(os.path.isfile(zipFn) == False):
                    prefix = self.config.getStr('Common_Section', 'DEM_Reference')
                    stdoutWrite('Trying to retrieve DEM from URL %s this may take some time ...\n' % prefix)
                    self.config.logger.info('Trying to retrieve DEM from URL: %s', prefix)
                    url = prefix + zipFn
                    webFile = urllib.urlopen(url)
                    localFile = open(url.split('/')[-1], 'wb')
                    localFile.write(webFile.read())
                    webFile.close()
                    localFile.close()

                zipf = zipfile.ZipFile(zipFn, mode='r')
                for subfile in zipf.namelist():
                    zipf.extract(subfile)
                zipf.close()

        command = 'gdalwarp '
        arguments = '-ot Int16 '
        srtmf_src = 'srtm_tmp_src.tif'

        if(lonMinId == lonMaxId) & (latMinId == latMaxId):
            srtmf_src = 'srtm_{:0>2d}_{:0>2d}.tif'.format(i,j)
        else:
            for i in range(lonMinId, lonMaxId+1):
                for j in range(latMinId, latMaxId+1):
                    arguments += 'srtm_{:0>2d}_{:0>2d}.tif '.format(i,j)

            callstr = command + arguments + srtmf_src + self._DEV0
            if(subprocess.call(callstr, shell=True) != 0):
                self.config.tracer.fatal('shell execution error using gdalwarp')
                self.config.exitError()

                return False

        hcsCode = xmlParser.root.Geometric_Info.Tile_Geocoding.HORIZONTAL_CS_CODE
        t_srs = '-t_srs ' + hcsCode

        te = ' -te %f %f %f %f' % (xy[0,0], xy[2,1], xy[2,0], xy[0,1])
        tr = ' -tr %d %d' % (self.config.resolution, self.config.resolution)
        t_warp = te + tr + ' -r cubicspline '

        if(os.path.isfile(self._TmpFile) == True):
            os.remove(self._TmpFile)

        arguments = '-ot Int16 '
        callstr = command + arguments + t_srs + t_warp + srtmf_src + ' ' + self._TmpFile + self._DEV0
        if(subprocess.call(callstr, shell=True) != 0):
            self.config.tracer.fatal('shell execution error using gdalwarp')
            os.remove(srtmf_src)
            self.config.exitError()

        self.ToH5_step2(self.DEM)
        os.rename(self._TmpFile, self._TmpDemFile)
        os.remove(srtmf_src)
        return True


    def gdalDEM_Shade(self):
        altitude = 90.0 - float32(mean(self.config.solze_arr))
        azimuth = float32(mean(self.config.solaz_arr))
        command = 'gdaldem hillshade '
        options = '-compute_edges -az ' + str(azimuth) + ' -alt ' + str(altitude)

        callstr = command + options + ' ' + self._TmpDemFile + ' ' + self._TmpFile + self._DEV0
        if(subprocess.call(callstr, shell=True) != 0):
            self.config.tracer.fatal('shell execution error using gdaldem')
            self.config.exitError()
            return False

        self.ToH5_step2(self.SDW)
        return True


    def gdalDEM_Slope(self):
        command = 'gdaldem slope '
        srtmf_src = 'srtm_tmp_src.tif'
        options = '-compute_edges'

        callstr = command + options + ' ' + self._TmpDemFile + ' ' + self._TmpFile + self._DEV0
        if(subprocess.call(callstr, shell=True) != 0):
            self.config.tracer.fatal('shell execution error using gdaldem')
            self.config.exitError()
            return False

        self.ToH5_step2(self.SLP)
        return True


    def gdalDEM_Aspect(self):
        command = 'gdaldem aspect '
        options = '-compute_edges'

        callstr = command + options + ' ' + self._TmpDemFile + ' ' + self._TmpFile + self._DEV0
        if(subprocess.call(callstr, shell=True) != 0):
            self.config.tracer.fatal('shell execution error using gdaldem')
            self.config.exitError()
            return False

        self.ToH5_step2(self.ASP)
        return True


    def ToH5_step2(self, index, filename = None):
        bandName = self.getBandNameFromIndex(index)
        if(filename == None):
            filename = self._TmpFile

        indataset = gdal.Open(filename, GA_ReadOnly)
        projectionRef = indataset.GetProjectionRef()
        cols = indataset.RasterXSize
        rows = indataset.RasterYSize

        h5file = openFile(self._ImageDataBase, mode='a', title =  str(self._resolution) + 'm bands')

        if(h5file.__contains__('/tmp') == False):
            h5file.createGroup('/', 'tmp', 'temporary data')

        if(h5file.__contains__('/metadata')):
            group = h5file.root.metadata
            table = h5file.root.metadata.META
        else:
            group = h5file.createGroup('/', 'metadata', 'metadata information')
            table = h5file.createTable(group, 'META', Particle, "Meta Data")

        if(h5file.__contains__('/arrays')):
            garrays = h5file.root.arrays
        else:
            garrays = h5file.createGroup('/', 'arrays', 'band arrays')

        if(h5file.__contains__('/arrays/' + bandName)):
            self.config.tracer.error( bandName + ' already exists in database ' + self._ImageDataBase)
            h5file.close()
            indataset = None
            return
        else:
            inband = indataset.GetRasterBand(1)
            dtOut = self.setDataType(inband.DataType)
            filters = Filters(complib="zlib", complevel=1)
            eArray = h5file.createEArray(garrays, bandName, dtOut, (0,inband.XSize), bandName, filters=filters)

        for i in range(inband.YSize):
            scanline = inband.ReadAsArray(0, i, inband.XSize, 1, inband.XSize, 1)
            scanline = choose( equal( scanline, None), (scanline, None) )
            eArray.append(scanline)

        particle = table.row
        particle['bandName'] = bandName
        particle['geoTransformation'] = self._geoTransformation
        particle['projectionRef'] = projectionRef
        particle['rasterXSize'] = cols
        particle['rasterYSize'] = rows
        particle['rasterCount'] = 1
        particle.append()
        table.flush()
        h5file.close()
        indataset = None
        self.config.tracer.info(bandName + ' imported')
        self.config.timestamp('L2A_Tables: band ' + bandName + ' imported')
        return


    def H5ToJ2k(self, scOnly=False):
        workDir = self._L2A_bandDir
        if(os.path.exists(workDir) == False):
            self.config.tracer.fatal('missing directory %s:' % workDir)
            self.config.exitError()
            return False

        os.chdir(workDir)
        tmpfile = self._TmpFile
        database = self._ImageDataBase

        self.config.timestamp('L2A_Tables: start export')
        if(self._resolution == 10):
            bandIndex = [1,2,3,7,13,17,18]

        elif(self._resolution == 20):
            bandIndex = [1,2,3,4,5,6,8,11,12,13,14,15,16,17,18,19]

        elif(self._resolution == 60):
            bandIndex = [0,1,2,3,4,5,6,8,9,11,12,13,14,15,16,17,18,19]

        #prepare the xml export
        Granules = objectify.Element('Granules')
        Granules.attrib['granuleIdentifier'] = self.config.L2A_TILE_ID
        Granules.attrib['datastripIdentifier'] = self.config.L2A_DS_ID
        Granules.attrib['imageFormat'] = 'JPEG2000'
 
        for index in bandIndex:
            bandName = self.getBandNameFromIndex(index)
            filename = self._L2A_Tile_BND_File
            tmpfile = os.path.basename(tmpfile)
            option2 = ' UInt16 '   
            filename = filename.replace('BXX', bandName)
            if (bandName == 'VIS'):
                filename = self._L2A_Tile_VIS_File                
                option2 = ' Byte '
            elif (bandName == 'SNW') :
                filename = self._L2A_Tile_SNW_File
                option2 = ' Byte '
            elif(bandName == 'CLD'):
                filename = self._L2A_Tile_CLD_File
                option2 = ' Byte '
            elif(bandName == 'SCL'):
                option2 = ' Byte '
                filename = self._L2A_Tile_SCL_File
            elif(bandName == 'AOT'):
                filename = self._L2A_Tile_AOT_File
            elif(bandName == 'WVP'):
                filename = self._L2A_Tile_WVP_File
            elif(bandName == 'DEM'):
                filename = self._L2A_Tile_DEM_File
                demDir = self.config.getStr('Common_Section', 'DEM_Directory')
                if demDir == 'false':
                    continue

            self.H5ToJ2k_step2(index)

            if os.name == 'posix':
                callstr = 'gdal_translate -of JPEG2000 -ot' + option2 + tmpfile + ' ' + filename + self._DEV0
            else: # windows
                callstr = 'geojasper -f ' + tmpfile + ' -F ' + filename + ' -T jp2 > NUL'
            if(subprocess.call(callstr, shell=True) != 0):
                self.config.tracer.fatal('shell execution error using gdal_translate')
                self.config.exitError()
                return False

            self.config.tracer.debug('Band: ' + bandName + ' converted')
            filename = os.path.basename(filename.strip('.jp2'))
            if (bandName != 'VIS'):
                imageId = etree.Element('IMAGE_ID_2A')
                imageId.text = filename
                Granules.append(imageId)
        '''
        if(self.config.traceLevel == 'DEBUG'):
            if(scOnly == False):
                b = self.getBand(self.AOT, uint16)
                self.saveArray('S2L2APP_OUT_AOT', b)
            if(self.config.resolution > 10):
                b = self.getBand(self.SCL, uint8)
                self.saveArray('S2L2APP_OUT_SCL', b)
                b = self.getBand(self.SNW, uint8)
                self.saveArray('S2L2APP_OUT_SNW', b)
                b = self.getBand(self.CLD, uint8)
                self.saveArray('S2L2APP_OUT_CLD', b)
        '''
        xp = L2A_XmlParser(self.config, 'UP2A')
        pi = xp.getTree('General_Info', 'L2A_Product_Info')
        gl = objectify.Element('Granule_List')        
        gl.append(Granules)
        po = pi.L2A_Product_Organisation
        po.append(gl)
        xp.export()

        # update on tile level
        if(self._resolution == 60):
            xp = L2A_XmlParser(self.config, 'T2A')
            gi = xp.getRoot('General_Info')
            tiOld = xp.getTree('General_Info', 'TILE_ID_2A')
            tiNew = etree.Element('TILE_ID_2A')
            tiNew.text = self.config.L2A_TILE_ID
            gi.replace(tiOld, tiNew)
            dsOld = xp.getTree('General_Info', 'DATASTRIP_ID_2A')
            dsNew = etree.Element('DATASTRIP_ID_2A')
            dsNew.text = self.config.L2A_DS_ID
            gi.replace(dsOld, dsNew)
            pxlqi2a = objectify.Element('L2A_Pixel_Level_QI')
            fn = os.path.basename(self._L2A_Tile_CLD_File)
            fn = fn.replace('_60m.jp2', '')
            pxlqi2a.CLOUD_CONFIDENCE_MASK = fn
            fn = os.path.basename(self._L2A_Tile_SNW_File) 
            fn = fn.replace('_60m.jp2', '')            
            pxlqi2a.SNOW_ICE_CONFIDENCE_MASK = fn		
            qii = xp.getRoot('Quality_Indicators_Info')
            qii.insert(3, pxlqi2a)
            pviOld = xp.getTree('Quality_Indicators_Info', 'PVI_FILENAME')
            pviNew = etree.Element('PVI_FILENAME')
            if(scOnly == False):
                self.createPreviewImage()
                self.config.timestamp('L2A_Tables: preview image exported')
                fn = os.path.basename(self._L2A_Tile_PVI_File)
                fn = fn.replace('.jp2', '')  
                pviNew.text = fn
            else:
                pviNew.text = 'NO'
            qii.replace(pviOld, pviNew)
            xp.export()
        
        # cleanup:
        if(os.path.isfile(self._TmpFile)):
            os.remove(self._TmpFile)
        if(os.path.isfile(database)):
            os.remove(database)
        globdir = self.config.L2A_UP_DIR + '/GRANULE/' + self.config.L2A_TILE_ID + '/*/*.jp2.aux.xml'
        for filename in glob.glob(globdir):
            os.remove(filename)
        globdir = self.config.L2A_UP_DIR + '/GRANULE/' + self.config.L2A_TILE_ID + '/*/*/*.jp2.aux.xml'
        for filename in glob.glob(globdir):
            os.remove(filename)

        self.config.timestamp('L2A_Tables: stop export')
        return True


    def H5ToJ2k_step2(self, index):
        bandName = self.getBandNameFromIndex(index)
        tmpfile = self._TmpFile
        database = self._ImageDataBase

        h5file = openFile(database, mode='r')
        if(h5file.__contains__('/arrays/' +  bandName)):
            node = h5file.getNode('/arrays', bandName)
        else:
            h5file.close()
            return

        array = node.read()
        out_driver = gdal.GetDriverByName('GTiff')
        table = h5file.root.metadata.META
        for x in table.iterrows():
            if(x['bandName'] == bandName):
                rasterXSize = x['rasterXSize']
                rasterYSize = x['rasterYSize']
                rasterCount = x['rasterCount']
                break

        table.flush()
        outdataset = out_driver.Create(tmpfile, rasterXSize, rasterYSize, rasterCount, GDT_UInt16)
        gt = self._geoTransformation
        if(gt is not None):
            outdataset.SetGeoTransform(gt)
        pr = self._projectionRef
        if(pr is not None):
            outdataset.SetProjection(pr.ExportToWkt())
        outband = outdataset.GetRasterBand(1)
        outband.WriteArray(array)
        h5file.close()
        outdataset = None
        self.config.tracer.info('Band ' + bandName + ' exported')
        self.config.timestamp('L2A_Tables: band ' + bandName + ' exported')
        return


    def createPreviewImage(self):
        self.config.tracer.debug('Creating Preview Image')
        workDir = self._L2A_QualityDataDir
        os.chdir(workDir)

        if(self._resolution != 60):
            self.config.tracer.fatal('wrong resolution for this procedure, must be 60m')
            self.config.exitError()
            return False

        filename = self._L2A_Tile_PVI_File
        acMode = self.acMode
        self.acMode = False
        b = self.getBand(self.B02)
        g = self.getBand(self.B03)
        r = self.getBand(self.B04)
        self.acMode = acMode

        b = self.scaleImgArray(b)
        g = self.scaleImgArray(g)
        r = self.scaleImgArray(r)

        b = Image.fromarray(b)
        g = Image.fromarray(g)
        r = Image.fromarray(r)

        try:
            out = Image.merge('RGB', (r,g,b))
            out.save(filename, 'PNG')
            self.config.tracer.debug('Preview Image created')
            return True
        except:
            self.config.tracer.fatal('Preview Image creation failed')
            self.config.exitError()
            return False


    def scaleImgArray(self, arr):
        if(arr.ndim) != 2:
            self.config.tracer.fatal('must be a 2 dimensional array')
            self.config.exitError()
            return False

        arrclip = arr.copy()
        min = 0.0
        max = 0.125
        scale = 255.0
        arr = clip(arrclip, min, max)
        #SIITBX-50: wrong scale was used: 
        scaledArr = uint8(arr*scale/max)
        return scaledArr


    def testDb(self):
        result = False
        try:
            h5file = openFile(self._ImageDataBase, mode='r')
            h5file.getNode('/arrays', 'B02')
            status = 'Database ' + self._ImageDataBase + ' exists and can be used'
            result = True
        except:
            status = 'Database  ' + self._ImageDataBase + ' will be removed due to corruption'
            os.remove(self._ImageDataBase)
            result = False

        h5file.close()
        self.config.tracer.info(status)
        return result


    def hasBand(self, index):
        result = False
        bandName = self.getBandNameFromIndex(index)
        try:
            h5file = openFile(self._ImageDataBase, mode='r')
            h5file.getNode('/arrays', bandName)
            self.config.tracer.debug('Channel %s is present', self.getBandNameFromIndex(index))
            result = True
        except:
            self.config.tracer.debug('Channel %s is not available', self.getBandNameFromIndex(index))
            result = False

        h5file.close()
        return result


    def getBandSize(self, index):
        h5file = openFile(self._ImageDataBase, mode='r')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META
        rasterXSize = -1
        rasterYSize = -1
        rasterCount = 0
        for x in table.iterrows():
            if(x['bandName'] == bandName):
                nrows = x['rasterYSize']
                ncols = x['rasterXSize']
                count = x['rasterCount']
        table.flush()
        h5file.close()
        return(nrows, ncols, count)


    def getBand(self, index, dataType=uint16):
        bandName = self.getBandNameFromIndex(index)
        if (bandName == 'SCL') | (bandName == 'CLD') | \
           (bandName == 'SNW') | (bandName == 'PRV') | (bandName == 'VIS'):
            dataType = uint8
        
        # the output is context sensitive
        # it will return TOA_reflectance (0:1) if self.acMode = False (this is the scene classification mode)
        # it will return the unmodified value for all channels > 12, these are all generated products
        # it will return the radiance if self.acMode = True (this is the atmospheric correction mode)

        h5file = openFile(self._ImageDataBase, mode='r')
        node = h5file.getNode('/arrays', bandName)

        if (node.dtype != dataType):
            self.config.tracer.fatal('Wrong data type, must be: ' + str(node.dtype))
            h5file.close()
            self.config.exitError()

        nrows, ncols, count = self.getBandSize(index)
        if (count < 1):
            self.config.tracer.fatal('Insufficient band size: ' + count)
            h5file.close()
            self.config.exitError()

        array = node.read()
        h5file.close()

        if(index > 12):
            return array # no further modification
        elif(self.acMode == True):
            return self.TOA_refl2rad(index, array) # return radiance
        else: # return reflectance value:
            return array / self.config.dnScale # scaling from 0:1


    def TOA_refl2rad(self, index, array):
        # this converts TOA reflectance to radiance:
        (nrows, ncols, count) = self.getBandSize(index)
        DN = array
        if(self.config.resolution == 10):
            validBand = [None,0,1,2,None,None,None,3,None,None,None,None]
        else:
            validBand = [0,1,2,3,4,5,6,None,7,8,9,10,11]
        bandIndex = validBand[index]
        if(bandIndex == None):
            self.config.tracer.debug('Wrong band index %02d for selected resolution %02d', index, self.config.resolution)
            self.config.exitError()

        #if index > 7:  # band 8 is not included !!!
        #    index -= 1
        c0 = self.config.c0[bandIndex]
        c1 = self.config.c1[bandIndex]
        e0 = self.config.e0[bandIndex]
        dnScale = self.config.dnScale

        rtoa = float32(c0 + DN / dnScale)
        d2 = self.config.d2
        x = arange(nrows, dtype=float32) / (nrows-1) * self.config.solze_arr.shape[0]
        y = arange(ncols, dtype=float32) / (ncols-1) * self.config.solze_arr.shape[1]
        szi = rectBivariateSpline(x,y,self.config.solze_arr)
        rad_szi = radians(szi)
        sza = cos(rad_szi)
        L = float32(rtoa * e0 * sza / (pi * d2) / c1)
        return L


    def TOA_rad2refl(self, index):
        # this converts TOA radiance to reflectance:
        # a helper function for converting ATCOR Input to L2A_SceneClass
        # only for testing purposes - not needed in normal operation
        DN = self.getBand(index)
        nrows, ncols, count = self.getBandSize(index)
        d2 = self.config.d2
        #if(index > 7): #band 8 is not included !!!
        #    index = index-1
        if(self.config.resolution == 10):
            validBand = [False,0,1,2,False,False,False,3,False,False,False,False]
        else:
            validBand = [0,1,2,3,4,5,6,False,7,8,9,10,11]
        bandIndex = validBand[index]
        if(bandIndex == False):
            self.config.tracer.debug('Wrong band index %02d for selected resolution %02d', index, self.config.resolution)
            self.config.exitError()

        e0 = self.config.e0[bandIndex]
        c0 = self.config.c0[bandIndex]
        c1 = self.config.c1[bandIndex]
        x = arange(nrows, dtype=float32) / (nrows-1) * self.config.solze_arr.shape[0]
        y = arange(ncols, dtype=float32) / (ncols-1) * self.config.solze_arr.shape[1]
        szi = rectBivariateSpline(x,y,self.config.solze_arr)
        sza = cos(radians(szi))
        L = c0 + c1 * DN
        rtoa = pi * d2 * L / (e0 * sza) * c1
        return rtoa


    def getDataType(self, index):
        h5file = openFile(self._ImageDataBase, mode='r')
        bandName = self.getBandNameFromIndex(index)
        node = h5file.getNode('/arrays', bandName)
        dt = node.dtype
        h5file.close()
        return(dt)


    def setBand(self, index, array):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)

        if(h5file.__contains__('/arrays/' + bandName)):
            node = h5file.getNode('/arrays', bandName)
            node.remove()

        arr = h5file.root.arrays
        dtIn = self.setDataType(array.dtype)
        filters = Filters(complib="zlib", complevel=1)
        node = h5file.createEArray(arr, bandName, dtIn, (0,array.shape[1]), bandName, filters=filters)
        self.config.tracer.debug('Channel %02d %s added to table', index, self.getBandNameFromIndex(index))
        node.append(array)

        table = h5file.root.metadata.META
        update = False
        # if row exists, change it:
        for row in table.iterrows():
            if(row['bandName'] == bandName):
                row['rasterXSize'] = array.shape[1]
                row['rasterYSize'] = array.shape[0]
                row['rasterCount'] = 1
                row.update()
                update = True
        # else append it:
        if(update == False):
            row = table.row
            row['bandName'] = bandName
            row['rasterXSize'] = array.shape[1]
            row['rasterYSize'] = array.shape[0]
            row['rasterCount'] = 1
            row.append()

        table.flush()
        h5file.close()
        return


    def removeBand(self, index):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META
        if(h5file.__contains__('/arrays/' + bandName)):
            node = h5file.getNode('/arrays', bandName)
            node.remove()
            table.flush()
            self.config.tracer.debug('Channel %02d %s removed from table', index, self.getBandNameFromIndex(index))


        h5file.close()
        return


    def removeAllBands(self):
        h5file = openFile(self._ImageDataBase, mode='a')
        if(h5file.__contains__('/arrays/')):
            node = h5file.getNode('/', 'arrays')
            del node
            self.config.tracer.debug('All channels removed from table')

        h5file.close()
        self.removeAllTmpBands()
        return


    def getTmpBand(self, index, dataType=int16):
        h5file = openFile(self._ImageDataBase, mode='r')
        bandName = self.getBandNameFromIndex(index)
        node = h5file.getNode('/tmp', bandName)

        if (node.dtype != dataType):
            self.config.tracer.error('wrong data type, must be: ' + str(node.dtype))
            h5file.close()
            self.config.exitError()

        nrows, ncols, count = self.getBandSize(index)
        if (count < 1):
            self.config.tracer.error('insufficient band size: ' + count)
            h5file.close()
            self.config.exitError()

        array = node.read()
        h5file.close()
        return array


    def setTmpBand(self, index, array):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META

        if(h5file.__contains__('/tmp/' + bandName)):
            node = h5file.getNode('/tmp', bandName)
            node.remove()
            table.flush()
        tmp = h5file.root.tmp
        dtIn = self.setDataType(array.dtype)
        filters = Filters(complib="zlib", complevel=1)
        node = h5file.createEArray(tmp, bandName, dtIn, (0,array.shape[1]), bandName, filters=filters)
        self.config.tracer.debug('Temporary channel ' + str(index) + ' added to table')
        node.append(array)
        table.flush()
        h5file.close()
        return


    def removeTmpBand(self, index):
        h5file = openFile(self._ImageDataBase, mode='a')
        bandName = self.getBandNameFromIndex(index)
        table = h5file.root.metadata.META

        if(h5file.__contains__('/tmp/' + bandName)):
            node = h5file.getNode('/tmp', bandName)
            node.remove()
            table.flush()
            self.config.tracer.debug('Temporary channel ' + str(index) + ' removed from table')

        h5file.close()
        return


    def removeAllTmpBands(self):
        h5file = openFile(self._ImageDataBase, mode='a')
        table = h5file.root.metadata.META
        if(h5file.__contains__('/tmp/' +'')):
            node = h5file.getNode('/tmp')
            del node
            table.flush()
            self.config.tracer.debug('All temporary channels removed from table')
        h5file.close()
        return


    def setDataType(self, dtIn):
        if(dtIn == uint8):
            dtOut = UInt8Atom()
        elif(dtIn == uint16):
            dtOut = UInt16Atom()
        elif(dtIn == int16):
            dtOut = Int16Atom()
        elif(dtIn == uint32):
            dtOut = UInt32Atom()
        elif(dtIn == int32):
            dtOut = Int32Atom()
        elif(dtIn == float32):
            dtOut = Float32Atom()
        elif(dtIn == float64):
            dtOut = Float64Atom()
        elif(dtIn == GDT_Byte):
            dtOut = UInt8Atom()
        elif(dtIn == GDT_UInt16):
            dtOut = UInt16Atom()
        elif(dtIn == GDT_Int16):
            dtOut = Int16Atom()
        elif(dtIn == GDT_UInt32):
            dtOut = UInt32Atom()
        elif(dtIn == GDT_Int32):
            dtOut = Int32Atom()
        elif(dtIn == GDT_Float32):
            dtOut = Float32Atom()
        elif(dtIn == GDT_Float64):
            dtOut = Float64Atom()

        return dtOut


    def getArray(self, filename):
        cfg = L2A_Config()
        filename = self._testdir + filename + '.npy'
        if((os.path.isfile(filename)) == False):
            cfg.tracer.critical('File ' + filename + ' not present')
            return False

        return load(filename)


    def getDiffFromArrays(self, filename1, filename2):
        cfg = L2A_Config()
        filename1 = self._testdir + filename1 + '.npy'
        filename2 = self._testdir + filename2 + '.npy'
        if((os.path.isfile(filename1)) == False):
            cfg.tracer.critical('File ' + filename1 + ' not present')
            return False

        if((os.path.isfile(filename2)) == False):
            cfg.tracer.critical('File ' + filename2 + ' not present')
            return False

        arr1 = load(filename1)
        arr2 = load(filename2)
        return (arr1-arr2)


    def saveArray(self, filename, arr):
        cfg = L2A_Config()
        filename = self._testdir + filename + '.npy'
        save(filename, arr)

        if(os.path.exists(self._L2A_bandDir) == False):
                os.makedirs(self._L2A_bandDir)

        cfg.tracer.info('File ' + filename + ' saved to disk')
        return


    def sceneCouldHaveSnow(self):
        globalSnowMapFn = self.config.getStr('Scene_Classification', 'Snow_Map_Reference')
        globalSnowMapFn = self.config.libDir + globalSnowMapFn
        if((os.path.isfile(globalSnowMapFn)) == False):
            self.config.tracer.error('global snow map not present, snow detection will be performed')
            return True

        img = Image.open(globalSnowMapFn)
        globalSnowMap = array(img)
        xy = self.cornerCoordinates
        xp = L2A_XmlParser(self.config, 'T2A')
        tg = xp.getTree('Geometric_Info', 'Tile_Geocoding')
        hcsName = tg.HORIZONTAL_CS_NAME.text
        zone = hcsName.split()[4]
        zone1 = int(zone[:-1])
        zone2 = zone[-1:].upper()
        lonMin, latMin, dummy = self.transform_utm_to_wgs84(xy[1,0], xy[1,1], zone1, zone2)
        lonMax, latMax, dummy = self.transform_utm_to_wgs84(xy[3,0], xy[3,1], zone1, zone2)

        # Snow map should have a dimension of 7200 x 3600, 20 pixels per degree:
        xMin = int((lonMin + 180.0) * 20.0 + 0.5)
        xMax = int((lonMax + 180.0) * 20.0 + 0.5)
        yMin = 3600 - int((latMax + 90.0) * 20.0 + 0.5) # Inverted by intention
        yMax = 3600 - int((latMin + 90.0) * 20.0 + 0.5) # Inverted by intention
        aoi = globalSnowMap[yMin:yMax,xMin:xMax]
        if(aoi.max() > 0):
            return True
        return False

