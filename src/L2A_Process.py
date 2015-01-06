#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from numpy import *
from tables import *
import sys
import os
import fnmatch
from time import time

from L2A_Config import L2A_Config
from L2A_Tables import L2A_Tables
from L2A_SceneClass import L2A_SceneClass
from L2A_AtmCorr import L2A_AtmCorr
from L2A_XmlParser import L2A_XmlParser
from L2A_Library import stdoutWrite, stderrWrite
from lxml import etree, objectify


class L2A_Process(object):
    def __init__(self, workdir):
        self._config = L2A_Config(workdir)
        self._tables = False
        self._processed60 = False
        self._processed20 = False
        self._processed10 = False
        self._scOnly = False


    def get_sc_only(self):
        return self._scOnly


    def set_sc_only(self, value):
        self._scOnly = value


    def del_sc_only(self):
        del self._scOnly


    def get_tables(self):
        return self._tables


    def set_tables(self, value):
        self._tables = value


    def del_tables(self):
        del self._tables


    def get_config(self):
        return self._config


    def set_config(self, value):
        self._config = value


    def del_config(self):
        del self._config


    def __exit__(self):
            sys.exit(-1)


    config = property(get_config, set_config, del_config, "config's docstring")
    tables = property(get_tables, set_tables, del_tables, "tables's docstring")
    scOnly = property(get_sc_only, set_sc_only, del_sc_only, "scOnly's docstring")


    def selectAndProcess(self, tile):
        if(self.config.resolution == 10):
            self.config.tracer.info('selected resolution is 10m')
            self.config.logger.info('selected resolution is 10m')
            if(self._processed20 == False):
                self.config.resolution = 20
                stdoutWrite('20m resolution must be processed first ...\n')
                self.config.tracer.info('20m resolution must be processed first')
                self.config.logger.info('20m resolution must be processed first')
                self.selectAndProcess(tile)

            self.config.resolution = 10
            self.config.readPreferences()
            self.tables = L2A_Tables(self.config, tile)
            self._processed10 = self.process()
            if(self._processed10 == False):
                return False

        elif(self.config.resolution == 20):
            self.config.tracer.info('selected resolution is 20m')
            self.config.logger.info('selected resolution is 20m')
            if(self._processed60 == False):
                self.config.resolution = 60
                stdoutWrite('60m resolution must be processed first ...\n')
                self.config.tracer.info('60m resolution must be processed first')
                self.config.logger.info('60m resolution must be processed first')
                self.selectAndProcess(tile)

            self.config.resolution = 20
            self.config.readPreferences()
            self.tables = L2A_Tables(self.config, tile)
            self._processed20 = self.process()
            if(self._processed20 == False):
                return False

        elif(self.config.resolution == 60):
            self.config.tracer.info('selected resolution is 60m')
            self.config.logger.info('selected resolution is 60m')
            self.config.readPreferences()
            self.tables = L2A_Tables(self.config, tile)
            self.config.readTileMetadata()
            self._processed60 = self.process()
            if(self._processed60 == False):
                return False
        else:
            self.config.logger.debug('wrong resolution for processing configured: ', str(self.config.resolution))
            return False


    def process(self):
        astr = 'L2A_Process: processing with resolution ' + str(self.config.resolution) + ' m'
        self.config.timestamp(astr)
        self.config.timestamp('L2A_Process: start of pre processing')
        if(self.preprocess() == False):
            return False

        if(self.config.resolution > 10):
            self.config.timestamp('L2A_Process: start of Scene Classification')
            sc = L2A_SceneClass(self.config, self.tables)
            self.config.tracer.info('Performing Scene Classification with resolution %d m', self.config.resolution)
            self.config.logger.info('Performing Scene Classification with resolution %d m', self.config.resolution)
            if(sc.process() == False):
                return False

        if(self.scOnly == False):
            self.config.timestamp('L2A_Process: start of Atmospheric Correction')
            self.config.tracer.info('Performing Atmospheric Correction with resolution %d m', self.config.resolution)
            self.config.logger.info('Performing Atmospheric Correction with resolution %d m', self.config.resolution)
            ac = L2A_AtmCorr(self.config, self.tables)
            if(ac.process() == False):
                return False

        self.config.timestamp('L2A_Process: start of post processing')
        if(self.postprocess() == False):
            return False

        return True


    def preprocess(self):
        self.config.tracer.info('Pre-processing with resolution %d m', self.config.resolution)
        self.config.logger.info('Pre-processing with resolution %d m', self.config.resolution)
        # this is to check the config for the L2A_AtmCorr in ahead.
        # This has historical reasons due to ATCOR porting.
        # Should be moved to the L2A_Config for better design:
        if(self.scOnly == False):
            dummy = L2A_AtmCorr(self.config, None)
            dummy.checkConfiguration()

        # validate the meta data:
        xp = L2A_XmlParser(self.config, 'UP1C')
        xp.validate()
        xp.export()
        xp = L2A_XmlParser(self.config, 'T1C')
        xp.validate()
        xp.export()
        xp = L2A_XmlParser(self.config, 'DS1C')
        xp.validate()
        xp.export()

        if(self.tables.J2kToH5() == False):
            return False   
        return True


    def postprocess(self):
        self.config.tracer.info('Post-processing with resolution %d m', self.config.resolution)
        self.config.logger.info('Post-processing with resolution %d m', self.config.resolution)
        
        res = self.tables.H5ToJ2k(self.scOnly)
        if(self.config.resolution == 60):
            self.config.postprocess()

        # validate the meta data:
        xp = L2A_XmlParser(self.config, 'UP2A')
        xp.validate()
        xp = L2A_XmlParser(self.config, 'T2A')
        xp.validate()
        xp = L2A_XmlParser(self.config, 'DS2A')
        xp.validate()
        return res


    def resetProcessingStatus(self):
        self._processed60 = False
        self._processed20 = False
        self._processed10 = False
        return


def main(args, config):

    if os.path.exists(args.directory) == False:
        stderrWrite('directory "%s" does not exist\n.' % args.directory)
        return False

    processor = L2A_Process(args.directory)
    processor.scOnly = args.sc_only

    HelloWorld = config.processorName +', '+ config.processorVersion +', created: '+ config.processorDate
    stdoutWrite('\n%s started ...\n' % HelloWorld)
    tStart = time()
    S2A_mask = 'S2A_*'

    # next statement creates L2A product Structure:
    tiles = config.createL2A_UserProduct()
    for tile in tiles:
        if(fnmatch.fnmatch(tile, S2A_mask) == False):
            continue

        processor.resetProcessingStatus()
        config.initLogAndTrace()
        config.tracer.info(HelloWorld)
        config.logger.info(HelloWorld)
        config.calcEarthSunDistance2(tile)
        if args.resolution == None:
            resolution = 60.0
        else:
            resolution = args.resolution

        config.resolution = resolution
        config.setTimeEstimation(resolution)
        config.logger.debug('Module L2A_Process initialized')

        result = processor.selectAndProcess(tile)
        if(result == False):
            stderrWrite('Application terminated with errors, see log file and traces.\n')
            return False

        tMeasure = time() - tStart
        config.writeTimeEstimation(resolution, tMeasure)
    
    stdoutWrite('\nApplication terminated successfully.\n')
    return True


if __name__ == "__main__":
    # Someone is launching this directly
    import argparse
    config = L2A_Config()
    descr = config.processorName +'. Version: '+ config.processorVersion + ', created: '+ config.processorDate + \
    ', supporting Level-1C product version: ' + config.productVersion + '.'
     
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('directory', help='Directory where the Level-1C input files are located')
    parser.add_argument('--resolution', type=int, choices=[10, 20, 60], help='Target resolution, must be 10, 20 or 60 [m]')
    parser.add_argument('--sc_only', action='store_true', help='Performs only the scene classification at 60m resolution')
    parser.add_argument('--profile', action='store_true', help='Performs a processor performance profile and displays the results')
    args = parser.parse_args()

    if(args.profile == True):
        import cProfile
        import pstats
        logdir = os.environ['S2L2APPHOME'] + '/log'
        profile = logdir + '/profile'
        cProfile.run('main(args, config)', profile)
        p = pstats.Stats(profile)
        p.strip_dirs().sort_stats('cumulative').print_stats(.25, 'L2A_')
    else:
        main(args, config)
        
