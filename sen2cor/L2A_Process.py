#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

### This module creates the central structure of the L2A_Product and calls the L2A_Schedule module

from numpy import *
from tables import *
import sys, os, logging, fnmatch, warnings, platform, multiprocessing
import cPickle as pickle
from shutil import copyfile
from time import time
from L2A_Logger import getLevel
from L2A_Schedule import L2A_Schedule
from L2A_Config import L2A_Config, getScriptDir
from L2A_XmlParser import L2A_XmlParser
from L2A_Library import stdoutWrite, stderrWrite
from L2A_Manifest import L2A_Manifest
from L2A_ProcessTile import SUCCESS, FAILURE

warnings.filterwarnings("ignore")
formatter = logging.Formatter('<check>\n<inspection execution=\"%(asctime)s\" level=\"%(levelname)s\" process=\"%(process)d\" module=\"%(module)s\" function=\"%(funcName)s\" line=\"%(lineno)d\"/>\n<message contentType=\"Text\">%(message)s</message>\n</check>')


def updateTiles(config):
    dirname, basename = os.path.split(config.workDir)    
    L2A_UP_ID = basename[:4] + 'USER' + basename[8:]
    L2A_UP_ID = L2A_UP_ID.replace('1C_', '2A_')    
    targetDir = config.targetDirectory
    if targetDir != 'DEFAULT':
        dirname = targetDir

    config.L2A_UP_DIR = os.path.join(dirname, L2A_UP_ID)
    L1C_TILES = config.createOrUpdateL2A_UserProduct()
    if L1C_TILES == False:
        return False
    filemask = 'S2A_OPER_*'
    picFn = 'configPic.p'
    L2A_TILES = []
    for tile in L1C_TILES:      
        if fnmatch.fnmatch(tile, filemask) == False:
            continue
        L2A_TILE_ID = config.create_L2A_Tile(tile)
        L2A_TILES.append(L2A_TILE_ID)
        picFnTile = os.path.join(L2A_TILE_ID, picFn)
        logger = config.logger
        config.logger = None
        if(os.path.isfile(picFnTile) == False or config.refresh == True):
            # next statement creates or overwrites the permanent config object:
            # we must remove the logger first ...
            try:
                src = open(picFn, 'wb')
                pickle.dump(config, src, 2)
                src.close()
                copyfile(picFn, picFnTile)
                os.remove(picFn)
            except:
                config.logger = logger
                logger.fatal('cannot create the config object %s' % picFn)
        else:
            tStart = config.tStart
            tEstimation = config.tEstimation
            try:
                f = open(picFnTile, 'rb')
                config = pickle.load(f)
                f.close()             
                config.tStart = tStart   
                config.tEstimation = tEstimation
                f = open(picFnTile, 'wb')                
                pickle.dump(config, f, 2)
                f.close()
            except:
                config.logger = logger
                logger.fatal('cannot update the config object %s with new time estimation' % f)
    
    return L2A_TILES


def postprocess(config):
    HTML = 'HTML'
    SEN2COR = 'SEN2COR'
    REPORT_XML = '_report.xml'
    
    basename = os.path.basename(config.L2A_UP_DIR)
    fileID = basename 
    fnOut =  os.path.join(config.L2A_UP_DIR, HTML, SEN2COR + REPORT_XML)
    filelist = sorted(os.listdir(config.logDir))
    for filename in filelist:
        if((fileID in filename) == False):
            continue
        try:
            fnIn = os.path.join(config.logDir, filename)
            f = open(fnIn, 'a')
            f.write('</Sen2Cor_Level-2A_Report_File>')
            f.close()
            copyfile(fnIn, fnOut)
            return True
        except:
            config.logger.error('cannot copy report file: %s' % fnIn)
            return False
        
    config.logger.error('report file not present: %s' % fnOut)
    return False


def main(args=None):
    import argparse
        
    config = L2A_Config(None)
    descr = config.processorName +'. Version: '+ config.processorVersion + ', created: '+ config.processorDate + \
    ', supporting Level-1C product version: ' + config.productVersion + '.'
     
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('directory', help='Directory where the Level-1C input files are located')
    parser.add_argument('--resolution', type=int, choices=[10, 20, 60], help='Target resolution, can be 10, 20 or 60m. If omitted, all resolutions will be processed')
    parser.add_argument('--sc_only', action='store_true', help='Performs only the scene classification at 60 or 20m resolution')
    parser.add_argument('--cr_only', action='store_true', help='Performs only the creation of the L2A product tree, no processing')
#     parser.add_argument('--profile', action='store_true', help='Profiles the processor\'s performance')
    parser.add_argument('--refresh', action='store_true', help='Performs a refresh of the persistent configuration before start')
    parser.add_argument('--GIP_L2A', help='Select the user GIPP')
    parser.add_argument('--GIP_L2A_SC', help='Select the scene classification GIPP')
    parser.add_argument('--GIP_L2A_AC', help='Select the atmospheric correction GIPP')
    args = parser.parse_args()
    
    # SIITBX-49: directory should not end with '/':
    directory = args.directory
    if directory[-1] == '/' or directory[-1] == '\\':
        directory = directory[:-1]

    # check if directory argument starts with a relative path. If not, expand: 
    if(os.path.isabs(directory)) == False:
        cwd = os.getcwd()
        directory = os.path.join(cwd, directory)
        
    directory = os.path.normpath(directory)
    if os.path.exists(directory) == False:
        stderrWrite('directory "%s" does not exist\n.' % directory)
        return FAILURE

    # check if directory argument contains a tile. If yes, split the tile from path,
    # put the tile in the config object created below as selected tile,
    # create the new path for the user directory.
    selectedTile = None
    if 'GRANULE' in directory:
        dirname, selectedTile = os.path.split(directory)
        directory = os.path.dirname(dirname)
    
    test = os.path.basename(directory)
    S2A_L1C_mask = 'S2A_????_???_???L1C*'
    if(fnmatch.fnmatch(test, S2A_L1C_mask) == False):
        stderrWrite('L1C user product directory must match the following mask: %s\n' % S2A_L1C_mask)
        stderrWrite('but is: %s\n' % test)
        return FAILURE

    config = L2A_Config(None, directory)
    HelloWorld = config.processorName +', '+ config.processorVersion +', created: '+ config.processorDate
    stdoutWrite('\n%s started ...\n' % HelloWorld)

#     if(args.profile == True):    
#         import cProfile, pstats, StringIO
#         pr = cProfile.Profile()
#         pr.enable()
        
    if args.resolution == None:
        resolution = 0
    else:
        resolution = args.resolution

    # create and initialize the base log system:
    dirname, basename = os.path.split(directory)
    L2A_UP_ID = basename[:4] + 'USER' + basename[8:]
    L2A_UP_ID = L2A_UP_ID.replace('1C_', '2A_')
    
    logName = L2A_UP_ID + '_report.xml'
    logDir = config.logDir
    logLevel = config.logLevel
    fnLog = os.path.join(logDir, logName)
    if not os.path.exists(logDir):
        os.mkdir(logDir)
        
    try:
        f = open(config.processingStatusFn, 'w')
        f.write('0.0\n')
        f.close()  
    except:
        stderrWrite('cannot create process status file: %s\n' % config.processingStatusFn)        
        return FAILURE
    try:
        f = open(fnLog, 'w')
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<Sen2Cor_Level-2A_Report_File>\n')
        f.close()
    except:    
        stderrWrite('cannot update the report file: %s\n' % fnLog)        
        return FAILURE

    # Just a normal logger
    logger = logging.getLogger('sen2cor')
    handler = logging.FileHandler(fnLog)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.level = getLevel(logLevel)
    logger.info('logging for the main process initialized')
    config.logger = logger
    
    CFG = 'cfg'
    if args.GIP_L2A != None:
        config._configFn =  os.path.join(config.home, CFG, args.GIP_L2A)
     
    if args.GIP_L2A_SC != None:
        config.configSC =  os.path.join(config.home, CFG, args.GIP_L2A_SC)

    if args.GIP_L2A_AC != None:
        config.configAC =  os.path.join(config.home , CFG, args.GIP_L2A_AC)

    config.workDir = directory
    config.resolution = resolution
    config.scOnly = args.sc_only
    config.crOnly = args.cr_only
    config.refresh = args.refresh
    config.selectedTile = selectedTile
    result = config.readPreferences()
    if result == False:
        return FAILURE

    config.tStart = time()
    config.setTimeEstimation(resolution)

    L2A_TILES = updateTiles(config)
    if L2A_TILES == False:
        return FAILURE
    
    result = SUCCESS
    
    if config.crOnly == False:
        scheduler = L2A_Schedule(config, L2A_TILES)
        result = scheduler.sync()
        config.logger = logger
        # validate the meta data on user product level:
        try:
            xp = L2A_XmlParser(config, 'UP2A')
            xp.validate()
        except:
            logger.error('parsing error for user product')
            result = FAILURE
    
#     if(args.profile == True):    
#         pr.disable()
#         s = StringIO.StringIO()
#         sortby = 'cumulative'
#         ps = pstats.Stats(pr, stream=s).sort_stats(sortby).print_stats(.25, 'L2A_')
#         ps.print_stats()
#         profile = s.getvalue()            
#         s.close()
#         with open(os.path.join(getScriptDir(), 'log', 'profile'), 'w') as textFile:
#             textFile.write(profile)
#             textFile.close()   
    else:
        config.logger = logger        
    #Create the manifest.safe (L2A)
    mn = L2A_Manifest(config)
    mn.generate(config.L2A_UP_DIR, config.L2A_MANIFEST_SAFE)
    try:
        xp = L2A_XmlParser(config, 'Manifest')
        xp.validate()
    except:
        logger.error('parsing error for manifest')
        result = FAILURE
    
    if postprocess(config) == False:
        result = FAILURE
        
    if result == FAILURE:
        stdoutWrite('Progress[%]: 100.00 : Application terminated with at least one error.\n')
    else:
        stdoutWrite('Progress[%]: 100.00 : Application terminated successfully.\n')
    
    return result
    

if __name__ == "__main__":
    if platform.system() == 'Windows':
        multiprocessing.freeze_support()
    sys.exit(main())
