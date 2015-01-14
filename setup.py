#!/usr/bin/env python

import os, sys, platform, subprocess
from shutil import copyfile
from distutils.core import setup

desription = 'SEN2COR is a Prototype Processor for processing \
Sentinel-2 Top of Atmosphere reflectance (Level 1C) data into Bottom of \
Atmoshperic corrected (Level 2A) data. It additionally performs a Scene \
Classification of the corresponding input. For details, read the Software \
User Manual which is attached to this python module.\
\
This Software contains an IPR of DyLR. Thus is is not Open Source and the \
module for the Atmospheric Correction is attached to this software in form \
of a binary library only.',

_name = 'SEN2COR'

setup (
    name = _name,
    version = '2.0.0',
    description = 'SEN2COR: Sentinel 2 Level 2A Prototype Processor',
    author = 'Telespazio VEGA Deutschland GmbH',
    author_email = 'info@telespazio-vega.de',
    url = 'www.telespazio-vega.de',
    download_url = 'www.telespazio-vega.de',
    platforms=['linux-x86_64', 'macosx-10.5-x86_64', 'win-amd64'],
)

user_input = raw_input('\nThis script will set up and test the configuration of the SEN2COR processor.\nOK to continue? [y/n]: ')

if user_input == 'n':
    sys.exit(0)

try:
    path = os.environ['PATH']
    if 'Anaconda' in path == False:
        raise KeyError('ANACONDA')
except KeyError as err:
    if err.args[0] == 'PATH':
        sys.stderr.write('No PATH variable configured.\n\
        Please follow the installation procedure of the user manual\n\
        before continuing ...')
    elif err.args[0] == 'ANACONDA':
        sys.stderr.write('Anaconda seems to be not installed or properly configured.\n\
        Please follow the installation procedure of the user manual\n\
        before continuing ...')    
    sys.exit(-1)

home = os.getcwd()+'/'+_name +'/'
build = home+'build/'
bin = home+'bin/'
src = home+'src/'

# Unlock files for editing:
for root, dirs, files in os.walk(home, topdown=False):
    for fn in files:
        fn = root + '/' + fn
        os.chmod(fn, 0666)

system = platform.system()
if system == 'Windows':
    fn = 'L2A_AtmCorr.pyd'
    libf = build+'lib.win-amd64-2.7/'+fn
    if os.path.exists(libf) == False:
        sys.stderr.write('File % does not exist!' % libf)

    srcf = src + fn
    copyfile(libf, srcf)
    copyfile(build+'L2A_process.bat', bin+'.')
    copyfile(build+'geojasper.exe', bin+'.')

elif (system == 'Linux') or (system == 'Darwin'):
    if os.path.exists(src) == False:
        sys.stderr.write('Directory % does not exist!' % src)
        sys.exit(-1)

    gdalBinConfigured = False
    gdalPckConfigured = False
    gdalResConfigured = False
    
    while True:
        # setting the environments for the application into L2A_Bashrc:
        s2l2apphome = home
        s2l2appcfg = home+'cfg'
        os.environ['S2L2APPHOME'] = s2l2apphome
        os.environ['S2L2APPCFG'] = s2l2appcfg
        L2A_Bashrc = 'export S2L2APPHOME='+s2l2apphome+'\n'
        L2A_Bashrc += 'export S2L2APPCFG='+s2l2appcfg+'\n'
        L2A_Bashrc += 'export PATH=$S2L2APPHOME/bin:$PATH\n'
    
        fn = 'L2A_AtmCorr.so'
        if system == 'Linux':
            libf = build+'lib.linux-x86_64-2.7/'+fn
            gdalBinaries = '/usr/bin'     
            gdalPackages = '/usr/lib64/python2.7/site-packages'
            gdalResources = '/usr/local/share/gdal'        
        elif system == 'Darwin':
            L2A_Bashrc += 'export LC_ALL=en_US.UTF-8\n' 
            L2A_Bashrc += 'export LANG=en_US.UTF-8\n'
            libf = build+'lib.macosx-10.5-x86_64-2.7/'+fn
            gdalBinaries = '/Library/Frameworks/GDAL.framework/Programs'
            gdalPackages = '/Library/Frameworks/GDAL.framework/Versions/Current/Python/2.7/site-packages'
            gdalResources = '/Library/Frameworks/GDAL.framework/Resources'

        break_condition = True
        sys.stdout.write('\nCheck if GDAL binaries are installed and configured:\n')
        try:
            command = gdalBinaries + '/gdalinfo --version'
            if subprocess.call(command, shell=True) != 0:
                raise
        except:
            sys.stdout.write('\nNo GDAL binaries found.\n')
            sys.stdout.write('The GDAL binaries are expected in the following directory:\n')
            sys.stdout.write(gdalBinaries+'\n')
            user_input = raw_input('Is this OK? [y/n]: ')
            if user_input == 'n':
                gdalBinaries = raw_input('New Path: ')
            break_condition = False
        L2A_Bashrc += 'export PATH='+ gdalBinaries + ':$PATH\n'
                
        sys.stdout.write('\nCheck if GDAL python bindings are installed and configured:\n')   
        try:
            pythonpath = os.environ['PYTHONPATH']
            if(pythonpath.find('GDAL')>=0):
                sys.stdout.write('GDAL is part of $PYTONPATH.\n')       
            else:
                raise KeyError
        except KeyError:
            sys.stdout.write('No GDAL python bindings configured via $PYTHONPATH.\n')
            sys.stdout.write('The GDAL python bindings are expected in the following directory:\n')
            sys.stdout.write(gdalPackages+'\n')
            user_input = raw_input('Is this OK? [y/n]: ')
            if user_input == 'n':
                gdalPackages = raw_input('New Path: ')
            os.environ['PYTHONPATH'] = gdalPackages + ':' + pythonpath
            break_condition = False
        L2A_Bashrc += 'export PYTHONPATH=' + gdalPackages + ':$PYTHONPATH\n'                
        
        sys.stdout.write('\nCheck if GDAL resources are installed and configured:\n')         
        try:
            gdalData = os.environ['GDAL_DATA']
            if os.path.exists(gdalData) == True:
                sys.stdout.write('GDAL resources are configured via $GDAL_DATA.\n')   
            else:
                raise
        except:
            sys.stdout.write('No GDAL resources configured via $GDAL_DATA.\n')
            sys.stdout.write('The GDAL resources are expected in the following directory:\n')
            sys.stdout.write(gdalResources+'\n')
            user_input = raw_input('Is this OK? [y/n]: ')
            if user_input == 'n':
                gdalResources = raw_input('New Path: ')
            os.environ['GDAL_DATA'] = gdalPackages
            break_condition = False                
        L2A_Bashrc += 'export GDAL_DATA='+ gdalResources + '\n'            
        
        if break_condition == True:
            break

    textFile = open(home+'L2A_Bashrc', "w")
    textFile.write(L2A_Bashrc)
    textFile.close()

    if os.path.exists(libf) == False:
        sys.stderr.write('File % does not exist!' % libf)
        sys.exit(-1)
    srcf = src + fn
    copyfile(libf, srcf)

    if((os.path.isfile(bin+'L2A_Process')) == True):
        os.remove(bin+'L2A_Process')
    os.symlink(src+'L2A_Process.py', bin+'L2A_Process')

    # Lock files for write protection:
    for root, dirs, files in os.walk(home, topdown=False):
        for fn in files:
            fn = root + '/' + fn
            os.chmod(fn, 0444)
    os.chmod(home+'L2A_Bashrc', 0644)
    os.chmod(src+'L2A_Process.py', 0554)

else:
    sys.stderr.wite('Unsupported Operating System,\n\
                     installation not performed!')
    sys.exit(-1)

logd = home+'log'
if not os.path.exists(logd):
    os.makedirs(logd)

sys.stdout.write('\nInstallation performed in '+home+'.\n') 
sys.stdout.write('All configuration is provided in the L2A_Bashrc shell script.\n')
sys.stdout.write('The script must be sourced either manually \"source L2A_Bashrc\",\n')
sys.stdout.write('or from the users .profile or from the Toolbox start script.\n')

sys.exit(0)

if __name__ == '__main__':
    pass