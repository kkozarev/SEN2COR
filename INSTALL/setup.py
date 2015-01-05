from distutils.core import setup
import os, sys

setup(
    name = 'SEN2COR_SCD',
    version = '1.0.1',
    description = 'Sentinel-2 Level 2A Prototype Processor',
    author = 'Telespazio VEGA Deutschland GmbH',
    author_email = 'info@telespazio-vega.de',
    url = 'www.telespazio-vega.de'
)

home = os.environ['PWD']
path = home+'/SEN2COR/src'

if(os.path.exists(path) == False):
    print 'Directory does not exist.'
    sys.exit()

for root, dirs, files in os.walk(path, topdown=False):
    for file in files:
        file = path + '/' + file 
        os.chmod(file, 0444)

os.chmod(home+'/SEN2COR/src/L2A_Process.py', 0554)
if((os.path.isfile(home+'/SEN2COR/bin/L2A_Process')) == True):
    os.remove(home+'/SEN2COR/bin/L2A_Process')
os.symlink(home+'/SEN2COR/src/L2A_Process.py',home+'/SEN2COR/bin/L2A_Process')
print 'Installation performed.'
