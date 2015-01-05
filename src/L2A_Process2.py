'''
Created on Jan 2, 2015

@author: umwilm
'''

import time
from sys import stdout

if __name__ == '__main__':
    for i in range(1,11):
        stdout.write('Progress[%%]: %03.2f\n' % (i*10))
        stdout.flush()
        time.sleep(1)
        