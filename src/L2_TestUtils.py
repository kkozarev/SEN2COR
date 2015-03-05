#!/usr/bin/env python
'''
Created on Feb 24, 2012

@author: umuellerwilm
'''

from numpy import *
#from spectral import *
import Image
import sys
import os
import fnmatch
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage.filters import uniform_filter

try:
    from osgeo import gdal,ogr,osr
    from osgeo.gdalconst import *
    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal
    from gdalconst import *


def showImage(arr):
    if(arr.ndim) != 2:
        print('Must be a two dimensional array')
        return False

    arrmin = arr.mean() - 3*arr.std()
    arrmax = arr.mean() + 3*arr.std()
    arrlen = arrmax-arrmin
    arr = clip(arr, arrmin, arrmax)
    scale = 255.0
    scaledArr = (arr-arrmin).astype(float32) / float32(arrlen) * scale
    arr = (scaledArr.astype(uint8))
    #plt.imshow(arr, interpolation='nearest')
    #plt.show()
    img = Image.fromarray(arr)
    img.show()
    return True


class L2A_TestUtils(object):
    print('Module L2A_TestUtils initialized')
    def __init__(self, rootDir, srcDir, resolution):
        self._tstDir = srcDir
        if resolution == 10:
            self._srcDir = rootDir + '/' + srcDir + '/TESTS_10'
        elif resolution == 20:
            self._srcDir = rootDir + '/' + srcDir + '/TESTS_20'
        else:
            self._srcDir = rootDir + '/' + srcDir + '/TESTS_60'
        self._tgtDir = rootDir + '/' + srcDir + '_results'
        if(os.path.exists(self._tgtDir) == False):
            os.mkdir(self._tgtDir)

        self._noData = 0
        self._saturatedDefective = 1
        self._darkFeatures = 2
        self._cloudShadows = 3
        self._vegetation = 4
        self._bareSoils = 5
        self._water = 6
        self._lowProbaClouds = 7
        self._medProbaClouds = 8
        self._highProbaClouds = 9
        self._thinCirrus = 10
        self._snowIce = 11

        self._SC = self.getArray('S2L2APP_OUT_SCL')

        self._reflArr = zeros([9,19], float32)
        self._reflArr[0,0] = 0.492225
        self._reflArr[1,0] = 0.560310
        self._reflArr[2,0] = 0.663085
        self._reflArr[3,0] = 0.703959
        self._reflArr[4,0] = 0.742381
        self._reflArr[5,0] = 0.781725
        self._reflArr[6,0] = 0.865816
        self._reflArr[7,0] = 1.609431
        self._reflArr[8,0] = 2.193888

    def __exit__(self):
        self._pp.close()


    def statistics(self, arr, min, max, comment = ''):
        maskInp = '*_INP_0*'
        maskOut = '*_OUT_0*'
 
        if max > 3000:
            max = 4000
        elif max > 2000:
            max = 3000
        elif max > 1000:
            max = 2000
        elif max > 500:
            max = 1000
        elif max > 255:
            max = 500
        else:
            max = 260

        if(fnmatch.fnmatch(comment, maskInp) == True):
            arr = arr.astype(float32) / 4000.0
            max = float32(max) / 4000.0
        elif(fnmatch.fnmatch(comment, maskOut) == True):
            arr = arr.astype(float32) / 2000.0
            max = float32(max) / 2000.0

        mu = arr[self._nodata > 0].mean() # mean of distribution
        sigma = arr[self._nodata > 0].std() # standard deviation of distribution
        #mu = arr.mean() # mean of distribution
        #sigma = arr.std() # standard deviation of distribution
        min = 0
        max = 1.0
        num_bins = 50

        dummy, basename = os.path.split(self._tstDir)
        x = arr[self._nodata > 0]
        #x = arr[arr>0]
        n, bins, patches = plt.hist(x, num_bins, range=[min,max], normed=True, facecolor='green')
        # add a 'best fit' line
        y = mlab.normpdf(bins, mu, sigma)
        plt.plot(bins, y, 'r--')
        plt.xlabel('Reflectance')
        plt.ylabel('Relative Frequency')
        plt.title(comment)
        # Tweak spacing to prevent clipping of ylabel
        plt.subplots_adjust(left=0.15)
        plt.savefig(self._tgtDir + '/' + comment + '_hist.jpg', dpi=100)
        plt.close()

        f = open(self._tgtDir + '/' + basename + '_statistics.txt', 'a')
        f.write('Image:' + str(comment) + '\n')
        f.write('--------------------' + '\n')
        f.write('shape: %s\n' % str(shape(arr)))
        f.write('mean : %05.3f\n' % mu)
        f.write('std  : %05.3f\n' % sigma)
        f.write('min  : %05.3f\n' % arr[self._nodata > 0].min())
        f.write('max  : %05.3f\n' % arr[self._nodata > 0].max())
        f.write('Histogram, bins = 50, range =' + str(min) + ':' + str(max) +'\n')
        h = histogram(arr[self._nodata > 0], num_bins, range=(min,max))
        '''
        f.write('min  : %05.3f\n' % arr.min())
        f.write('max  : %05.3f\n' % arr.max())
        f.write('Histogram, bins = 50, range =' + str(min) + ':' + str(max) +'\n')
        h = histogram(arr, num_bins, range=(min,max))
        '''
        f.write(str(h[0].astype(float32)/h[0].sum()) + '\n\n')
        f.close()

        return


    def getArray(self, filename):
        filename = self._srcDir + '/' + filename + '.npy'
        if((os.path.isfile(filename)) == False):
            print('File ' + filename + ' not present')
            return False

        return load(filename)


    def setNoData(self):
        B03 = self.getArray('S2L2APP_INP_003')
        B8A = self.getArray('S2L2APP_INP_008')
        self._nodata = ones(B03.shape, uint8)
        self._nodata[(B03==0) & (B8A==0)] = 0
        return


    def getDiffFromArrays(self, filename1, filename2):
        filename1 = self._srcDir + '/' + filename1 + '.npy'
        filename2 = self._srcDir + '/' + filename2 + '.npy'
        if((os.path.isfile(filename1)) == False):
            print('File ' + filename1 + ' not present')
            return False

        if((os.path.isfile(filename2)) == False):
            print('File ' + filename2 + ' not present')
            return False

        arr1 = load(filename1)
        arr2 = load(filename2)
        arr = (arr1 - arr2)

        # treatment of outliers: mean + 3 * std dev is maximum for display
        mean_arr = arr.mean()
        std_arr = arr.std()
        min_arr = mean_arr - 3 * std_arr
        max_arr = mean_arr + 3 * std_arr
        arr = clip(arr, min_arr, max_arr)

        return arr


    def fft2(self, filename1, filename2):
        filename1 = self._srcDir + '/' + filename1 + '.npy'
        filename2 = self._srcDir + '/' + filename2 + '.npy'
        if((os.path.isfile(filename1)) == False):
            print('File ' + filename1 + ' not present')
            return False

        if((os.path.isfile(filename2)) == False):
            print('File ' + filename2 + ' not present')
            return False

        arr1 = load(filename1)
        arr2 = load(filename2)
        m1 = fft.fft2(arr1)
        m2 = fft.fft2(arr2)
        s = sum(abs(m1-m2))
        return s


    def saveArray(self, filename, arr):
        filename = self._srcDir + '/' + filename + '.npy'
        save(filename, arr)
        return


    def saveScaledTestImage(self, arr, min_arr, max_arr, filename):
        if(arr.ndim) != 2:
            print('Must be a two dimensional array')
            return False

        self.statistics(arr, min_arr, max_arr, filename)
        filename = self._tgtDir + '/' + filename + '.jpg'
        arrlen = (max_arr-min_arr)
        arr = clip(arr, min_arr, max_arr)
        scale = 255.0
        scaledArr = (arr-min_arr).astype(float32) / arrlen * scale
        arr = (scaledArr.astype(uint8))
        img = Image.fromarray(arr)
        img.save(filename, "JPEG")
        return True


    def showPreview(self):
        a = Image.open(self._srcDir + '/S2L2APP_INP_RGB.jpg')
        b = Image.open(self._srcDir + '/S2L2APP_OUT_RGB.jpg')
        c = Image.open(self._srcDir + '/ATCOR_OUT_RGB.jpg')
        a.show()
        b.show()
        c.show()
        return


    def toSC(self):
        myColors = zeros([12, 3], uint)
        myColors[0,:] = [0,0,0]      # 0: SHADOW, Black
        myColors[1,:] = [255,0,0]    # 1: DEFECTIVE, Red
        myColors[2,:] = [50,50,50]   # 2: SHADOW, Very Dark Gray
        myColors[3,:] = [128,64,0]    # 3: CLOUD_SHADOW, Brown
        myColors[4,:] = [0,128,0]     # 4: VEGETATION, Green
        myColors[5,:] = [255,255,0]   # 5: BARE_SOILS, Yellow
        myColors[6,:] = [0,0,255]     # 6: WATER  Blue
        myColors[7,:] = [100,100,100] # 7: CLOUD_LOW_PROBA, Dark Gray
        myColors[8,:] = [150,150,150] # 8: CLOUD_MED_PROBA, Light Gray
        myColors[9,:] = [255,255,255] # 9: CLOUD_HI_PROBA, White
        myColors[10,:] = [128,255,255] # 10: THIN_CIRRUS, Light Blue
        myColors[11,:] = [255,128,255] # 11: SNOW, Pink

        filename = 'S2L2APP_OUT_SCL'
        arr = self.getArray(filename)
        #save_rgb(self._tgtDir + '/S2L2APP_OUT_SCL.jpg', arr, colors=myColors)
        return


    def toRGB(self):
        b = Image.open(self._tgtDir + '/S2L2APP_INP_001.jpg')
        g = Image.open(self._tgtDir + '/S2L2APP_INP_002.jpg')
        r = Image.open(self._tgtDir + '/S2L2APP_INP_003.jpg')

        out = Image.merge('RGB', (r,g,b))
        out.save(self._tgtDir + '/S2L2APP_INP_RGB.jpg')
        self.statistics(array(out), 0, 255, 'S2L2APP_INP_RGB')

        b = Image.open(self._tgtDir + '/S2L2APP_OUT_001.jpg')
        g = Image.open(self._tgtDir + '/S2L2APP_OUT_002.jpg')
        r = Image.open(self._tgtDir + '/S2L2APP_OUT_003.jpg')

        out = Image.merge('RGB', (r,g,b))
        out.save(self._tgtDir + '/S2L2APP_OUT_RGB.jpg')
        self.statistics(array(out), 0, 255, 'S2L2APP_OUT_RGB')
        return


    def calcReflectanceSpectra(self, arr, i, type):
        #wl, mw, sdw, mv, sdv, ms, sds,
        if(type == 'INP'):
            arr = arr.astype(float32) / 4000.0
            self._reflArr[i,1] = arr[self._SC == self._bareSoils].mean()
            self._reflArr[i,2] = arr[self._SC == self._bareSoils].std()
            self._reflArr[i,3] = arr[self._SC == self._bareSoils].size

            self._reflArr[i,4] = arr[self._SC == self._vegetation].mean()
            self._reflArr[i,5] = arr[self._SC == self._vegetation].std()
            self._reflArr[i,6] = arr[self._SC == self._vegetation].size

            self._reflArr[i,7] = arr[self._SC == self._water].mean()
            self._reflArr[i,8] = arr[self._SC == self._water].std()
            self._reflArr[i,9] = arr[self._SC == self._water].size

        elif(type == 'OUT'):
            arr = arr.astype(float32) / 2000.0
            self._reflArr[i,10] = arr[self._SC == self._bareSoils].mean()
            self._reflArr[i,11] = arr[self._SC == self._bareSoils].std()
            self._reflArr[i,12] = arr[self._SC == self._bareSoils].size

            self._reflArr[i,13] = arr[self._SC == self._vegetation].mean()
            self._reflArr[i,14] = arr[self._SC == self._vegetation].std()
            self._reflArr[i,15] = arr[self._SC == self._vegetation].size

            self._reflArr[i,16] = arr[self._SC == self._water].mean()
            self._reflArr[i,17] = arr[self._SC == self._water].std()
            self._reflArr[i,18] = arr[self._SC == self._water].size


    def plotReflectanceSpectra(self):
        # example data
        plt.xlim(0.4,2.3)
        plt.ylim(0,0.4)
        x = self._reflArr[:,0]
        plt.title('Reflectance Spectra')
        plt.xlabel('Wavelength [$\\mu$m]')
        plt.ylabel('Reflectance')

        y = self._reflArr[:,10]
        e = self._reflArr[:,11]
        plt.errorbar(x, y, yerr=e, linewidth=3, color='yellow')
        y = self._reflArr[:,13]
        e = self._reflArr[:,14]
        plt.errorbar(x, y, yerr=e, linewidth=3, color='green')
        y = self._reflArr[:,16]
        e = self._reflArr[:,17]
        plt.errorbar(x, y, yerr=e, linewidth=3, color='blue')                
        y = self._reflArr[:,1]
        e = self._reflArr[:,2]
        plt.errorbar(x, y, yerr=e, linewidth=1, color='black')
        y = self._reflArr[:,4]
        e = self._reflArr[:,5]
        plt.errorbar(x, y, yerr=e, linewidth=1, color='black')
        y = self._reflArr[:,7]
        e = self._reflArr[:,8]
        plt.errorbar(x, y, yerr=e, linewidth=1, color='black')

        ## legend
        plt.legend((r'Bare Soils', r'Vegetation', r'Water', r'TOA'), loc = (0.8, 0.8))
        ltext = plt.gca().get_legend().get_texts()
        plt.setp(ltext[0], fontsize = 12, color = 'y')
        plt.setp(ltext[1], fontsize = 12, color = 'g')
        plt.setp(ltext[2], fontsize = 12, color = 'b')
        plt.setp(ltext[3], fontsize = 12, color = 'black')
        plt.savefig(self._tgtDir + '/Reflectance_Spectra.jpg', dpi=100)
        plt.close()


def main():
    import argparse
    descr = 'S2L2APP test Utility'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('directory', help='Working directory for module')
    parser.add_argument('--resolution', type=int, choices=[10, 20, 60], help='Target resolution, must be 10, 20 or 60 [m]')

    args = parser.parse_args()
    set_printoptions(precision=3, suppress=True)
    rootDir = args.directory
    bandIndex = [1,2,3,4,5,6,8,11,12]

    if args.resolution == None:
        resolution = 60.0
    else:
        resolution = args.resolution

    filemask = 'S2A_*'
    dirs = sorted(os.listdir(rootDir))
    for dir in dirs:
        if(fnmatch.fnmatch(dir, filemask) == False):
            continue
        tutil = L2A_TestUtils(rootDir, dir, resolution)
        tutil.setNoData()
        '''
        f8 = 'S2L2APP_INP_008'
        f9 = 'S2L2APP_INP_009'
        a8 = tutil.getArray(f8)
        a9 = tutil.getArray(f9)
        #showImage(clip(a8/maximum(a9,0.001),0,10))
        #a0 = clip(a8[0:4998,0:1755]/maximum(a9[1:4999,1:1756],0.001),0,10)
        #a1 = clip(a9[1:4999,1:1756]/maximum(a8[0:4998,0:1755],0.001),0,10)
        #a0 = clip(a8/maximum(a9,0.001),0,10)
        a1 = clip(a9/maximum(a8,0.000001),0,1)
        #showImage(a8)
        showImage(a9-a8)
        showImage(a9[1:4999,1:1756]-a8[0:4998,0:1755])
        showImage(a8-a9)
        showImage(a8[0:4998,0:1755]-a9[1:4999,1:1756])
        #showImage(a0)
        #showImage(a1)
        #showImage(uniform_filter(a1,size=3))
        return
        '''
        # create S2L2A input Images from numpy arrays:
        count = 0
        for i in arange(0, 13):
            filename = 'S2L2APP_INP_%03d' % i
            arr = tutil.getArray(filename)
            if(isinstance(arr, ndarray)):
                min_arr = 0
                max_arr = 1000
                tutil.saveScaledTestImage(arr, min_arr, max_arr, filename)
                if(i in bandIndex):
                    tutil.calcReflectanceSpectra(arr, count, 'INP')
                    count += 1
        # create S2L2A output Images from numpy arrays:
        count = 0
        for i in arange(0, 13):
            filename = 'S2L2APP_OUT_%03d' % i
            arr = tutil.getArray(filename)
            if(isinstance(arr, ndarray)):
                min_arr = 0
                max_arr = 500
                tutil.saveScaledTestImage(arr, min_arr, max_arr, filename)
                if(i in bandIndex):
                    tutil.calcReflectanceSpectra(arr, count, 'OUT')
                    count += 1
        tutil.plotReflectanceSpectra()
        print tutil._reflArr

        # do the same for the generated output:
        #for s in ['AOT', 'VIS', 'WVP']:
        for s in ['AOT', 'WVP', 'DEM', 'ASP', 'SDW', 'SLP']:
            filename = 'S2L2APP_OUT_%s' % s
            arr = tutil.getArray(filename)
            if(isinstance(arr, ndarray)):
                min_arr = arr[tutil._nodata > 0].mean() - 3 * arr[tutil._nodata > 0].std()
                #min_arr = 0
                max_arr = arr[tutil._nodata > 0].mean() + 3 * arr[tutil._nodata > 0].std()
                tutil.saveScaledTestImage(arr, min_arr, max_arr, filename)

        tutil.toSC()
        tutil.toRGB()

    print 'all done!'
    return


if __name__ == "__main__":
    # Someone is launching this directly
    main()
