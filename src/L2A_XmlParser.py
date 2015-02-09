from numpy import *
import os, sys
from lxml import etree, objectify
from L2A_Library import stdoutWrite, stderrWrite
from L2A_Borg import Borg

class L2A_XmlParser(Borg):
    def __init__(self, config, product):
        self._config = config
        self._product = product
        self._xmlFn = None
        self._xmlName = None
        self._root = None
        self._tree = None
        self._scheme = None
        
        try:
            doc = objectify.parse(config.configFn)
            root = doc.getroot()
            cs = root.Common_Section
            upScheme1c = cs.UP_Scheme_1C.text
            upScheme2a = cs.UP_Scheme_2A.text
            tileScheme1c = cs.Tile_Scheme_1C.text
            tileScheme2a = cs.Tile_Scheme_2A.text
            dsScheme1c = cs.DS_Scheme_1C.text
            dsScheme2a = cs.DS_Scheme_2A.text
        except:
            config.tracer.fatal('Error in parsing config file.')
            config.exitError();

        if(product == 'UP1C'):
            self._xmlFn = config.L1C_UP_MTD_XML
            self._scheme = upScheme1c
        elif(product == 'UP2A'):
            self._xmlFn = config.L2A_UP_MTD_XML
            self._scheme = upScheme2a
        elif(product == 'DS1C'):
            self._xmlFn = config.L1C_DS_MTD_XML
            self._scheme = dsScheme1c
        elif(product == 'DS2A'):
            self._xmlFn = config.L2A_DS_MTD_XML
            self._scheme = dsScheme2a
        elif(product == 'T1C'):
            self._xmlFn = config.L1C_TILE_MTD_XML
            self._scheme = tileScheme1c
        elif(product == 'T2A'):
            self._xmlFn = config.L2A_TILE_MTD_XML
            self._scheme = tileScheme2a
        elif(product == 'GIPP'):
            self._xmlFn = config.configFn
            self._scheme = None
        else:
            config.tracer.fatal('wrong identifier for xml structure: ' + product)
            config.exitError()
        
        self.setRoot();
        return


    def getRoot(self, key=None):
        try:
            if key == None:
                return self._root
            else:
                root = self._root[key]
                return root
        except:
            return False


    def setRoot(self):
        if self._root != None:
            return True
        try:
            doc = objectify.parse(self._xmlFn)
            self._root = doc.getroot()
            return True
        except:
            return False


    def getTree(self, key, subkey):
        try:
            tree = self._root[key]    
            return tree['{}' + subkey]
        except:
            return False


    def setTree(self, key, subkey):
        try:
            root = self._root[key]
        except:
            return False
        try:
            self._tree = root['{}' + subkey]
            return True
        except:
            self._tree = root
            if(self.append(subkey, '') == True):
                try:
                    self._tree = root['{}' + subkey]
                    self.export()
                    return True        
                except:
                    return False
        return False

    
    def validate(self):
        dummy, fn = os.path.split(self._xmlFn)
        stdoutWrite('Validating metadata %s against scheme ...\n' % fn)       
        try:
            schema = etree.XMLSchema(file = self._config.configDir + '/' + self._scheme)
            parser = etree.XMLParser(schema = schema)
            objectify.parse(self._xmlFn, parser)
            stdoutWrite('Metadata is valid.\n')                
            return True
        except etree.XMLSyntaxError, err:
            stderrWrite('Error in validation:\n')
            stderrWrite('- Schema file: %s\n' % self._scheme)
            stderrWrite('- Details: %s\n\n' % str(err))
            stderrWrite('Application will be forced to terminate,\n')
            stderrWrite('Please correct the errors and restart.\n')
            sys.exit(-1)
            return False


    def append(self, key, value):
        try:
            e = etree.Element(key)
            e.text = value
            self._tree.append(e)
            return True
        except:
            return False


    def export(self):
        import codecs
        outfile = codecs.open(self._xmlFn, 'w', 'utf-8')
        outfile.write('<?xml version="1.0"  encoding="UTF-8"?>\n')
        objectify.deannotate(self._root, xsi_nil=True, cleanup_namespaces=True)
        outstr = etree.tostring(self._root, pretty_print=True)
        outfile.write(outstr)        
        outfile.close()
        return self.setRoot()

    
    def convert(self):
        import codecs
        outfile = codecs.open(self._xmlFn, 'w', 'utf-8')
        outfile.write('<?xml version="1.0"  encoding="UTF-8"?>\n')
        objectify.deannotate(self._root, xsi_nil=True, cleanup_namespaces=True)
        outstr = etree.tostring(self._root, pretty_print=True)
        outstr = outstr.replace('-1C', '-2A')
        outstr = outstr.replace('Product_Info>', 'L2A_Product_Info>')
        outstr = outstr.replace('Product_Organisation>', 'L2A_Product_Organisation>')
        outstr = outstr.replace('Product_Image_Characteristics>', 'L2A_Product_Image_Characteristics>')
        outstr = outstr.replace('Pixel_Level_QI', 'L1C_Pixel_Level_QI')
        outstr = outstr.replace('TILE_ID>', 'TILE_ID_2A>')
        outstr = outstr.replace('DATASTRIP_ID>', 'DATASTRIP_ID_2A>')
        if self._product == 'T2A':
            outstr = outstr.replace('Image_Content_QI>', 'L1C_Image_Content_QI>')
        if self._product == 'UP2A':
            outstr = outstr.replace('QUANTIFICATION_VALUE', 'L1C_L2A_Quantification_Values_List')
            outstr = outstr.replace('</n1:Auxiliary_Data_Info>', '</n1:Auxiliary_Data_Info>\n'\
                                '<n1:L2A_Auxiliary_Data_Info/>')
            outstr = outstr.replace('</n1:Quality_Indicators_Info>', '</n1:Quality_Indicators_Info>\n'\
                                '<n1:L2A_Quality_Indicators_Info/>')
        outfile.write(outstr)
        outfile.close()
        return self.setRoot()


    def getIntArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False

        ncols = len(node[0].text.split())
        a = zeros([nrows,ncols],dtype=int)        

        for i in range(nrows):
            a[i,:] = nan_to_num(array(node[i].text.split(),dtype(int)))
            
        return a


    def getUintArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False
        
        ncols = len(node[0].text.split())
        a = zeros([nrows,ncols],dtype=uint)
        
        for i in range(nrows):
            a[i,:] = nan_to_num(array(node[i].text.split(),dtype(uint)))
            
        return a

    
    def getFloatArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False
        
        ncols = len(node[0].text.split())
        a = zeros([nrows,ncols],dtype=float32)
        
        for i in range(nrows):
            a[i,:] = nan_to_num(array(node[i].text.split(),dtype(float32)))
            
        return a


    def getStringArray(self, node):
        nrows = len(node)
        if nrows < 0:
            return False
        
        ncols = len(node[0].text.split())
        a = zeros([nrows,ncols],dtype=str)
        
        for i in range(nrows):
            a[i,:] = array(node[i].text.split(),dtype(str))

        return a


    def setArrayAsStr(self, node, a):
        set_printoptions(precision=6)
        if a.ndim == 1:
            nrows = a.shape[0]
            for i in nrows:
                node[i] = a[i],dtype=str
                        
        elif a.ndim == 2:
            nrows = a.shape[0]
            ncols = a.shape[1]
            for i in range(nrows):
                aStr = array_str(a[i,:]).strip('[]')
                node[i] = aStr
            return True
        else:
            return False

    
    def getViewingIncidenceAnglesArray(self, node, bandId, detectorId, _type='Zenith'):
        nrows = len(node)
        for i in range(nrows):
            if((int(node[i].bandId) == bandId) and (int(node[i].detectorId) == detectorId)):
                if _type == 'Zenith':
                    a = self.getFloatArray(node[i].Zenith.Values_List.VALUES)
                elif _type == 'Azimuth':
                    a = self.getFloatArray(node[i].Azimuth.Values_List.VALUES)
                return a
        return False


    def setViewingIncidenceAnglesArray(self, node, arr, bandId, detectorId, _type='Zenith'):
        nrows = len(node)
        for i in range(nrows):
            if((int(node[i].bandId) == bandId) and (int(node[i].detectorId) == detectorId)):
                if _type == 'Zenith':
                    return self.setArrayAsStr(node[i].Zenith.Values_List.VALUES, arr)
                elif _type == 'Azimuth':
                    return self.setArrayAsStr(node[i].Azimuth.Values_List.VALUES, arr)

        return False

    