from numpy import *
import os, sys
from lxml import etree, objectify
from L2A_Library import stdoutWrite, stderrWrite
from multiprocessing import Lock
l = Lock()

class L2A_XmlParser(object):
    def __init__(self, config, product):
        self._config = config
        self._logger = config.logger
        self._product = product
        self._xmlFn = None
        self._xmlName = None
        self._root = None
        self._tree = None
        self._scheme = None
        
        l.acquire()
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
            # gippScheme = cs.GIPP_Scheme.text
            # scScheme = cs.SC_Scheme.text
            # acScheme = cs.AC_Scheme.text
    
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
                self._scheme = "L2A_GIPP.xsd"
            elif(product == 'SC_GIPP'):
                self._xmlFn = config.configSC
                self._scheme = 'L2A_CAL_SC_GIPP.xsd'
            elif(product == 'AC_GIPP'):
                self._xmlFn = config.configAC
                self._scheme = "L2A_CAL_AC_GIPP.xsd"   
            elif(product == 'Manifest'):
                self._xmlFn = config.L2A_MANIFEST_SAFE
                self._scheme = 'S2-PDGS-TAS-DI-PSD-V13.1_SAFE/resources/xsd/int/esa/safe/sentinel/1.1/sentinel-2/msi/archive_l2a_user_product/xfdu.xsd'     
            else:
                self._logger.fatal('wrong identifier for xml structure: ' + product)
        finally:
            l.release()
             
        self.setRoot();


    def setRoot(self):
        l.acquire()
        try:
            doc = objectify.parse(self._xmlFn)
            self._root = doc.getroot()
            ret = True
        except:
            ret = False
        finally:
            l.release()
            
        return ret


    def getRoot(self, key=None):
        l.acquire()
        try:
            if key == None:
                root = self._root
            else:
                root = self._root[key]
        except:
            root = False
        finally:
            l.release()
            
        return root


    def getTree(self, key, subkey):
        l.acquire()
        try:
            tree = self._root[key]    
            ret = tree['{}' + subkey]
        except:
            ret = False
        finally:
            l.release()
            
        return ret
    
    
    def validate(self):
        fn = os.path.basename(self._xmlFn)
        self._logger.info('validating metadata file %s against scheme' % fn)
        l.acquire()
        try:
            schema = etree.XMLSchema(file = os.path.join(self._config.configDir, self._scheme))
            parser = etree.XMLParser(schema = schema)
            objectify.parse(self._xmlFn, parser)
            self._logger.info('metadata file is valid')                
            ret = True
        except etree.XMLSyntaxError, err:
            stdoutWrite('Metadata file is invalid, see report file for details.\n')
            self._logger.error('Schema file: %s' % self._scheme)
            self._logger.error('Details: %s' % str(err))
            ret = False
        except:
            stdoutWrite('Unspecific Error in metadata.\n')
            self._logger.error('unspecific error in metadata')
            ret = False
        finally:                  
            l.release()
            
        return ret


    def append(self, key, value):
        ret = False
        l.acquire()
        try:
            e = etree.Element(key)
            e.text = value
            self._tree.append(e)
            ret = True
        except:
            ret = False
        finally:
            l.release()
            
        return ret
        

    def export(self):
        import codecs
        l.acquire()
        try:
            outfile = codecs.open(self._xmlFn, 'w', 'utf-8')         
            outfile.write('<?xml version="1.0"  encoding="UTF-8"?>\n')
            objectify.deannotate(self._root, xsi_nil=True, cleanup_namespaces=True)
            outstr = etree.tostring(self._root, pretty_print=True)
            outfile.write(outstr)        
            outfile.close()
        finally:
            l.release()

        return self.setRoot()            

    

    def convert(self):
        import codecs
        objectify.deannotate(self._root, xsi_nil=True, cleanup_namespaces=True)
        outstr = etree.tostring(self._root, pretty_print=True)
        outstr = outstr.replace('-1C', '-2A')
        outstr = outstr.replace('<TILE_ID',                        '<TILE_ID_2A')
        outstr = outstr.replace('<DATASTRIP_ID',                   '<DATASTRIP_ID_2A')
        outstr = outstr.replace('<Product_Info',                   '<L2A_Product_Info')
        outstr = outstr.replace('<Product_Organisation',           '<L2A_Product_Organisation')
        outstr = outstr.replace('<Product_Image_Characteristics',  '<L2A_Product_Image_Characteristics')
        outstr = outstr.replace('<Pixel_Level_QI',                 '<L1C_Pixel_Level_QI')

        outstr = outstr.replace('</TILE_ID',                       '</TILE_ID_2A')
        outstr = outstr.replace('</DATASTRIP_ID',                  '</DATASTRIP_ID_2A')
        outstr = outstr.replace('</Product_Info',                  '</L2A_Product_Info')
        outstr = outstr.replace('</Product_Organisation',          '</L2A_Product_Organisation')
        outstr = outstr.replace('</Product_Image_Characteristics', '</L2A_Product_Image_Characteristics')
        outstr = outstr.replace('</Pixel_Level_QI',                '</L1C_Pixel_Level_QI')
        
        if self._product == 'T2A':
            outstr = outstr.replace('Image_Content_QI>', 'L1C_Image_Content_QI>')
        if self._product == 'UP2A':
            outstr = outstr.replace('QUANTIFICATION_VALUE', 'L1C_L2A_Quantification_Values_List')
            outstr = outstr.replace('</n1:Auxiliary_Data_Info>', '</n1:Auxiliary_Data_Info>\n'\
                                '<n1:L2A_Auxiliary_Data_Info/>')
            outstr = outstr.replace('</n1:Quality_Indicators_Info>', '</n1:Quality_Indicators_Info>\n'\
                                '<n1:L2A_Quality_Indicators_Info/>')
            
        l.acquire()
        try:
            outfile = codecs.open(self._xmlFn, 'w', 'utf-8')
            outfile.write('<?xml version="1.0"  encoding="UTF-8"?>\n')
            outfile.write(outstr)
            outfile.close()
        except:
            stdoutWrite('Error in writing file: %s\n' % self._xmlFn)
            self._logger.error('error in writing file: %s\n' % self._xmlFn)
        finally:
            l.release()
 
        return self.setRoot()
