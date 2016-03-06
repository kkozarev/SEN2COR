import os, sys
import re
import datetime as dt
from hashlib import md5
import lxml.etree as ET


class L2A_Manifest():
    def __init__(self, config):
        self.configDir = config.configDir
        self.DICO_IDS = { 
        # PRODUCT ROOT -->
            "S2_Level-2A_Product_Metadata": "S2A_USER_MTD_SAFL2A_PDMC",
            "INSPIRE_Metadata": "INSPIRE.xml",
            "HTML_Presentation": "UserProduct_index.html",
            "HTML_Presentation_Stylesheet": "UserProduct_index.xsl",
        # DATASTRIP -->
            "S2_Level-2A_Datastrip([1-9]?[0-9]{1,6})_Metadata": "S2A_USER_MTD_L2A_DS_",
        # DATASTRIP/QI_DATA -->  
            "Sensor_OLQC_Report_Datastrip([1-9]?[0-9]{1,4})_InformationData":   "S2A_OPER_MSI_L1C_DS_.*SENSOR_QUALITY",
            "Geometric_OLQC_Report_Datastrip([1-9]?[0-9]{1,4})_InformationData": "S2A_OPER_MSI_L1C_DS_.*GEOMETRIC_QUALITY",
            "General_OLQC_Report_Datastrip([1-9]?[0-9]{1,4})_InformationData": "S2A_OPER_MSI_L1C_DS_.*GENERAL_QUALITY",
            "Format_OLQC_Report_Datastrip([1-9]?[0-9]{1,4})_InformationData": "S2A_OPER_MSI_L1C_DS_.*FORMAT_CORRECTNESS",
            "Radiometric_OLQC_Report_Datastrip([1-9]?[0-9]{1,4})_InformationData": "S2A_OPER_MSI_L1C_DS_.*RADIOMETRIC_QUALITY",
        # GRANULE-->
            "S2_Level-2A_Tile([1-9]?[0-9]{1,6})_Metadata": "S2A_USER_MTD_L2A_TL__........T......_......._......\.xml",
            #"S2_Level-2A_Tile([1-9]?[0-9]{1,6})_Metadata": "S2A_USER_MTD_L2A_TL__20151021T153950_A001515_T53LNA.xml",
        # GRANULE/IMG_DATA -->
            "IMG_DATA_Band_10m_[1-4]_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_MSI_L2A_TL_.*B[0-1][0-9]_10m.jp2",
            "IMG_DATA_Band_20m_[1-9]_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_MSI_L2A_TL_.*B[0-1][0-9]_20m.jp2",
            "IMG_DATA_Band_60m_([1-9]|1[0-2])_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_MSI_L2A_TL_.*B[0-1][0-9]_60m.jp2",
            "AOT_DATA_Band_10m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_AOT_L2A_TL_.*_10m.jp2",
            "AOT_DATA_Band_20m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_AOT_L2A_TL_.*_20m.jp2",
            "AOT_DATA_Band_60m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_AOT_L2A_TL_.*_60m.jp2",
            "WVP_DATA_Band_10m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_WVP_L2A_TL_.*_10m.jp2",
            "WVP_DATA_Band_20m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_WVP_L2A_TL_.*_20m.jp2",
            "WVP_DATA_Band_60m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_WVP_L2A_TL_.*_60m.jp2",
            "VIS_DATA_Band_20m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_VIS_L2A_TL_.*_20m.jp2",
            "VIS_DATA_Band_60m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_VIS_L2A_TL_.*_60m.jp2",
            "SCL_DATA_Band_20m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_SCL_L2A_TL_.*_20m.jp2",
            "SCL_DATA_Band_60m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_SCL_L2A_TL_.*_60m.jp2",          
        # GRANULE/QI_DATA -->
            "Sensor_OLQC_Report_Tile([1-9]?[0-9]{1,4})_InformationData": "S2A_USER_MSI_L2A_TL_.*SENSOR_QUALITY", 
            "Geometric_OLQC_Report_Tile([1-9]?[0-9]{1,4})_InformationData": "S2A_USER_MSI_L2A_TL_.*GEOMETRIC_QUALITY",
            "General_OLQC_Report_Tile([1-9]?[0-9]{1,4})_InformationData": "S2A_USER_MSI_L2A_TL_.*GENERAL_QUALITY",
            "Format_OLQC_Report_Tile([1-9]?[0-9]{1,4})_InformationData": "S2A_USER_MSI_L2A_TL_.*FORMAT_CORRECTNESS",  
            "DefectivePixelsMask_Band_([1-9]|1[0-3])_Tile([1-9]?[0-9]{1,6})_Data": "S2A_OPER_MSK_DEFECT_",
            "DetectorFootprintMask_Band_([1-9]|1[0-3])_Tile([1-9]?[0-9]{1,6})_Data": "S2A_OPER_MSK_DETFOO_",
            "NodataPixelsMask_Band_([1-9]|1[0-3])_Tile([1-9]?[0-9]{1,6})_Data": "S2A_OPER_MSK_NODATA_",
            "SaturatedPixelsMask_Band_([1-9]|1[0-3])_Tile([1-9]?[0-9]{1,6})_Data": "S2A_OPER_MSK_SATURA_",
            "TechQualityMask_Band_([1-9]|1[0-3])_Tile([1-9]?[0-9]{1,6})_Data": "S2A_OPER_MSK_TECQUA_",
            "S2_Level-1C_Preview_Tile([1-9]?[0-9]{1,6})_Data": "S2A_OPER_PVI_L1C_TL_.*.jp2",
            "CloudMask_Level-1C_Tile([1-9]?[0-9]{1,6})_Data": "S2A_OPER_MSK_CLOUDS_.*.gml", # discard .gfs
            "CloudMask_Band_20m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_CLD_L2A_TL_.*_20m.jp2",
            "CloudMask_Band_60m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_CLD_L2A_TL_.*_60m.jp2",
            "Preview_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_PVI_L2A_TL_.*.jp2",
            "SnowMask_20m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_SNW_L2A_TL_.*_20m.jp2",
            "SnowMask_60m_Tile([1-9]?[0-9]{1,6})_Data": "S2A_USER_SNW_L2A_TL_.*_60m.jp2",
            "Processing_Report_Tile([1-9]?[0-9]{1,6})_InformationData": "S2A_USER_MTD_L2A_TL_.*Report",
        # GRANULE/AUX_DATA -->
            "S2_GIP_L2A_Tile([1-9]?[0-9]{1,6})_Metadata": "S2A_USER_GIP_L2A_TL_" }
        
        
        self.NSMAP={'xfdu' :"urn:ccsds:schema:xfdu:1",
                'safe':"http://www.esa.int/safe/sentinel/1.1",
                'gml' :"http://www.opengis.net/gml"}
    

    
    def match(self, pattern, stringlist):
        """
        test a list of string against a regex pattern
        """
        results=[]
        prog = re.compile(pattern)
        for strin in stringlist:
            #print pattern, string
            if prog.search(strin) is not None:
                results.append(strin)
                
        return results
        
    def listToIds(self, dataObjectSection, s2alist):
        """
        for each id in DICO_IDS, retrieve which filepaths are matching
        """
        ids={}
        for key, pattern in self.DICO_IDS.iteritems():
            results = self.match(pattern, s2alist)
            
            iTile=0
            iBand=1
            currenttile=None
            for result in results:
                if 'GRANULE' in result or 'DATASTRIP' in result:
                    tilename=result.split('GRANULE/')[-1].split('DATASTRIP/')[-1].split('/')[0]
                    if currenttile != tilename:
                        iTile += 1
                        currenttile = tilename
                        iBand = 1

                a = key.find('Band')
                id = key
                if a!= -1:
                    id = id.replace(key[a:], 'Band_' + str(iBand) + '_Tile' + str(iTile) + '_Data')            
                
                ids[result] = id
                iBand+=1
        return ids
    
    
    
    def getMimeType(self, href):
        """
        Returns mime type function of file extension
        """
        dicoMimeType = { 'xml': "text/xml",
                        'gml': "application/xml",
                        'html': "octet",
                        'xsl': "octet",
                        'jp2': "application/octet-stream" }
        return dicoMimeType[href.split(".")[-1]]


    def getChecksum(self, href, l2productpath):
        """
        Computes MD5 checksum
        """
        filename = os.path.join(l2productpath, href)
        if os.path.exists(filename):
            if href.split(".")[-1] == 'jp2':
                return md5(open(filename, 'rb').read()).hexdigest()
            else:
                return md5(open(filename, 'r').read()).hexdigest()
        else:
            return None
    
    def addDataObject(self, parent, href, id, l2productpath):
        """
        add a data object element with its subelements
        """
        dataObject = ET.SubElement(parent, "dataObject")
        byteStream = ET.SubElement(dataObject, "byteStream")
        fileLocation = ET.SubElement(byteStream, "fileLocation")
        checksum = ET.SubElement(byteStream, "checksum")
        dataObject.attrib["ID"]  = id
        byteStream.attrib["mimeType"] = self.getMimeType(href)

        l2productpath_ = unicode(os.path.join(l2productpath, href))
        if os.path.exists(l2productpath_):
            byteStream.attrib["size"] = str(l2productpath_)
        else:
            byteStream.attrib["size"] = str(None)            
        
        fileLocation.attrib["locatorType"] = "URL"
        fileLocation.attrib["href"] = href
        checksum.attrib["checksumName"] = "MD5"
        checksum.text = self.getChecksum(href, l2productpath_)

    
    def addContentUnit(self, parent, ID=None, unitType=None, textInfo=None, pdiID=None, dmdID=None, dataObjectID=None):
        """
        Add content unit element to given parent element
        """
        contentUnit = ET.SubElement(parent, '{'+self.NSMAP['xfdu']+'}'+'contentUnit')
        if ID is not None:
            contentUnit.attrib['ID'] = ID
        if unitType is not None:
            contentUnit.attrib['unitType'] = unitType
        if textInfo is not None:
            contentUnit.attrib['textInfo'] = textInfo
        if pdiID is not None:
            contentUnit.attrib['pdiID'] = pdiID
        if dataObjectID is not None:
            dataObjectPointer = ET.SubElement(contentUnit, "dataObjectPointer")
            dataObjectPointer.attrib['dataObjectID'] = dataObjectID
        return contentUnit
    
    def addMetadataObject(self, parent, ID=None, classification=None, category=None, dataObjectID=None):
        """
        Add metadata object element to given parent element
        """
        metadataObject = ET.SubElement(parent, "metadataObject")
        if ID is not None:
            metadataObject.attrib['ID'] = ID
        if classification is not None:
            metadataObject.attrib['classification'] = classification
        if category is not None:
            metadataObject.attrib['category'] = category
        if dataObjectID is not None:
            dataObjectPointer = ET.SubElement(metadataObject, "dataObjectPointer")
            dataObjectPointer.attrib['dataObjectID'] = dataObjectID
        return metadataObject
    
    def splitPath(self, filepath):
        """
        Returns datastrip, granule and sub directory names
        from a L2A product relative file path
        """
        tilename = None
        stripname = None
        subdir = None
        splitted = os.path.dirname(filepath).split('/')
        if len(splitted) > 1:
            # tile or datastrip
            if splitted[0] == 'GRANULE':
                tilename = splitted[1]
            elif splitted[0] == 'DATASTRIP':
                stripname = splitted[1]
            # subdir (IMG_DATA, QI_DATA ...)
            if len(splitted) > 2:
                subdir = splitted[2]
        return (tilename, stripname, subdir)
    
    def updateProcessingelement(self, metadataObject):
        """
        Update processing element from L1C values to L2A
        """
        processingelement = metadataObject.find("metadataWrap").find("xmlData").find("{"+self.NSMAP["safe"]+"}processing")
        processingelement.attrib["name"] = processingelement.attrib["name"].replace("1C","2A")
        processingelement.attrib["start"] = dt.datetime.utcnow().isoformat()+"Z"
    
    
    def append(self, parent, child):
        """
        append a child by recopy
        """
        newchild = ET.SubElement(parent, child.tag)
        try:
            child.text.strip()
            newchild.text = child.text
        except:
            pass 
        for key in child.attrib:
            newchild.attrib[key] = child.attrib[key]
        for childOfChild in child.getchildren():
            self.append(newchild, childOfChild)
        return newchild
    
    def generate(self, l2productpath, safeL1cFn):
        """
        generates the L2A safe manifest
        by parsing the LA product 
        and copying some info from the L1C safe manifest
        """
    
        # == INITIALIZE XML == #
        # namespaces and root tree
        for key in self.NSMAP.keys():
            ET.register_namespace(key, self.NSMAP[key]) #Registers a namespace prefix. The registry is global, and any existing mapping for either the given prefix or the namespace URI will be removed. prefix is a namespace prefix. uri is a namespace uri. Tags and attributes in this namespace will be serialized with the given prefix, if at all possible.
        root = ET.Element('{'+self.NSMAP['xfdu']+'}'+'XFDU', nsmap=self.NSMAP) #Element class. This class defines the Element interface, and provides a reference implementation of this interface.
        root.attrib['version']="esa/safe/sentinel/1.1/sentinel-2/msi/archive_l2a_user_product"
    
        # Create the 3 main sections
        informationPackageMap = ET.SubElement(root, "informationPackageMap")
        metadataSection = ET.SubElement(root, "metadataSection")
        dataObjectSection = ET.SubElement(root, "dataObjectSection")
        
    
        # ======== DATA OBJECT SECTION ===========
        s2alist = []
        # 1) parse product for all files
        for rootdir, dirnames, filenames in os.walk(l2productpath): #check all the dirs and filenames contained in the rootdir
            reldir = os.path.relpath(rootdir, l2productpath) #return the relative path from l2productpath to rootdir
            for filename in filenames:
                s2alist.append(os.path.join(reldir, filename))
        s2alist.sort() # s2list is a list with all the files (with relative paths)
        
        # 2) parse list of files and guess a data object id 
        ids = self.listToIds(dataObjectSection, s2alist)
        
        # 3) creates the data objects element
        for filepath in s2alist:
            if filepath in ids:
                self.addDataObject(dataObjectSection, filepath, ids[filepath], l2productpath)
        
        # ======== INFORMATION PACKAGE MAP SECTION ===========
        
        # 1) Product root (generic) + tiles + strips main elements
        productroot = self.addContentUnit(informationPackageMap, unitType="Product_Level-2A", textInfo="SENTINEL-2 MSI Level-2A User Product", pdiID="processing", dmdID="acquisitionPeriod platform")
        tileselement = self.addContentUnit(productroot, ID="Tiles", textInfo="Tiles Container")
        stripselement = self.addContentUnit(productroot, ID="Datastrips", textInfo="Datastrips Container")
    
        # 2) Create elements for each granule, each datastrip, including their sub folders (IMG_DATA, QI_DATA...)
        tilemap={}
        stripmap={}
        for filepath in s2alist:
            if filepath in ids:
                (tilename, stripname, subdir) = self.splitPath(filepath)
                ## GRANULES
                if tilename is not None:
                    # create elements for new tile
                    if tilename not in tilemap:
                        tileid = len(tilemap.keys())+1
                        tilemap[tilename] = self.addContentUnit(tileselement, ID="Tile%d"%(tileid), textInfo="Granule Container")
                        self.addContentUnit(tilemap[tilename], ID="QI_DATA_Tile%d"%(tileid), textInfo="Quality Control Container")
                        self.addContentUnit(tilemap[tilename], ID="IMG_DATA_Tile%d"%(tileid), textInfo="Image Data Container")
                        self.addContentUnit(tilemap[tilename], ID="AUX_DATA_Tile%d"%(tileid), textInfo="Auxiliary Data container")
                ## DATASTRIPS
                elif stripname is not None:
                    if stripname not in stripmap:
                        # create elements for new datastrip
                        stripid = len(stripmap.keys())+1
                        stripmap[stripname] = self.addContentUnit(stripselement, ID="Datastrip%d"%(stripid), textInfo="Granule Container")
                        self.addContentUnit(stripmap[stripname], ID="QI_DATA_Datastrip%d"%(stripid), textInfo="Quality Control Container")
    
        # 3) Then create a element for each element below the correct parent
        for filepath in s2alist:
            if filepath in ids:
                id = ids[filepath]
                (tilename, stripname, subdir) = self.splitPath(filepath)
                ## GRANULES
                if tilename is not None:
                    # then add new content unit
                    if subdir is None:
                        self.addContentUnit(tilemap[tilename], ID=id+"_Unit", unitType="Metadata Unit", dataObjectID=id)
                    else:
                        if subdir == "QI_DATA":
                            self.addContentUnit(tilemap[tilename].getchildren()[0], ID=id+"_Unit", unitType="Measurement Data Unit", textInfo="Quality Control Container", dataObjectID=id)
                        elif subdir == "IMG_DATA":
                            self.addContentUnit(tilemap[tilename].getchildren()[1], ID=id+"_Unit", unitType="Measurement Data Unit", textInfo="Image Data Container", dataObjectID=id)
                        elif subdir == "AUX_DATA":
                            self.addContentUnit(tilemap[tilename].getchildren()[2], ID=id+"_Unit", unitType="Measurement Data Unit", textInfo="Auxiliary Data Container", dataObjectID=id)
                ## DATASTRIPS
                elif stripname is not None:
                    if subdir is None:
                        self.addContentUnit(stripmap[stripname], ID=id+"_Unit", unitType="Metadata Unit", dataObjectID=id)
                    else:
                        if subdir == "QI_DATA":
                            self.addContentUnit(stripmap[stripname].getchildren()[0], ID=id+"_Unit", unitType="Metadata Unit", dataObjectID=id)
                else:
                    self.addContentUnit(productroot, ID=id+"_Unit", unitType="Metadata Unit", dataObjectID=id)
        
        
        # ======== METADATA OBJECT SECTION ===========
        
        # 1) Recopy 1st part from L1C manifest 
            ### root ###
        treeL1C = ET.parse(safeL1cFn)
        manifestL1C = treeL1C.getroot()
        # metadataSection
        metadataSectionL1C = manifestL1C.find('metadataSection')
        L1C_KEEP_LIST = ["acquisitionPeriod", "measurementOrbitReference", "platform", "processing" , "measurementFrameSet"]
        for childL1C in  metadataSectionL1C.getchildren():
            if "ID" in childL1C.attrib.keys() and childL1C.attrib["ID"] in L1C_KEEP_LIST:
                #metadataSection.append(childL1C)
                childL2A = self.append(metadataSection, childL1C)
                if childL2A.attrib["ID"] == "processing":
                    self.updateProcessingelement(childL2A) 
        
        # 2) Add metadata elements (only) in 2d part
        for filepath in s2alist:
            if filepath in ids:
                id = ids[filepath]
                if "Metadata" in id:
                    self.addMetadataObject(metadataSection, ID=id+"_Information", classification="DESCRIPTION", category="DMD", dataObjectID=id)
        
        # == WRITE XML == #
        #print ET.tostring(root, pretty_print=True)
        fn = os.path.join(l2productpath, 'L2A_Manifest.xml')
        with open(fn,"w") as o:
            o.write(ET.tostring(root, encoding='UTF-8', pretty_print=True, xml_declaration=True))
        #os.system("xmllint -format toto.tmp > out.xml")
    

    
