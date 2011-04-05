__all__=['SimulationSystem']

import re
import string
from xml.dom import minidom


def toDomQV(dom, v):
    elem = None
    if isinstance(v, str):
        elem = dom.createElement('qvString')
        elem.appendChild( dom.createTextNode(v) )
    elif isinstance(v, int):
        elem = dom.createElement('qvInt')
        s = '%ld' % v
        elem.appendChild( dom.createTextNode(s) )
    elif isinstance(v, float):
        elem = dom.createElement('qvFloat')
        s = '%.16g' % v
        elem.appendChild( dom.createTextNode(s) )
    elif isinstance(v, list):
        elem = dom.createElement('qvList')
        for vv in v:
            elem.appendChild( toDomQV(dom, vv) )
    elif isinstance(v, dict):
        elem = dom.createElement('qvMap')
        for kk,vv in v.items():
            ePair = dom.createElement('qvPair')
            ePair.appendChild( toDomQV(dom, kk) )
            ePair.appendChild( toDomQV(dom, vv) )
            elem.appendChild(ePair)
    else:
        elem = dom.createElement('qvInvalid')
    return elem


class SimulationSystem(object):
    def __init__(self, sys):
        super(SimulationSystem, self).__init__()
        self.sys = sys
    
    def toText(self):
        result = ''
        nRegion = self.sys.n_regions()
        for i in xrange(0,nRegion):
            region = self.sys.region(i)
            result += region.name() + '\n'
            
        nBC = self.sys.n_bcs()
        for i in xrange(0,nBC):
            bc = self.sys.get_bc(i)
            result += bc.label() + '\n'
        
        return result  
        
    def toXML(self, syntax='VisualTCAD'):
        dom = minidom.getDOMImplementation().createDocument(None, "structure", None)
        eRoot = dom.documentElement
        
        eRegions = dom.createElement('regions')
        nRegion = self.sys.n_regions()
        for i in xrange(0,nRegion):
            region = self.sys.region(i)
            
            eRegion = dom.createElement('region')
            
            eLabel = dom.createElement('label')
            eLabel.appendChild( toDomQV(dom, region.name()) )
            eRegion.appendChild(eLabel)

            eType = dom.createElement('type')
            eType.appendChild( toDomQV(dom, region.type_name() ) )
            eRegion.appendChild(eType)
            
            eMaterial = dom.createElement('material')
            eMaterial.appendChild( toDomQV(dom, region.material()) )
            eRegion.appendChild(eMaterial)

            eRegions.appendChild(eRegion)
        
        eRoot.appendChild(eRegions)
        
        # { <electrode label>: [<bc label>, <bc label>]}
        electrodes = {}
        
        eBnds = dom.createElement('boundaries')
        nBC = self.sys.n_bcs()
        for i in xrange(0, nBC):
            bc = self.sys.get_bc(i)
            
            eBnd = dom.createElement('boundary')
            
            eLabel = dom.createElement('label')
            eLabel.appendChild( toDomQV(dom, bc.label() ) )
            eBnd.appendChild(eLabel)
            
            eType = dom.createElement('type')
            eType.appendChild( toDomQV(dom, bc.bc_type_name() ) )
            eBnd.appendChild(eType)
            
            if bc.is_electrode():
                elecLabel = bc.electrode_label()
                if not len(elecLabel)>0:
                    elecLabel = bc.label()
                    
                # if electrode not existing yet, create it
                if not electrodes.has_key(elecLabel):
                    electrodes[elecLabel] = []
                # append bc to electrode
                electrodes[elecLabel].append(bc.label())
                
            eBnds.appendChild(eBnd)
        eRoot.appendChild(eBnds)

        eElecs = dom.createElement('default-contacts')
        for elecLabel, bcs in electrodes.items():
            eElec = dom.createElement('contact')
            
            eLabel = dom.createElement('label')
            eLabel.appendChild( toDomQV(dom, elecLabel) )
            eElec.appendChild(eLabel)
            
            eBnds = dom.createElement('boundaries')
            eBnds.appendChild( toDomQV(dom, bcs) )
            eElec.appendChild(eBnds)
            
            eElecs.appendChild(eElec)
        eRoot.appendChild(eElecs)

        #return dom.toprettyxml(indent='  ')
        return dom.toxml()

