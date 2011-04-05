#!/usr/bin/env python

import getopt
import os
import sys

class GeniusError(Exception):
    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return str(self.msg)
    
def addDLPath(path):
    if os.sysname == 'Windows':
        pass
    else:
        if os.environ.has_key(os.varDLPath):
            os.environ[os.varDLPath] += path
        else:
            os.environ[os.varDLPath] = path
    
def loadGenius():
    if sys.platform=='win32':
        modulePath = '/lib/genius.pyd'
    else:
        modulePath = '/lib/genius.so'

    searchPath = [ os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + '/..'),
                  '.', '..']
    if os.environ.has_key('GENIUS_DIR'):
        searchPath.insert(0, os.environ['GENIUS_DIR'])
        
    GeniusDir = None
    for path in searchPath:
        if os.path.exists(path+modulePath):
            GeniusDir = os.path.abspath(path)
            break
    
    if GeniusDir==None:
        raise GeniusError('Genius library binaries are not found.')
            
    sys.path.insert(0, GeniusDir+'/lib')
    
    global genius
    try:
        flagsave = sys.getdlopenflags()
        if sys.platform=='linux2':
            import DLFCN
            sys.setdlopenflags(DLFCN.RTLD_NOW|DLFCN.RTLD_GLOBAL)
        import genius
        sys.setdlopenflags(flagsave)
    except ImportError:
        raise GeniusError('Failed loading Genius binaries.')

    genius.Genius.set_genius_dir(GeniusDir)
    genius.baseDir = GeniusDir

def usage():
    sys.stderr.write('Usage:')
    sys.stderr.write('  genius -i <input file>\n')

def print_info(simuSys):    
    nRegion = simuSys.n_regions()
    for i in xrange(0,nRegion):
        region = simuSys.region(i)
        print i, region.name()
        
    nBC = simuSys.n_bcs() 
    for i in xrange(0, nBC):
        bc = simuSys.get_bc(i)
        print i, bc.label(), bc.bc_type()

def main():
    try:
        optsGenius, args = getopt.getopt(sys.argv[1:], 'i:')
    except getopt.GetoptError, err:
        print(err)
        usage()
        exit(-1)
        
    inputFile=None
    for opt,val in optsGenius:
        if opt=='-i':
            inputFile = val
    if inputFile==None:
        sys.stderr.write('Error: Input file not specified.\n')
        usage()
        sys.exit(-1)

    try:
        loadGenius()
    except GeniusError, e:
        sys.stderr.write('Error: '+str(e)+'\n')
        sys.exit(-1)
        
    logfilename = 'genius.log.%d' % genius.Genius.processor_id()
    genius.Genius.add_log_filename('file', logfilename);
    if genius.Genius.processor_id()==0:
      genius.Genius.add_log_file('console', sys.stdout);
    LOG = genius.Genius.log

    materialFile = genius.baseDir + '/lib/material.def'
    genius.Material.init_material_define(materialFile)

    patternFile = genius.baseDir + "/lib/GeniusSyntax.xml"

    petscOpts = ['dummy', '-on_error_attach_debugger', 'noxterm']
    genius.Genius.init_processors(petscOpts)
    LOG("This is process #%d of %d\n" % (genius.Genius.processor_id(), genius.Genius.n_processors()))

    LOG(patternFile + "\n")
    pattern = genius.Parser.Pattern()
    if not os.path.exists(patternFile):
        LOG("Error: Pattern file %s does not exist." % patternFile)
        sys.exit(-1)
    if not pattern.get_from_XML(patternFile)==0:
        LOG("error reading pattern\n")
    LOG("pattern loaded\n")
    
    text = None
    if genius.Genius.processor_id()==0:
        if not os.path.exists(inputFile):
            LOG("Error: Input file %s does not exist." % inputFile)
            sys.exit(-1)
        ppFile = genius.Parser.FilePreProcess(inputFile).output()
        tmpfile = open(ppFile)
        text = tmpfile.read()
        text = genius.Parallel.broadcast(text, 0)
        tmpfile.close()
        os.remove(ppFile)
    else:
        text = genius.Parallel.broadcast(text, 0)
    localFile = inputFile+".tmp"+str(genius.Genius.processor_id())
    tmpfile = open(localFile, "w")
    tmpfile.write(text)
    tmpfile.close()

    inputDeck = genius.Parser.InputParser(pattern)
    if not inputDeck.read_card_file(localFile)==0:
        LOG("Input file parsing error.\n")
        sys.exit(-1)
    os.remove(localFile)
    
    control = genius.SolverControl()
    control.setDecks(inputDeck)
    control.mainloop()
    
    print_info(control.system())
 
    if genius.Genius.processor_id()==0:
      genius.Genius.remove_log_file('console');
    genius.Genius.remove_log_file('file');
    genius.Genius.clean_processors()

#if Genius.processor_id()==0 :
#    msg = 13.4
#
#msg = Parallel.broadcast(msg, 0)
#Genius.log(str(msg)+"\n")    

if __name__=='__main__':
    main()
    sys.exit(0)
