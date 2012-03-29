__all__=['GeniusError', 'CmdQueue', 'GeniusCommand', 'CmdRunDeck', 'CmdQueryStruct', 'GeniusCtrl']

import os
import sys
import threading
import time
import re
import string
import GeniusLib
import sip

class GeniusError(Exception):
    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return str(self.msg)

class GeniusCommand(object):
    opcode = 'nop'
    def __init__(self):
        self.finished = threading.Event()
    
    def __str__(self):
        return self.toStr()
    
    def toStr(self):
        return 'nop'
    
    def fromStr(self, str):
        if not re.match('^nop\s*$', str, re.I):
            raise ValueError 
        return self
        
    def do(self):
        return True


def _mkCmdFromStr(cmdStr):
    cmd=None
    for cc in ValidCommands:
        if cmdStr.startswith(cc.opcode):
            cmd = cc().fromStr(cmdStr)
    if cmd==None:
        cmd=GeniusCommand() # nop
        
    return cmd


class CmdQuit(GeniusCommand):
    opcode = 'quit'
    def __init__(self):
        super(CmdQuit,self).__init__()
    
    def toStr(self):
        return 'quit'
    
    def fromStr(self, str):
        if not re.match('^quit\s*$', str, re.I):
            raise ValueError
        return self
        
    def do(self):
        return True

class CmdReset(GeniusCommand):
    opcode = 'reset'
    def __init__(self):
        super(CmdReset, self).__init__()
        
    def toStr(self):
        return 'reset'
    
    def fromStr(self, str):
        if not re.match('^reset\s*$', str, re.I):
            raise ValueError
        return self
    
    def do(self):
        GeniusCtrl.Genius.log("Simulation System Reset.\n")
        GeniusCtrl.SolverControl.reset_simulation_system();

class CmdChdir(GeniusCommand):
    opcode = 'chdir'
    def __init__(self, dirName=None):
        super(CmdChdir, self).__init__()
        self.dirName = dirName

    def toStr(self):
        return 'chdir %s' % self.dirName

    def fromStr(self, str):
        m = re.match('chdir\s+(.*)', str, re.I)
        if m==None:
            raise ValueError

        self.dirName = m.group(1).rstrip()
        return self

    def do(self):
        GeniusCtrl.Genius.log('Changing working directory to %s\n' % self.dirName)
        os.chdir(self.dirName)

class CmdQueryStruct(GeniusCommand):
    opcode = 'query-struct'
    def __init__(self, queryStr=None):
        super(CmdQueryStruct, self).__init__()
        self.queryStr = queryStr
        self.result = None;

    def toStr(self):
        return 'query-struct %s' % self.queryStr

    def fromStr(self, str):
        m = re.match('query-struct\s+(.*)', str, re.I)
        if m==None:
            raise ValueError

        self.queryStr = m.group(1).rstrip()
        self.result = None;
        return self

    def do(self):
        GeniusCtrl.Genius.log('Query device structure.')
        control = GeniusCtrl.SolverControl
        _simuSys = control.system()
        sip.transferto(_simuSys, None) ## temporary workaround
        simuSys = GeniusLib.SimulationSystem(_simuSys)
        if simuSys:
            if re.match('xml', self.queryStr, re.I):
                self.result = simuSys.toXML()
            else:
                self.result = simuSys.toText()

class CmdRunPython(GeniusCommand):
    opcode = 'run-python'
    def __init__(self, code=None):
        super(CmdRunPython,self).__init__()
        self.code = code
        self.result = None
        
    def toStr(self):
        return 'run-python \n%s' % self.code
    
    def fromStr(self, str):
        ilb = str.find('\n')
        line1 = str[0:ilb]
        
        if not line1.startswith('run-python'):
            return ValueError
        if not ilb+1<len(str):
            return ValueError
        
        self.code = str[ilb+1:]
        self.result = None
        
    def do(self):
        exec self.code
    
class CmdRunDeck(GeniusCommand):
    opcode = 'run-deck'
    def __init__(self, fname=None):
        super(CmdRunDeck,self).__init__()
        self.fname = fname
        
    def toStr(self):
        return 'run-deck %s' % self.fname
    
    def fromStr(self, str):
        m = re.match('^run-deck\s+(.*)$', str, re.I )
        if m==None:
            raise ValueError
            
        self.fname = m.group(1).rstrip()
        return self
        
    def do(self):
        if self.fname==None:
            return False
            
        text = None
        fname=self.fname
        
        if GeniusCtrl.Genius.processor_id()==0:
            if not os.path.exists(fname):
                raise GeniusError("Error: Input file %s does not exist." % fname)
            ppFile = GeniusCtrl.Parser.FilePreProcess(fname).output()
            tmpfile = open(ppFile)
            text = tmpfile.read()
            text = GeniusCtrl.Parallel.broadcast(text, 0)
            tmpfile.close()
            os.remove(ppFile)
        else:
            text = GeniusCtrl.Parallel.broadcast(text, 0)
        localFile = self.fname+".tmp"+str(GeniusCtrl.Genius.processor_id())
        tmpfile = open(localFile, "w")
        tmpfile.write(text)
        tmpfile.close()
    
        inputDeck = GeniusCtrl.Parser.InputParser(GeniusCtrl.pattern)
        if not inputDeck.read_card_file(localFile)==0:
            os.remove(localFile)
            raise GeniusError("Input file parsing error.")
        os.remove(localFile)
        
        control = GeniusCtrl.SolverControl
        control.setDecks(inputDeck)
        
        sol_fname = fname+'.sol'
        if os.path.exists(sol_fname):
          os.remove(sol_fname);
        control.setSolutionFile(sol_fname)
        control.mainloop()

        return True

        
class CmdQueue:
    ''' Command queue. Only the first processor maintains this queue, 
        which accepts commands from the (HTTP) front-end.
        An internal lock is maintained to avoid race condition
    '''
    def __init__(self):
        self.queue = []
        self.lock = threading.Lock()
        
    def getCmd(self):
        ''' Remove and return the first item in queue, returns None if queue is empty
        '''
        self.lock.acquire()
        tmp = None
        if len(self.queue)>0:
            tmp = self.queue[0]
            self.queue = self.queue[1:]
        self.lock.release()
        return tmp
    
    def addCmd(self, cmd_or_str):
        ''' Append a command to queue
        @param: cmd_or_str     GeniusCommand object or string representation of the command
        @return: a reference to the added command, useful when input is a string
        '''
        self.lock.acquire()
        if isinstance(cmd_or_str, GeniusCommand):
            cmd = cmd_or_str
        else:
            cmd = _mkCmdFromStr(cmd_or_str)
        self.queue.append(cmd)
        self.lock.release()
        return cmd

class GeniusCtrlClass:
    def __init__(self):
        self._g = {}
        self.baseDir = None
        self.inputDeck = None
        self.ready = False
        self.error = False
        self.running = False
        self._cmdQueue = None
        
        self.patternFile = None
        
    def __getattr__(self, name):
        if name=='SolverControl':
            sKey = 'SolverControlSingleton'
            if not self._g.__dict__.has_key(sKey):
                self._g.__dict__[sKey] = self._g.SolverControl()
            return self._g.__dict__[sKey]
        if self._g and self._g.__dict__.has_key(name):
            return self._g.__dict__[name]
        else:
            raise AttributeError
            
    def _loadGenius(self):
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
        
        try:
            if sys.platform=='win32':
                import genius
            else:
                flagsave = sys.getdlopenflags()
                if sys.platform=='linux2':
                    import DLFCN
                    if hasattr(DLFCN, 'RTLD_DEEPBIND'):
                      sys.setdlopenflags(DLFCN.RTLD_NOW|DLFCN.RTLD_GLOBAL|DLFCN.RTLD_DEEPBIND)
                    else:
                      sys.setdlopenflags(DLFCN.RTLD_NOW|DLFCN.RTLD_GLOBAL|0x00008)
                import genius
                sys.setdlopenflags(flagsave)
        except ImportError, e:
            raise GeniusError('Failed loading Genius binaries. %s' % str(e))
    
        self.baseDir = GeniusDir
        genius.Genius.set_genius_dir(GeniusDir)
        self._g = genius

    def initialize(self, jobName='genius'):
        try:
            self._loadGenius()
        except GeniusError, e:
            sys.stderr.write('Error: '+str(e)+'\n')
            sys.exit(-1)
        
        self.jobName = jobName
        petscOpts = ['dummy', '-on_error_attach_debugger', 'noxterm']
        self.Genius.init_processors(petscOpts)
       
        self.logFilename = os.path.abspath('%s.log.%d' % (self.jobName, self.Genius.processor_id()))
        self.Genius.add_log_filename('file', self.logFilename)
        self.log = self.Genius.log
        if self.Genius.processor_id()==0 and sys.stdout.isatty():
          self.Genius.add_log_file('console', sys.stdout);
        
        self.log("Genius initialized successfully.\n")
        self.log("This is process #%d of %d\n" % (self.Genius.processor_id(), self.Genius.n_processors()))
                
        patternFile = GeniusCtrl.baseDir + "/lib/GeniusSyntax.xml"
    
        self.pattern = GeniusCtrl.Parser.Pattern()
        if not os.path.exists(patternFile):
            raise GeniusError("Pattern file %s does not exist." % patternFile)
            sys.exit(-1)
        if not self.pattern.get_from_XML(patternFile)==0:
            raise GeniusError("Failed reading syntax definition.\n")

        materialFile = GeniusCtrl.baseDir + '/lib/material.def'
        GeniusCtrl.Material.init_material_define(materialFile)

        self.ready = True

    def finalize(self):
        self.log("Shutting down Genius...\n")
         
        if self.Genius.processor_id()==0:
            self.Genius.remove_log_file('console');
        self.Genius.remove_log_file('file');
        
        self.Genius.clean_processors()
        self.inputDeck = None
        self.ready = False
        self.running = False

    def run(self):
        nopCnt=10
        while(True):
            time.sleep(0.5)
            cmd = None
            if self.Genius.processor_id()==0:
                # processor 0 gets the original command object
                if self._cmdQueue==None:
                    raise ValueError

                cmd = self._cmdQueue.getCmd()
                
                if cmd==None:
                    # broadcast an NOP command after 10 cycles, to keep everybody alive.
                    if nopCnt==0:
                        nopCnt=10
                        cmd = GeniusCommand() # NOP
                    else:
                        nopCnt -= 1
                        continue
                elif not isinstance(cmd, GeniusCommand):
                    raise TypeError

                self.Parallel.broadcast(cmd.toStr(),0)

            else:
                # other processors make a copy of the command from the string 
                # some properties (e.g. finished event) are lost
                cmdStr=None
                cmdStr=self.Parallel.broadcast(cmdStr,0)
                cmd = _mkCmdFromStr(cmdStr)

            if isinstance(cmd, CmdQuit):
                cmd.finished.set()
                break;
            else:
                try:
                    cmd.do()
                except Exception, e:
                    self.log(str(e)+'\n')
                finally:
                    cmd.finished.set()

    def setCmdQueue(self, queue):
        self._cmdQueue=queue
        
# module initialization
GeniusCtrl = GeniusCtrlClass()
ValidCommands = [GeniusCommand, CmdQuit, CmdRunDeck, CmdReset, CmdChdir, CmdQueryStruct]
