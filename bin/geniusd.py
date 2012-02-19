#!/usr/bin/env python
'''
Genius Daemon with HTTP frontend
'''

import os
import sys
import threading
import time
import getopt
import random

from GeniusCtrl import *
from HTTPFE import *

class NullStream:
  def close(self): pass
  def flush(self): pass
  def write(self, str): pass
  def writelines(self, seq): pass

def usage():
    print '''geniusd.py -a <server name>'''

def main():

    # parse options
    try:
        cmdOpts, args = getopt.getopt(sys.argv[1:], 'ha:')
    except getopt.GetoptError, err:
        print(err)
        usage()
        exit(-1)
        
    serverName=None
    for opt,val in cmdOpts:
        if opt=='-a':
            serverName = val
        elif opt=='-h':
            usage()
            exit(0)
    if not serverName:
        serverName = 'Genius-' + str(random.randint(100000,999999))

    # initialize Genius
    try:
        GeniusCtrl.initialize(serverName)
    except GeniusError, e:
        sys.stderr.write(e.message+'\n')
        sys.exit(-1)

    # if we don't have a console (in pythonw.exe)
    # don't write to stdout
    if not sys.stdout.isatty():  
      sys.stdout = NullStream()
      sys.stderr = sys.stdout

    # Genius command queue and worker thread
    cmdQueue = CmdQueue()
    GeniusCtrl.setCmdQueue(cmdQueue)

    worker = threading.Thread(target=GeniusCtrl.run)
    worker.start()

    if GeniusCtrl.Genius.processor_id()==0:
        # start HTTP front-end at process 0
        fe = None
        try:
            fe = HTTPFrontEnd(GeniusCtrl, cmdQueue, serverName)
        except Exception, e:
            GeniusCtrl.log(str(e)+"\n")
            GeniusCtrl.log("Failed to start HTTP front end, shutting down\n")
            cmdQueue.addCmd('quit')
        if fe:
            fe.run()

    if not GeniusCtrl.Genius.error():
        # HTTP front-end finished, there is no error
        # a quit command should have been sent to the Genius command queue
        # wait for the Genius worker thread to finish
        worker.join()
        GeniusCtrl.finalize()
    else:
        exit(-1) # some error occurred, simply quit.

if __name__=='__main__':
    if sys.platform=='win32':
        # append self path to PATH
        selfdir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
        envpath = os.environ['PATH']
        envpath+=';'+selfdir
        os.environ['PATH']=envpath
    else:
        sys.stdin.close()
        try: os.close(0)
        except OSError: pass

    main()

