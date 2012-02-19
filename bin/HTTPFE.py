__all__=['HTTPFrontEnd']

from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer
import SocketServer
import re
import os
import socket
import sys
import threading
import traceback
import tempfile
import errno

def findFreePort(begin, end):
    '''
    Find an available TCP port between begin and end
        @param: begin   integer. lower bound for port number
        @param: end     integer. upper bound for port number
        @return: integer. an available tcp port number
    '''
    for port in xrange(begin, end+1):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.bind(('', port))
        except socket.error, (errno, msg):
            if errno==98 or errno==48 or errno==10048: # address already in use
                continue
            raise socket.error(errno, msg)
        
        # bind is successful
        s.close()
        return port


class GeniusRequest(BaseHTTPRequestHandler):
    server_version='GeniusHTTPService/0.1'
    def do_GET(self):
        cur_thread = threading.current_thread()
        self.server.dispatcher.dispatch(self)
    
    def do_POST(self):
        cur_thread = threading.current_thread()
        self.server.dispatcher.dispatch(self)

        
class Dispatcher:
    def __init__(self, map, errorHandler=None):
        '''
            @param: map [(re_pattern, handler)]
            @param: errHandler 
        '''
        self.map = map
        self.errorHandler = errorHandler

    def addURL(self, pattern, handler):
        self.map.append((pattern, handler))
        
    def dispatch(self, request):
        if not isinstance(request, GeniusRequest):
            raise TypeError
        
        found = False
        for pattern,handler in self.map:
            m = re.match(pattern, request.path, re.I)
            if not m==None:
                found = True
                try:
                    handler(m.groupdict(), request)
                except Exception, e:
                    if self.errorHandler:
                        self.errorHandler(e)
                finally:
                    break;
        
        if not found:
            request.send_error(404)


class HTTPFrontEnd(SocketServer.ThreadingMixIn, HTTPServer):
    def __init__(self, GeniusCtrl, cmdQueue, serverName):
        self.GeniusCtrl = GeniusCtrl
        self.cmdQueue = cmdQueue
        self.serverName = serverName

        port = findFreePort(11211,11220)
        HTTPServer.__init__(self, ('', port), GeniusRequest)
        self.GeniusCtrl.log('HTTP service %s started on port %d.\n' % (self.serverName, port))
        self.running = True
        self.quitEvent = None
        self.quitTimer = None

        self.socket_errors_to_ignore = []
        # Not all of these names will be defined for every platform.
        for _ in ("EPIPE", "ETIMEDOUT", "ECONNREFUSED", "ECONNRESET",
                  "EHOSTDOWN", "EHOSTUNREACH",
                  "WSAECONNABORTED", "WSAECONNREFUSED", "WSAECONNRESET",
                  "WSAENETRESET", "WSAETIMEDOUT"):
            if _ in dir(errno):
                self.socket_errors_to_ignore.append(getattr(errno, _))
        # de-dupe the list
        self.socket_errors_to_ignore = dict.fromkeys(self.socket_errors_to_ignore).keys()

        map = [
               ('/status/log', self.doLog),
               ('/status(/(?P<cmd>quit|exit|reset|chdir|tmpdir))', self.doStatusCmd),
               ('/status(/(?P<cmd>name|pwd|wait))?(\?(?P<arg>.*))?', self.doStatusQuery),
               ('/deck(/(?P<type>file|text))', self.doRunDeck),
               ('/solution', self.doSolution),
               ('/struct(/(?P<format>xml|text))', self.doStructQuery),
               ('/file/(?P<path>.*)', self.doFile)
               ]
        self.dispatcher = Dispatcher(map, self.errorHandler)

    def set_quit(self):
        self.running = False

    def get_request(self):
        while True:
            if not self.running:
                raise socket.error

            if self.quitEvent and self.quitEvent.isSet():
                if not self.quitTimer:
                    self.quitTimer = threading.Timer(10.0, self.set_quit)
                    self.quitTimer.start()

            self.socket.settimeout(1.0)
            try:
                sock, addr = self.socket.accept()
                sock.settimeout(None)
                return (sock, addr)
            except socket.timeout:
                self.socket.settimeout(None)

    def errorHandler(self, e):
        self.GeniusCtrl.log(str(e)+'\n')
        et, ev, tb = sys.exc_info()                                                 
        while tb :
            co = tb.tb_frame.f_code
            filename = str(co.co_filename)
            line_no =  str(traceback.tb_lineno(tb))
            self.GeniusCtrl.log(filename+':'+line_no+'\n')
            tb = tb.tb_next
        self.GeniusCtrl.log('Error: '+et+':'+ev+'\n')
        self.running=False
        
    def run(self):
        try:
            while(self.running):
                self.handle_request()
                  
        except KeyboardInterrupt:
            self.GeniusCtrl.log('Interrupted by user, shutting down HTTP service.\n')
            self.socket.close()
            self.cmdQueue.addCmd('quit')
    
    def _checkError(self, request):
        if self.GeniusCtrl.Genius.error():
            request.send_error(500)
            raise SystemError
        
    def doRunDeck(self, param, request):
        self._checkError(request)
        if not request.command=='POST':
            request.send_error(405) # method not allowed
            return

        type = param['type']
        if type=='file':
            fname = request.rfile.readline().strip()
            if not os.path.exists(fname):
                request.send_error(501)  # file not found
                return
            
            cmd = self.cmdQueue.addCmd('run-deck %s' % fname)
            request.send_response(201)
            request.send_header('Content-Type', 'text/plain')
            request.end_headers()
            request.wfile.write('Input file %s queued for execution.\n' % fname)
        else:
            request.send_error(501)
            return
                        
    
    def doStatusCmd(self, param, request):
        '''
            POST /status/quit   quit
            POST /status/exit   quit
            POST /status/reset  reset
            POST /status/chdir  change working dir
            POST /status/tmpdir create a temporary dir and use it as working dir
        '''
        self._checkError(request)
        if not request.command=='POST':
            request.send_error(405) # method not allowed
            return

        if param.has_key('cmd'):
            cmd = param['cmd']
            if cmd=='quit' or cmd=='exit':
                gCmd = self.cmdQueue.addCmd('quit')
                if self.quitEvent == None:
                    self.quitEvent = gCmd.finished

                request.send_response(200)
                request.send_header('Content-Type', 'text/plain')
                request.end_headers()
                request.wfile.write('Genius shut down is scheduled.\n')
                return
            elif cmd=='reset':
                gCmd = self.cmdQueue.addCmd('reset')
                #gCmd.finished.wait()
                request.send_response(200)
                request.send_header('Content-Type', 'text/plain')
                request.end_headers()
                request.wfile.write('Genius reset.\n')
                return
            elif cmd=='chdir':
                dirname = request.rfile.readline().strip()
                if not os.path.exists(dirname):
                    request.send_error(501)  # dir not found
                    return
                gCmd = self.cmdQueue.addCmd('chdir %s' % dirname)
                gCmd.finished.wait()
                request.send_response(200)
                request.send_header('Content-Type', 'text/plain')
                request.end_headers()
                request.wfile.write('ok.\n')
            elif cmd=='tmpdir':
                dirname = tempfile.mkdtemp()
                gCmd = self.cmdQueue.addCmd('chdir %s' % dirname)
                gCmd.finished.wait()
                request.send_response(200)
                request.send_header('Content-Type', 'text/plain')
                request.end_headers()
                request.wfile.write(dirname+'\n')
        else:
            request.send_error(400) # we shouldn't reach here
            return

    def doStatusQuery(self, param, request):
        '''
            GET /status/name    return the server name, which corresponds to the mpd job alias
            GET /status/pwd     return the current working directory
        '''
        self._checkError(request)

        if not request.command=='GET':
            request.send_error(405) # method not allowed
            return

        if param.has_key('cmd'):
            cmd = param['cmd']
            arg = param.get('arg', None)

            if cmd=='name':
                request.send_response(200)
                request.send_header('Content-Type', 'text/plain')
                request.end_headers()
                request.wfile.write(self.serverName)
                return
            elif cmd=='pwd':
                request.send_response(200)
                request.send_header('Content-Type', 'text/plain')
                request.end_headers()
                request.wfile.write(os.getcwd())
                return
            elif cmd=='wait':
                # /status/wait          timeout=0
                # /status/wait?5        timeout=5
                try:
                    timeout = int(arg)
                except:
                    timeout = None

                gCmd = self.cmdQueue.addCmd('nop')
                gCmd.finished.wait(timeout)
                while not gCmd.finished.is_set():
                    request.send_response(100)
                    gCmd.finished.wait(timeout)
                request.send_response(200)
                request.send_header('Content-Type', 'text/plain')
                request.end_headers()
                request.wfile.write('ok')
                return

        request.send_response(200)
        request.wfile.write('status\n')

    def readFile(self, fname, request, blocksize=1024):
        if not (os.path.exists(fname) and os.path.isfile(fname)):
            request.send_error(404) # file not found
            return
        if not os.access(fname, os.R_OK):
            request.send_error(403) # not access right
            return

        filesize = os.path.getsize(fname)

        start=0
        end=-1
        if request.headers.has_key('Range'):
            range = request.headers['Range']
            m = re.match('bytes\s*=\s*(\d+)?\s*-\s*(\d+)?\s*', range)
            if m:
                tmp = m.group(1)
                if tmp:
                    start = int(tmp)
                tmp = m.group(2)
                if tmp:
                    end = int(tmp)
        else:
            start=0
            end=filesize-1
        
        if end<0:
            readsize=blocksize
        else:
            readsize=end-start+1

        if start < filesize:
            f = open(fname, 'rb')
            f.seek(start)
            buf=f.read(readsize)
            f.close()
            end=start+len(buf)-1
           
            if start>0 or end<filesize-1:
                request.send_response(206) # partial content
                request.send_header('Content-Range', 'bytes %d-%d/*' % (start, end))
            else:
                request.send_response(200) # full content
            request.send_header('Content-Length', len(buf))
            request.end_headers()
            request.wfile.write(buf)
            return
        else:
            request.send_error(416) # requested range not satisfiable
            return

    def writeFile(self, fname, request, blocksize=1024):
        if os.path.exists(fname):
            if os.path.isdir(fname):
                request.send_error(403) # can not overwrite a dir
                return
            if not os.access(fname, os.W_OK):
                request.send_error(403) # can not overwrite a dir
                return
        else:
            dir = os.path.dirname(fname)
            if not os.path.exists(dir):
                request.send_error(404) # dir does not exist
                return
            if not os.access(dir, os.W_OK):
                request.send_error(403) # can not overwrite a dir
                return

        length=1024*1024*100 # max file size
        if request.headers.has_key('Content-Length'):
            length = int(request.headers['Content-Length'])

        f = open(fname, 'wb')
        try:
        
            if blocksize>length:
                blocksize = length
            buf = request.rfile.read(blocksize)
            while buf:
                f.write(buf)

                length -= len(buf)

                if blocksize>length:
                    blocksize = length
                buf = request.rfile.read(blocksize)
            request.send_response(200)
            request.end_headers()
            request.wfile.write('ok')
        except socket.error, (errno, msg):
            if not errno in self.socket_errors_to_ignore:
                self.GeniusCtrl.log('Socket error %s.\n' % msg)
        finally:
            f.close()


    def doLog(self, param, request):
        #self._checkError(request)
        if not request.command=='GET':
            request.send_error(405) # method not allowed
            return
            
        fname = self.GeniusCtrl.logFilename
        return self.readFile(fname, request)

    def doStructQuery(self, param, request):
        self._checkError(request)
        
        format = param['format']

        if format=='xml':
            gCmd = self.cmdQueue.addCmd( 'query-struct xml')
        elif format=='text':
            gCmd = self.cmdQueue.addCmd( 'query-struct text')

        
        if gCmd:
            request.send_response(200)
            if format=='xml':
                request.send_header('Content-Type', 'text/xml')
            elif format=='text':
                request.send_header('Content-Type', 'text/plain')
            request.end_headers()
            
            gCmd.finished.wait()
            data = gCmd.result

            request.wfile.write(data)
        else:
            request.send_error(404) # content no found
            
        return

    def doFile(self, param, request):
        self._checkError(request)

        path = param['path']
        path = os.path.join(os.getcwd(), path)

        if os.path.isdir(path):
            pass
        else:
            if request.command=='GET':
                return self.readFile(path, request, 32*1024)
            if request.command=='POST':
                return self.writeFile(path, request, 32*1024)
        
        request.send_error(400)


    def doSolution(self, param, request):
        self._checkError(request)
        
        gCmd = self.cmdQueue.addCmd('query-solution '+param['cmd'])
        gCmd.finished.wait()
        data = gCmd.result
        
        if not data==None:
            request.send_response(200)
            request.wfile.write(str(data))
        else:
            request.send_error(404)

    def doSolution(self, param, request):
        self._checkError(request)
        if not request.command=='GET':
            request.send_error(405) # method not allowed
            return
            
        fname = self.GeniusCtrl.SolverControl.getSolutionFile();
        return self.readFile(fname, request)

