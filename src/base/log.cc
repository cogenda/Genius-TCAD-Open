/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/

#include "log.h"

GENIUS_LOG_STREAM genius_log;

GENIUS_LOG_STREAM::GENIUS_LOG_STREAM()
{
}

GENIUS_LOG_STREAM::~GENIUS_LOG_STREAM()
{
  for (std::map<std::string, std::ostream*>::iterator it = _streams.begin();
       it != _streams.end(); it++)
  {
    delete it->second;
  }
}

void GENIUS_LOG_STREAM::record()
{
  for (std::map<std::string, std::ostream*>::iterator it = _streams.begin();
       it != _streams.end(); it++)
  {
    std::string msg(_sstream.str());
    (*it->second) << msg;
    it->second->flush();
  }
  _sstream.str("");
}

void GENIUS_LOG_STREAM::addStream(const std::string &name, std::streambuf* buf)
{
  if (buf)
  {
    std::ostream *os = new std::ostream(buf);
    os->setf(std::ios::scientific);
    _streams.insert(std::pair<std::string, std::ostream*>(name, os) );
  }
}

void GENIUS_LOG_STREAM::addStream(const std::string &name, const std::string &fname)
{
  std::filebuf *buf = new std::filebuf;
  buf->open(fname.c_str(), std::ios::out);

  std::ostream *os = new std::ostream(buf);

  _streams.insert(std::pair<std::string, std::ostream*>(name, os) );
  _bufs   .insert(std::pair<std::string, std::filebuf*>(name, buf));
}

void GENIUS_LOG_STREAM::removeStream(const std::string &name)
{
  {
    std::map<std::string, std::ostream*>::iterator it = _streams.find(name);
    if (it != _streams.end())
    {
      _streams.erase(it);
      delete it->second;
    }
  }

  {
    std::map<std::string, std::filebuf*>::iterator it = _bufs.find(name);
    if (it != _bufs.end())
    {
      _bufs.erase(it);
      it->second->close();
      delete it->second;
    }
  }
}
