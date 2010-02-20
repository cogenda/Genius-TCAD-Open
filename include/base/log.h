#ifndef __log_h__
#define __log_h__

#include <fstream>
#include <sstream>
#include <string>
#include <map>

/**
 *  We need two streams for message output: screen and log file.
 */
class GENIUS_LOG_STREAM
{
private:
  static GENIUS_LOG_STREAM *_instance;

  std::map<std::string, std::ostream*> _streams;  // output streams.
  std::map<std::string, std::filebuf*> _bufs;  // file buffers opened in genius.
  std::ostringstream _sstream;

public:
  /**
   * constructor
   */
  GENIUS_LOG_STREAM();

  /**
   * destructor
   */
  ~GENIUS_LOG_STREAM();

  /**
   * attach a stream buffer to which log message is written.
   * Typical choice of buf can be std::cerr.rdbuf()
   * or __gnu_cxx::stdio_filebuf(...)
   */
  void addStream(const std::string &name, std::streambuf* buf);

  /**
   * open a file and attach it for log message
   */
  void addStream(const std::string &name, const std::string &fname);

  /**
   * remove a stream buffer
   */
  void removeStream(const std::string &name);

  /**
   * stream to add log message
   */
  std::ostringstream& log_stream() { return _sstream; }

  /**
   * flush the buffer
   */
  void record();
};

extern GENIUS_LOG_STREAM genius_log;

#define   MESSAGE   genius_log.log_stream()
#define   RECORD()  genius_log.record()


#endif
