#ifndef __memory_log_h__
#define __memory_log_h__

#include <memory>


/**
 *  this class measure the memory useage of the process
 */ 
class MMU
{

private:
  
  MMU() {}
  ~MMU() {}
  
  static std::auto_ptr<MMU> _instance;
  friend class std::auto_ptr<MMU>;
   
  static int _vmsize;
  static int _vmpeak;
  static int _vmrss;
  static int _vmhwm;
  
public:

  static MMU * instance();

  /**
   * call this function to do the memory usage measurement 
   */ 
  void measure();
  
  /// Current virtual memory usage
  int vmsize() { return _vmsize; }
  
  /// Peak virtual memory usage
  int vmpeak() { return _vmpeak; }
   
  /// Current resident set size
  int vmrss() { return _vmrss; }
     
  /// Peak resident set size
  int vmhwm() { return _vmhwm; }

private:
  
  /**
   * do the statistic, @return 0 for success
   */
  int _statistic(char *pidstatus);

};




#endif //__memory_log_h__

