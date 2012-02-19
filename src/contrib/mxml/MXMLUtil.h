#ifndef _MXMLUTIL_H_
#define _MXMLUTIL_H_

#include <string>
#include "mxml.h"

#if defined(WINDOWS) || defined(WIN32)
inline int strncasecmp(const char* s1, const char* s2, size_t n)
{
  return _strnicmp(s1, s2, n);
}
#endif

namespace DomMXMLParser
{
  extern std::string getElementAttr(mxml_node_t* node, const char* name);

  extern std::string getElementText(mxml_node_t* node);

  extern std::string getElementCData(mxml_node_t* node);

  extern void dumpTree(mxml_node_t* node, int level=0);
}

namespace MXMLQVariant
{
  extern mxml_node_t* makeQVFloat(double v);

  extern mxml_node_t* makeQVInt(int value);

  extern mxml_node_t* makeQVInt(long value);

  extern mxml_node_t* makeQVString(const std::string &v);
}

#endif
