#include "MXMLUtil.h"
#include <iostream>
#include <sstream>

namespace DomMXMLParser
{

std::string getElementAttr(mxml_node_t* node, const char* name)
{
  std::string str;

  if (node==NULL || name==NULL)
    return str;

  const char *s = mxmlElementGetAttr(node, name);
  if (s)
    str = s;

  return str;
}

std::string getElementText(mxml_node_t* node)
{
  std::stringstream str;

  if (node && node->type==MXML_ELEMENT && node->child)
  {
    for (mxml_node_t* tn = mxmlWalkNext(node, node, MXML_DESCEND);
         tn != NULL; tn = mxmlWalkNext(tn, node, MXML_NO_DESCEND))
    {
      if (tn->type != MXML_TEXT)
        str<<" ";
      else
      {
        if (tn->value.text.whitespace)
          str<<" ";
        str<< tn->value.text.string;
      }
    }
  }
  return str.str();
}

std::string getElementCData(mxml_node_t *node)
{
  if (node && node->type==MXML_ELEMENT && node->child)
  {
    for (mxml_node_t* n = mxmlWalkNext(node, node, MXML_DESCEND);
         n != NULL; n = mxmlWalkNext(n, node, MXML_NO_DESCEND))
    {
      char* data = n->value.element.name;

      if (n->type == MXML_ELEMENT && strncasecmp(data, "![CDATA[", 8)==0)
      {
        std::string str(data+8, strlen(data)-8-2);
        return str;
      }
    }
  }
  return std::string();
}

void dumpTree(mxml_node_t *node, int level)
{
  if (node==NULL)
    return;

  for (int i=0; i<level; i++)
    std::cout << " ";

  switch (node->type)
  {
  case MXML_ELEMENT:
    std::cout << "Element: " << node->value.element.name << std::endl;
    break;
  case MXML_TEXT:
    std::cout << "Text: " << node->value.text.string << std::endl;
    break;
  default:
    std::cout << "other" << std::endl;
  }
  for (mxml_node_t* tn = mxmlWalkNext(node, node, MXML_DESCEND);
       tn != NULL; tn = mxmlWalkNext(tn, node, MXML_NO_DESCEND))
  {
    dumpTree(tn, level+1);
  }
}


}


namespace MXMLQVariant
{

mxml_node_t* makeQVFloat(double value)
{
  mxml_node_t *elem = mxmlNewElement(NULL, "qvFloat");

  std::stringstream ss;
  ss.precision(15);
  ss << value;
  mxmlNewText(elem, 0, ss.str().c_str());

  return elem;
}

mxml_node_t* makeQVInt(int value)
{
  mxml_node_t *elem = mxmlNewElement(NULL, "qvInt");

  std::stringstream ss;
  ss << value;
  mxmlNewText(elem, 0, ss.str().c_str());

  return elem;
}

mxml_node_t* makeQVInt(long value)
{
  mxml_node_t *elem = mxmlNewElement(NULL, "qvInt");

  std::stringstream ss;
  ss << value;
  mxmlNewText(elem, 0, ss.str().c_str());

  return elem;
}

mxml_node_t* makeQVString(const std::string &value)
{
  mxml_node_t *elem = mxmlNewElement(NULL, "qvString");
  mxmlNewText(elem, 0, value.c_str());
  return elem;
}

}
