#include "config.h"
#include "pattern.h"
#include "mxml.h"
#include "MXMLUtil.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>

using namespace Parser;
using namespace DomMXMLParser;

namespace Parser
{

class DomPatternParser
{
public:
  Parameter parseParameter(mxml_node_t* node);
  PatternCard parseCommand(mxml_node_t* node);

  std::vector<std::string> _errors;
};

}

int Pattern::save_to_XML(const std::string &fname) const
{
  mxml_node_t *xml, *root;
  xml = mxmlNewXML("1.0");
  root = mxmlNewElement(xml, "Genius-syntax");

  for ( std::map< std::string, PatternCard>::const_iterator it = _pattern_card_map.begin();
        it != _pattern_card_map.end(); it++)
  {
    const PatternCard & command  = it->second;

    mxml_node_t *cmd = mxmlNewElement(root, "command");
    mxmlElementSetAttr(cmd, "name", command._key.c_str());

    mxml_node_t *desc = mxmlNewElement(cmd, "description");
    mxmlNewText(desc, 0, command._description.c_str());

    for ( std::map<std::string, Parameter>::const_iterator pp = command._parameter_map.begin();
          pp != command._parameter_map.end(); pp++ )
    {
      const Parameter& parameter = pp->second;

      mxml_node_t *param = mxmlNewElement(cmd, "parameter");
      mxmlElementSetAttr(param, "name", parameter.name().c_str());

      mxml_node_t *desc = mxmlNewElement(param, "description");
      mxmlNewText(desc, 0, parameter.description().c_str());

      switch (parameter.type())
      {
      case BOOL:
        mxmlElementSetAttr(param, "type", "bool");
        mxmlElementSetAttr(param, "default", parameter.get_bool() ? "yes": "no");
        break;
      case INTEGER:
      {
        mxmlElementSetAttr(param, "type", "int");
        std::stringstream str;
        str<<parameter.get_int();
        mxmlElementSetAttr(param, "default", str.str().c_str());
        break;
      }
      case REAL:
      {
        mxmlElementSetAttr(param, "type", "num");
        std::stringstream str;
        str<<parameter.get_real();
        mxmlElementSetAttr(param, "default", str.str().c_str());
        break;
      }
      case STRING:
        mxmlElementSetAttr(param, "type", "string");
        mxmlElementSetAttr(param, "default", parameter.get_string().c_str());
        break;
      case ENUM:
      {
        mxmlElementSetAttr(param, "type", "enum");
        for (Parameter::StringEnumIterator it = parameter.stringPatternBegin();
             it != parameter.stringPatternEnd(); it++)
        {
          mxml_node_t *e = mxmlNewElement(param, "enum");
          mxmlNewText(e, 0, it->c_str());
        }
        mxmlElementSetAttr(param, "default", parameter.get_string().c_str());
        break;
      }
      default:
        break;
      }
    }
  }

  FILE *fp;
  fp = fopen(fname.c_str(), "w");
  if (fp==NULL)
    return -1;

  int rc = mxmlSaveFile(xml, fp, MXML_NO_CALLBACK); // zero means successful
  fclose(fp);

  return rc;
}


int Pattern::get_from_XML(const std::string & fname)
{
  mxml_node_t *root;
  FILE* file;

  file = fopen(fname.c_str(), "r");
  root = mxmlLoadFile(NULL, file, MXML_TEXT_CALLBACK);

  if (root==NULL)
  {
    return -1;
  }

  DomPatternParser parser;
  for (mxml_node_t * nodeCmd = mxmlFindElement(root, root, "command", NULL, NULL, MXML_DESCEND);
       nodeCmd != NULL;
       nodeCmd = mxmlFindElement(nodeCmd, root, "command", NULL, NULL, MXML_NO_DESCEND))
  {
    PatternCard cmd = parser.parseCommand(nodeCmd);
    _pattern_card_map.insert(std::pair<std::string, PatternCard>(cmd._key, cmd));
  }

  mxmlRelease(root);
  fclose(file);

  if (parser._errors.empty())
    return 0;


  for (std::vector<std::string>::const_iterator it = parser._errors.begin();
       it != parser._errors.end(); it++)
    std::cout << *it << std::endl;
  return -1;

}

Parameter DomPatternParser::parseParameter(mxml_node_t* node)
{
  Parameter param;

  param.set_name( getElementAttr(node, "name") );
  if (param.name().length()<=0)
    _errors.push_back( "DOM Error: Parameter name missing");
  std::string type = getElementAttr(node, "type");
  if (type.length() <= 0)
    _errors.push_back( "DOM Error: Parameter type missing");

  if (mxml_node_t *nd = mxmlFindElement(node, node, "description", NULL, NULL, MXML_DESCEND_FIRST))
    param.set_description( getElementText(nd) );

  if (strncasecmp(type.c_str(), "enum", 10) == 0)
    param.set_type(ENUM);
  else if (strncasecmp(type.c_str(), "bool", 10) == 0)
    param.set_type(BOOL);
  else if (strncasecmp(type.c_str(), "num", 10) == 0)
    param.set_type(REAL);
  else if (strncasecmp(type.c_str(), "int", 10) == 0)
    param.set_type(INTEGER);
  else if (strncasecmp(type.c_str(), "string", 10) == 0)
    param.set_type(STRING);
  else
    param.set_type(REAL);

  switch (param.type())
  {
  case ENUM:
  {
    for (mxml_node_t *ne = mxmlFindElement(node, node, "enum", NULL, NULL, MXML_DESCEND_FIRST);
         ne != NULL; ne = mxmlFindElement(ne, node, "enum", NULL, NULL, MXML_NO_DESCEND))
    {
      param.add_string_pattern(getElementText(ne));
    }
    std::string strDefault = getElementAttr(node, "default");

    if (!strDefault.empty())
      if (param.string_pattern_match(strDefault))
        _errors.push_back("DOM Error: invalid enum default value "+strDefault);
    break;
  }
  case REAL:
  {
    double val;
    std::stringstream ss(std::string( getElementAttr(node, "default") ));
    ss >> val;
    if (ss.fail())
    {
      param.set_real(0.0);
      _errors.push_back("DOM Error: invalid numerical default value "+ss.str());
    }
    else
      param.set_real(val);
    break;
  }
  case INTEGER:
  {
    int val;
    std::stringstream ss(std::string( getElementAttr(node, "default") ));
    ss >> val;
    if (ss.fail())
    {
      param.set_int(0);
      _errors.push_back("DOM Error: invalid integer default value "+ss.str());
    }
    else
      param.set_int(val);
    break;
  }
  case BOOL:
  {
    std::string strBool = getElementAttr(node, "default");
    if (strncasecmp(strBool.c_str(), "yes", 10) ==0 ||
        strncasecmp(strBool.c_str(), "true", 10) ==0 ||
        strncasecmp(strBool.c_str(), "on", 10) ==0 )
      param.set_bool(true);
    else if (strncasecmp(strBool.c_str(), "no", 10) ==0 ||
        strncasecmp(strBool.c_str(), "false", 10) ==0 ||
        strncasecmp(strBool.c_str(), "off", 10) ==0 )
      param.set_bool(false);
    else
    {
      param.set_bool(false);
      _errors.push_back("DOM Error: invalid bool default value "+getElementAttr(node, "default"));
    }
    break;
  }
  case STRING:
  {
    param.set_string(getElementAttr(node, "default"));
    break;
  }
  default:
    _errors.push_back("DOM Error: invalid parameter type.");
    break;
  }

  return param;
}

PatternCard DomPatternParser::parseCommand(mxml_node_t* node)
{
  PatternCard cmd;
  cmd._key = getElementAttr(node, "name");
  if (cmd._key.length()<=0)
    _errors.push_back("DOM Error: Command name missing");

  for (mxml_node_t * n = mxmlFindElement(node, node, "parameter", NULL, NULL, MXML_DESCEND);
       n!=NULL; n = mxmlFindElement(n, node, "parameter", NULL, NULL, MXML_DESCEND))
  {
    Parameter param = parseParameter(n);
    cmd._parameter_map.insert(std::pair<const std::string, const Parser::Parameter>(param.name(),param));
  }
  return cmd;
}
