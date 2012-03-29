
// {{{ namespace Parser
namespace Parser
{
%TypeHeaderCode
#include "parser_card.h"
%End
  enum  ElemType {BOOL, INTEGER, REAL, STRING, ENUM, INVALID};

  // {{{ class Parameter
  class Parameter
  {

    public:


      Parameter(const std::string &name, const bool   v /Constrained/ );
      Parameter(const std::string &name, const int    v /Constrained/ );
      Parameter(const std::string &name, const double v /Constrained/ );
      Parameter(const std::string &name, const char * v /Constrained/ );

      void set_name (const std::string &name);
      const std::string name() const;

      void set_type (const Parser::ElemType type);
      Parser::ElemType type() const;

      void   set_bool(const bool v);
      void   set_int(const int v);
      void   set_real(const double v);
      void   set_string(const std::string & s);
      int    set_enum(const std::string & s);

      void clear ();

      bool    get_bool()   const;
      int     get_int()    const;
      double  get_real()   const;
      std::string  get_string() const;

      bool   is_user_defined();
  };
  // }}}

  // {{{ class Card
  class Card
  {
  public:
    Card();
    Card(const std::string &key);
    void insert(const Parser::Parameter &p);
    void insert(const std::string &name, const bool value /Constrained/ );
    void insert(const std::string &name, const int value /Constrained/ );
    void insert(const std::string &name, const double value /Constrained/ );
    void insert(const std::string &name, const char * value /Constrained/ );

    bool    get_bool(const char * name , bool default_value, const char * alias=NULL)   const;
    int     get_int(const char * name , int default_value, const char * alias=NULL)   const;
    double  get_real(const char * name , double default_value, const char * alias=NULL)   const;
    std::string  get_string(const char * name , std::string default_value, const char * alias=NULL)   const;
    std::string  get_string_lower_case(const char * name , std::string default_value, const char * alias=NULL)   const;

    bool  is_parameter_exist(const std::string & parameter_name) const;
    int parameter_count(const std::string & parameter_name) const;
    int parameter_size() const;

    const Parser::Parameter & get_parameter(unsigned int idx) const;
    void set_parameter(const Parser::Parameter &p, unsigned int idx);

    std::string key()  const;
    void set_key(const std::string key);

    void clear();

    const std::string & get_fileline() const;
    void set_fileline(const std::string & file_line);

    int get_lineno() const;
    void set_lineno(int line);
  };
  // }}}

  // {{{ class Pattern
  class Pattern
  {
%TypeHeaderCode
#include "pattern.h"
%End
  public:
    Pattern();

    int get_from_XML(const std::string &);
    int check_detail(Parser::Card & c) const;
  };
  // }}}

  // {{{ class FilePreProcess
  class FilePreProcess
  {
%TypeHeaderCode
#include "file_include.h"
%End

  public:
    FilePreProcess(const char * filename);
    std::string output();
  };
  // }}}

  // {{{ class InputParser
  class InputParser
  {
%TypeHeaderCode
#include "parser.h"
%End
  public:
    InputParser(Parser::Pattern & p);
    int read_card_file(const char * file);
  };
  // }}}
};
// }}}

// {{{ class SolverControl
class SolverControl  /DelayDtor/
{
%TypeHeaderCode
#include "control.h"
%End

public:
  SolverControl() /Transfer/;

  const SimulationSystem& system() const /NoCopy/;


  void setDecks(Parser::InputParser* p /Transfer/);
  void setSolutionFile(const std::string &fname_result);
  std::string getSolutionFile();

  int mainloop() /ReleaseGIL/;
  int reset_simulation_system() /ReleaseGIL/;
  int do_process() /ReleaseGIL/;
  int do_mesh () /ReleaseGIL/;

  int set_method ( const Parser::Card & c ) /ReleaseGIL/;
  int set_model ( const Parser::Card & c ) /ReleaseGIL/;
  int do_hook ( const Parser::Card & c) /ReleaseGIL/;
  int do_solve   ( const Parser::Card & c ) /ReleaseGIL/;
  int do_export  ( const Parser::Card & c ) /ReleaseGIL/;
  int do_import  ( const Parser::Card & c ) /ReleaseGIL/;
  int set_initial_node_voltage  ( const Parser::Card & c ) /ReleaseGIL/;
  int set_electrode_source  ( const Parser::Card & c ) /ReleaseGIL/;
  int set_physical_model  ( const Parser::Card & c ) /ReleaseGIL/;
  int do_refine_conform ( const Parser::Card & c ) /ReleaseGIL/;
  int do_refine_hierarchical ( const Parser::Card & c ) /ReleaseGIL/;
  int do_refine_uniform  ( const Parser::Card & c ) /ReleaseGIL/;
  int extend_to_3d ( const Parser::Card & c ) /ReleaseGIL/;
  int plot_mesh ( const Parser::Card & c ) /ReleaseGIL/;

  int get_dom_solution_count();

  SIP_PYOBJECT get_dom_solution();
%MethodCode
    sipRes = Py_None;

    mxml_node_t* dom = sipCpp->get_dom_solution();
    char *buf = mxmlSaveAllocString(dom, MXML_NO_CALLBACK);
    if (buf)
    {
      sipRes = PyString_FromString(buf);
      free(buf);
    }
%End

private:
  SolverControl(const SolverControl&);

};
// }}}
