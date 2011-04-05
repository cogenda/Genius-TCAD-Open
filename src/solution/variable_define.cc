#include "variable_define.h"
#include "expr_evaluate.h"

SimulationVariable::SimulationVariable(const std::string &name, DataType data_type, DataLocation data_location,
                     const std::string & unit_string, unsigned int index, bool valid, bool user_defined)
  :variable_name(name), variable_data_type(data_type), variable_data_location(data_location),
      variable_unit_string(unit_string), variable_index(index), variable_valid(valid),
      variable_user_defined(user_defined)
{
  if(variable_unit_string.empty())
    variable_unit = 1.0;
  else
    variable_unit = ConstanteExprEvalute(variable_unit_string)();
}
