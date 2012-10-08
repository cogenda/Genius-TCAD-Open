HERE=$(cd `dirname ${BASH_SOURCE[0]}` && pwd)
BASEDIR=`cd $HERE/../.. && pwd`

H_MATH="acosh.hpp adolc.h asinh.hpp atanh.hpp type_vector.h vector_value.h type_tensor.h tensor_value.h"
H_GEOM="point.h "
H_SOLUTION="data_storage.h data_object.h fvm_node_data.h variable_define.h"
H_MATERIAL="physical_unit.h PMI.h"
H_ENUM="enum_data_location.h enum_data_type.h enum_solution.h"
H_BASE="base/genius_common.h base/genius_dll.h parser/parser_parameter.h utils/compare_types.h"

for f in $H_MATH;     do cp $BASEDIR/include/math/$f $HERE/include/; done
for f in $H_GEOM;     do cp $BASEDIR/include/geom/$f $HERE/include/; done
for f in $H_SOLUTION; do cp $BASEDIR/include/solution/$f $HERE/include/; done
for f in $H_MATERIAL; do cp $BASEDIR/include/material/$f $HERE/include/; done
for f in $H_ENUM;     do cp $BASEDIR/include/enums/$f $HERE/include/; done
for f in $H_BASE;     do cp $BASEDIR/include/$f $HERE/include/; done

CC_BASE="math/adolc_init.cc parser/parser_parameter.cc material/PMI.cc"
for f in $CC_BASE;     do cp $BASEDIR/src/$f $HERE/; done

cp $BASEDIR/src/material/4H-SiC/*.cc $HERE/4H-SiC

