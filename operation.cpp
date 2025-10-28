#include "operation.h"

namespace snrk {

Operation::Operation(opFunc_t opFunc)
    : m_opFunc{opFunc}
{
}

mpf_class Operation::operator()(input_t ipnut)
{
    return m_opFunc(ipnut);
}

}
