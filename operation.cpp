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


#define OPHEADER [](Operation::input_t input) -> ValueType
#define OP(code) {OPHEADER code}

std::unordered_map<GateType_t, Operation> Operate =
{
    {Sum, OP({
        return input.a + input.b;
    })
    },
    {Product, OP({
         return input.a * input.b;
     })

    }
};

#undef OPHEADER
#undef OP

}
