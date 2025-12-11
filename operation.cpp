#include "operation.h"

namespace snrk {

Operation::Operation(opFunc_t opFunc)
    : m_opFunc{opFunc}
{
}

value_t Operation::operator()(input_t ipnut)
{
    return m_opFunc(ipnut);
}


#define OPHEADER [](Operation::input_t input) -> value_t
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
    },
    {Minus, OP({
         return input.a - input.b;
    })
    },
    {Devide, OP({
         return input.a / input.b;
    })
    }
};

#undef OPHEADER
#undef OP

}
