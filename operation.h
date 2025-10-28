#ifndef OPERATION_H
#define OPERATION_H

#include "circut.h"

namespace snrk {

class Operation
{
public:
    using input_t = struct{ValueType a; ValueType b;};
    using opFunc_t = std::function<ValueType(input_t)>;

    Operation(opFunc_t opFunc);

    ValueType operator()(input_t input);

private:
    opFunc_t m_opFunc;
};

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

#define FOROPS for(const auto &[operation, _] : Operate)

}

#endif // OPERATION_H
