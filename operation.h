#ifndef OPERATION_H
#define OPERATION_H

#include "circut.h"

namespace snrk {

class Operation
{
public:
    using input_t = struct{value_t a; value_t b;};
    using opFunc_t = std::function<value_t(input_t)>;

    Operation(opFunc_t opFunc);

    value_t operator()(input_t input);

private:
    opFunc_t m_opFunc;
};

extern std::unordered_map<GateType_t, Operation> Operate;

#define FOROPS for(const auto &[operation, _] : snrk::Operate)

}

#endif // OPERATION_H
