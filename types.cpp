#include "types.h"

namespace snrk
{

Value::operator int() const
{
    return static_cast<int>(get_si());
}

Value::operator const unsigned long() const
{
    return get_ui();
}

Value::operator double() const
{
    return get_d();
}

Value::operator GateType_t() const
{
    if (floor(*this) != *this) {
        return Unknown;
    }

    return static_cast<GateType_t>(get_ui());
}

Value::operator witness_t() const
{
    return static_cast<witness_t>(get_ui());
}

bool Value::operator==(int a) const
{
    return cmp(*this, a) == 0;
}

address_t Value::address() const
{
    return address_t(this);
}

bool operator<(const Value &a, const Value &b)
{
    return cmp(a, b) == -1;
}

//bool operator==(const Value &a, const Value &b)
//{
//    return a.m_value.get() == b.m_value.get();
//}


}
