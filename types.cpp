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

address_t Value::address() const
{
    return address_t(this);
}

bool operator<(const Value &a, const Value &b)
{
    return cmp(a, b) == -1;
}

bool operator==(const Value &a, const Value &b)
{
    return cmp(a, b) == 0;
}

json_t Jsonable::toJson() const
{
    return "{}";
}

bool Jsonable::fromJson(json_t json)
{
    return false;
}

}
