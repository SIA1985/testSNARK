#include "types.h"

namespace snrk
{

void to_json(json_t& j, const witnesses_t& ws)
{
    j["array"] = static_cast<const std::vector<witness_t>&>(ws);
}

void from_json(const snrk::json_t& j, witnesses_t& ws)
{
    ws = j.at("array").get<std::vector<witness_t>>();
}

Value::Value(mcl::Fr fr)
{
    *this = fr;
}

Value::operator int() const
{
    return static_cast<int>(getInt64());
}

Value::operator const unsigned long() const
{
    return getUint64();
}

Value::operator double() const
{
    return getInt64();
}

Value::operator GateType_t() const
{
//    if (floor(*this) != *this) {
//        return Unknown;
//    }

    return static_cast<GateType_t>(*getUnit());
}

//Value::operator witness_t() const
//{
//    return static_cast<witness_t>(getUint64());
//}

Value::operator mpz_class() const
{
    return getMpz();
}

//Value::operator mcl::Fr() const
//{
//    return get_ui();
//}

void to_json(json_t& j, const dot_t& dot) {
    j["x"] = dot.x;
    j["y"] = dot.y;
}

void from_json(const snrk::json_t& j, dot_t& dot) {
    dot.x = j.at("x").get<value_t>();
    dot.y = j.at("y").get<value_t>();
}

address_t Value::address() const
{
    return address_t(this);
}

//bool operator<(const Value &a, const Value &b)
//{
//    return cmp(a, b) == -1;
//}

//bool operator==(const Value &a, const Value &b)
//{
//    return cmp(a, b) == 0;
//}

//todo:
void to_json(json_t& j, const value_t& value)
{
    j["double"] = value.getUint64();
}

void from_json(const snrk::json_t& j, value_t& value)
{
    value = j.at("double").get<double>();
}

//todo:
void to_json(json_t& j, const GPK_t& gpk)
{
//    j["t"] = tG.t;
//    j["G"] = tG.G;
}

//todo:
void from_json(const snrk::json_t& j, GPK_t& gpk)
{
//    tG.t = j.at("t").get<value_t>();
//    tG.G = j.at("G");
}

json_t Jsonable::toJson() const
{
    return "{}";
}

bool Jsonable::fromJson(const json_t &json)
{
    return false;
}

void to_json(json_t& j, const Jsonable& jsonable)
{
    j["obj"] = jsonable.toJson();
}

void from_json(const snrk::json_t& j, Jsonable& jsonable)
{
    jsonable.fromJson(j.at("obj"));
}

}
