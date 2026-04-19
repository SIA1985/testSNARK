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

//Value::operator double() const
//{
//    return getInt64();
//}

//Value::operator mcl::Fr() const
//{
//    return get_ui();
//}

//todo:
void to_json(json_t& j, const dot_t& dot) {
//    j["x"] = dot.x;
//    j["y"] = dot.y;
}

//todo:
void from_json(const snrk::json_t& j, dot_t& dot) {
//    dot.x = j.at("x").get<value_t>();
//    dot.y = j.at("y").get<value_t>();
}

//todo:
void to_json(json_t& j, const value_t& value)
{
//    j["double"] = static_cast<double>(value);
}

//todo: получаем всегда уникальное значение, как быть с адресом?
void from_json(const snrk::json_t& j, value_t& value)
{
//    value = j.at("double").get<double>();
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
