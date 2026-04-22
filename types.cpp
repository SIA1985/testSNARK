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

//todo:
hash_t hash(std::string toHash)
{
    return toHash.substr(0, 3);
}

MerkleTree::MerkleTree(const std::vector<std::string> &data)
{
    assert(data.size() != 0);


    leafsCount = data.size();
    int layers = std::ceil(std::log2(leafsCount)) + 1;

    m_tree.reserve(layers);

    hashes_t lay, temp;
    for(const auto &d : data) {
        lay.push_back(hash(d));
    }

    for(int i = 0; i < layers; i++) {
        temp.clear();

        if (lay.size() % 2 == 1) {
            lay.push_back(lay.back());
        }

        for(std::size_t j = 0; j < lay.size(); j += 2) {
            auto toHash = lay[j] + lay[j + 1];
            temp.push_back(hash(toHash));
        }

        lay = temp;
        m_tree.push_back(lay);
    }

}

hash_t MerkleTree::root() const
{
    return m_tree.back().back();
}

hashes_t MerkleTree::path(hash_t leaf) const
{

}

}
