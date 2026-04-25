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
    if (toHash.length() <= 3) {
        return toHash;
    }

    return toHash.substr(0, 1) + toHash.substr(toHash.length() - 2);
}

MerkleTree::MerkleTree(const std::vector<std::string> &data)
{
    assert(!data.empty());


    leafsCount = data.size();

    //todo:
//    m_tree.reserve(layers);

    hashes_t lay;
    for(const auto &d : data) {
        lay.push_back(hash(d));
    }
    m_tree.push_back(lay);

    while(lay.size() > 1) {
        //todo: уязвимость с последним дубликатом
        if (lay.size() % 2 == 1) {
            lay.push_back(lay.back());
        }

        hashes_t temp;
        for(std::size_t i = 0; i < lay.size(); i += 2) {
            temp.push_back(hash(lay[i] + lay[i + 1]));
        }

        lay = temp;
        m_tree.push_back(lay);
    }

}

hash_t MerkleTree::root() const
{
    return m_tree.back().back();
}

MerkleTree::path_t MerkleTree::path(hash_t leaf) const
{
    if (m_tree.empty()) {
        return std::nullopt;
    }

    auto &leafs = m_tree.front();

    auto it = std::find(leafs.begin(), leafs.end(), leaf);
    if (it == leafs.end()) {
        return std::nullopt;
    }

    std::size_t index = std::distance(leafs.begin(), it);
    hashes_t proof;

    for (std::size_t i = 0; i < m_tree.size() - 1; ++i) {
        const auto &layer = m_tree[i];

        std::size_t pairIndex = (index % 2 == 0) ? index + 1 : index - 1;

        if (pairIndex >= layer.size()) {
            pairIndex = index;
        }

        proof.push_back(layer[pairIndex]);
        index /= 2;
    }

    return proof;
}

bool MerkleTree::verify(const hashes_t &path, hash_t leaf, bool isLeft, hash_t root)
{
    hash_t result = leaf;

    for(const auto &pair : path) {
        if (isLeft) {
            result = hash(result + pair);
        } else {
            result = hash(pair + result);
        }
    }

    return root == result;
}

}
