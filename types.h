#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <gmpxx.h>
#include <set>
#include <memory>
#include "nlohmann/json.hpp"
#include <string>
#include <vector>
#include <cmath>
#define MCL_USE_GMP 1
#include <mcl/bls12_381.hpp>

namespace snrk {

typedef std::size_t address_t;
typedef nlohmann::json json_t;

enum GateType_t : char {
    Unknown = 0,

    Sum,
    Product,
    Minus,
    Devide,
};

typedef unsigned long witness_t;

class witnesses_t: public std::vector<witness_t>
{
private:
    witnesses_t() = default;
    using std::vector<witness_t>::operator=;

    friend witnesses_t genWitnesses(witness_t start, std::size_t count, witness_t wStep);
    friend class GlobalParams;

    friend void to_json(json_t& j, const witnesses_t& ws);
    friend void from_json(const snrk::json_t& j, witnesses_t& ws);
};
void to_json(json_t& j, const witnesses_t& ws);
void from_json(const snrk::json_t& j, witnesses_t& ws);

template<typename T>
class CircutValue
{
    using ptr_t = std::shared_ptr<T>;

public:
    CircutValue() {
        m_data = ptr_t(new T);
    }
    CircutValue(T data) {
        m_data = ptr_t(new T(data));
    }
    CircutValue(const CircutValue &value) {
        m_data = value.m_data;
    }

    CircutValue& operator=(const CircutValue&) = delete;

    bool operator<(const CircutValue &value) {
        return *m_data < *value.m_data;
    }

    address_t address() const {
        return address_t(m_data.get());
    }

    T value() const {
        return *m_data;
    }

private:
    ptr_t m_data;
};

typedef mcl::Fr value_t;
typedef CircutValue<value_t> CircutValue_t;
void to_json(json_t& j, const value_t& value);
void from_json(const snrk::json_t& j, value_t& value);

typedef std::vector<value_t> values_t;
typedef std::vector<CircutValue_t> CircutValues_t;

typedef value_t X_t;
typedef value_t Y_t;

template<typename X = X_t, typename Y = Y_t>
struct DotType{X x; Y y;};

typedef DotType<X_t, Y_t> dot_t;
void to_json(json_t& j, const dot_t& dot);
void from_json(const snrk::json_t& j, dot_t& dot);

typedef std::vector<dot_t> dots_t;

typedef mcl::bn::G1 G1;
typedef mcl::bn::G2 G2;
typedef mcl::bn::GT GT;
typedef G1 commit_t;
typedef std::vector<G1> keys_t;
struct GPK_t {keys_t keys; G1 g1;};
void to_json(json_t& j, const GPK_t& gpk);
void from_json(const snrk::json_t& j, GPK_t& gpk);

typedef std::set<X_t> xs_t;

class InterpolationPolynom;
typedef InterpolationPolynom T_t;
typedef InterpolationPolynom S_t;
typedef InterpolationPolynom W_t;

#define Name(a) #a
#define FromJson(json, a) json[Name(a)].get_to(a)
#define ToJson(json, a) json[Name(a)] = a
class Jsonable
{
protected:
    virtual json_t toJson() const;
    virtual bool fromJson(const json_t &json);

    friend void to_json(json_t& j, const Jsonable& jsonable);
    friend void from_json(const snrk::json_t& j, Jsonable& jsonable);
};
void to_json(json_t& j, const Jsonable& jsonable);
void from_json(const snrk::json_t& j, Jsonable& jsonable);

typedef std::string hash_t; //todo: char[256]
typedef std::vector<hash_t> hashes_t;

hash_t hash(std::string toHash);

class MerkleTree {
public:
    MerkleTree(const std::vector<std::string> &data);

    hash_t root() const;

    hashes_t path(hash_t leaf) const;

private:
    std::vector<hashes_t> m_tree;
    std::size_t leafsCount;
};
}

#endif // TYPES_H
