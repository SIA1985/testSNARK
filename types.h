#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <gmpxx.h>
#include <set>
#include <memory>
#include "nlohmann/json.hpp"
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
    friend class ZeroTestProof;

    friend void to_json(json_t& j, const witnesses_t& ws);
    friend void from_json(const snrk::json_t& j, witnesses_t& ws);
};
void to_json(json_t& j, const witnesses_t& ws);
void from_json(const snrk::json_t& j, witnesses_t& ws);

class Value : public mpf_class
{
public:
    Value() = default;
    using mpf_class::mpf_class;

    operator int() const;
    operator const unsigned long() const;
    operator double() const;
    operator GateType_t() const;
    operator witness_t() const;
    operator mcl::Fr() const;
    using mpf_class::operator=;

    address_t address() const;

private:
    friend bool operator<(const Value &a, const Value &b);
//    friend bool operator==(const Value &a, const Value &b);
};

typedef Value value_t;
bool operator<(const value_t &a, const value_t &b);
bool operator==(const value_t &a, const value_t &b);
void to_json(json_t& j, const value_t& value);
void from_json(const snrk::json_t& j, value_t& value);

typedef std::vector<value_t> values_t;

typedef value_t X_t;
typedef value_t Y_t;

template<typename X = X_t, typename Y = Y_t>
struct DotType{X x; Y y;};

typedef DotType<X_t, Y_t> dot_t;
void to_json(json_t& j, const dot_t& dot);
void from_json(const snrk::json_t& j, dot_t& dot);

typedef std::vector<dot_t> dots_t;

typedef std::vector<mcl::G1> GPK_t;
void to_json(json_t& j, const GPK_t& gpk);
void from_json(const snrk::json_t& j, GPK_t& gpk);

typedef std::set<X_t> xs_t;

class InterpolationPolynom;
typedef InterpolationPolynom T_t;
typedef InterpolationPolynom S_t;
typedef InterpolationPolynom W_t;
typedef InterpolationPolynom WT_t;

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

typedef mcl::G1 commit_t;

}

#endif // TYPES_H
