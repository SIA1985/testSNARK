#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <gmpxx.h>
#include <set>
#include <memory>
#include "nlohmann/json.hpp"

namespace snrk {

typedef std::size_t address_t;

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

    friend witnesses_t genWitnesses(witness_t start, std::size_t count, witness_t wStep);
    friend class GlobalParams;
    friend class ZeroTestProof;
};

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
    using mpf_class::operator=;

    address_t address() const;

private:
    friend bool operator<(const Value &a, const Value &b);
//    friend bool operator==(const Value &a, const Value &b);
};

typedef Value value_t;
bool operator<(const value_t &a, const value_t &b);
bool operator==(const value_t &a, const value_t &b);
typedef std::vector<value_t> values_t;

typedef value_t X_t;
typedef value_t Y_t;

struct dot_t {X_t x; Y_t y;};
typedef std::vector<dot_t> dots_t;
struct TG_t {X_t t; int G;};
typedef std::set<X_t> xs_t;

class InterpolationPolynom;
typedef InterpolationPolynom T_t;
typedef InterpolationPolynom S_t;
typedef InterpolationPolynom W_t;
typedef InterpolationPolynom WT_t;

typedef nlohmann::json json_t;
class Jsonable
{
public:
    virtual json_t toJson() const;
    virtual bool fromJson(const json_t &json);
};

}

#endif // TYPES_H
