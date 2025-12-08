#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <gmpxx.h>
#include <set>
#include <complex>

namespace snrk {

//todo: заменяемый тип
#ifndef ValueType
    #define ValueType mpf_class
#endif

enum GateType_t : char {
    Unknown = 0,

    Sum,
    Product,
    Minus,
    Devide,
};

class Value;
typedef Value value_t;
bool operator<(const Value &a, const Value &b);
bool operator==(const Value &a, const Value &b);
typedef std::vector<value_t> values_t;

typedef ValueType X_t;
typedef ValueType Y_t;

struct dot_t {X_t x; Y_t y;};
typedef std::vector<dot_t> dots_t;
struct TG_t {X_t t; int G;};
typedef std::set<X_t> xs_t;
typedef std::complex<ValueType> complexValue_t;
typedef std::vector<complexValue_t> complexValues_t;


typedef unsigned long witness_t;
typedef std::vector<witness_t> witnesses_t;

class InterpolationPolynom;
typedef InterpolationPolynom T_t;
typedef InterpolationPolynom S_t;
typedef InterpolationPolynom W_t;

}

#endif // TYPES_H
