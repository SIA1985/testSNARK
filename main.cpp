#include "snark.h"
#include <iostream>

/*todo:*/
using X = float;
using Y = float;

X t = 1;
int G = 2;

int main(int argc, char *argv[])
{
    /*todo: пример из тетради проверить*/
    auto x1 = snrk::Value(5);
    auto x2 = snrk::Value(6);
    auto w1 = snrk::Value(1);

    snrk::Circut c({x1, x2}, {w1});

    auto out1 = snrk::Value(11);
    c.addGate({snrk::Gate::Sum, {x1, x2}, out1});
    auto out2 = snrk::Value(7);
    c.addGate({snrk::Gate::Sum, {x2, w1}, out2});
    auto out3 = snrk::Value(77);
    c.addGate({snrk::Gate::Product, {out1, out2}, {out3}});

    snrk::GlobalParams gp(c);

    auto funcT = snrk::Lagrange::generate({{1, 3}, {2, 5}, {4, 2}});
    auto funcQ = snrk::CustomPolynom::generate([&funcT](snrk::X_t x) -> snrk::Y_t
    {
        return (funcT(x) - 2) / (x - 4);
    });


    snrk::PolynomSubstitutionProof<X, Y> proof(funcT.commit(t, G), funcQ.commit(t, G), {.u = 4, .v = 2});

    proof.setGp(t, G);

    std::cout << proof.check() << std::endl;
}

/* ЭТАПЫ
 * [V]1. Получение С - скорее в табличном виде
 * [V]2. Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * [50/50]3. P строит T(x) и получает comt (Доразобраться с gp и сделать норм. commit)
 * []4. T корректно кодирует входы
 * []5. Каждый вентиль корректно посчитан
 * []6. Стрелки соответствуют С
 * []7. Выход последнего вентиля =0
 * []8. Оптимизации (сложностей О)
*/
