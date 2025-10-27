#include "snark.h"
#include <iostream>

int main(int argc, char *argv[])
{
    mpf_set_default_prec(256);

    /*todo: пример из тетради проверить*/
    auto x1 = snrk::Value(5);
    auto x2 = snrk::Value(6);
    auto w1 = snrk::Value(1);

    snrk::Circut c({x1, x2}, {w1});

    auto out1 = snrk::Value(11);
    c.addGate({snrk::Sum, {x1, x2}, out1});
    auto out2 = snrk::Value(7);
    c.addGate({snrk::Sum, {x2, w1}, out2});
    auto out3 = snrk::Value(77);
    c.addGate({snrk::Product, {out1, out2}, {out3}});

    snrk::GlobalParams gp(c);

    auto funcT = gp.PP().t.toCanonicPolynom();
    auto funcV = snrk::InterpolationPolynom::generate({{1, 5}, {2, 6}, {3, 1}, {4, 5}, {5, 6}, {6, 11}, {7, 6}, {8, 1}, {9, 7}, {10, 11}, {11, 7}, {12, 77}}).toCanonicPolynom();

//    auto funcT = snrk::InterpolationPolynom::generate({{1, 0}, {2, 2}, {3, 5}, {4, 70}}).toCanonicPolynom();
//    auto funcV = snrk::InterpolationPolynom::generate({{1, 2}, {2, 5}, {3, 4}}).toCanonicPolynom();

    auto proof = snrk::ZeroTestProof::forProver(funcT, funcV, gp.TG());

    std::cout << proof->check() << std::endl;
}

/*todo:
 * 1. Если t == x, тогда PolynomSubstitutionProof.check() выдаёт 0!
 *
*/
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
