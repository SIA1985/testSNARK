#include "snark.h"
#include <iostream>

int main(int argc, char *argv[])
{
    /*todo: пример из тетради проверить*/
//    auto x1 = snrk::Value(5);
//    auto x2 = snrk::Value(6);
//    auto w1 = snrk::Value(1);

//    snrk::Circut c({x1, x2}, {w1});

//    auto out1 = snrk::Value(11);
//    c.addGate({snrk::Sum, {x1, x2}, out1});
//    auto out2 = snrk::Value(7);
//    c.addGate({snrk::Sum, {x2, w1}, out2});
//    auto out3 = snrk::Value(77);
//    c.addGate({snrk::Product, {out1, out2}, {out3}});

//    snrk::GlobalParams gp(c);

//    auto funcT = gp.PP().t;
//    snrk::dot_t dot{1, 6};

//    auto proof = snrk::PolynomSubstitutionProof::forProver(funcT, dot, gp.TG());


//    std::cout << proof->check() << " " << funcT(dot.x) << std::endl;

    auto f = snrk::ClassicPolynom::generate({0, 1, 2, 3.5});
    auto v = snrk::ClassicPolynom::generate({5, 2, 2, 3, 4, 5, 6});

    std::cout << (f / v)(3) << std::endl;
}

/*todo:
 * 1. Если t == x, тогда PolynomSubstitutionProof.check() выдаёт 0!
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
