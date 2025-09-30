#include "snark.h"
#include <iostream>

/*todo:*/
int t = 1;
int G = 2;

int main(int argc, char *argv[])
{
    /*todo: пример из тетради проверить*/
    auto x1 = Value(5);
    auto x2 = Value(6);
    auto w1 = Value(1);

    snrk::Circut c({x1, x2}, {w1});

    auto out1 = Value(11);
    c.addGate({snrk::Gate::Sum, {x1, x2}, out1});
    auto out2 = Value(7);
    c.addGate({snrk::Gate::Sum, {x2, w1}, out2});
    auto out3 = Value(77);
    c.addGate({snrk::Gate::Product, {out1, out2}, {out3}});

    auto gp = snrk::setup(c);

    auto funcT = snrk::lagrange<int, int>::generate({{1, 3}, {2, 5}, {4, 2}});

    std::cout << funcT.commit(t, G) << std::endl;
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
