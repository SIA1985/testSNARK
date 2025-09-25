#include "snark.h"
#include <iostream>

int main(int argc, char *argv[])
{
    /*todo: пример из тетради проверить*/
    snrk::values_t inputX = {5, 6};
    snrk::values_t inputW = {1};

    snrk::Circut c(inputX, inputW);

    c.addGate({snrk::Gate::Sum, inputX, snrk::value_t(11)});

    auto gp = snrk::setup(c);
}

/* ЭТАПЫ
 * []1. Получение С - скорее в табличном виде
 * []2. (мб начать с этого?) Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * []3. P строит T(x) и получает comt
 * []4. T корректно кодирует входы
 * []5. Каждый вентиль корректно посчитан
 * []6. Стрелки соответствуют С
 * []7. Выход последнего вентиля =0
*/
