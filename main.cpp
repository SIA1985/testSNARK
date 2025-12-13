#include "snark.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <csignal>


void sigFpeHandler(int signum) {
    std::cout << "Ошибка в создании доказательства!" << std::endl;
    exit(signum);
}

void sigTermHandler(int signum) {
    std::cout << "Внутренняя ошибка!" << std::endl;
    exit(signum);
}


/*
 * Схема (С):
 *
 *----------------------
 * w1==1 | w2==2 | w3==3 - Входы
 * ---------------------
 * w4==4 | w5==5 | w6==6
 * --------------------- - Алгоримт схемы
 *  ...  |  ...  |  ...
 * ---------------------
 *  ...  |  ...  | wn==n - Выход
 * ---------------------
*/

int main(int argc, char *argv[])
{
    if (std::signal(SIGFPE, sigFpeHandler) == SIG_ERR) {
        exit(2);
    }

    if (std::signal(SIGTERM, sigTermHandler) == SIG_ERR) {
        exit(2);
    }

    mpf_set_default_prec(128);

    /*todo: пример из тетради проверить*/
    auto x1 = snrk::value_t(5);
    auto x2 = snrk::value_t(6);
    auto w1 = snrk::value_t(1);

    snrk::Circut c({x1, x2}, {w1});

    for(int i = 0; i < 5000; i++) {
    auto out1 = snrk::value_t(11);
    c.addGate({snrk::Sum, {x1, x2}, {out1}});
    auto out2 = snrk::value_t(7);
    c.addGate({snrk::Sum, {x2, {w1}}, {out2}});
    auto out3 = snrk::value_t(77);
    c.addGate({snrk::Product, {out1, out2}, {out3}});
    auto out4 = snrk::value_t(70);
    c.addGate({snrk::Minus, {out3, out2}, {out4}});
    auto out5 = snrk::value_t(7);
    c.addGate({snrk::Minus, {out3, out4}, {out5}});
    auto out6 = snrk::value_t(490);
    c.addGate({snrk::Product, {out4, out5}, {out6}});
    auto out7 = snrk::value_t(-4);
    c.addGate({snrk::Minus, {out2, out1}, {out7}});
    auto out8 = snrk::value_t(-28);
    c.addGate({snrk::Product, {out5, out7}, {out8}});
    auto out9 = snrk::value_t(42);
    c.addGate({snrk::Sum, {out8, out4}, {out9}});
    auto out10 = snrk::value_t(14);
    c.addGate({snrk::Sum, {out9, out8}, {out10}});
    }

    snrk::GlobalParams gp(c);

    snrk::ProverProof proof(gp, {x1, x2, {w1}}, {14});

    if (proof.check()) {
        std::cout << "Ok!" << std::endl;
    } else {
        std::cout << "Not ok!" << std::endl;
    }
}

/*todo:
 * ! Вопрос: проверить, правильно, ли что q(m_R) выдаёт иногда нуль + посмотреть на comQ!
 * ! Вопрос: почему Partition > 3 не работает
 * (Тут дело в если точек меньше, чем Partition(=4), то раз на входе 3 точки -> полином(2)/полином(3))
 *
 * 1. Доразобраться с gp и сделать норм. commit
 * !2. Распараллелить вычисления в сплайновый (при создании 0-полинома долго)
*/
/* ЭТАПЫ
 * [V]1. Получение С - скорее в табличном виде
 * [V]2. Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * [V]3. P строит T(x) и получает comt
 * [V]4. T корректно кодирует входы
 * [V]5. Каждый вентиль корректно посчитан (иначе не создает proof)
 * [V]6. Стрелки соответствуют С
 * [V]7. Выход последнего вентиля =0 (как-то через ZeroTest?, скорее через Substitution)
 * [V]8. Оптимизации (сложностей О)
 * []9. Распраллеливание в нагруженных функциях
*/
