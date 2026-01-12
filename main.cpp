#define THREAD_SPLINE 3

#include "snark.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <csignal>
#include "crypto.h"


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

//    auto x1 = snrk::value_t(5);
//    auto x2 = snrk::value_t(6);
//    auto w1 = snrk::value_t(1);

//    snrk::Circut c({x1, x2}, {w1});

//    for(int i = 0; i < 100; i++) {
//    auto out1 = snrk::value_t(5);
//    c.addGate({snrk::Devide, {x1, w1}, out1});
//    auto out2 = snrk::value_t(1);
//    c.addGate({snrk::Devide, {x1, out1}, out2});
//    auto out3 = snrk::value_t(5);
//    c.addGate({snrk::Product, {x1, w1}, out3});
//    auto out4 = snrk::value_t(4);
//    c.addGate({snrk::Minus, {x1, w1}, out4});
//    auto out5 = snrk::value_t(6);
//    c.addGate({snrk::Sum, {x1, w1}, out5});
//    auto out6 = snrk::value_t(5);
//    c.addGate({snrk::Devide, {x1, w1}, out6});
//    }

    //todo: GPK
    initPairing(mcl::BLS12_381); //временно

    mcl::G1 g;
    mcl::mapToG1(g, 1);
    mcl::Fr t = 3;
    snrk::GPK_t GPK;
    GPK.push_back(g);

    for (int i = 0; i < 2; i++) {
        mcl::G1::mul(g, g, t);
        GPK.push_back(g);
    }

    snrk::CanonicPolynom p({1, 2, 3});
    std::cout << p.commit(GPK) << std::endl;


//    snrk::GlobalParams gp(c, GPK);

//    snrk::ProverProof proof(gp, {x1, x2, {w1}}, {5});

//    if (proof.check()) {
//        std::cout << "Ok!" << std::endl;
//    } else {
//        std::cout << "Not ok!" << std::endl;
//    }
}

/*todo:
 * ! Вопрос: проверить, правильно, ли что q(m_R) выдаёт иногда нуль + посмотреть на comQ!
 * ! Вопрос: почему Partition > 3 не работает
 * (Тут дело в если точек меньше, чем Partition(=4), то раз на входе 3 точки -> полином(^2)/полином(^3))
 *
 * 1. Доразобраться с gp
 * 2. Commit для фукнции с 1м сплайном
 * 3. Сделать схему PLONK
 * 4. Агрегация коммитов
 * 5. Адаптация под PLONK
 *
 * Зависимости:
 * mcl, gmp++, nlohmanjson
*/
/* Спринт 1(Plonk):
 * ЭТАПЫ
 * [V]1. Получение С - скорее в табличном виде
 * [V]2. Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * [V]3. P строит T(x) и получает comt
 * [V]4. T корректно кодирует входы
 * [V]5. Каждый вентиль корректно посчитан (иначе не создает proof)
 * [V]6. Стрелки соответствуют С
 * [V]7. Выход последнего вентиля =0 (как-то через ZeroTest?, скорее через Substitution)
 * [V]8. Оптимизации (сложностей О)
 * [V]9. Распраллеливание в нагруженных функциях
 *
 * Спринт 2 (VM):
 * []1. Расширить круг операций (+ деление [V], логические операции, побитовые операии?)
 * []2. Оформить библиотеку plonk для виртуальной машины
 * []3. Разработка виртуальной машины (разбить на задачи для этого спринта)
 * (модули: [вход, доказательство] -> [алгоритм, схема] -> [выход, доказательство])
*/
