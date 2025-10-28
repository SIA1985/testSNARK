#include "snark.h"
#include <iostream>

bool correctInputs(const snrk::T_t &t, snrk::values_t inputs, snrk::TG_t tG)
{
    auto funcT = t.toCanonicPolynom();

    snrk::dots_t inputsW;
    for(std::size_t i = 0; i < inputs.size(); i++) {
        inputsW.push_back({i + 1, inputs[i]});
    }
    auto funcV = snrk::InterpolationPolynom::generate(inputsW).toCanonicPolynom();

    auto proof = snrk::ZeroTestProof::forProver(funcT, funcV, tG, snrk::wGeneratorDefault);

    return proof->check();
}

bool correctGates(const snrk::T_t &t, const snrk::S_t &s, snrk::TG_t tG)
{
    /* До этого t и s перевести в Canonic !
     * auto funcF = [s(x) == Sum]*(t(x) + t(x+1)) + [s(x) == Product]*(t(x) * t(x+1)) */;
    /*[s(x) == Sum] <=> InterpolationPolynom({{Sum, 1}, {Product, 0}}) либо отдельный тип полинома (по сути конструктор для операций)*/

    auto funcG = s.toCanonicPolynom();

    auto proof = snrk::ZeroTestProof::forProver(funcF, funcG, tG);

    return proof->check();
}

int main(int argc, char *argv[])
{
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

    if (!correctInputs(gp.PP().t, {x1, x2, w1}, gp.TG())) {
        std::cout << "Некорректные входы!" << std::endl;
        return 1;
    }

    if (!correctGates(gp.PP().t, gp.PP().s, gp.TG())) {
        std::cout << std::endl;
        return 1;
    }

    std::cout << "Ok!" << std::endl;
}

/*todo:
 * 1. Если t == x, тогда PolynomSubstitutionProof.check() выдаёт 0!
 * 2. Генерация свидетелей в ZeroTestProof с 1, если для 7 этапа, то нужно передать номер свидетеля, либо сделать отдельный класс Proof
 *
 * (мало ли mpf_set_default_prec(256);)
*/
/* ЭТАПЫ
 * [V]1. Получение С - скорее в табличном виде
 * [V]2. Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * [50/50]3. P строит T(x) и получает comt (Доразобраться с gp и сделать норм. commit)
 * [V]4. T корректно кодирует входы
 * []5. Каждый вентиль корректно посчитан
 * []6. Стрелки соответствуют С
 * []7. Выход последнего вентиля =0 (как-то через ZeroTest?)
 * []8. Оптимизации (сложностей О)
*/
