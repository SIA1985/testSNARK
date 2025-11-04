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

    auto proof = snrk::ZeroTestProof::forProver(funcT, funcV, tG);

    return proof->check();
}

bool correctGates(const snrk::T_t &t, const snrk::S_t &s, snrk::TG_t tG)
{
    auto tCanonic = t.toCanonicPolynom();
    auto tCanonic3wPlus1 = tCanonic(snrk::CanonicPolynom::generate({1, 3}));
    auto tCanonic3wPlus2 = tCanonic(snrk::CanonicPolynom::generate({2, 3}));
    auto tCanonic3wPlus3 = tCanonic(snrk::CanonicPolynom::generate({3, 3}));

    auto funcF = snrk::CanonicPolynom::generate({0});
    auto sCanonic = s.toCanonicPolynom();

    FOROPS {
        auto currentOperation = operation;

        snrk::dots_t dots;
        FOROPS {
            dots.push_back(currentOperation == operation ? snrk::dot_t{operation, 1} : snrk::dot_t{operation, 0});
        }

        auto isOperation = snrk::InterpolationPolynom::generate(dots).toCanonicPolynom();

        switch(currentOperation) {
        case snrk::Sum: {
            funcF += (tCanonic3wPlus1 + tCanonic3wPlus2) * isOperation(sCanonic);
            break;
        }
        case snrk::Product: {
            funcF += (tCanonic3wPlus1 * tCanonic3wPlus2) * isOperation(sCanonic);
            break;
        }
        default:
            break;
        }
    }

    auto a = 3;
    std::cout << tCanonic3wPlus1(a) + tCanonic3wPlus2(a) << " " << tCanonic3wPlus3(a) << std::endl;
    std::cout << funcF(a) << " " << tCanonic3wPlus3(a) << std::endl;

    auto proof = snrk::ZeroTestProof::forProver(funcF, tCanonic3wPlus3, tG);

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
        std::cout << "Некорректные переходы!" << std::endl;
        return 1;
    }

//    std::cout << (snrk::CanonicPolynom::generate({1, 2, 3}) + snrk::CanonicPolynom::generate({-1, -2, -3}))(1) << std::endl;

    std::cout << "Ok!" << std::endl;
}

/*todo:
 * 1. Если t == x, тогда PolynomSubstitutionProof.check() выдаёт 0!
 * 2. Генерация свидетелей в ZeroTestProof с 1, если для 7 этапа, то нужно передать номер свидетеля, либо сделать отдельный класс Proof
 * 3. Графически (в комментариях) представить таблицу (начиная с 1 и тп)
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
