#include "snark.h"
#include <iostream>
#include <iomanip>

bool correctInputs(const snrk::T_t &t, snrk::values_t inputs, snrk::TG_t tG)
{
    auto funcT = t.toCanonicPolynom();

    snrk::dots_t inputsW;
    for(std::size_t i = 0; i < inputs.size(); i++) {
        inputsW.push_back({i + 1, inputs[i]});
    }
    auto funcV = snrk::InterpolationPolynom::generate(inputsW).toCanonicPolynom();

    snrk::xs_t witness;
    for(std::size_t i = 1; i <= funcV.degree() + 1; i++) {
        witness.insert(i);
    }

    auto proof = snrk::ZeroTestProof::forProver(funcT, funcV, tG, witness);

    return proof->check();
}

bool correctGates(const snrk::T_t &t, const snrk::S_t &s, snrk::TG_t tG)
{
    /*Работало с 3мя точками
    auto tCanonic = t.toCanonicPolynom();

    auto sCanonic = s.toCanonicPolynom();

    snrk::dots_t dots3wPlus1, dots3wPlus2, dots3wPlus3;
    for(std::size_t i = 1; i <= sCanonic.degree() + 1; i++) {
        dots3wPlus1.push_back({i, t(3*i + 1)});
        dots3wPlus2.push_back({i, t(3*i + 2)});
        dots3wPlus3.push_back({i, t(3*i + 3)});
    }

    auto tCanonic3wPlus1 = snrk::InterpolationPolynom::generate(dots3wPlus1).toCanonicPolynom();
    auto tCanonic3wPlus2 = snrk::InterpolationPolynom::generate(dots3wPlus2).toCanonicPolynom();
    auto tCanonic3wPlus3 = snrk::InterpolationPolynom::generate(dots3wPlus3).toCanonicPolynom();

    + аналогичное с isOperation (проходимся по S и где операторы совпадают -> 1, иначе -> 0)
    */

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
        case snrk::Minus: {
            funcF += (tCanonic3wPlus1 - tCanonic3wPlus2) * isOperation(sCanonic);
            break;
        }
            //todo: после деления имеем CustomPolynom
//        case snrk::Devide: {
//            funcF += (tCanonic3wPlus1 / tCanonic3wPlus2) * isOperation(sCanonic);
//            break;
//        }
        default:
            break;
        }
    }

    snrk::xs_t witness;
    for(std::size_t i = 1; i <= sCanonic.degree() + 1; i++) {
        witness.insert(i);
        std::cout << std::setprecision(20) << funcF(i) << " " << tCanonic3wPlus3(i) << std::endl;
    }

    //todo: дело определённо в полиноме funcF, так как при вычисление с большей точностью все работает (долнжо по к.м.)
    auto proof = snrk::ZeroTestProof::forProver(funcF, tCanonic3wPlus3, tG, witness);

    return proof->check();
}

bool currentOutput(const snrk::T_t &t, snrk::value_t output, snrk::TG_t tG) {
    auto tCanonic = t.toCanonicPolynom();
    auto outputDot = snrk::dot_t{tCanonic.degree() + 1, output};

    auto proof = snrk::PolynomSubstitutionProof::forProver(tCanonic, outputDot, tG);

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

    //todo: 6.

    if (!currentOutput(gp.PP().t, out3, gp.TG())) {
        std::cout << "Некорректный выход!" << std::endl;
        return 1;
    }

    std::cout << "Ok!" << std::endl;
}

/*todo:
 * 1. Если t == x, тогда PolynomSubstitutionProof.check() выдаёт 0!
 * 2. Генерация свидетелей в ZeroTestProof с 1, если для 7 этапа, то нужно передать номер свидетеля, либо сделать отдельный класс Proof
 * 3. Графически (в комментариях) представить таблицу (начиная с 1 и тп)
 * 4. Доразобраться с gp и сделать норм. commit
 * !5. Заменить (но не удалить, а на его основе) CanonicPolynom на кубические сплайны (сначала попробовать имеющимися средствами, потом из примера в ИИ)!
 * (мало ли mpf_set_default_prec(256);)
 * 6. generate -> constructor
*/
/* ЭТАПЫ
 * [V]1. Получение С - скорее в табличном виде
 * [V]2. Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * [V]3. P строит T(x) и получает comt
 * [V]4. T корректно кодирует входы
 * []5. Каждый вентиль корректно посчитан
 * []6. Стрелки соответствуют С
 * [V]7. Выход последнего вентиля =0 (как-то через ZeroTest?, скорее через Substitution)
 * []8. Оптимизации (сложностей О)
*/
