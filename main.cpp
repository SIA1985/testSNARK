#include "snark.h"
#include <iostream>
#include <iomanip>

bool correctInputs(const snrk::T_t &t, snrk::values_t inputs, snrk::TG_t tG)
{
    auto funcT = t.toPartedCanonicPolynom();

    snrk::dots_t inputsW;
    for(std::size_t i = 0; i < inputs.size(); i++) {
        inputsW.push_back({i + 1, inputs[i]});
    }
    auto funcV = snrk::InterpolationPolynom::generate(inputsW).toPartedCanonicPolynom();

    snrk::xs_t witness;
    for(std::size_t i = 1; i <= inputs.size(); i++) {
        witness.insert(i);
    }

    auto proof = snrk::ZeroTestProof::forProver(funcT, funcV, tG, witness);

    return proof->check();
}

bool correctGates(const snrk::T_t &t, const snrk::S_t &s, snrk::TG_t tG)
{
    auto tCanonic = t.toPartedCanonicPolynom();
    auto tCanonic3wPlus1 = tCanonic(snrk::CanonicPolynom::generate({1, 3}));
    auto tCanonic3wPlus2 = tCanonic(snrk::CanonicPolynom::generate({2, 3}));
    auto tCanonic3wPlus3 = tCanonic(snrk::CanonicPolynom::generate({3, 3}));

    auto funcF = snrk::PartedCanonicPolynom::generate(snrk::PartedCanonicPolynom::map{});
    auto sCanonic = s.toPartedCanonicPolynom();

    FOROPS {
        //todo: вынести в генерацию S
        snrk::dots_t dots;
        for(std::size_t i = 1; i <= s.toCanonicPolynom().degree() + 1; i++) {
            dots.push_back(operation == sCanonic(i) ? snrk::dot_t{i, 1} : snrk::dot_t{i, 0});
        }

        auto isOperation = snrk::InterpolationPolynom::generate(dots).toPartedCanonicPolynom();

        switch(operation) {
        case snrk::Sum: {
            funcF += (tCanonic3wPlus1 + tCanonic3wPlus2) * isOperation;
            break;
        }
        case snrk::Product: {
            funcF += (tCanonic3wPlus1 * tCanonic3wPlus2) * isOperation;
            break;
        }
        case snrk::Minus: {
            funcF += (tCanonic3wPlus1 - tCanonic3wPlus2) * isOperation;
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
    for(std::size_t i = 1; i <= s.toCanonicPolynom().degree() + 1; i++) {
        witness.insert(i);
    }

    auto proof = snrk::ZeroTestProof::forProver(funcF, tCanonic3wPlus3, tG, witness);

    return proof->check();
}

bool currentOutput(const snrk::T_t &t, snrk::value_t output, snrk::TG_t tG) {
    auto tCanonic = t.toPartedCanonicPolynom();
    auto outputDot = snrk::dot_t{t.toCanonicPolynom().degree() + 1, output};

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
    auto out4 = snrk::Value(70);
    c.addGate({snrk::Minus, {out3, out2}, {out4}});

    snrk::GlobalParams gp(c);

//    if (!correctInputs(gp.PP().t, {x1, x2, {w1}}, gp.TG())) {
//        std::cout << "Некорректные входы!" << std::endl;
//        return 1;
//    }

    //вернуть предыдущую isOperation и добавить решение уравнение 2й степени
    if (!correctGates(gp.PP().t, gp.PP().s, gp.TG())) {
        std::cout << "Некорректные переходы!" << std::endl;
        return 1;
    }

    //todo: 6.

    if (!currentOutput(gp.PP().t, out4, gp.TG())) {
        std::cout << "Некорректный выход!" << std::endl;
        return 1;
    }

    //todo: нет assert на размер диапазона
//    snrk::PartedCanonicPolynom::map m;
//    m.insert(snrk::Range{1, 3.5}, snrk::CanonicPolynom::generate({1}));
//    m.insert(snrk::Range{3.5, 4.5}, snrk::CanonicPolynom::generate({2}));
//    m.insert(snrk::Range{4.5, 5}, snrk::CanonicPolynom::generate({3}));

//    auto a = snrk::PartedCanonicPolynom::generate(m);

//    snrk::PartedCanonicPolynom::map m2;
//    m2.insert(snrk::Range{3, 5}, snrk::CanonicPolynom::generate({1})); //4*4 + 4*6 == 40
//    m2.insert(snrk::Range{5, 7}, snrk::CanonicPolynom::generate({4}));

//    auto b = snrk::PartedCanonicPolynom::generate(m2);

//    std::cout << (a / b)(4) << " " << a(4) / b(4) << std::endl;

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
 * 7. Вынести типы в types.h
 * 8. assert на непрерывные диапазоны и их длинну из Range и RangeMap в классы, что в таких нуждаются
 * 9. Опять 2й проход по некоторым диапазонам
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
