#include "snark.h"
#include <iostream>
#include <iomanip>
#include <cmath>

bool correctInputs(const snrk::T_t &t, snrk::values_t inputs, const snrk::witnesses_t &ws, snrk::TG_t tG)
{
    auto funcT = t.toPartedCanonicPolynom();

    snrk::dots_t inputsW;
    snrk::xs_t witness;
    for(std::size_t i = 0; i < inputs.size(); i++) {
        witness.insert(ws[i]);
        inputsW.push_back({ws[i], inputs[i]});
    }
    auto funcV = snrk::InterpolationPolynom::generate(inputsW).toPartedCanonicPolynom();

    auto proof = snrk::ZeroTestProof::forProver(funcT, funcV, tG, witness);

    return proof->check();
}

bool correctGates(const snrk::SplittedT_t &t, const snrk::S_t &s, const snrk::witnesses_t &ws, snrk::TG_t tG)
{
    auto left = t.left.toPartedCanonicPolynom();
    auto right = t.right.toPartedCanonicPolynom();
    auto result = t.result.toPartedCanonicPolynom();

    auto funcF = snrk::PartedCanonicPolynom::generate(snrk::PartedCanonicPolynom::map{});
    auto sCanonic = s.toPartedCanonicPolynom();

    snrk::xs_t witness;
    for(const auto &w : ws) {
        witness.insert(w);
    }

    FOROPS {
        //todo: вынести в генерацию S
        snrk::dots_t dots;
        for(const auto &w : witness) {
            dots.push_back(operation == (snrk::GateType_t)std::round(sCanonic(w).get_d()) ? snrk::dot_t{w, 1} : snrk::dot_t{w, 0});
        }

        auto isOperation = snrk::InterpolationPolynom::generate(dots).toPartedCanonicPolynom();

        switch(operation) {
        case snrk::Sum: {
            funcF += (left + right) * isOperation;
            break;
        }
        case snrk::Product: {
            funcF += (left * right) * isOperation;
            break;
        }
        case snrk::Minus: {
            funcF += (left - right) * isOperation;
            break;
        }
            //todo: после деления имеем CustomPolynom (если кол-во опреаций деления != 0 -> выполнять)
        case snrk::Devide: {
//            funcF += (snrk::InterpolationPolynom::generate((left / right).dots(witness))).toPartedCanonicPolynom() * isOperation(sCanonic);
            break;
        }
        default:
            break;
        }
    }

    auto proof = snrk::ZeroTestProof::forProver(funcF, result, tG, witness);

    return proof->check();
}

bool currentOutput(const snrk::T_t &t, snrk::value_t output, std::size_t lastWNum, snrk::TG_t tG) {
    auto tCanonic = t.toPartedCanonicPolynom();
    //todo: изменить получение степени
    //todo: последний свидетель
    auto outputDot = snrk::dot_t{lastWNum, output};

    auto proof = snrk::PolynomSubstitutionProof::forProver(tCanonic, outputDot, tG);

    return proof->check();
}

/*
 * Если все переходы и входы корректны, то в любом свидетеле будет пройден ZeroTest,
 * если какой-то переход неверен, то в его свидетеле rem != 0 -> деление на 0 -> ошибка.
 * Выходит, нужно проверять по 1му свидетелю для каждого диапазона. (при partition = 2 <=> SubTest)
 */

int main(int argc, char *argv[])
{
    /*todo: пример из тетради проверить*/
    auto x1 = snrk::Value(5);
    auto x2 = snrk::Value(6);
    auto w1 = snrk::Value(1);

    snrk::Circut c({x1, x2}, {w1});

    auto out1 = snrk::Value(11);
    c.addGate({snrk::Sum, {x1, x2}, {out1}});
    auto out2 = snrk::Value(7);
    auto out3 = snrk::Value(77);
    for(int i = 0; i < 100; i++) {
        c.addGate({snrk::Sum, {x2, {w1}}, {out2}});

        c.addGate({snrk::Product, {out1, out2}, {out3}});
    }
    auto out4 = snrk::Value(70);
    c.addGate({snrk::Minus, {out3, out2}, {out4}});

//    c.addGate({snrk::Sum, {{1}, {0}}, {1}});

    snrk::GlobalParams gp(c);

    if (!correctInputs(gp.PP().t, {x1, x2, {w1}}, gp.witnesses(), gp.TG())) {
        std::cout << "Некорректные входы!" << std::endl;
        return 1;
    }

    if (!correctGates(gp.PP().splittedT, gp.PP().s, gp.SWitnesses(), gp.TG())) {
        std::cout << "Некорректные переходы!" << std::endl;
        return 1;
    }

    //todo: 6.

    if (!currentOutput(gp.PP().t, out4, gp.witnesses().size(), gp.TG())) {
        std::cout << "Некорректный выход!" << std::endl;
        return 1;
    }

    //todo: нет assert на размер диапазона
//    snrk::PartedCanonicPolynom::map m;
//    m.insert(snrk::Range{1, 3}, snrk::CanonicPolynom::generate({1}));
//    m.insert(snrk::Range{3, 5}, snrk::CanonicPolynom::generate({2}));

//    auto a = snrk::PartedCanonicPolynom::generate(m);

//    snrk::PartedCanonicPolynom::map m2;
//    m2.insert(snrk::Range{3, 5}, snrk::CanonicPolynom::generate({1})); //4*4 + 4*6 == 40
//    m2.insert(snrk::Range{5, 7}, snrk::CanonicPolynom::generate({4}));

//    auto b = snrk::PartedCanonicPolynom::generate(m2);

//    std::cout << (a / b)(3) << " " << a(3) / b(3) << std::endl;

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
