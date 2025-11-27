#include "snark.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <csignal>

bool correctInputs(const snrk::T_t &t, snrk::values_t inputs, const snrk::witnesses_t &ws, snrk::TG_t tG)
{
    auto funcT = t.toPartedCanonicPolynom();

    snrk::dots_t inputsW;
    snrk::xs_t witness;
    for(std::size_t i = 0; i < inputs.size(); i++) {
        witness.insert(ws[i]);
        inputsW.push_back({ws[i], inputs[i]});
    }
    auto funcV = snrk::InterpolationPolynom(inputsW).toPartedCanonicPolynom();

    auto proof = snrk::ZeroTestProof::forProver(funcT, funcV, tG, witness);

    return proof->check();
}

bool correctGates(const snrk::SplittedT_t &t, const snrk::GlobalParams::SParams_t SParams, const snrk::witnesses_t &ws, snrk::TG_t tG)
{
    auto left = t.left.toPartedCanonicPolynom();
    auto right = t.right.toPartedCanonicPolynom();
    auto result = t.result.toPartedCanonicPolynom();

    auto funcF = snrk::PartedCanonicPolynom(snrk::PartedCanonicPolynom::map{});

    snrk::xs_t witness;
    for(const auto &w : ws) {
        witness.insert(w);
    }

    for(const auto &[operation, dots] : SParams.opsFromS) {
        auto isOperation = snrk::InterpolationPolynom(dots).toPartedCanonicPolynom();

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
            //todo:
        case snrk::Devide: {
//            funcF += (snrk::InterpolationPolynom((left / right).dots(witness))).toPartedCanonicPolynom() * isOperation;
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

    auto outputDot = snrk::dot_t{lastWNum, output};

    auto proof = snrk::PolynomSubstitutionProof::forProver(tCanonic, outputDot, tG);

    return proof->check();
}

void sigfpeHandler(int signum) {
    std::cout << "Ошибка в создании доказательства!" << std::endl;
    exit(signum);
}


int main(int argc, char *argv[])
{
    if (std::signal(SIGFPE, sigfpeHandler) == SIG_ERR) {
        exit(2);
    }

    /*todo: пример из тетради проверить*/
    auto x1 = snrk::Value(5);
    auto x2 = snrk::Value(6);
    auto w1 = snrk::Value(1);

    snrk::Circut c({x1, x2}, {w1});

    auto out1 = snrk::Value(11);
    c.addGate({snrk::Sum, {x1, x2}, {out1}});
    auto out2 = snrk::Value(7);
    c.addGate({snrk::Sum, {x2, {w1}}, {out2}});
    auto out3 = snrk::Value(77);
    c.addGate({snrk::Product, {out1, out2}, {out3}});
    auto out4 = snrk::Value(70);
    c.addGate({snrk::Minus, {out3, out2}, {out4}});

    snrk::GlobalParams gp(c);
    auto TParams = gp.PP().TParams;

    if (!correctInputs(TParams.t, {x1, x2, {w1}}, gp.witnesses(), gp.TG())) {
        std::cout << "Некорректные входы!" << std::endl;
        exit(1);
    }

    //todo: если падает -> не построить полином -> не получить доказательство
    if (!correctGates(TParams.splittedT, gp.PP().SParams, gp.SWitnesses(), gp.TG())) {
        std::cout << "Некорректные переходы!" << std::endl;
        exit(1);
    }

    //todo: 6.

    if (!currentOutput(TParams.t, out4, gp.witnesses().size(), gp.TG())) {
        std::cout << "Некорректный выход!" << std::endl;
        exit(1);
    }

    std::cout << "Ok!" << std::endl;
}

/*todo:
 * 1. Если t == x, тогда PolynomSubstitutionProof.check() выдаёт 0!
 * 2. Графически (в комментариях) представить таблицу (начиная с 1 и тп)
 * 3. Доразобраться с gp и сделать норм. commit
 * 4. assert на непрерывные диапазоны и их длинну из Range и RangeMap в классы, что в таких нуждаются
 * 5. Разобраться с делением
*/
/* ЭТАПЫ
 * [V]1. Получение С - скорее в табличном виде
 * [V]2. Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * [V]3. P строит T(x) и получает comt
 * [V]4. T корректно кодирует входы
 * [V]5. Каждый вентиль корректно посчитан (иначе не создает proof)
 * []6. Стрелки соответствуют С
 * [V]7. Выход последнего вентиля =0 (как-то через ZeroTest?, скорее через Substitution)
 * []8. Оптимизации (сложностей О)
*/
