#include "snark.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <csignal>
#include <chrono>

#define printDur(text, end, start)     std::cout << text << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

//Подготовка увеличивается согласно О(n^2)
bool correctInputs(const snrk::T_t &t, snrk::values_t inputs, const snrk::witnesses_t &ws, snrk::TG_t tG)
{
    auto start = std::chrono::steady_clock::now();

    snrk::xs_t witness;
    snrk::dots_t inputsW;
    inputsW.reserve(inputs.size());
    for(std::size_t i = 0; i < inputs.size(); i++) {
        witness.insert(ws[i]);
        inputsW.push_back({ws[i], inputs[i]});
    }
    auto funcV = snrk::InterpolationPolynom(inputsW).toPartedCanonicPolynom();
    auto funcTCut = t.toPartedCanonicPolynom().cut(funcV.distance());

    auto end = std::chrono::steady_clock::now();
    printDur("Подготовка correctInputs: ", end, start);

    start = std::chrono::steady_clock::now();
    auto proof = snrk::ZeroTestProof::forProver(funcTCut, funcV, tG, witness);
    auto result = proof->check();
    end = std::chrono::steady_clock::now();
    printDur("Проверка correctInputs: ", end, start);

    return result;
}

//Подготовка увеличивается согласно О(n^2)
bool correctGates(const snrk::SplittedT_t &t, const snrk::GlobalParams::SParams_t SParams, const snrk::witnesses_t &ws, snrk::TG_t tG)
{
    auto start = std::chrono::steady_clock::now();

    auto left = t.left.toPartedCanonicPolynom();
    auto right = t.right.toPartedCanonicPolynom();
    auto result = t.result.toPartedCanonicPolynom();

    snrk::xs_t witnesses(ws.begin(), ws.end());
    //400ms - 10k свидетелей

    auto funcF = snrk::PartedCanonicPolynom(snrk::PartedCanonicPolynom::map{});
    for(const auto &[operation, dots] : SParams.opsFromS) {

        //
        auto isOperation = snrk::InterpolationPolynom(dots).toPartedCanonicPolynom();
        //150ms - 10k свидетелей!

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
        //60ms - 10k свидетелей

    }
    //1100ms - 10k свидетелей

    auto end = std::chrono::steady_clock::now();
    printDur("Подготовка correctGates: ", end, start);

    //
    auto proof = snrk::ZeroTestProof::forProver(funcF, result, tG, witnesses);
    //150ms - 10k свидетелей

    return proof->check();
    //150ms - 10k свидетелей
}

//мб дело в том, что надо проверять не t, а такое t, что выводит адреса
bool currentVars(const snrk::WT_t &wt, const snrk::T_t &t, const snrk::witnesses_t ws, snrk::TG_t tG) {
    auto start = std::chrono::steady_clock::now();

    auto tCanonic = t.toPartedCanonicPolynom();
    auto wtCanonic = wt.toPartedCanonicPolynom();

    snrk::xs_t witnesses(ws.begin(), ws.end());


    auto end = std::chrono::steady_clock::now();
    printDur("Подготовка currentVars: ", end, start);

    start = std::chrono::steady_clock::now();
    auto proof = snrk::ZeroTestProof::forProver(tCanonic, wtCanonic, tG, witnesses);
    auto result = proof->check();
    end = std::chrono::steady_clock::now();
    printDur("Проверка currentVars: ", end, start);

    return result;
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

    mpf_set_default_prec(128);

    /*todo: пример из тетради проверить*/
    auto x1 = snrk::Value(5);
    auto x2 = snrk::Value(6);
    auto w1 = snrk::Value(1);

    snrk::Circut c({x1, x2}, {w1});

    //15000ms - 30k свидетелей (< 1000)
    for(int i = 0; i < 100; i++) {
    auto out1 = snrk::Value(11);
    c.addGate({snrk::Sum, {x1, x2}, {out1}});
    auto out2 = snrk::Value(7);
    c.addGate({snrk::Sum, {x2, {w1}}, {out2}});
    auto out3 = snrk::Value(77);
    c.addGate({snrk::Product, {out1, out2}, {out3}});
    auto out4 = snrk::Value(70);
    c.addGate({snrk::Minus, {out3, out2}, {out4}});
    auto out5 = snrk::Value(7);
    c.addGate({snrk::Minus, {out3, out4}, {out5}});
    auto out6 = snrk::Value(490);
    c.addGate({snrk::Product, {out4, out5}, {out6}});
    auto out7 = snrk::Value(-4);
    c.addGate({snrk::Minus, {out2, out1}, {out7}});
    auto out8 = snrk::Value(-28);
    c.addGate({snrk::Product, {out5, out7}, {out8}});
    auto out9 = snrk::Value(42);
    c.addGate({snrk::Sum, {out8, out4}, {out9}});
    auto out10 = snrk::Value(14);
    c.addGate({snrk::Sum, {out9, out8}, {out10}});
    }

    snrk::GlobalParams gp(c);
    auto TParams = gp.PP().TParams;
    auto witnesses = gp.witnesses();
    auto tG = gp.TG();

    if (!correctInputs(TParams.t, {x1, x2, {w1}}, witnesses, tG)) {
        std::cout << "Некорректные входы!" << std::endl;
        exit(1);
    }

    if (!correctGates(TParams.splittedT, gp.PP().SParams, gp.SWitnesses(), tG)) {
        std::cout << "Некорректные переходы!" << std::endl;
        exit(1);
    }

    if (!currentVars(gp.PP().WParams.wt, TParams.t, witnesses, tG)) {
        std::cout << "Некорректные переменные!" << std::endl;
        exit(1);
    }

    if (!currentOutput(TParams.t, {14}, witnesses.size(), tG)) {
        std::cout << "Некорректный выход!" << std::endl;
        exit(1);
    }

    std::cout << "Ok!" << std::endl;
}

/*todo:
 * ! Баг: на 1318 свидетеле эффект Рунге, хотя используем сплайны ! (вроде пофиксил повышением точности)
 * ! Вопрос: проверить, правильно, ли что q(m_R) выдаёт иногда нуль + посмотреть на comQ!
 *
 * 1. Если t == x, тогда PolynomSubstitutionProof.check() выдаёт 0! (проверка при генерации r, что r != t?)
 * 2. Графически (в комментариях) представить таблицу (начиная с 1 и тп)
 * 3. Доразобраться с gp и сделать норм. commit
 * 4. assert на непрерывные диапазоны и их длинну из Range и RangeMap в классы, что в таких нуждаются
 * 5. Перевод proof в json и обратно
 * !6. Распараллелить вычисления в сплайновый (при создании 0-полинома долго)
 * 7. Убрать свидетелей из доказательства (оставить только количество)
 * 8. Объединить построение T, S, W в один цикл (хотя и так довольно быстро, ибо линейно)
*/
/* ЭТАПЫ
 * [V]1. Получение С - скорее в табличном виде
 * [V]2. Setup(C): генерация S(x) - селекторного полинома и W(o) - полином ограничений на конкретную перестановку
 * [V]3. P строит T(x) и получает comt
 * [V]4. T корректно кодирует входы
 * [V]5. Каждый вентиль корректно посчитан (иначе не создает proof)
 * [V]6. Стрелки соответствуют С
 * [V]7. Выход последнего вентиля =0 (как-то через ZeroTest?, скорее через Substitution)
 * []8. Оптимизации (сложностей О)
*/
