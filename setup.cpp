#include "setup.h"
#include "operation.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <thread>

namespace snrk {

witnesses_t genWitnesses(witness_t start, std::size_t count, witness_t wStep)
{
    witnesses_t witnesses;
    for(witness_t i = 0; i < count; i++) {
        witnesses.push_back(start);
        start += wStep;
    }

    return witnesses;
}

GlobalParams::GlobalParams(const Circut &circut, const GPK_t &GPK)
    : m_GPK{GPK}
    , m_circutSize{circut.size()}
{
    m_witnesses = genWitnesses(wStart, circut.degree(), wStep);

    std::thread tS(&GlobalParams::generateS, this, std::ref(circut));
    std::thread tW(&GlobalParams::generateW, this, std::ref(circut));

    tS.join();
    tW.join();
}

witnesses_t GlobalParams::witnesses() const
{
    return m_witnesses;
}

witnesses_t GlobalParams::SWitnesses() const
{
    return m_SWitnesses;
}

std::size_t GlobalParams::witnessesCount() const
{
    return m_witnesses.size();
}

GPK_t GlobalParams::GPK() const
{
    return m_GPK;
}

GlobalParams::ProverParams_t GlobalParams::PP() const
{
    return {.TParams = {m_T, m_splittedT}, .SParams = {m_S, m_opsFromS}, .WParams = {m_WT, m_WI}};
}

std::size_t GlobalParams::circutSize() const
{
    return m_circutSize;
}

void GlobalParams::generateT(const Circut &circut)
{
    dots_t dots, leftDots, rightDots, resultDots;

    auto cw = m_witnesses.cbegin();
    auto fillMap = [&cw, &dots](const values_t &row)
    {
        for(const auto& element : row) {
            dots.push_back({*cw, element});
            cw++;
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    m_SWitnesses = genWitnesses(wStart, circut.m_gates.size(), wStep);
    auto i = m_SWitnesses.begin();

    for(const auto& gate : circut.m_gates) {
        dots.push_back({*cw, gate.m_input.a});
        leftDots.push_back({*i, gate.m_input.a});
        cw++;

        dots.push_back({*cw, gate.m_input.b});
        rightDots.push_back({*i, gate.m_input.b});
        cw++;

        dots.push_back({*cw, gate.m_output});
        resultDots.push_back({*i, gate.m_output});
        cw++;

        i++;
    }

    m_T = T_t(dots);
    m_splittedT = {
                    .left = InterpolationPolynom(leftDots),
                    .right = InterpolationPolynom(rightDots),
                    .result = InterpolationPolynom(resultDots)
                  };
}

void GlobalParams::generateS(const Circut &circut)
{
    dots_t dots;

    for(std::size_t i = 0; i < circut.size(); i++) {
        auto currentOperation = circut.m_gates[i].m_type;
        auto sw = m_SWitnesses[i];

        dots.push_back({sw, currentOperation});

        for(auto operation : {Sum, Product}) {
            m_opsFromS[operation].push_back(
                        operation == currentOperation ? dot_t{sw, 1} : dot_t{sw, 0}
            );
        }
    }

    m_S = S_t(dots);
}

void GlobalParams::generateW(const Circut &circut)
{
    /*addr -> свидетели*/
    std::unordered_map<std::size_t, std::shared_ptr<cond_t>> duplicates;
    auto insert = [&duplicates](value_t key, witness_t value)
    {
        auto address = (std::size_t)key.address();
        if (!duplicates[address]) {
            duplicates[address] = std::make_shared<cond_t>();
        }

        duplicates[address]->insert(value);
    };

    auto cw = m_witnesses.cbegin();

    auto fillMap = [&cw, &insert](const values_t &row)
    {
        for(const auto &elment : row) {
            insert(elment, *cw++);
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    for(const auto& gate : circut.m_gates) {
        insert(gate.m_input.a, *cw++);
        insert(gate.m_input.b, *cw++);
        insert(gate.m_output, *cw++);
    }


    auto circleIterator = [](cond_t::iterator begin, cond_t::iterator it, cond_t::iterator end)
    {
        if (it == end) {
            return begin;
        }

        return it;
    };

    dots_t dotsWI;
    dots_t dotsWT;
    dotsWI.reserve(duplicates.size());
    dotsWT.reserve(duplicates.size());

    for(const auto &[_, condition] : duplicates) {
        auto begin = condition->begin();
        auto end = condition->end();

        for(auto it = begin; it != end;) {
             dotsWT.push_back({*it, *it});
            dotsWI.push_back({*it, *circleIterator(begin, ++it, end)});
        }
    }

    m_WT = W_t(dotsWT); //witness -> circut value
    m_WI = W_t(dotsWI); //witness -> next_this_value_wintess

    /*// 1. Собираем все дубликаты (твой исходный блок без изменений)
std::unordered_map<std::size_t, std::shared_ptr<cond_t>> duplicates;
auto insert = [&duplicates](value_t key, witness_t value) {
    auto address = (std::size_t)key.address();
    if (!duplicates[address]) {
        duplicates[address] = std::make_shared<cond_t>();
    }
    duplicates[address]->insert(value);
};

auto cw = m_witnesses.cbegin();
auto fillMap = [&cw, &insert](const values_t &row) {
    for(const auto &element : row) {
        insert(element, *cw++);
    }
};

fillMap(circut.m_inputX);
fillMap(circut.m_inputW);

for(const auto& gate : circut.m_gates) {
    insert(gate.m_input.a, *cw++);
    insert(gate.m_input.b, *cw++);
    insert(gate.m_output, *cw++);
}

// 2. Инициализация ПОЛНЫХ векторов перестановок
// n_total — это общее количество всех witness ячеек в схеме
size_t n_total = m_witnesses.size();
std::vector<value_t> sigma_values(n_total);
std::vector<value_t> identity_values(n_total);

for (size_t i = 0; i begin();
    auto end = condition->end();

    for(auto it = begin; it != end; ) {
        size_t current_idx = *it;
        // Находим следующий индекс в цикле
        size_t next_idx = *circleIterator(begin, ++it, end);

        // Записываем перестановку: ячейка current_idx теперь указывает на next_idx
        sigma_values[current_idx] = (value_t)next_idx;
    }
}

// 4. Формирование точек (dots) для сплайнов на базе ПОЛНОГО диапазона
dots_t dotsWI; // Sigma
dots_t dotsWT; // Identity
dotsWI.reserve(n_total);
dotsWT.reserve(n_total);

for(size_t i = 0; i < n_total; ++i) {
    // В сплайне X - это координата (0, 1, 2...), Y - это значение индекса
    dotsWT.push_back({(value_t)i, identity_values[i]});
    dotsWI.push_back({(value_t)i, sigma_values[i]});
}

// 5. Создание итоговых полиномов-сплайнов
m_WI = W_t(dotsWI); // Полином перестановки (Sigma)
m_WT = WT_t(dotsWT); // Полином идентичности (ID)


1. Подготовка «Карты дорог» (WI и WT)
Представь, что твоя таблица — это город. Каждая ячейка (адрес в памяти) — это дом.

    WT (Identity): Это карта с «пропиской». Она говорит: «В доме №1 живет житель №1, в доме №2 — житель №2». Это просто нумерация.
    WI (Sigma): Это карта «связей». Если данные из дома №1 должны быть такими же, как в доме №5, то на этой карте в доме №1 будет написано «Иди в дом №5».
    Главное правило: Мы берем все номера домов и просто меняем их местами. Набор номеров остается тем же, меняется только их распределение по адресам.

2. Заполнение таблицы жителями (Witness T)
Ты берешь свои реальные вычисления и рассаживаешь их по домам (в сплайн
).

    Если твоя карта (WI) говорит, что дома №1 и №5 связаны, ты обязан поселить туда одинаковых людей (одинаковые числа). Если числа будут разными, «замок» не закроется.

3. Запуск «Бухгалтера» (Аккумулятор Z)
Теперь мы создаем специальный полином-счетчик
. Он идет по таблице от начала до конца:

    В каждом доме он берет жителя и его «прописку» (из WT), а потом делит это на того же жителя с его «картой связей» (из WI).
    Магия: Если житель в доме №1 и №5 один и тот же, то когда бухгалтер дойдет до конца, все числители сократятся со знаменателями.
    Итог: Если в самом конце пути (в последней ячейке) бухгалтер получил ровно 1, значит, все связи в городе соблюдены честно.

4. Создание «Слепка» (Merkle Tree)
Поскольку город огромный, ты не можешь показать его целиком. Ты разбиваешь его на районы (сегменты сплайна).

    Для каждого района ты делаешь «фотографию» (KZG-коммитмент), которая включает всё сразу: жителей, карту и отчет бухгалтера.
    Все эти фото ты складываешь в пирамиду (Дерево Меркла), на вершине которой остается один короткий хеш — Root. Это и есть твое компактное доказательство всей таблицы.

5. Проверка (Verification)
Когда приходит проверяющий:

    Он тыкает в случайный адрес (точку
    ).
    Ты показываешь ему «фотографию» только этого района и доказываешь через Merkle-путь, что это фото из настоящей пирамиды.
    Ты даешь ему значения жителей и бухгалтера в этой точке.
    Финальный тест: Проверяющий подставляет эти числа в формулу бухгалтера. Если
    соотносится с
    так, как велит карта (WI/WT), и всё это подтверждается твоим «фото» (паринг на кривых) — значит, ты не соврал ни в одной ячейке из миллионов.

Итог: Мы превратили огромную проверку связей в проверку одной математической формулы в одной случайной точке.*/
}

}
