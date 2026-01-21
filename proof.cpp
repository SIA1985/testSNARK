#include "proof.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <nlohmann/json.hpp>
#include <thread>
#include <csignal>

namespace snrk {

PlonkProof::ptr_t PlonkProof::forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, GPK_t &GPK, const witnesses_t &witness, Y_t u)
{
    assert(g.distance() == p.distance());


    auto ptr = ptr_t(new PlonkProof);

    ptr->m_GPK = GPK;

    auto f = g - p;
    ptr->m_comF = f.commit(GPK);

    auto z = ZeroWitnessPolynom(witness).toPartedCanonicPolynom();

    auto q = f.mustDevide(z);
    ptr->m_comQ = q.commit(GPK);

    ptr->m_toProve.x = u;
    ptr->m_toProve.y = q(u);

    return ptr;
}

PlonkProof::ptr_t PlonkProof::forVerifier(commit_t comF, commit_t comQ, dot_t toProve, GPK_t GPK)
{
    auto ptr = ptr_t(new PlonkProof);

    ptr->m_comF = comF;
    ptr->m_comQ = comQ;
    ptr->m_toProve = toProve;
    ptr->m_GPK = GPK;

    return ptr;
}

bool PlonkProof::check(G2 tG2, G2 g2)
{
    G2 zG2, forE1;
    G2::mul(zG2, g2, m_toProve.x);
    G2::sub(forE1, tG2, zG2);

    GT e1;
    mcl::pairing(e1, m_comQ, forE1);

    G1 yG, forE2;
    G1::mul(yG, m_GPK.g, m_toProve.y);
    G1::sub(forE2, m_comF, yG);

    GT e2;
    mcl::pairing(e2, forE2, g2);

    return e1 == e2;
}

json_t PlonkProof::toJson() const
{
    //todo:
    return {};
}

bool PlonkProof::fromJson(const json_t &json)
{
    //todo:
    return false;
}

PolynomSubstitutionProof::ptr_t PolynomSubstitutionProof::forProver(CanonicPolynom &f, dot_t toProve, GPK_t GPK)
{
    //todo: мб потребуется взять за основу PlonkProof
    auto ptr = ptr_t(new PolynomSubstitutionProof);
    auto u = toProve.x;
    auto v = f(u);
    auto xSubU = CanonicPolynom({-u, 1});

    auto q = (f - CanonicPolynom({v, 0})).mustDevide(xSubU);

    ptr->m_comF = f.commit(GPK);
    ptr->m_comQ = q.commit(GPK);
    ptr->m_toProve = toProve;
    ptr->m_GPK = GPK;

    return ptr;
}

PolynomSubstitutionProof::ptr_t PolynomSubstitutionProof::forVerifier(commit_t comF, commit_t comQ, dot_t toProve, GPK_t GPK)
{
    auto ptr = ptr_t(new PolynomSubstitutionProof);

    ptr->m_comF = comF;
    ptr->m_comQ = comQ;
    ptr->m_toProve = toProve;
    ptr->m_GPK = GPK;

    return ptr;
}

bool PolynomSubstitutionProof::check(G2 tG2, G2 g2)
{
    G2 zG2, forE1;
    G2::mul(zG2, g2, m_toProve.x);
    G2::sub(forE1, tG2, zG2);

    GT e1;
    mcl::pairing(e1, m_comQ, forE1);

    G1 yG, forE2;
    G1::mul(yG, m_GPK.g, m_toProve.y);
    G1::sub(forE2, m_comF, yG);

    GT e2;
    mcl::pairing(e2, forE2, g2);

    return e1 == e2;
}

json_t PolynomSubstitutionProof::toJson() const
{
    json_t json;

    //todo:

//    ToJson(json, m_comF);
//    ToJson(json, m_comQ);

    ToJson(json, m_toProve);
//    ToJson(json, m_tG);

    return json;
}

bool PolynomSubstitutionProof::fromJson(const json_t &json)
{
    //todo:

//    FromJson(json, m_comF);
//    FromJson(json, m_comQ);
    FromJson(json, m_toProve);
//    FromJson(json, m_tG);

    return true;
}

/* ПРИМЕР ТРАНСКРИПТА
 * class Transcript {
private:
    std::string state; // Внутреннее состояние "губки"
    const std::string label;

    // Вспомогательная функция для обновления хеша
    void update_state(const void* data, size_t size) {
        const char* p = reinterpret_cast<const char*>(data);
        state.append(p, size);

        // В реальных системах 2026 года здесь используется Keccak256 или Poseidon.
        // Для примера используем встроенный в mcl метод хеширования строки.
        mcl::bn::Fr hash_value;
        hash_value.setHashOf(state);
        state = hash_value.getStr(16); // Обновляем состояние хешем от самого себя
    }

public:
    explicit Transcript(std::string protocol_label) : label(std::move(protocol_label)) {
        state = label;
    }

    // Поглощение точки на кривой G1 (Commitment)
    void append_point(const std::string& label, const mcl::bn::G1& point) {
        state += label;
        // Сериализуем точку в байты (сжатый формат 48 байт)
        char buf[1024];
        size_t n = point.serialize(buf, sizeof(buf));
        update_state(buf, n);
    }

    // Поглощение числа из поля Fr (Evaluation)
    void append_scalar(const std::string& label, const mcl::bn::Fr& scalar) {
        state += label;
        char buf[1024];
        size_t n = scalar.serialize(buf, sizeof(buf));
        update_state(buf, n);
    }

    // Выжимание случайного вызова (Challenge)
    mcl::bn::Fr get_challenge(const std::string& label) {
        state += label;

        mcl::bn::Fr challenge;
        challenge.setHashOf(state);

        // Обновляем состояние, чтобы следующий вызов был другим
        update_state(state.data(), state.size());

        return challenge;
    }
};
*/

ProverProof::ProverProof(const GlobalParams &gp, const values_t &input, value_t output)
{
    auto TParams = gp.PP().TParams;
    auto GPK = gp.GPK();

    //1. Инициализация Транскрипта

    //2. Получаем обязательства splittedT и добавляем в T

    auto Ca = TParams.splittedT.left.commit(GPK);
    auto Cb = TParams.splittedT.right.commit(GPK);
    auto Cc = TParams.splittedT.result.commit(GPK);


    //3. Какой-то полином-аккумулятор, что у меня является W(x) и WT(x)

    //4. Деление на Zero-полином (comT -> T)

    //5. Получение точки раскрытия из Т

    //6. Получение значений всех полиномов выше (a, b, c, w, q) (для последующей проверки в SubstitutionProof,
    //   точнее он будет фигурировать в проверке в "разобранном" виде)

    auto witnesses = gp.witnesses();
    auto tCanonic = TParams.t.toPartedCanonicPolynom();

    correctInputs(tCanonic, input, witnesses, GPK);
//    currentOutput(tCanonic, output, witnesses.size(), GPK);

    #define r(a) std::ref(a)

    std::thread th1(&ProverProof::correctGates, this, r(TParams.splittedT), gp.PP().SParams, gp.SWitnesses(), GPK);
    std::thread th2(&ProverProof::currentVars, this, gp.PP().WParams.wt, r(tCanonic), r(witnesses), GPK);

    #undef r

    th1.join();
    th2.join();
}

bool ProverProof::check(G2 tG2, G2 g2)
{

    return true;
}

json_t ProverProof::toJson() const
{
    //todo:
    json_t json;

    return json;
}

bool ProverProof::fromJson(const json_t &json)
{
    //todo:
    return false;
}

void ProverProof::correctInputs(const PartedCanonicPolynom &tCanonic, values_t inputs, const witnesses_t &ws, GPK_t GPK)
{
    dots_t inputsW;
    inputsW.reserve(inputs.size());
    for(std::size_t i = 0; i < inputs.size(); i++) {
        inputsW.push_back({ws[i], inputs[i]});
    }

    auto funcV = InterpolationPolynom(inputsW).toPartedCanonicPolynom();
    auto funcTCut = tCanonic.cut(funcV.distance());

    auto wStep = *(++ws.begin()) - ws.front();
    witnesses_t witness = genWitnesses(ws.front(), inputs.size(), wStep);
}

void ProverProof::correctGates(const SplittedT_t &t, const GlobalParams::SParams_t SParams, const witnesses_t &ws, GPK_t GPK)
{
    auto left = t.left.toPartedCanonicPolynom();
    auto right = t.right.toPartedCanonicPolynom();
    auto result = t.result.toPartedCanonicPolynom();

    auto funcF = PartedCanonicPolynom(PartedCanonicPolynom::map{});
    for(const auto &[operation, dots] : SParams.opsFromS) {

        auto isOperation = InterpolationPolynom(dots).toPartedCanonicPolynom();

        switch(operation) {
        case Sum: {
            funcF += (left + right) * isOperation;
            break;
        }
        case Product: {
            funcF += (left * right) * isOperation;
            break;
        }
        default:
            assert(false);
        }

    }

//    auto wStep = *(++ws.begin()) - ws.front();
}

void ProverProof::currentVars(const WT_t &wt, PartedCanonicPolynom &tCanonic, const witnesses_t &ws, GPK_t GPK)
{
//    auto wtCanonic = wt.toPartedCanonicPolynom();

//    auto wStep = *(++ws.begin()) - ws.front();
}

void ProverProof::currentOutput(PartedCanonicPolynom &tCanonic, value_t output, std::size_t lastWNum, GPK_t GPK)
{
//    auto outputDot = dot_t{lastWNum, output};

}

}
