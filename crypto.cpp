#include "crypto.h"

#include "types.h"
#define MCL_USE_GMP 1
#include <mcl/bls12_381.hpp>

namespace snrk {

void x()
{
    using namespace mcl::bn;

    // 1. Инициализация параметров кривой BLS12-381
    initPairing(mcl::BLS12_381);

    // 2. Генераторы групп G1 и G2
    G1 P;
    G2 Q;
    mapToG1(P, 1); // Базовая точка G1
    mapToG2(Q, 1); // Базовая точка G2

    // 3. Создаем два скаляра a и b
    Fr a, b;
    a.setByCSPRNG(); // Случайное число
    b.setByCSPRNG(); // Случайное число

    // 4. Вычисляем точки: A = a*P, B = b*Q
    G1 A;
    G2 B;
    G1::mul(A, P, a); // A = aP
    G2::mul(B, Q, b); // B = bQ

    // 5. Вычисляем левую часть: e(aP, bQ)
    GT e1;
    pairing(e1, A, B);

    // 6. Вычисляем правую часть: e(P, Q)^(a*b)
    GT e2;
    pairing(e2, P, Q);         // Сначала само спаривание
    GT::pow(e2, e2, a * b);    // Затем возведение в степень (a*b)

    // 7. Проверка равенства
    if (e1 == e2) {
//        std::cout << "Успех: Билинейное свойство подтверждено!" << std::endl;
//        std::cout << "e(aP, bQ) == e(P, Q)^(ab)" << std::endl;
    } else {
//        std::cout << "Ошибка: спаривания не совпали." << std::endl;
    }

    // 8. Типичный пример для проверки подписи BLS:
    // e(aP, Q) == e(P, aQ)
    G1 aP;
    G2 aQ;
    G1::mul(aP, P, a);
    G2::mul(aQ, Q, a);

    GT e3, e4;
    pairing(e3, aP, Q);
    pairing(e4, P, aQ);

    if (e3 == e4) {
//        std::cout << "Проверка e(aP, Q) == e(P, aQ) пройдена." << std::endl;
    }
}

}
