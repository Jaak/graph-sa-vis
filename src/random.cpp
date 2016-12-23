#include "random.h"

#include <random>

std::random_device rd;
std::mt19937 gen{rd()};

float randFloat(float low, float high) {
    std::uniform_real_distribution<float> dist(low, high);
    return dist(gen);
}

vec2f sampleCircle() {
    for (;;) {
        const auto x1 = randFloat(-1, 1);
        const auto x2 = randFloat(-1, 1);
        const auto x1Sqr = x1 * x1;
        const auto x2Sqr = x2 * x2;
        const auto sqrSum = x1Sqr + x2Sqr;
        if (sqrSum < 1.0f) {
            const auto invSqrSum = 1.0f / sqrSum;
            return {(x1Sqr - x2Sqr) * invSqrSum, 2*x1*x2*invSqrSum };
        }
    }
}

int randInt(int low, int high) {
    std::uniform_int_distribution<int> dist(low, high);
    return dist(gen);
}