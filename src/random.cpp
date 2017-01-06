#include "random.h"

#include <random>
#include <iostream>

namespace /* anonymous */ {

uint64_t s[2];

struct SeedInitializer {
    SeedInitializer() {
        std::random_device rd;
        s[0] = rd();
        s[1] = rd();
    }
};

SeedInitializer __useless;

inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

inline uint64_t next(void) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;

	s1 ^= s0;
	s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
	s[1] = rotl(s1, 36); // c

	return result;
}

inline double u64ToDouble(uint64_t x) {
    const uint64_t u = UINT64_C(0x3FF) << 52 | x >> 12;
    double result;
    memcpy(&result, &u, sizeof(uint64_t));
    return result - 1.0;
}

inline double randomDouble01() {
    return u64ToDouble(next());
}

} // namespace anonymous

float randFloat(float low, float high) {
    return randomDouble01()*(high - low) + low;
}

vec2f sampleCircle() {
    for (;;) {
        const auto x1 = randFloat(-1.0f, 1.0f);
        const auto x2 = randFloat(-1.0f, 1.0f);
        const auto x1Sqr = x1 * x1;
        const auto x2Sqr = x2 * x2;
        const auto sqrSum = x1Sqr + x2Sqr;
        if (sqrSum < 1.0f) {
            const auto invSqrSum = 1.0f / sqrSum;
            return {(x1Sqr - x2Sqr) * invSqrSum, 2.0f*x1*x2*invSqrSum };
        }
    }
}

int randInt(int low, int high) {
    const uint64_t maxVal = high - low;
    uint64_t mask = maxVal;
    for (int i = 1; i < 64; i *= 2)
        mask = mask | (mask >> i);

    std::cout << maxVal << "    " << mask << std::endl;

    for (;;) {
        const uint64_t val = next() & mask;
        if (val <= maxVal) {
            return val + low;
        }
    }
}