#include "annealing.h"
#include "options.h"
#include "random.h"

#include <iostream>
#include <future>
#include <thread>

using namespace std::chrono;

namespace /* anonymous */ {

inline float acceptanceP(float e0, float e1, float T) {
    if (e1 < e0)
        return 1;
    else
        return std::exp(-(e1 - e0)/T);
}

std::vector<event_type<float>> trackEvents(float& e, std::atomic_bool& running) {
    const size_t samplesPerSecond = 4;
    const auto delayInMilliseconds = 100 / samplesPerSecond;

    std::vector<event_type<float>> energies;
    while (running) {
        std::this_thread::sleep_for(milliseconds(delayInMilliseconds));
        const auto energy = e;
        const auto timeNow = steady_clock::now();
        energies.emplace_back(timeNow, energy);
    }

    return energies;
}

} // namespace anonymous {

anneal_result anneal(GraphInstance s) {
    auto e = s.totalEnergy;
    auto T = InitialTemperature;
    const auto n = s.gr->num_vertices;
    float R = 2.0f * std::min(BoundingBoxHeight, BoundingBoxWidth) / std::sqrt(n);

    std::cout << "InitialEnergy: " << e << std::endl;

    std::atomic_bool running {true};

    std::thread printerThread {
        [&e, &running]() {
            while (running) {
                std::this_thread::sleep_for(milliseconds(50));
                std::cout << "Energy:        " << e << '\r';
                std::cout.flush();
            }
        }
    };

    auto energiesHandle = std::async(std::launch::async, trackEvents, std::ref(e), std::ref(running));

    const auto limit = n * StepFactor;
    size_t stageCount = 0;
    const auto startTime = steady_clock::now();

    for (;;) {
        uint32_t n = 0;
        float mean = 0.0f;
        float M2 = 0.0f;
        for (uint step = 0; step < limit; ++ step) {
            const auto s1 = neighbour(s, R);
            const auto e1 = s1.totalEnergy;
            if (acceptanceP(e, e1, T) >= randFloat(0, 1)) {
                s = s1;
                e = e1;
            }

            n += 1;
            const auto delta = e - mean;
            mean += delta / n;
            const auto delta2 = e - mean;
            M2 += delta*delta2;
        }

        if (M2 < (n - 1) * 1e-5f) {
            break;
        }

        T = T * GammaFactor;
        R = std::max(R - 0.5f, MinimumMoveRadius);
        ++ stageCount;
    }

    running = false;
    const auto endTime = steady_clock::now();
    const duration<double> timeDelta = endTime - startTime;

    printerThread.join();
    const auto energies = energiesHandle.get();

    std::cout << "StageCount:    " << stageCount << "        " << std::endl;
    std::cout << "HistorySize:   " << energies.size() << std::endl;
    std::cout << "FinalEnergy:   " << e << std::endl;
    std::cout << "Time:          " << timeDelta.count() << "s" << std::endl;
    return {s, energies};
}

anneal_result anneal(const Graph& gr) {
    return anneal(GraphInstance::randomised(gr));
}