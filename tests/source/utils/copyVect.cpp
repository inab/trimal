#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-flp30-c"

#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("copyVect", "[utils]") {
    GIVEN("Two int vectors") {
        std::unique_ptr<int[]> a(nullptr), b(nullptr);

        WHEN("Copying partially") {
            a.reset(new int[4]{1, 2, 3, 4});
            b.reset(new int[4]{0, 0, 4, 5});

            utils::copyVect(a.get(), b.get(), 2);

            CHECK(b[0] == a[0]);
            CHECK(b[1] == a[1]);
            CHECK(b[2] == 4);
            CHECK(b[3] == 5);
        }

        WHEN("Copying totally") {
            a.reset(new int[4]{1, 2, 3, 4});
            b.reset(new int[4]{0, 0, 0, 0});

            utils::copyVect(a.get(), b.get(), 4);
            for (int i = 0; i < 4; i++)
                CHECK(a[i] == b[i]);
        }

    }

    GIVEN("Two float vectors") {
        std::unique_ptr<float[]> a(nullptr), b(nullptr);

        WHEN("Copying partially") {
            a.reset(new float[4]{1, 2, 3, 4});
            b.reset(new float[4]{0, 0, 4, 5});
            utils::copyVect(a.get(), b.get(), 2);

            CHECK(b[0] == a[0]);
            CHECK(b[1] == a[1]);
            CHECK(b[2] == 4);
            CHECK(b[3] == 5);
        }

        WHEN("Copying totally") {
            a.reset(new float[4]{1, 2, 3, 4});
            b.reset(new float[4]{0, 0, 0, 0});
            utils::copyVect(a.get(), b.get(), 4);
            for (int i = 0; i < 4; i++)
                CHECK(a[i] == b[i]);
        }
    }
}

TEST_CASE("[B] copyVect", "[!benchmark][utils]") {
    GIVEN("Two int vectors") {
        std::unique_ptr<int[]> a(nullptr), b(nullptr);

        for (int i = 1; i < 10; i++) {
            int size = pow(2, i);

            a.reset(new int[size]);
            b.reset(new int[size]);
            BENCHMARK("Copying int vector with " + std::to_string(size) + " elements") {
                utils::copyVect(a.get(), b.get(), size);
            }
        }
    }

    GIVEN("Two float vectors") {
        std::unique_ptr<float[]> a(nullptr), b(nullptr);

        for (int i = 1; i < 10; i++) {
            int size = pow(2, i);

            a.reset(new float[size]);
            b.reset(new float[size]);
            BENCHMARK("Copying int vector with " + std::to_string(size) + " elements") {
                utils::copyVect(a.get(), b.get(), size);
            }
        }
    }
}


#pragma clang diagnostic pop