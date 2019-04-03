#include "../../catch.hpp"
#include "utils.h"
#include "random"

#ifndef PotencyOfTwoLevels
#define PotencyOfTwoLevels 12
#endif

TEST_CASE("quicksort", "[utils]") {

    WHEN("Float array") {
        std::unique_ptr<float[]> a(nullptr), b(nullptr), c(nullptr);

        GIVEN("Nullptr, init 0, fin 0") {
            utils::quicksort(a.get(), 0, 0);
        }

        GIVEN("Ordered float array") {
            a.reset(new float[10]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
            b.reset(new float[10]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
            c.reset(new float[10]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

            utils::quicksort(a.get(), 0, 9);
            std::sort(b.get(), b.get() + 10);
            std::stable_sort(c.get(), c.get() + 10);

            for (int i = 0; i < 10; i++) {
                CHECK(a[i] == i);
                CHECK(a[i] == b[i]);
                CHECK(a[i] == c[i]);
            }
        }

        GIVEN("Unordered float array") {
            a.reset(new float[10]{6, 7, 8, 9, 0, 1, 2, 3, 4, 5});
            b.reset(new float[10]{6, 7, 8, 9, 0, 1, 2, 3, 4, 5});
            c.reset(new float[10]{6, 7, 8, 9, 0, 1, 2, 3, 4, 5});

            utils::quicksort(a.get(), 0, 9);
            std::sort(b.get(), b.get() + 10);
            std::stable_sort(c.get(), c.get() + 10);

            for (int i = 0; i < 10; i++) {
                CHECK(a[i] == i);
                CHECK(a[i] == b[i]);
                CHECK(a[i] == c[i]);
            }
        }
    }

    WHEN("Int array") {
        std::unique_ptr<int[]> a(nullptr), b(nullptr), c(nullptr);

        GIVEN("Nullptr, init 0, fin 0") {
            utils::quicksort(a.get(), 0, 0);
        }

        GIVEN("Ordered int array") {
            a.reset(new int[10]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
            b.reset(new int[10]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
            c.reset(new int[10]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

            utils::quicksort(a.get(), 0, 9);
            std::sort(b.get(), b.get() + 10);
            std::stable_sort(c.get(), c.get() + 10);

            for (int i = 0; i < 10; i++) {
                CHECK(a[i] == i);
                CHECK(a[i] == b[i]);
                CHECK(a[i] == c[i]);
            }
        }

        GIVEN("Unordered int array") {
            a.reset(new int[10]{6, 7, 8, 9, 0, 1, 2, 3, 4, 5});
            b.reset(new int[10]{6, 7, 8, 9, 0, 1, 2, 3, 4, 5});
            c.reset(new int[10]{6, 7, 8, 9, 0, 1, 2, 3, 4, 5});

            utils::quicksort(a.get(), 0, 9);
            std::sort(b.get(), b.get() + 10);
            std::stable_sort(c.get(), c.get() + 10);

            for (int i = 0; i < 10; i++) {
                CHECK(a[i] == i);
                CHECK(a[i] == b[i]);
                CHECK(a[i] == c[i]);
            }
        }
    }
}

TEST_CASE("doublequicksort", "[utils]")
{
    int** a;

    GIVEN("Ordered int array") {
        a = new int*[4] {
                new int[4]{0, 1, 2, 3},
                new int[4]{2, 3, 4, 5},
                new int[4]{4, 5, 6, 7},
                new int[4]{6, 7, 8, 9}
        };

        utils::quicksort(a, 0, 3);

        CAPTURE(a[0][0], a[0][1], a[0][2], a[0][3]);
        CHECK(a[0][0] == 0);
        CHECK(a[0][1] == 1);
        CHECK(a[0][2] == 2);
        CHECK(a[0][3] == 3);

        CAPTURE(a[1][0], a[1][1], a[1][2], a[1][3]);
        CHECK(a[1][0] == 2);
        CHECK(a[1][1] == 3);
        CHECK(a[1][2] == 4);
        CHECK(a[1][3] == 5);

        CAPTURE(a[2][0], a[2][1], a[2][2], a[2][3]);
        CHECK(a[2][0] == 4);
        CHECK(a[2][1] == 5);
        CHECK(a[2][2] == 6);
        CHECK(a[2][3] == 7);

        CAPTURE(a[3][0], a[3][1], a[3][2], a[3][3]);
        CHECK(a[3][0] == 6);
        CHECK(a[3][1] == 7);
        CHECK(a[3][2] == 8);
        CHECK(a[3][3] == 9);

        for (int i = 0; i < 4; i++)
            delete [] a[i];
        delete [] a;
    }

    GIVEN("Unordered int double array") {
        a = new int*[4] {
                new int[4]{0, 1, 2, 3},
                new int[4]{4, 5, 6, 7},
                new int[4]{2, 3, 4, 5},
                new int[4]{6, 7, 8, 9}
        };

        utils::quicksort(a, 0, 3);

        CAPTURE(a[0][0], a[0][1], a[0][2], a[0][3]);
        CHECK(a[0][0] == 0);
        CHECK(a[0][1] == 1);
        CHECK(a[0][2] == 2);
        CHECK(a[0][3] == 3);

        CAPTURE(a[1][0], a[1][1], a[1][2], a[1][3]);
        CHECK(a[1][0] == 2);
        CHECK(a[1][1] == 3);
        CHECK(a[1][2] == 4);
        CHECK(a[1][3] == 5);

        CAPTURE(a[2][0], a[2][1], a[2][2], a[2][3]);
        CHECK(a[2][0] == 4);
        CHECK(a[2][1] == 5);
        CHECK(a[2][2] == 6);
        CHECK(a[2][3] == 7);

        CAPTURE(a[3][0], a[3][1], a[3][2], a[3][3]);
        CHECK(a[3][0] == 6);
        CHECK(a[3][1] == 7);
        CHECK(a[3][2] == 8);
        CHECK(a[3][3] == 9);

        for (int i = 0; i < 4; i++)
            delete [] a[i];
        delete [] a;
    }

    GIVEN("Unordered int array") {
        a = new int*[4] {
                new int(0),
                new int(4),
                new int(5),
                new int(2)
            };

        utils::quicksort(a, 0, 3);

        CAPTURE(a[0], a[1], a[2], a[3]);
        CHECK(*a[0] == 0);
        CHECK(*a[1] == 2);
        CHECK(*a[2] == 4);
        CHECK(*a[3] == 5);

        for (int i = 0; i < 4; i++)
            delete a[i];
        delete [] a;
    }
}

TEST_CASE("[B] quicksort", "[!benchmark][utils]") {

    static std::default_random_engine generator;
    static std::uniform_real_distribution<float> distribution;
    static auto getRandomFloat = std::bind ( distribution, generator );

    WHEN("Float array") {
        std::unique_ptr<float[]> a(nullptr), b(nullptr), c(nullptr);

        WHEN("Benchmarking") {
            GIVEN("Random arrays") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new float[size]);
                    b.reset(new float[size]);
                    c.reset(new float[size]);

                    for (int _ = 0; _ < size; _++)
                        a[_] = getRandomFloat();

                    std::copy(a.get(), a.get() + size, b.get());
                    std::copy(a.get(), a.get() + size, c.get());

                    BENCHMARK("std::sort        (float) " + std::to_string(size) + " random elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (float) " + std::to_string(size) + " random elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (float) " + std::to_string(size) + " random elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }

            GIVEN("Ordered arrays") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new float[size]);
                    b.reset(new float[size]);
                    c.reset(new float[size]);

                    for (int x = 0; x < size; x++)
                        a[x] = b[x] = c[x] = x;

                    BENCHMARK("std::sort        (float) " + std::to_string(size) + " ordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (float) " + std::to_string(size) + " ordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (float) " + std::to_string(size) + " ordered elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }

            GIVEN("Unordered arrays") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new float[size]);
                    b.reset(new float[size]);
                    c.reset(new float[size]);

                    for (int x = 0; x < size; x++)
                        a[x] = c[x] = b[x] = size - x;

                    BENCHMARK("std::sort        (float) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (float) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (float) " + std::to_string(size) + " unordered elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }

            GIVEN("Unordered arrays type 2") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new float[size]);
                    b.reset(new float[size]);
                    c.reset(new float[size]);

                    for (int x = 0; x < size; x += 2)
                        a[x] = b[x] = c[x] = size - x;
                    for (int x = 1; x < size; x += 2)
                        a[x] = b[x] = c[x] = x;

                    BENCHMARK("std::sort        (float) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (float) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (float) " + std::to_string(size) + " unordered elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }
        }
    }

    WHEN("Int array") {
        std::unique_ptr<int[]> a(nullptr), b(nullptr), c(nullptr);

        WHEN("Benchmarking") {
            GIVEN("Random arrays") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new int[size]);
                    b.reset(new int[size]);
                    c.reset(new int[size]);

                    for (int _ = 0; _ < size; _++)
                        a[_] = getRandomFloat();

                    std::copy(a.get(), a.get() + size, b.get());
                    std::copy(a.get(), a.get() + size, c.get());

                    BENCHMARK("std::sort        (int) " + std::to_string(size) + " random elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (int) " + std::to_string(size) + " random elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (int) " + std::to_string(size) + " random elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }

            GIVEN("Ordered arrays") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new int[size]);
                    b.reset(new int[size]);
                    c.reset(new int[size]);

                    for (int x = 0; x < size; x++)
                        a[x] = b[x] = c[x] = x;

                    BENCHMARK("std::sort        (int) " + std::to_string(size) + " ordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (int) " + std::to_string(size) + " ordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (int) " + std::to_string(size) + " ordered elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }

            GIVEN("Unordered arrays") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new int[size]);
                    b.reset(new int[size]);
                    c.reset(new int[size]);

                    for (int x = 0; x < size; x++)
                        a[x] = c[x] = b[x] = size - x;

                    BENCHMARK("std::sort        (int) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (int) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (int) " + std::to_string(size) + " unordered elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }

            GIVEN("Unordered arrays type 2") {
                for (int i = 1; i < PotencyOfTwoLevels; i++) {
                    int size = (int) pow(2, i);

                    a.reset(new int[size]);
                    b.reset(new int[size]);
                    c.reset(new int[size]);

                    for (int x = 0; x < size; x += 2)
                        a[x] = b[x] = c[x] = size - x;
                    for (int x = 1; x < size; x += 2)
                        a[x] = b[x] = c[x] = x;

                    BENCHMARK("std::sort        (int) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("std::stable_sort (int) " + std::to_string(size) + " unordered elements")std::sort(b.get(), b.get() + size - 1);

                    BENCHMARK("utils::quicksort (int) " + std::to_string(size) + " unordered elements")utils::quicksort(a.get(), 0, size - 1);

                    std::cerr << "\n";
                }
            }
        }
    }
}
#undef PotencyOfTwoLevels