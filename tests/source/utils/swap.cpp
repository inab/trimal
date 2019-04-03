#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-flp30-c"

#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("swap", "[utils]") {
    WHEN("float numbers") {
        GIVEN("Two numbers") {
            float a = 0, b = 10;
            utils::swap(&a, &b);
            CHECK(a == 10);
            CHECK(b == 0);
        }

        GIVEN("Two non-null pointers") {
            float *a = new float(0), *b = new float(10);
            utils::swap(a, b);
            CHECK(*a == 10);
            CHECK(*b == 0);

            delete a;
            delete b;
        }

        GIVEN("Two pointers, a nullptr") {
            WARN("Fatal error");
            return;
            float *a = nullptr, *b = new float(10);
            utils::swap(a, b);
            CHECK(*a == 10);
            CHECK(b == nullptr);

            delete a;
            delete b;
        }

        GIVEN("Two pointers, b nullptr") {
            WARN("Fatal error");
            return;
            float *a = new float(0), *b = nullptr;
            utils::swap(a, b);
            CHECK(a == nullptr);
            CHECK(*b == 0);

            delete a;
            delete b;
        }
    }

    WHEN("int numbers") {
        GIVEN("Two numbers") {
            int a = 0, b = 10;
            utils::swap(&a, &b);
            CHECK(a == 10);
            CHECK(b == 0);
        }

        GIVEN("Two non-null pointers") {
            int *a = new int(0), *b = new int(10);
            utils::swap(a, b);
            CHECK(*a == 10);
            CHECK(*b == 0);

            delete a;
            delete b;
        }

        GIVEN("Two pointers, a nullptr") {
            WARN("Fatal error");
            return;
            int *a = nullptr, *b = new int(10);
            utils::swap(a, b);
            CHECK(*a == 10);
            CHECK(b == nullptr);

            delete a;
            delete b;
        }

        GIVEN("Two pointers, b nullptr") {
            WARN("Fatal error");
            return;
            int *a = new int(0), *b = nullptr;
            utils::swap(a, b);
            CHECK(a == nullptr);
            CHECK(*b == 0);

            delete a;
            delete b;
        }
    }

    WHEN("int **") {

        int **a, **b;

        GIVEN("Two non-null pointers") {
            a = new int *[4]{
                    new int(1),
                    new int(2),
                    new int(3),
                    new int(4)
            };

            b = new int *[4]{
                    new int(10),
                    new int(20),
                    new int(30),
                    new int(40)
            };

            utils::swap(&a[0], &b[0]);
            THEN("First elements are swapped")
            {
                CHECK(a[0][0] == 10);
                CHECK(b[0][0] == 1);

            }

            THEN("Other elements remain untouched")
            {
                CHECK(a[1][0] == 2);
                CHECK(b[1][0] == 20);

                CHECK(a[2][0] == 3);
                CHECK(b[2][0] == 30);

                CHECK(a[3][0] == 4);
                CHECK(b[3][0] == 40);
            }

            for (int i = 0; i < 4; i++)
            {
                delete a[i];
                delete b[i];
            }

            delete[] a;
            delete[] b;
        }

        GIVEN("Two pointers, a nullptr") {
            WARN("Fatal error");
        }

        GIVEN("Two pointers, b nullptr") {
            WARN("Fatal error");
        }
    }

    WHEN("int ** array") {

        int **a, **b;

        GIVEN("Two non-null pointers") {
            a = new int *[4]{
                    new int[2]{1, 2},
                    new int[2]{2, 2},
                    new int[2]{3, 2},
                    new int[2]{4, 2},
            };

            b = new int *[4]{
                    new int[2]{10, 20},
                    new int[2]{20, 20},
                    new int[2]{30, 20},
                    new int[2]{40, 20},
            };

            utils::swap(a, b);

            THEN("First elements are swapped")
            {
                CHECK(a[0][0] == 10);
                CHECK(b[0][0] == 1);

                CHECK(a[0][1] == 20);
                CHECK(b[0][1] == 2);
            }

            THEN("Other elements remain untouched")
            {
                CHECK(a[1][0] == 2);
                CHECK(b[1][0] == 20);

                CHECK(a[1][1] == 2);
                CHECK(b[1][1] == 20);

                CHECK(a[2][0] == 3);
                CHECK(b[2][0] == 30);

                CHECK(a[2][1] == 2);
                CHECK(b[2][1] == 20);

                CHECK(a[3][0] == 4);
                CHECK(b[3][0] == 40);

                CHECK(a[3][1] == 2);
                CHECK(b[3][1] == 20);
            }

            for (int i = 0; i < 4; i++)
            {
                delete[] a[i];
                delete[] b[i];
            }

            delete[] a;
            delete[] b;
        }

        GIVEN("Two pointers, a nullptr") {
            WARN("Fatal error");
        }

        GIVEN("Two pointers, b nullptr") {
            WARN("Fatal error");
        }
    }
}

#pragma clang diagnostic pop