//
// Created by vfernandez on 05/03/19.
//

#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("sandbox","[.][sandbox][utils]")
{

    int **a, **b;

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

    WHEN("Straightforward")
    {
        utils::swap(a, b);

        CAPTURE(a[0][0], a[1][0], a[2][0], a[3][0]);
        CAPTURE(a[0][1], a[1][1], a[2][1], a[3][1]);
        CAPTURE(nullptr);
        CAPTURE(b[0][0], b[1][0], b[2][0], b[3][0]);
        CAPTURE(b[0][1], b[1][1], b[2][1], b[3][1]);
        FAIL("Check 1");
    }

    WHEN("Indirect")
    {
        utils::swap(&a[0], &b[0]);

        CAPTURE(a[0][0], a[1][0], a[2][0], a[3][0]);
        CAPTURE(a[0][1], a[1][1], a[2][1], a[3][1]);
        CAPTURE(nullptr);
        CAPTURE(b[0][0], b[1][0], b[2][0], b[3][0]);
        CAPTURE(b[0][1], b[1][1], b[2][1], b[3][1]);
        FAIL("Check 1");
    }


    for (int i = 0; i < 4; i++)
    {
        delete[] a[i];
        delete[] b[i];
    }

    delete[] a;
    delete[] b;
}