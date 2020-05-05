#include "rotate.h"

#include <stdio.h>

#include <cassert>
#include <cmath>

void rotate(struct image_t &out, const struct image_t &in, int angle) {
    out.size = in.size;

    int m = in.size;
    int n = in.size;

    if (angle == 0) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                out[i * n + j] = in[i * n + j];
            }
        }
        return;
    }

    if (angle == 90) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                // first row has to be last column
                // (too be honest it was trial and error)
                out[j * n + (m - i - 1)] = in[i * n + j];
            }
        }
        return;
    }

    if (angle == 180) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                // first row has to be last row reversed
                out[(m - i - 1) * n + (n - j - 1)] = in[i * n + j];
            }
        }
        return;
    }

    if (angle == 270) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                // first row has to be first column reversed
                out[(m - j - 1) * n + i] = in[i * n + j];
            }
        }
        return;
    }

    throw std::logic_error("given angle is not supported");
}

