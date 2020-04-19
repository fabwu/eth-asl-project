#include "rotate.h"

#include <stdio.h>

#include <cassert>
#include <cmath>

void rotate(struct image_t *in, struct image_t *out, int angle) {
  out->size = in->size;

  int m = in->size;
  int n = in->size;

  if (angle == 0) {
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        out->data[i * n + j] = in->data[i * n + j];
      }
    }
    return;
  }

  if (angle == 90) {
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        // first row has to be last column
        // (too be honest it was trial and error)
        out->data[j * n + (m - i - 1)] = in->data[i * n + j];
      }
    }
    return;
  }

  if (angle == 180) {
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        // first row has to be last row reversed
        out->data[(m - i - 1) * n + (n - j - 1)] = in->data[i * n + j];
      }
    }
    return;
  }

  if (angle == 270) {
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        // first row has to be first column reversed
        out->data[(m - j - 1) * n + i] = in->data[i * n + j];
      }
    }
    return;
  }

  throw std::logic_error("given angle is not supported");
}

