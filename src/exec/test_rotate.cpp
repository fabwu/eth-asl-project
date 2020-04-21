#include <stdio.h>

#include <cassert>
#include <cmath>

#include "../rotate.h"

static void compare(struct image_t &a, struct image_t &b) {
    bool is_equal = true;

    for (int i = 0; i < a.size; ++i) {
        for (int j = 0; j < a.size; ++j) {
            double val_a = a[i * a.size + j];
            double val_b = b[i * b.size + j];
            if (fabs(val_a - val_b) > 0.001) {
                printf("images different at pos (%d/%d) values %.2f/%.2f\n", i, j,
                       val_a, val_b);
                is_equal = false;
            }
        }
    }

    if (!is_equal) {
        printf("images are NOT equal\n");
        assert(0);
    } else {
        printf("images are equal\n");
    }
}

int main(int argc, char const *argv[]) {
    printf("Test rotate image\n");

    double empty[] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    double in_data[] = {
            1, 2, 3, 4, 5, 6, 7, 8, 9,
    };

    struct image_t in = image_t(in_data, 3);

    struct image_t out(empty, -1);

    struct image_t comp(empty, -1);

    printf("BEGIN rotate image by 0 degree\n");
    // reset out img
    out.size = -1;
    out.data = empty;

    rotate(in, out, 0);

    compare(in, out);

    printf("END rotate image by 0 degree\n");

    printf("BEGIN rotate image by 90 degree\n");
    // reset out img
    out.size = -1;
    out.data = empty;

    rotate(in, out, 90);

    double data90[] = {
            7, 4, 1, 8, 5, 2, 9, 6, 3,
    };
    comp.data = data90;

    compare(comp, out);

    printf("END rotate image by 90 degree\n");

    printf("BEGIN rotate image by 180 degree\n");
    // reset out img
    out.size = -1;
    out.data = empty;

    rotate(in, out, 180);

    double data180[] = {
            9, 8, 7, 6, 5, 4, 3, 2, 1,
    };
    comp.data = data180;

    compare(comp, out);

    printf("END rotate image by 180 degree\n");

    printf("BEGIN rotate image by 270 degree\n");
    // reset out img
    out.size = -1;
    out.data = empty;

    rotate(in, out, 270);

    double data270[] = {
            3, 6, 9, 2, 5, 8, 1, 4, 7,
    };
    comp.data = data270;

    compare(comp, out);

    printf("END rotate image by 270 degree\n");
#if 0
    for(int i = 0; i < out.height; ++i) {
      for(int j = 0; j < out.width; ++j) {
        printf("%f ", out.data[i*out.width + j]);
      }
      printf("\n");
    }
#endif
}
