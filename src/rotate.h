#include "common.h"

enum angle { deg90, deg180, deg270 };

/**
 * Rotates an image by the given angle.
 *
 * The data array of out needs to be initialized to zero otw. 
 * you get a seg fault.
 */
void rotate(struct image_t *in, struct image_t *out, enum angle angle);
