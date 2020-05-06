#ifndef ROTATE_H
#define ROTATE_H

#include "types.h"

/**
 * Rotates an image by the given angle.
 *
 * The data array of out needs to be initialized to zero otw.
 * you get a seg fault.
 */
void rotate(struct image_t *out, const struct image_t *in, int angle);

#endif  // ROTATE_H