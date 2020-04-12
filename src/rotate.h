enum angle { deg90, deg180, deg270 };

struct image {
  int height;
  int width;
  double* data;
};

/**
 * Rotates an image by the given angle.
 *
 * The data array of out needs to be initialized to zero otw. 
 * you get a seg fault.
 */
void rotate(struct image *in, struct image *out, enum angle angle);
