#include <stdio.h>
#include <cassert>
#include <cmath>
#include "../rotate.h"

static void compare(struct image *a, struct image *b) {
  assert(a->width == b->width);
  assert(a->height == b->height);
  bool is_equal = true;

  for(int i = 0; i < a->height; ++i) {
    for(int j = 0; j < a->width; ++j) {
      double val_a = a->data[i*a->width + j];
      double val_b = b->data[i*b->width + j];
      if(fabs(val_a - val_b) > 0.001) {
        printf("images different at pos (%d/%d) values %.2f/%.2f\n", i,j,val_a,val_b);
        is_equal = false;
      }
    }
  }

  if(!is_equal) {
    printf("images are NOT equal\n");
    assert(0);
  } else {
    printf("images are equal\n");
  }
}

void rotate(struct image *in, struct image *out, enum angle angle) {
  assert(in->height == in->width); // only square images supported
  out->height = in->height;
  out->width = in->width;

  int m = in->height;
  int n = in->width;

  if(angle == deg90) {
    for(int i = 0; i < m; ++i) {
      for(int j = 0; j < n; ++j) {
        // first row has to be last column 
        // (too be honest it was trial and error)
        out->data[j*n + (m - i - 1)] = in->data[i*n + j]; 
      }
    } 
  }
  
  if(angle == deg180) {
    for(int i = 0; i < m; ++i) {
      for(int j = 0; j < n; ++j) {
        // first row has to be last row reversed 
        out->data[(m - i - 1)*n + (n - j - 1)] = in->data[i*n + j]; 
      }
    } 
  }
  
  if(angle == deg270) {
    for(int i = 0; i < m; ++i) {
      for(int j = 0; j < n; ++j) {
        // first row has to be first column reversed 
        out->data[(m - j - 1)*n + i] = in->data[i*n + j]; 
      }
    } 
  }
}

int main(int argc, char const *argv[]) {
  printf("Test rotate image\n");

  double empty[] = {
    0,0,0,
    0,0,0,
    0,0,0,
  };

  double in_data[] = {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9,
  };

  struct image in = {
    .height = 3,
    .width = 3,
    .data = in_data
  }; 

  struct image out;

  struct image comp = {
    .height = in.height,
    .width = in.width,
  };

  printf("BEGIN rotate image by 90 degree\n");
  // reset out img
  out.height = -1;
  out.width = -1;
  out.data = empty;

  rotate(&in, &out, deg90);

  double data90[] = {
    7, 4, 1,
    8, 5, 2,
    9, 6, 3,
  };
  comp.data = data90;

  compare(&comp, &out);

  printf("END rotate image by 90 degree\n");

  printf("BEGIN rotate image by 180 degree\n");
  // reset out img
  out.height = -1;
  out.width = -1;
  out.data = empty;

  rotate(&in, &out, deg180);

  double data180[] = {
    9, 8, 7,
    6, 5, 4,
    3, 2, 1,
  };
  comp.data = data180;

  compare(&comp, &out);

  printf("END rotate image by 180 degree\n");

  printf("BEGIN rotate image by 270 degree\n");
  // reset out img
  out.height = -1;
  out.width = -1;
  out.data = empty;

  rotate(&in, &out, deg270);

  double data270[] = {
    3, 6, 9,
    2, 5, 8,
    1, 4, 7,
  };
  comp.data = data270;

  compare(&comp, &out);

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
