/**
 * \file Common C types and functions.
 */

#ifndef TYPES_H
#define TYPES_H

#include "queue.h"

#define ALLOCATE(size) (aligned_alloc(32, size));

#define BLOCK_CORD_REL_X(block_id, block_size, image_size) \
    (block_id % (image_size / block_size)) * block_size
#define BLOCK_CORD_REL_Y(block_id, block_size, image_size) \
    ((int)(block_id / (image_size / block_size))) * block_size
#define BLOCK_CORD_X(block_id, block_size, image_size) \
    (block_id % (image_size / block_size))
#define BLOCK_CORD_Y(block_id, block_size, image_size) \
    ((int)(block_id / (image_size / block_size)))

struct image_t {
    double *data;
    int size;
};

struct image_t make_image(int size, int randomize_data);

void free_image_data(const struct image_t *image);

struct block_t {
    int rel_x, rel_y;
    int width, height;
    double sum, sum_squared;
};

struct block_t make_block(int x, int y, int width, int height);

int get_index_in_image(const struct block_t *block, const int y_rel_block,
                       const int x_rel_block, const struct image_t *full_image);

void print_block(const struct block_t *block, const struct image_t *image);

void quad(const struct block_t *block, struct queue *q);

struct transformation_t {
    struct block_t domain_block, range_block;
    double contrast, brightness;
    int angle;
};

typedef struct queue *(*compress_func_type)(const struct image_t *image,
                                            const int error_threshold);

typedef void (*decompress_func_type)(struct image_t *image,
                                     const struct queue *transformations,
                                     const int iterations);

struct func_suite_t {
    compress_func_type compress_func;
    decompress_func_type decompress_func;
};

#endif  // TYPES_H
