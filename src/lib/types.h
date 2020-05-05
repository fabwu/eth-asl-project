/**
 * \file Common C types and functions.
 */

#ifndef TYPES_H
#define TYPES_H

#include <stdio.h>
#include <stdlib.h>
#include "queue.h"

struct image_t {
    double *data;
    int size;
};

struct image_t make_image(int size, bool randomize_data) {
    struct image_t image;
    image.data = (double *)malloc(size * size * sizeof(double));
    if (randomize_data) {
        for (int i = 0; i < size * size; ++i) image.data[i] = rand() % 256;
    }
    return image;
}

void free_image(const image_t *image) { free(image->data); }

struct block_t {
    int rel_x, rel_y;
    int width, height;
};

struct block_t make_block(int x, int y, int width, int height) {
    struct block_t block;
    block.rel_x = x;
    block.rel_y = y;
    block.width = width;
    block.height = height;
    return block;
}

int get_index_in_image(const struct block_t *block, const int y_rel_block,
                       const int x_rel_block,
                       const struct image_t *full_image) {
    return (block->rel_y + y_rel_block) * full_image->size + block->rel_x +
           x_rel_block;
}

void print_block(const struct block_t *block, const struct image_t *image) {
    printf("Block rel_x: %d, rel_y: %d, width: %d, height: %d\n", block->rel_x,
           block->rel_y, block->width, block->height);
    for (int i = 0; i < block->height; ++i) {
        for (int j = 0; j < block->width; ++j) {
            int index = get_index_in_image(block, i, j, image);
            printf("%.1f, ", image->data[index]);
        }
        printf("\n");
    }
}

void quad(const struct block_t *block, struct queue *q) {
    assert(block->width % 2 == 0);
    assert(block->height % 2 == 0);
    assert(block->width >= 2);
    assert(block->height >= 2);

    const int quad_width = block->width / 2;
    const int quad_height = block->height / 2;

    struct block_t *blocks =
        (struct block_t *)malloc(4 * sizeof(struct block_t));
    blocks[0] = make_block(block->rel_x, block->rel_y, quad_width, quad_height);
    blocks[1] = make_block(block->rel_x + quad_width, block->rel_y, quad_width,
                           quad_height);
    blocks[2] = make_block(block->rel_x, block->rel_y + quad_height, quad_width,
                           quad_height);
    blocks[3] = make_block(block->rel_x + quad_width,
                           block->rel_y + quad_height, quad_width, quad_height);

    for (int i = 0; i < 4; ++i) {
        enqueue(q, blocks + i);
    }
}

struct transformation_t {
    block_t domain_block, range_block;
    double contrast, brightness;
    int angle;
};

#endif  // TYPES_H
