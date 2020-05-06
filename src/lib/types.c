#include "types.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

struct image_t make_image(int size, int randomize_data) {
    struct image_t image;
    image.size = size;
    image.data = (double *)malloc(size * size * sizeof(double));
    if (randomize_data) {
        for (int i = 0; i < size * size; ++i) image.data[i] = rand() % 256;
    }
    return image;
}

void free_image_data(const struct image_t *image) { free(image->data); }

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

    // Has to be 4 different malloc calls, because we want to free them
    // independently later on
    struct block_t *block_0 = (struct block_t *)malloc(sizeof(struct block_t));
    *block_0 = make_block(block->rel_x, block->rel_y, quad_width, quad_height);
    enqueue(q, block_0);

    struct block_t *block_1 = (struct block_t *)malloc(sizeof(struct block_t));
    *block_1 = make_block(block->rel_x + quad_width, block->rel_y, quad_width,
                          quad_height);
    enqueue(q, block_1);

    struct block_t *block_2 = (struct block_t *)malloc(sizeof(struct block_t));
    *block_2 = make_block(block->rel_x, block->rel_y + quad_height, quad_width,
                          quad_height);
    enqueue(q, block_2);

    struct block_t *block_3 = (struct block_t *)malloc(sizeof(struct block_t));
    *block_3 = make_block(block->rel_x + quad_width, block->rel_y + quad_height,
                          quad_width, quad_height);
    enqueue(q, block_3);
}

typedef struct queue *(*compress_func_type)(const struct image_t *image,
                                            const int block_size_domain);

typedef void (*decompress_func_type)(struct image_t *image,
                                     const struct queue *transformations,
                                     const int iterations);
