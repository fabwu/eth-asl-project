// uncomment to disable assert()
// #define NDEBUG
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>

#include "lib/types.h"
#include "lib/performance.h"
#include "lib/rotate.h"
#include "lib/queue.h"

#define ALL_ANGLES_LENGTH 4
static const int ALL_ANGLES[] = {0, 90, 180, 270};

static struct block_t *create_squared_blocks(const int image_size,
                                             const int block_size,
                                             const int x_offset,
                                             const int y_offset) {
    assert(image_size % block_size == 0);

    const int num_blocks =
        (image_size / block_size) * (image_size / block_size);
    struct block_t *blocks =
        (struct block_t *)malloc(num_blocks * sizeof(struct block_t));
    int index = 0;
    for (int i = 0; i < image_size; i += block_size) {
        for (int j = 0; j < image_size; j += block_size) {
            blocks[index] = make_block(j + x_offset, i + y_offset, block_size, block_size);
            ++index;
        }
    }

    return blocks;
}

struct image_t *scale_block(const struct image_t *image, const struct block_t *block, int width, int height) {
    assert(block->width >= width);
    assert(block->height >= height);
    assert(block->width == block->height);  // just for simplicity
    assert(block->width % width == 0);      // just for simplicity
    assert(block->height % height == 0);    // just for simplicity
    assert(width == height);                // just for simplicity
    assert(block->width == 2 * width);
    assert(block->height == 2 * height);

    struct image_t *scaled_image = (struct image_t *)malloc(sizeof(struct image_t));
    *scaled_image = make_image(width, 0);

    int idx=0;
    for (int y = 0; y < block->height; y+=2) {
        for (int x = 0; x < block->width; x+=2) {
            double val = 0.0;
            val += image->data[get_index_in_image(block, y, x, image)];
            val += image->data[get_index_in_image(block, y, x + 1, image)];
            val += image->data[get_index_in_image(block, y + 1, x, image)];
            val += image->data[get_index_in_image(block, y + 1, x + 1, image)];

            scaled_image->data[idx] = val / 4.0;
            idx++;
            __record_double_flops(5);
        }
    }

    return scaled_image;
}

/**
 * Returns brightness (o) and contrast (s), such that the RMS between the domain
 * image block and the range block is minmal. Mathematically, this is a least
 * squares problem. The Python implementation we based this naive code on uses
 * such an approach.
 *
 * The Fractal Image Compression book uses a bit a different, more analytical
 * approach. It formulates the error as formula with respect to s and o, and
 * then does partial differentiation with respect to the two variables to get a
 * "direct" solution. Details can be found in the book on page 20 and 21.
 *
 */
double compute_brightness_and_contrast_with_error(
    const struct image_t *image, const struct image_t *domain_block_image,
    const struct block_t *range_block, double *ret_brightness,
    double *ret_contrast) {
    assert(ret_brightness != NULL);
    assert(ret_contrast != NULL);
    assert(domain_block_image->size == range_block->height);
    assert(domain_block_image->size == range_block->width);
    assert(range_block->height == range_block->width);

    const int n = range_block->width;
    const int num_pixels = n * n;

    double sum_domain = 0.0;
    double sum_range = 0.0;
    double sum_range_times_domain = 0.0;
    double sum_domain_squared = 0.0;
    double sum_range_squared = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double di = domain_block_image->data[i * n + j];
            double ri =
                image->data[get_index_in_image(range_block, i, j, image)];
            sum_domain += di;
            sum_range += ri;
            sum_range_times_domain += ri * di;
            sum_domain_squared += di * di;
            sum_range_squared += ri * ri;
        }
    }

    __record_double_flops(n * n * 8);

    double denominator =
        (num_pixels * sum_domain_squared - (sum_domain * sum_domain));
    __record_double_flops(3);
    if (denominator == 0) {
        *ret_contrast = 0.0;
    } else {
        *ret_contrast =
            (num_pixels * sum_range_times_domain - sum_domain * sum_range) /
            denominator;
        __record_double_flops(4);
    }
    *ret_brightness = (sum_range - (*ret_contrast) * sum_domain) / num_pixels;
    __record_double_flops(3);

    // Directly compute the error
    double error =
        (sum_range_squared +
         (*ret_contrast) *
             ((*ret_contrast) * sum_domain_squared -
              2 * sum_range_times_domain + 2 * (*ret_brightness) * sum_domain) +
         (*ret_brightness) * (num_pixels * (*ret_brightness) - 2 * sum_range)) /
        num_pixels;
    __record_double_flops(14);

    if(*ret_contrast>1.0 || *ret_contrast < -1.0) error=DBL_MAX;

    return error;
}

static void rotate_domain_blocks(const struct image_t *domain_block,
                                 struct image_t *result) {
    result[0] = *domain_block;

    for (size_t i = 1; i < ALL_ANGLES_LENGTH; ++i) {
        result[i] = make_image(domain_block->size, 0);
        rotate(&result[i], domain_block, ALL_ANGLES[i]);
    }
}

struct prepared_block_t {
    struct block_t *domain_block;
    struct image_t angles[ALL_ANGLES_LENGTH];
};

void free_prepared_blocks(struct prepared_block_t *prepared_domain_blocks, int size) {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < ALL_ANGLES_LENGTH; ++j) {
            free_image_data(prepared_domain_blocks[i].angles + j);
        }
    }
    free(prepared_domain_blocks);
}

void prepare_domain_blocks(struct prepared_block_t *prepared_domain_blocks,
                           const struct image_t *image,
                           const struct block_t *domain_blocks,
                           const int domain_blocks_length,
                           const int range_block_size) {
    for (size_t i = 0; i < domain_blocks_length; ++i) {
        const struct image_t *scaled_domain_block =
            scale_block(image, domain_blocks + i, range_block_size,
                        range_block_size);

        prepared_domain_blocks[i].domain_block = domain_blocks + i;
        rotate_domain_blocks(scaled_domain_block, prepared_domain_blocks[i].angles);
    }
}

struct queue *compress(const struct image_t *image, const int error_threshold) {
    int current_domain_block_size = image->size;
    const int initial_range_block_size = current_domain_block_size / 2;

    struct block_t *initial_range_blocks =
        create_squared_blocks(image->size, initial_range_block_size, 0, 0);
    const size_t initial_range_blocks_length = 4;
    struct block_t *domain_blocks =
        create_squared_blocks(image->size, current_domain_block_size, 0, 0);
    size_t domain_blocks_length = 1;

    // Need to compress domain_block, such that size(domain_block) ==
    // size(range_block) in order to compare difference! That means that a
    // square of n pixels need to be compressed to one pixel
    struct prepared_block_t *prepared_domain_blocks = (struct prepared_block_t *)malloc(domain_blocks_length * sizeof(struct prepared_block_t));
    prepare_domain_blocks(prepared_domain_blocks, image, domain_blocks, domain_blocks_length, initial_range_block_size);

    // The range blocks we start with
    struct queue remaining_range_blocks = make_queue();
    for (size_t i = 0; i < initial_range_blocks_length; ++i) {
        enqueue(&remaining_range_blocks, initial_range_blocks + i);
    }

    // Learn mappings from domain blocks to range blocks
    // That is, find a domain block for every range block, such that their
    // difference is minimal
    struct queue *transformations = (struct queue *)malloc(sizeof(struct queue));
    *transformations = make_queue();
    int current_range_block_size = initial_range_block_size;
    while (!queue_empty(&remaining_range_blocks)) {
        struct block_t *range_block = (struct block_t *)dequeue(&remaining_range_blocks);

        // Should hold because the queue is FIFO and handles all
        // range blocks of a size before reaching smaller sizes
        assert(range_block->width == current_range_block_size ||
               range_block->width == current_range_block_size / 2);
        assert(range_block->width == range_block->height);
        if (range_block->width < current_range_block_size) {
            free(domain_blocks);
            free_prepared_blocks(prepared_domain_blocks, domain_blocks_length);

            current_domain_block_size = 2*range_block->width;
            domain_blocks =
                create_squared_blocks(image->size, current_domain_block_size, 0, 0);
            domain_blocks_length = (image->size / current_domain_block_size) * (image->size / current_domain_block_size);

            prepared_domain_blocks = (struct prepared_block_t *)malloc(domain_blocks_length * sizeof(struct prepared_block_t));
            prepare_domain_blocks(prepared_domain_blocks, image, domain_blocks, domain_blocks_length, range_block->width);

            current_range_block_size = range_block->width;
        }

        double best_error = DBL_MAX;
        struct transformation_t *best_transformation = (struct transformation_t *)malloc(sizeof(struct transformation_t));

        for (size_t i = 0; i < domain_blocks_length; ++i) {
            struct prepared_block_t *prepared_domain_block =
                prepared_domain_blocks + i;

            assert(prepared_domain_block->domain_block->width == 2 * range_block->width);

            for (size_t j = 0; j < ALL_ANGLES_LENGTH; ++j) {
                const struct image_t *rotated_domain_block =
                    prepared_domain_block->angles + j;
                const int angle = ALL_ANGLES[j];

                double brightness, contrast, error;
                error = compute_brightness_and_contrast_with_error(
                    image, rotated_domain_block, range_block, &brightness,
                    &contrast);

                if (error < best_error) {
                    best_error = error;
                    best_transformation->domain_block =
                        *prepared_domain_block->domain_block;
                    best_transformation->range_block = *range_block;
                    best_transformation->contrast = contrast;
                    best_transformation->brightness = brightness;
                    best_transformation->angle = angle;
                }
            }
        }

        if (best_error > error_threshold) {
            assert(range_block->width >= 2);
            assert(range_block->height >= 2);

            quad(range_block, &remaining_range_blocks);
        } else {
            enqueue(transformations, best_transformation);
        }

        // If range block is not one of the initial range blocks, the range
        // block was created by the quad function and can be freed
        if (!(initial_range_blocks <= range_block &&
              range_block < initial_range_blocks + initial_range_blocks_length)) {
            free(range_block);
        }
    }

    // Free all intermediate values
    free(initial_range_blocks);
    free(domain_blocks);
    free_prepared_blocks(prepared_domain_blocks, domain_blocks_length);
    free_queue(&remaining_range_blocks);

    return transformations;
}

void apply_transformation(struct image_t *image,
                          const struct transformation_t *t) {
    assert(t->domain_block.width == t->domain_block.height);
    assert(t->range_block.width == t->range_block.height);

    struct image_t *scaled_domain_block = scale_block(
        image, &t->domain_block, t->range_block.width, t->range_block.height);
    struct image_t rotated_domain_block = make_image(t->range_block.height, 0);
    rotate(&rotated_domain_block, scaled_domain_block, t->angle);

    for (int i = 0; i < t->range_block.height; ++i) {
        for (int j = 0; j < t->range_block.width; ++j) {
            double value =
                rotated_domain_block.data[i * rotated_domain_block.size + j];
            int idx = get_index_in_image(&t->range_block, i, j, image);
            int new_pixel_value = value * t->contrast + t->brightness;
            if (new_pixel_value < 0) new_pixel_value = 0;
            if (new_pixel_value > 255) new_pixel_value = 255;
            image->data[idx] = new_pixel_value;
            __record_double_flops(2);
        }
    }

    free_image_data(scaled_domain_block);
    free(scaled_domain_block);
    free_image_data(&rotated_domain_block);
}

void decompress(struct image_t *decompressed_image,
                const struct queue *transformations, const int num_iterations) {
    for (int iter = 0; iter < num_iterations; ++iter) {
        const struct queue_node *current = transformations->front;
        while (current != transformations->back) {
            apply_transformation(
                decompressed_image,
                (const struct transformation_t *)current->data);
            current = current->next;
        }
    }
}

struct func_suite_t register_suite(void) {
    struct func_suite_t suite = {.compress_func = &compress,
                                 .decompress_func = &decompress};
    return suite;
}
