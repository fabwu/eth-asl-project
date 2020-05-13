/**
 * IMPORTANT CHANGES MADE ("algorithmic")
 * - MAX_QUADTREE_DEPTH
 * - MIN_QUADTREE_DEPTH
 * - disallow domain block sizes of 1 pixel
 *
 *
 */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "lib/performance.h"
#include "lib/queue.h"
#include "lib/rotate.h"
#include "lib/types.h"

/**
 * The inital domain block size is image_size / 2^MIN_QUADTREE_DEPTH
 * Should be >=1, because we can then can assume that the amount of domain
 * blocks is divisible by 2.
 */
#define MIN_QUADTREE_DEPTH 1
#define MAX_QUADTREE_DEPTH 7

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
            blocks[index] =
                make_block(j + x_offset, i + y_offset, block_size, block_size);
            ++index;
        }
    }

    return blocks;
}

struct image_t scale_block(const struct image_t *image,
                           const struct block_t *block, int width, int height) {
    assert(block->width >= width);
    assert(block->height >= height);
    assert(block->width == block->height);  // just for simplicity
    assert(block->width % width == 0);      // just for simplicity
    assert(block->height % height == 0);    // just for simplicity
    assert(width == height);                // just for simplicity
    assert(block->width == 2 * width);
    assert(block->height == 2 * height);

    struct image_t scaled_image = make_image(width, 0);

    size_t scaled_image_idx = 0;
    size_t original_image_idx = block->rel_y * image->size + block->rel_x;
    for (int y = 0; y < block->height; y += 2) {
        for (int x = 0; x < block->width; x += 2) {
            double val = 0.0;
            val += image->data[original_image_idx];
            val += image->data[original_image_idx + 1];
            val += image->data[original_image_idx + image->size];
            val += image->data[original_image_idx + image->size + 1];
            scaled_image.data[scaled_image_idx] = val / 4.0;
            scaled_image_idx++;
            original_image_idx += 2;
            __record_double_flops(5);
        }
        original_image_idx += image->size;
    }

    return scaled_image;
}

/**
 * Returns brightness (o) and contrast (s), such that the RMS between the domain
 * image block and the range block is minimal. Mathematically, this is a least
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
    const int range_block_size, double *ret_brightness, double *ret_contrast,
    const double sum_range, const double sum_range_squared,
    const double sum_domain, const double sum_domain_squared,
    const double sum_range_times_domain) {
    assert(ret_brightness != NULL);
    assert(ret_contrast != NULL);

    const int num_pixels = range_block_size * range_block_size;

    double contrast;
    double brightness;

    double denominator =
        (num_pixels * sum_domain_squared - (sum_domain * sum_domain));
    __record_double_flops(3);
    if (denominator == 0) {
        contrast = 0.0;
    } else {
        contrast =
            (num_pixels * sum_range_times_domain - sum_domain * sum_range) /
            denominator;
        __record_double_flops(4);
    }
    brightness = (sum_range - contrast * sum_domain) / num_pixels;
    __record_double_flops(3);

    // Directly compute the error
    double error =
        (sum_range_squared +
         contrast * (contrast * sum_domain_squared -
                     2 * sum_range_times_domain + 2 * brightness * sum_domain) +
         brightness * (num_pixels * brightness - 2 * sum_range)) /
        num_pixels;
    __record_double_flops(14);

    if (contrast > 1.0 || contrast < -1.0) error = DBL_MAX;

    *ret_contrast = contrast;
    *ret_brightness = brightness;
    return error;
}

void free_prepared_blocks(struct image_t *prepared_domain_blocks,
                          const int length) {
    for (size_t i = 0; i < length; ++i) {
        free((prepared_domain_blocks + i)->data);
    }
    free(prepared_domain_blocks);
}

void prepare_domain_blocks_norotation(struct image_t *prepared_domain_blocks,
                                      double *sums, double *sums_squared,
                                      const struct image_t *image,
                                      const struct block_t *domain_blocks,
                                      const int domain_blocks_length,
                                      const int range_block_size) {
    for (size_t i = 0; i < domain_blocks_length; ++i) {
        const struct image_t scaled_domain_block = scale_block(
            image, domain_blocks + i, range_block_size, range_block_size);

        double sum = 0;
        double sum_squared = 0;
        for (size_t j = 0;
             j < scaled_domain_block.size * scaled_domain_block.size; j++) {
            sum += scaled_domain_block.data[j];
            sum_squared +=
                scaled_domain_block.data[j] * scaled_domain_block.data[j];
        }
        __record_double_flops(scaled_domain_block.size *
                              scaled_domain_block.size * 3);

        sums[i] = sum;
        sums_squared[i] = sum_squared;

        prepared_domain_blocks[i] = scaled_domain_block;
    }
}

typedef double (*rtd_func_type)(const struct image_t *image,
                                const struct image_t *rotated_domain_block,
                                const struct block_t *range_block);

typedef size_t (*rotation_index_transformer_type)(const size_t i,
                                                  const size_t j, const int n,
                                                  const int m);
size_t index_rot0(const size_t i, const size_t j, const int n, const int m) {
    return i * n + j;
}
size_t index_rot90(const size_t i, const size_t j, const int n, const int m) {
    return (m - j - 1) * n + i;
}
size_t index_rot180(const size_t i, const size_t j, const int n, const int m) {
    return (m - i - 1) * n + (n - j - 1);
}
size_t index_rot270(const size_t i, const size_t j, const int n, const int m) {
    return j * n + (m - i - 1);
}

//double rtd_generic_with_rot(const struct image_t *image,
//                            const struct image_t *domain_block,
//                            const struct block_t *range_block,
//                            const rotation_index_transformer_type rot) {
//    const int dbs = domain_block->size;
//
//    double rtd_sum1 = 0;
//    double rtd_sum2 = 0;
//
//    size_t idx_rb1 = range_block->rel_y * image->size + range_block->rel_x;
//
//    for (size_t i = 0; i < dbs; i++) {
//        for (size_t j = 0; j < dbs; j += 2) {
//            const size_t idx_db1 = rot(i, j, dbs, dbs);
//
//            double ri1 = image->data[idx_rb1];
//            double di1 = domain_block->data[idx_db1];
//            rtd_sum1 = fma(ri1, di1, rtd_sum1);
//            idx_rb1++;
//
//            const size_t idx_db2 = rot(i, j + 1, dbs, dbs);
//            double ri2 = image->data[idx_rb1];
//            double di2 = domain_block->data[idx_db2];
//            rtd_sum2 = fma(ri2, di2, rtd_sum2);
//            idx_rb1++;
//        }
//        idx_rb1 += image->size - dbs;
//    }
//
//    __record_double_flops(dbs * dbs * 2);
//
//    return rtd_sum1 + rtd_sum2;
//}

double rtd_generic_with_rot(const struct image_t *image,
                            const struct image_t *domain_block,
                            const struct block_t *range_block,
                            const rotation_index_transformer_type rot) {
    const int dbs = domain_block->size;

    double rtd_sum1 = 0;
    double rtd_sum2 = 0;

    size_t idx_rb1 = range_block->rel_y * image->size + range_block->rel_x;

    for (size_t i = 0; i < dbs; i++) {
        for (size_t j = 0; j < dbs; j += 2) {
            const size_t idx_db1 = rot(i, j, dbs, dbs);

            double ri1 = image->data[idx_rb1];
            double di1 = domain_block->data[idx_db1];
            rtd_sum1 = fma(ri1, di1, rtd_sum1);
            idx_rb1++;

            const size_t idx_db2 = rot(i, j + 1, dbs, dbs);
            double ri2 = image->data[idx_rb1];
            double di2 = domain_block->data[idx_db2];
            rtd_sum2 = fma(ri2, di2, rtd_sum2);
            idx_rb1++;
        }
        idx_rb1 += image->size - dbs;
    }

    __record_double_flops(dbs * dbs * 2);

    return rtd_sum1 + rtd_sum2;
}


void precompute_sums(double *sums, double *sums_squared,
                     const struct block_t *blocks, const int blocks_length,
                     const struct image_t *image) {
    for (int i = 0; i < blocks_length; ++i) {
        struct block_t *rb = blocks + i;
        double sum = 0;
        double sum_squared = 0;
        size_t ind = rb->rel_y * image->size + rb->rel_x;
        for (size_t l = 0; l < rb->width; l++) {
            for (size_t m = 0; m < rb->height; m++) {
                sum += image->data[ind];
                sum_squared += image->data[ind] * image->data[ind];
                ind++;
            }
            ind += image->size - rb->height;
        }

        sums[i] = sum;
        sums_squared[i] = sum_squared;
        __record_double_flops(rb->width * rb->height * 3);
    }
}

void precompute_rtd_with_rotation0(
    double *rtd_sums, const struct image_t *image, struct block_t *range_blocks,
    const int range_blocks_length,
    const struct image_t *downsampled_domain_blocks,
    const int domain_blocks_length) {
    assert(domain_blocks_length % 2 == 0);

    for (size_t idx_db = 0; idx_db < domain_blocks_length; idx_db++) {
        const struct image_t *domain_block = downsampled_domain_blocks + idx_db;
        const int dbs = domain_block->size;
        for (int idx_rb = 0; idx_rb < range_blocks_length; ++idx_rb) {
            struct block_t *range_block = range_blocks + idx_rb;

            // Compute the sum
            double rtd_sum = 0;
            double rtd_sum2 = 0;

            size_t idx1_rb =
                range_block->rel_y * image->size + range_block->rel_x;
            size_t idx1_db = 0;
            for (size_t i = 0; i < dbs; i++) {
                for (size_t j = 0; j < dbs; j++) {
                    double ri = image->data[idx1_rb];
                    double di = domain_block->data[idx1_db];

                    rtd_sum += ri * di;
                    idx1_rb++;
                    idx1_db++;
                }
                idx1_rb += image->size - dbs;
            }
            __record_double_flops(dbs * dbs * 2);

            rtd_sums[idx_rb * (domain_blocks_length) + idx_db] =
                rtd_sum + rtd_sum2;
        }
    }
}

void precompute_rtd_with_rotation(
    double *rtd_sum, const struct image_t *image, struct block_t *range_blocks,
    const int range_blocks_length,
    const struct image_t *downsampled_domain_blocks,
    const int domain_blocks_length, rotation_index_transformer_type rot) {
    for (size_t idx_db = 0; idx_db < domain_blocks_length; ++idx_db) {
        const struct image_t *domain_block = downsampled_domain_blocks + idx_db;

        for (int idx_rb = 0; idx_rb < range_blocks_length; ++idx_rb) {
            struct block_t *range_block = range_blocks + idx_rb;

            double rtd =
                rtd_generic_with_rot(image, domain_block, range_block, rot);

            rtd_sum[idx_rb * (domain_blocks_length) + idx_db] = rtd;
        }
    }
}

void quad2(const struct block_t *block, struct block_t *list) {
    assert(block->width % 2 == 0);
    assert(block->height % 2 == 0);
    assert(block->width >= 2);
    assert(block->height >= 2);

    const int quad_width = block->width / 2;
    const int quad_height = block->height / 2;

    *list = make_block(block->rel_x, block->rel_y, quad_width, quad_height);
    *(list + 1) = make_block(block->rel_x + quad_width, block->rel_y,
                             quad_width, quad_height);
    *(list + 2) = make_block(block->rel_x, block->rel_y + quad_height,
                             quad_width, quad_height);
    *(list + 3) =
        make_block(block->rel_x + quad_width, block->rel_y + quad_height,
                   quad_width, quad_height);
}

struct queue *compress(const struct image_t *image, const int error_threshold) {
    const int initial_domain_block_size =
        image->size / (int)pow(2.0, (double)MIN_QUADTREE_DEPTH);
    const int initial_range_block_size = initial_domain_block_size / 2;

    struct block_t *initial_range_blocks =
        create_squared_blocks(image->size, initial_range_block_size, 0, 0);
    const size_t initial_range_blocks_length =
        (image->size / initial_range_block_size) *
        (image->size / initial_range_block_size);

    // Forward declarations
    struct block_t *domain_blocks;
    size_t domain_blocks_length = -1;
    int domain_block_size_current_iteration = -1;
    struct image_t *downsampled_domain_blocks;
    double *rtd_sum_rot0;
    double *rtd_sum_rot90;
    double *rtd_sum_rot180;
    double *rtd_sum_rot270;
    double *range_block_sums;
    double *range_block_sums_squared;
    double *domain_block_sums;
    double *domain_block_sums_squared;

    // Queue for saving transformations
    struct queue *transformations =
        (struct queue *)malloc(sizeof(struct queue));
    *transformations = make_queue();

    int current_quadtree_depth = MIN_QUADTREE_DEPTH - 1;
    struct block_t *range_blocks_curr_iteration = NULL;
    int range_blocks_length_current_iteration = -1;
    int range_blocks_size_current_iteration = -1;

    int range_blocks_length_next_iteration = initial_range_blocks_length;
    int range_blocks_size_next_iteration = initial_range_block_size;
    struct block_t *range_blocks_next_iteration = initial_range_blocks;
    bool has_remaining_range_blocks = true;

    while (has_remaining_range_blocks) {
        has_remaining_range_blocks = false;
        current_quadtree_depth++;

        // Init current iteration over range blocks
        // Prepare for next run
        {
            range_blocks_length_current_iteration =
                range_blocks_length_next_iteration;
            range_blocks_curr_iteration = range_blocks_next_iteration;
            range_blocks_size_current_iteration =
                range_blocks_size_next_iteration;

            range_blocks_length_next_iteration = 0;
            range_blocks_size_next_iteration /= 2;
            range_blocks_next_iteration =
                malloc(4 * range_blocks_length_current_iteration *
                       sizeof(struct block_t));

            __record_double_flops(1);
        }

        // Prepare the domain blocks
        {
            domain_block_size_current_iteration =
                2 * range_blocks_size_current_iteration;
            assert(domain_block_size_current_iteration >= 2);
            domain_blocks = create_squared_blocks(
                image->size, domain_block_size_current_iteration, 0, 0);
            domain_blocks_length =
                (image->size / domain_block_size_current_iteration) *
                (image->size / domain_block_size_current_iteration);

            downsampled_domain_blocks =
                malloc(domain_blocks_length * sizeof(struct image_t));

            domain_block_sums = malloc(sizeof(double) * domain_blocks_length);
            domain_block_sums_squared =
                malloc(sizeof(double) * domain_blocks_length);

            prepare_domain_blocks_norotation(
                downsampled_domain_blocks, domain_block_sums,
                domain_block_sums_squared, image, domain_blocks,
                domain_blocks_length, range_blocks_size_current_iteration);

            __record_double_flops(4);
        }

        // Precomputations on domain / range blocks
        {
            range_block_sums =
                malloc(sizeof(double) * range_blocks_length_current_iteration);
            assert(range_block_sums != NULL);
            range_block_sums_squared =
                malloc(sizeof(double) * range_blocks_length_current_iteration);
            assert(range_block_sums_squared != NULL);
            precompute_sums(range_block_sums, range_block_sums_squared,
                            range_blocks_curr_iteration,
                            range_blocks_length_current_iteration, image);

            rtd_sum_rot0 = malloc(sizeof(double) * domain_blocks_length *
                                  range_blocks_length_current_iteration);
            assert(rtd_sum_rot0 != NULL);
            precompute_rtd_with_rotation(
                rtd_sum_rot0, image, range_blocks_curr_iteration,
                range_blocks_length_current_iteration,
                downsampled_domain_blocks, domain_blocks_length, &index_rot0);

            rtd_sum_rot90 = malloc(sizeof(double) * domain_blocks_length *
                                   range_blocks_length_current_iteration);
            assert(rtd_sum_rot90 != NULL);
            precompute_rtd_with_rotation(
                rtd_sum_rot90, image, range_blocks_curr_iteration,
                range_blocks_length_current_iteration,
                downsampled_domain_blocks, domain_blocks_length, &index_rot90);

            rtd_sum_rot180 = malloc(sizeof(double) * domain_blocks_length *
                                    range_blocks_length_current_iteration);
            assert(rtd_sum_rot180 != NULL);
            precompute_rtd_with_rotation(
                rtd_sum_rot180, image, range_blocks_curr_iteration,
                range_blocks_length_current_iteration,
                downsampled_domain_blocks, domain_blocks_length, &index_rot180);

            rtd_sum_rot270 = malloc(sizeof(double) * domain_blocks_length *
                                    range_blocks_length_current_iteration);
            assert(rtd_sum_rot270 != NULL);
            precompute_rtd_with_rotation(
                rtd_sum_rot270, image, range_blocks_curr_iteration,
                range_blocks_length_current_iteration,
                downsampled_domain_blocks, domain_blocks_length, &index_rot270);
        }

        // Process each range block
        for (size_t idx_rb = 0; idx_rb < range_blocks_length_current_iteration;
             ++idx_rb) {
            struct block_t *range_block = range_blocks_curr_iteration + idx_rb;

            double best_error = DBL_MAX;

            int best_range_block_rel_x;
            int best_range_block_rel_y;
            int best_domain_block_idx;
            double best_contrast;
            double best_brightness;
            int best_angle;

            const double range_sum = range_block_sums[idx_rb];
            const double range_sum_squared = range_block_sums_squared[idx_rb];

            for (size_t idx_db = 0; idx_db < domain_blocks_length; ++idx_db) {
                assert(domain_blocks[idx_db].width == 2 * range_block->width);

                const double domain_sum = domain_block_sums[idx_db];
                const double domain_sum_squared =
                    domain_block_sums_squared[idx_db];

                // ROTATION 0
                double brightness_rot0, contrast_rot0, error_rot0;
                error_rot0 = compute_brightness_and_contrast_with_error(
                    range_blocks_size_current_iteration, &brightness_rot0,
                    &contrast_rot0, range_sum, range_sum_squared, domain_sum,
                    domain_sum_squared,
                    rtd_sum_rot0[idx_rb * domain_blocks_length + idx_db]);
                if (error_rot0 < best_error) {
                    best_error = error_rot0;
                    best_domain_block_idx = idx_db;
                    best_range_block_rel_x = range_block->rel_x;
                    best_range_block_rel_y = range_block->rel_y;
                    best_contrast = contrast_rot0;
                    best_brightness = brightness_rot0;
                    best_angle = 0;
                }

                // ROTATION 90
                double brightness_rot90, contrast_rot90, error_rot90;
                error_rot90 = compute_brightness_and_contrast_with_error(
                    range_blocks_size_current_iteration, &brightness_rot90,
                    &contrast_rot90, range_sum, range_sum_squared, domain_sum,
                    domain_sum_squared,
                    rtd_sum_rot90[idx_rb * domain_blocks_length + idx_db]);
                if (error_rot90 < best_error) {
                    best_error = error_rot90;
                    best_domain_block_idx = idx_db;
                    best_range_block_rel_x = range_block->rel_x;
                    best_range_block_rel_y = range_block->rel_y;
                    best_contrast = contrast_rot90;
                    best_brightness = brightness_rot90;
                    best_angle = 90;
                }

                // ROTATION 180
                double brightness_rot180, contrast_rot180, error_rot180;
                error_rot180 = compute_brightness_and_contrast_with_error(
                    range_blocks_size_current_iteration, &brightness_rot180,
                    &contrast_rot180, range_sum, range_sum_squared, domain_sum,
                    domain_sum_squared,
                    rtd_sum_rot180[idx_rb * domain_blocks_length + idx_db]);
                if (error_rot180 < best_error) {
                    best_error = error_rot180;
                    best_domain_block_idx = idx_db;
                    best_range_block_rel_x = range_block->rel_x;
                    best_range_block_rel_y = range_block->rel_y;
                    best_contrast = contrast_rot180;
                    best_brightness = brightness_rot180;
                    best_angle = 180;
                }

                // ROTATION 270
                double brightness_rot270, contrast_rot270, error_rot270;
                error_rot270 = compute_brightness_and_contrast_with_error(
                    range_blocks_size_current_iteration, &brightness_rot270,
                    &contrast_rot270, range_sum, range_sum_squared, domain_sum,
                    domain_sum_squared,
                    rtd_sum_rot270[idx_rb * domain_blocks_length + idx_db]);
                if (error_rot270 < best_error) {
                    best_error = error_rot270;
                    best_domain_block_idx = idx_db;
                    best_range_block_rel_x = range_block->rel_x;
                    best_range_block_rel_y = range_block->rel_y;
                    best_contrast = contrast_rot270;
                    best_brightness = brightness_rot270;
                    best_angle = 270;
                }
            }

            if (best_error > error_threshold &&
                current_quadtree_depth < MAX_QUADTREE_DEPTH) {
                assert(range_block->width >= 2);
                assert(range_block->height >= 2);

                quad2(range_block, range_blocks_next_iteration +
                                   range_blocks_length_next_iteration);
                range_blocks_length_next_iteration += 4;
                has_remaining_range_blocks = true;
            } else {
                struct transformation_t *best_transformation =
                    malloc(sizeof(struct transformation_t));
                best_transformation->range_block =
                    make_block(best_range_block_rel_x, best_range_block_rel_y,
                               range_blocks_size_current_iteration,
                               range_blocks_size_current_iteration);
                best_transformation->domain_block =
                    domain_blocks[best_domain_block_idx];
                best_transformation->brightness = best_brightness;
                best_transformation->contrast = best_contrast;
                best_transformation->angle = best_angle;
                enqueue(transformations, best_transformation);
            }
        }

        // Cleanups
        free(range_blocks_curr_iteration);
        free(rtd_sum_rot0);
        free(rtd_sum_rot90);
        free(rtd_sum_rot180);
        free(rtd_sum_rot270);
        free(range_block_sums);
        free(range_block_sums_squared);
        free(domain_block_sums);
        free(domain_block_sums_squared);
        free(domain_blocks);
        free_prepared_blocks(downsampled_domain_blocks, domain_blocks_length);
    }

    return transformations;
}

void apply_transformation(struct image_t *image,
                          const struct transformation_t *t) {
    assert(t->domain_block.width == t->domain_block.height);
    assert(t->range_block.width == t->range_block.height);

    struct image_t scaled_domain_block = scale_block(
        image, &t->domain_block, t->range_block.width, t->range_block.height);
    struct image_t rotated_domain_block = make_image(t->range_block.height, 0);
    rotate(&rotated_domain_block, &scaled_domain_block, t->angle);

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

    free_image_data(&scaled_domain_block);
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
