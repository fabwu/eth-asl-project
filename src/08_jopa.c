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
#include <stdio.h>
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
            scaled_image.data[scaled_image_idx] = val * 0.25;
            scaled_image_idx++;
            original_image_idx += 2;
            __record_double_flops(5);
        }
        original_image_idx += image->size;
    }

    return scaled_image;
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

double rtd_generic_with_rot0(const struct image_t *image,
                             const struct image_t *domain_block,
                             const struct block_t *range_block) {
    const int dbs = domain_block->size;

    double rtd_sum1 = 0;
    double rtd_sum2 = 0;

    size_t idx_rb = range_block->rel_y * image->size + range_block->rel_x;

    for (size_t i = 0; i < dbs; i++) {
        for (size_t j = 0; j < dbs; j += 2) {
            const size_t idx_db1 = i * dbs + j;

            double ri1 = image->data[idx_rb];
            double di1 = domain_block->data[idx_db1];
            rtd_sum1 = fma(ri1, di1, rtd_sum1);
            idx_rb++;

            const size_t idx_db2 = i * dbs + j + 1;
            double ri2 = image->data[idx_rb];
            double di2 = domain_block->data[idx_db2];
            rtd_sum2 = fma(ri2, di2, rtd_sum2);
            idx_rb++;
        }
        idx_rb += image->size - dbs;
    }

    __record_double_flops(dbs * dbs * 2);

    return rtd_sum1 + rtd_sum2;
}

double rtd_generic_with_rot90(const struct image_t *image,
                              const struct image_t *domain_block,
                              const struct block_t *range_block) {
    const int dbs = domain_block->size;

    double rtd_sum1 = 0;
    double rtd_sum2 = 0;

    size_t idx_rb1 = range_block->rel_y * image->size + range_block->rel_x;

    for (size_t i = 0; i < dbs; i++) {
        for (size_t j = 0; j < dbs; j += 2) {
            const size_t idx_db1 = (dbs - j - 1) * dbs + i;

            double ri1 = image->data[idx_rb1];
            double di1 = domain_block->data[idx_db1];
            rtd_sum1 = fma(ri1, di1, rtd_sum1);
            idx_rb1++;

            const size_t idx_db2 = (dbs - (j + 1) - 1) * dbs + i;
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

double rtd_generic_with_rot180(const struct image_t *image,
                               const struct image_t *domain_block,
                               const struct block_t *range_block) {
    const int dbs = domain_block->size;

    double rtd_sum1 = 0;
    double rtd_sum2 = 0;

    size_t idx_rb1 = range_block->rel_y * image->size + range_block->rel_x;

    for (size_t i = 0; i < dbs; i++) {
        for (size_t j = 0; j < dbs; j += 2) {
            const size_t idx_db1 = (dbs - i - 1) * dbs + (dbs - j - 1);

            double ri1 = image->data[idx_rb1];
            double di1 = domain_block->data[idx_db1];
            rtd_sum1 = fma(ri1, di1, rtd_sum1);
            idx_rb1++;

            const size_t idx_db2 = (dbs - i - 1) * dbs + (dbs - (j + 1) - 1);
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

double rtd_generic_with_rot270(const struct image_t *image,
                               const struct image_t *domain_block,
                               const struct block_t *range_block) {
    const int dbs = domain_block->size;

    double rtd_sum1 = 0;
    double rtd_sum2 = 0;

    size_t idx_rb1 = range_block->rel_y * image->size + range_block->rel_x;

    for (size_t i = 0; i < dbs; i++) {
        for (size_t j = 0; j < dbs; j += 2) {
            const size_t idx_db1 = j * dbs + (dbs - i - 1);

            double ri1 = image->data[idx_rb1];
            double di1 = domain_block->data[idx_db1];
            rtd_sum1 = fma(ri1, di1, rtd_sum1);
            idx_rb1++;

            const size_t idx_db2 = (j + 1) * dbs + (dbs - i - 1);
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
    double *range_block_sums;
    double *range_block_sums_squared;
    double *domain_block_sums;
    double *domain_block_sums_squared;
    int num_pixels;
    double num_pixels_of_blocks_inv;

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
            num_pixels = range_blocks_size_current_iteration *
                         range_blocks_size_current_iteration;
            num_pixels_of_blocks_inv = 1.0 / num_pixels;

            range_blocks_length_next_iteration = 0;
            range_blocks_size_next_iteration /= 2;
            range_blocks_next_iteration =
                malloc(4 * range_blocks_length_current_iteration *
                       sizeof(struct block_t));

            __record_double_flops(2);
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
            const double sr_x_2 = 2 * range_sum;
            __record_double_flops(1);

            for (size_t idx_db = 0; idx_db < domain_blocks_length; ++idx_db) {
                assert(domain_blocks[idx_db].width == 2 * range_block->width);
                struct image_t *downsampled_domain_block =
                    downsampled_domain_blocks + idx_db;

                const double domain_sum = domain_block_sums[idx_db];
                const double domain_sum_squared =
                    domain_block_sums_squared[idx_db];

                const double ds_x_ds = domain_sum * domain_sum;
                const double num_pixels_x_dss = num_pixels * domain_sum_squared;
                const double denominator = num_pixels_x_dss - ds_x_ds;
                __record_double_flops(3);

                if (denominator == 0) {
                    double brightness = range_sum * num_pixels_of_blocks_inv;
                    __record_double_flops(1);
                    double error;
                    if (num_pixels == 1) {
                        error = 0.0;
                    } else {
                        error =
                            (range_sum_squared +
                             brightness * (num_pixels * brightness - sr_x_2)) *
                            num_pixels_of_blocks_inv;
                        __record_double_flops(5);
                    }

                    if (error < best_error) {
                        best_error = error;
                        best_domain_block_idx = idx_db;
                        best_range_block_rel_x = range_block->rel_x;
                        best_range_block_rel_y = range_block->rel_y;
                        best_contrast = 0.0;
                        best_brightness = brightness;
                        best_angle = 0;
                    }
                }

                const double sd_x_sr = range_sum * domain_sum;
                const double denominator_inv = 1.0 / denominator;
                const double sd_x_2 = 2 * domain_sum;
                __record_double_flops(3);

                // ROTATION 0
                {
                    double sum_range_times_domain = rtd_generic_with_rot0(
                        image, downsampled_domain_block, range_block);

                    double contrast =
                        fma(num_pixels, sum_range_times_domain, -sd_x_sr) *
                        denominator_inv;
                    double brightness = fma(contrast, -domain_sum, range_sum) *
                                        num_pixels_of_blocks_inv;
                    double a1 = fma(num_pixels, brightness, -sr_x_2);
                    double a2 = fma(brightness, a1, range_sum_squared);
                    double a3 =
                        fma(brightness, sd_x_2, -sum_range_times_domain);
                    double a4 = fma(contrast, domain_sum_squared,
                                    -sum_range_times_domain);
                    double a5 = a3 + a4;
                    double error = fma(a5, contrast, a2);
                    error *= num_pixels_of_blocks_inv;
                    __record_double_flops(18);

                    /* TODO: should fp comparissons count as flop? */
                    if (contrast > 1.0 || contrast < -1.0) error = DBL_MAX;

                    if (error < best_error) {
                        best_error = error;
                        best_domain_block_idx = idx_db;
                        best_range_block_rel_x = range_block->rel_x;
                        best_range_block_rel_y = range_block->rel_y;
                        best_contrast = contrast;
                        best_brightness = brightness;
                        best_angle = 0;
                    }
                }

                // ROTATION 90
                {
                    double sum_range_times_domain = rtd_generic_with_rot90(
                        image, downsampled_domain_block, range_block);

                    double contrast =
                        fma(num_pixels, sum_range_times_domain, -sd_x_sr) *
                        denominator_inv;
                    double brightness = fma(contrast, -domain_sum, range_sum) *
                                        num_pixels_of_blocks_inv;
                    double a1 = fma(num_pixels, brightness, -sr_x_2);
                    double a2 = fma(brightness, a1, range_sum_squared);
                    double a3 =
                        fma(brightness, sd_x_2, -sum_range_times_domain);
                    double a4 = fma(contrast, domain_sum_squared,
                                    -sum_range_times_domain);
                    double a5 = a3 + a4;
                    double error = fma(a5, contrast, a2);
                    error *= num_pixels_of_blocks_inv;
                    __record_double_flops(18);

                    /* TODO: should fp comparissons count as flop? */
                    if (contrast > 1.0 || contrast < -1.0) error = DBL_MAX;

                    if (error < best_error) {
                        best_error = error;
                        best_domain_block_idx = idx_db;
                        best_range_block_rel_x = range_block->rel_x;
                        best_range_block_rel_y = range_block->rel_y;
                        best_contrast = contrast;
                        best_brightness = brightness;
                        best_angle = 90;
                    }
                }

                // ROTATION 180
                {
                    double sum_range_times_domain = rtd_generic_with_rot180(
                        image, downsampled_domain_block, range_block);

                    double contrast =
                        fma(num_pixels, sum_range_times_domain, -sd_x_sr) *
                        denominator_inv;
                    double brightness = fma(contrast, -domain_sum, range_sum) *
                                        num_pixels_of_blocks_inv;
                    double a1 = fma(num_pixels, brightness, -sr_x_2);
                    double a2 = fma(brightness, a1, range_sum_squared);
                    double a3 =
                        fma(brightness, sd_x_2, -sum_range_times_domain);
                    double a4 = fma(contrast, domain_sum_squared,
                                    -sum_range_times_domain);
                    double a5 = a3 + a4;
                    double error = fma(a5, contrast, a2);
                    error *= num_pixels_of_blocks_inv;
                    __record_double_flops(18);

                    /* TODO: should fp comparissons count as flop? */
                    if (contrast > 1.0 || contrast < -1.0) error = DBL_MAX;

                    if (error < best_error) {
                        best_error = error;
                        best_domain_block_idx = idx_db;
                        best_range_block_rel_x = range_block->rel_x;
                        best_range_block_rel_y = range_block->rel_y;
                        best_contrast = contrast;
                        best_brightness = brightness;
                        best_angle = 180;
                    }
                }

                // ROTATION 270
                {
                    double sum_range_times_domain = rtd_generic_with_rot270(
                        image, downsampled_domain_block, range_block);

                    double contrast =
                        fma(num_pixels, sum_range_times_domain, -sd_x_sr) *
                        denominator_inv;
                    double brightness = fma(contrast, -domain_sum, range_sum) *
                                        num_pixels_of_blocks_inv;
                    double a1 = fma(num_pixels, brightness, -sr_x_2);
                    double a2 = fma(brightness, a1, range_sum_squared);
                    double a3 =
                        fma(brightness, sd_x_2, -sum_range_times_domain);
                    double a4 = fma(contrast, domain_sum_squared,
                                    -sum_range_times_domain);
                    double a5 = a3 + a4;
                    double error = fma(a5, contrast, a2);
                    error *= num_pixels_of_blocks_inv;
                    __record_double_flops(18);

                    /* TODO: should fp comparissons count as flop? */
                    if (contrast > 1.0 || contrast < -1.0) error = DBL_MAX;

                    if (error < best_error) {
                        best_error = error;
                        best_domain_block_idx = idx_db;
                        best_range_block_rel_x = range_block->rel_x;
                        best_range_block_rel_y = range_block->rel_y;
                        best_contrast = contrast;
                        best_brightness = brightness;
                        best_angle = 270;
                    }
                }
            }

            if (best_error > error_threshold &&
                current_quadtree_depth < MAX_QUADTREE_DEPTH &&
                range_blocks_size_next_iteration % 2 == 0) {
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
