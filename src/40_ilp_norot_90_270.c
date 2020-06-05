#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "lib/performance.h"
#include "lib/queue.h"
#include "lib/types.h"

/**
 * The inital domain block size is image_size / 2^MIN_QUADTREE_DEPTH
 * Should be >=1, because we can then can assume that the amount of domain
 * blocks is divisible by 2.
 */
#define MIN_QUADTREE_DEPTH 1
#define MAX_QUADTREE_DEPTH 7
#define MIN_RANGE_BLOCK_SIZE 4

static inline void rotate_raw_0(double *out, const double *in, int size) {
    int m = size;
    int n = size;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            out[i * n + j] = in[i * n + j];
        }
    }
}

static inline void rotate_raw_90(double *out, const double *in, int size) {
    int m = size;
    int n = size;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // first row has to be last column
            // (too be honest it was trial and error)
            out[j * n + (m - i - 1)] = in[i * n + j];
        }
    }
}

static inline void rotate_raw_180(double *out, const double *in, int size) {
    int m = size;
    int n = size;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // first row has to be last row reversed
            out[(m - i - 1) * n + (n - j - 1)] = in[i * n + j];
        }
    }
}

static inline void rotate_raw_270(double *out, const double *in, int size) {
    int m = size;
    int n = size;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // first row has to be first column reversed
            out[(m - j - 1) * n + i] = in[i * n + j];
        }
    }
}

static void load_block(double *ret_out, double *ret_sum, double *ret_sum_squared, const int id,
                       const int block_size, const double *image, const int image_size) {
    int block_rel_x = BLOCK_CORD_REL_X(id, block_size, image_size);
    int block_rel_y = BLOCK_CORD_REL_Y(id, block_size, image_size);

    double sum = 0;
    double sum_squared = 0;

    int idx = 0;
    int idx_in_image = block_rel_y * image_size + block_rel_x;
    for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
            double val = image[idx_in_image];
            sum += val;
            sum_squared = fma(val, val, sum_squared);
            ret_out[idx] = val;
            idx++;
            idx_in_image++;
        }
        idx_in_image += image_size - block_size;
    }

    __record_double_flops(block_size * block_size * 3);

    *ret_sum = sum;
    *ret_sum_squared = sum_squared;
}

static void scale_block(double *out, const double *image, const int image_size, const int block_rel_x,
                        const int block_rel_y, const int block_size) {
    size_t scaled_image_idx = 0;
    size_t original_image_idx = block_rel_y * image_size + block_rel_x;

    for (int y = 0; y < block_size; y += 2) {
        for (int x = 0; x < block_size; x += 2) {
            double val = 0.0;
            val += image[original_image_idx];
            val += image[original_image_idx + 1];
            val += image[original_image_idx + image_size];
            val += image[original_image_idx + image_size + 1];
            out[scaled_image_idx] = val * 0.25;
            scaled_image_idx++;
            original_image_idx += 2;
            __record_double_flops(5);
        }
        original_image_idx += image_size;
    }
}

static void prepare_domain_blocks_norotation(double *prepared_domain_blocks, double *sums,
                                             double *sums_squared, const double *image, const int image_size,
                                             const int domain_block_size, const int domain_blocks_length,
                                             const int range_block_size) {
    for (size_t i = 0; i < domain_blocks_length; ++i) {
        const int db_x = BLOCK_CORD_REL_X(i, domain_block_size, image_size);
        const int db_y = BLOCK_CORD_REL_Y(i, domain_block_size, image_size);
        scale_block(prepared_domain_blocks + i * range_block_size * range_block_size, image, image_size, db_x,
                    db_y, domain_block_size);

        double sum = 0;
        double sum_squared = 0;
        for (size_t j = 0; j < range_block_size * range_block_size; j++) {
            double val = *(prepared_domain_blocks + i * range_block_size * range_block_size + j);
            sum += val;
            sum_squared += val * val;
        }
        __record_double_flops(range_block_size * range_block_size * 3);

        sums[i] = sum;
        sums_squared[i] = sum_squared;
    }
}

static void quad3(int *list, const int id, const int curr_block_size, const int image_size) {
    const int blocks_per_row = image_size / curr_block_size;
    int rel_x = BLOCK_CORD_X(id, curr_block_size, image_size);
    int rel_y = BLOCK_CORD_Y(id, curr_block_size, image_size);

    int next_id_1 = 4 * rel_y * blocks_per_row + 2 * rel_x;
    int next_id_2 = next_id_1 + 1;
    int next_id_3 = next_id_1 + 2 * blocks_per_row;
    int next_id_4 = next_id_3 + 1;
    *(list + 0) = next_id_1;
    *(list + 1) = next_id_2;
    *(list + 2) = next_id_3;
    *(list + 3) = next_id_4;
}

static struct queue *compress(const struct image_t *image, const int error_threshold) {
    const int initial_domain_block_size = image->size / (int)pow(2.0, (double)MIN_QUADTREE_DEPTH);
    const int initial_range_block_size = initial_domain_block_size / 2;

    const size_t initial_range_blocks_length =
        (image->size / initial_range_block_size) * (image->size / initial_range_block_size);

    // Forward declarations
    size_t domain_blocks_length;
    int domain_block_size_current_iteration;
    int num_pixels;
    double num_pixels_of_blocks_inv;

    // Queue for saving transformations
    struct queue *transformations = (struct queue *)ALLOCATE(sizeof(struct queue));
    *transformations = make_queue();

    int range_blocks_length_current_iteration = -1;
    int range_blocks_size_current_iteration = -1;

    int range_blocks_length_next_iteration = initial_range_blocks_length;
    int range_blocks_size_next_iteration = initial_range_block_size;

    // Allocate memory
    const int upper_bound_domain_blocks = (int)pow(4.0, MAX_QUADTREE_DEPTH);
    const int upper_bound_range_blocks = (int)pow(4.0, MAX_QUADTREE_DEPTH + 1);

    double *prep_domain_blocks = ALLOCATE(sizeof(double) * image->size * image->size / 4);

    double *domain_block_sums = ALLOCATE(sizeof(double) * upper_bound_domain_blocks);
    double *domain_block_sums_squared = ALLOCATE(sizeof(double) * upper_bound_domain_blocks);
    int *range_blocks_idx_curr_iteration = ALLOCATE(upper_bound_range_blocks * sizeof(int));
    int *range_blocks_idx_next_iteration = ALLOCATE(upper_bound_range_blocks * sizeof(int));
    for (int i = 0; i < initial_range_blocks_length; ++i) {
        range_blocks_idx_next_iteration[i] = i;
    }

    double *prepared_range_block = ALLOCATE(sizeof(double) * image->size * image->size / 4);

    for (int current_quadtree_depth = MIN_QUADTREE_DEPTH; current_quadtree_depth <= MAX_QUADTREE_DEPTH;
         ++current_quadtree_depth) {
        if (range_blocks_length_next_iteration == 0) break;

        range_blocks_length_current_iteration = range_blocks_length_next_iteration;
        memcpy(range_blocks_idx_curr_iteration, range_blocks_idx_next_iteration,
               sizeof(int) * range_blocks_length_next_iteration);
        range_blocks_size_current_iteration = range_blocks_size_next_iteration;
        num_pixels = range_blocks_size_current_iteration * range_blocks_size_current_iteration;
        num_pixels_of_blocks_inv = 1.0 / num_pixels;

        range_blocks_length_next_iteration = 0;
        range_blocks_size_next_iteration /= 2;

        __record_double_flops(2);

        domain_block_size_current_iteration = 2 * range_blocks_size_current_iteration;
        assert(domain_block_size_current_iteration >= 2);
        domain_blocks_length = (image->size / domain_block_size_current_iteration) *
                               (image->size / domain_block_size_current_iteration);
        prepare_domain_blocks_norotation(prep_domain_blocks, domain_block_sums, domain_block_sums_squared,
                                         image->data, image->size, domain_block_size_current_iteration,
                                         domain_blocks_length, range_blocks_size_current_iteration);
        __record_double_flops(4);

        // Process each range block
        for (size_t idx_rb = 0; idx_rb < range_blocks_length_current_iteration; ++idx_rb) {
            int curr_relative_rb_idx = range_blocks_idx_curr_iteration[idx_rb];
            double range_sum, range_sum_squared;
            load_block(prepared_range_block, &range_sum, &range_sum_squared, curr_relative_rb_idx,
                       range_blocks_size_current_iteration, image->data, image->size);

            double best_error = DBL_MAX;

            int best_range_block_idx = -1;
            int best_domain_block_idx = -1;
            double best_contrast = -1;
            double best_brightness = -1;
            int best_angle = -1;

            const double sr_x_2 = 2 * range_sum;
            __record_double_flops(1);

            for (size_t idx_db = 0; idx_db < domain_blocks_length; ++idx_db) {
                double *prep_domain_block = prep_domain_blocks + idx_db *
                                                                 range_blocks_size_current_iteration *
                                                                 range_blocks_size_current_iteration;

                const double domain_sum = domain_block_sums[idx_db];
                const double domain_sum_squared = domain_block_sums_squared[idx_db];

                const double ds_x_ds = domain_sum * domain_sum;
                const double num_pixels_x_dss = num_pixels * domain_sum_squared;
                const double denominator = num_pixels_x_dss - ds_x_ds;
                __record_double_flops(3);

                const double sd_x_sr = range_sum * domain_sum;
                const double denominator_inv = 1.0 / denominator;
                const double sd_x_2 = 2 * domain_sum;
                __record_double_flops(3);

                // BEGIN precompute rtd
                int rtd_idx_rb = 0;
                int dbs = range_blocks_size_current_iteration;
                int dbs_dbs = dbs * dbs;

                double rtd_sum_0_1 = 0;
                double rtd_sum_0_2 = 0;
                double rtd_sum_180_1 = 0;
                double rtd_sum_180_2 = 0;

                int dbs_i = 0;
                for (int i = 0; i < dbs; i++) {
                    int dbs_j = 0;
                    for (int j = 0; j < dbs; j += 2) {
                        int idx_0_db1 = dbs_i + j;
                        int idx_0_db2 = idx_0_db1 + 1;
                        int idx_180_db1 = dbs_dbs - dbs_i - j - 1;
                        int idx_180_db2 = idx_180_db1 - 1;

                        double ri1 = prepared_range_block[rtd_idx_rb];
                        double ri2 = prepared_range_block[rtd_idx_rb + 1];

                        double di_0_1 = prep_domain_block[idx_0_db1];
                        double di_0_2 = prep_domain_block[idx_0_db2];
                        double di_180_1 = prep_domain_block[idx_180_db1];
                        double di_180_2 = prep_domain_block[idx_180_db2];

                        rtd_sum_0_1 += ri1 * di_0_1;
                        rtd_sum_0_2 += ri2 * di_0_2;

                        rtd_sum_180_1 += ri1 * di_180_1;
                        rtd_sum_180_2 += ri2 * di_180_2;

                        rtd_idx_rb += 2;
                        dbs_j = dbs_j + dbs + dbs;
                    }
//                    rtd_idx_rb += dbs;
                    dbs_i += dbs;
                }

                __record_double_flops(dbs * dbs * 4);

                double rtd_0 = rtd_sum_0_1 + rtd_sum_0_2;
                double rtd_180 = rtd_sum_180_1 + rtd_sum_180_2;

                __record_double_flops(4);
                // END precompute rtd

                // ROTATION 0
                {
                    double contrast = fma(num_pixels, rtd_0, -sd_x_sr) * denominator_inv;
                    double brightness = fma(contrast, -domain_sum, range_sum) * num_pixels_of_blocks_inv;
                    double a1 = fma(num_pixels, brightness, -sr_x_2);
                    double a2 = fma(brightness, a1, range_sum_squared);
                    double a3 = fma(brightness, sd_x_2, -rtd_0);
                    double a4 = fma(contrast, domain_sum_squared, -rtd_0);
                    double a5 = a3 + a4;
                    double error = fma(a5, contrast, a2);
                    error *= num_pixels_of_blocks_inv;
                    __record_double_flops(18);

                    if (error < best_error && contrast < 1.0 && contrast > -1.0) {
                        best_error = error;
                        best_domain_block_idx = idx_db;
                        best_range_block_idx = curr_relative_rb_idx;
                        best_contrast = contrast;
                        best_brightness = brightness;
                        best_angle = 0;
                    }
                }

                // ROTATION 90
//                {
//                    double contrast = fma(num_pixels, rtd_90, -sd_x_sr) * denominator_inv;
//                    double brightness = fma(contrast, -domain_sum, range_sum) * num_pixels_of_blocks_inv;
//                    double a1 = fma(num_pixels, brightness, -sr_x_2);
//                    double a2 = fma(brightness, a1, range_sum_squared);
//                    double a3 = fma(brightness, sd_x_2, -rtd_90);
//                    double a4 = fma(contrast, domain_sum_squared, -rtd_90);
//                    double a5 = a3 + a4;
//                    double error = fma(a5, contrast, a2);
//                    error *= num_pixels_of_blocks_inv;
//                    __record_double_flops(18);
//
//                    if (error < best_error && contrast < 1.0 && contrast > -1.0) {
//                        best_error = error;
//                        best_domain_block_idx = idx_db;
//                        best_range_block_idx = curr_relative_rb_idx;
//                        best_contrast = contrast;
//                        best_brightness = brightness;
//                        best_angle = 90;
//                    }
//                }

                // ROTATION 180
                {
                    double contrast = fma(num_pixels, rtd_180, -sd_x_sr) * denominator_inv;
                    double brightness = fma(contrast, -domain_sum, range_sum) * num_pixels_of_blocks_inv;
                    double a1 = fma(num_pixels, brightness, -sr_x_2);
                    double a2 = fma(brightness, a1, range_sum_squared);
                    double a3 = fma(brightness, sd_x_2, -rtd_180);
                    double a4 = fma(contrast, domain_sum_squared, -rtd_180);
                    double a5 = a3 + a4;
                    double error = fma(a5, contrast, a2);
                    error *= num_pixels_of_blocks_inv;
                    __record_double_flops(18);

                    if (error < best_error && contrast < 1.0 && contrast > -1.0) {
                        best_error = error;
                        best_domain_block_idx = idx_db;
                        best_range_block_idx = curr_relative_rb_idx;
                        best_contrast = contrast;
                        best_brightness = brightness;
                        best_angle = 180;
                    }
                }

                // ROTATION 270
//                {
//                    double contrast = fma(num_pixels, rtd_270, -sd_x_sr) * denominator_inv;
//                    double brightness = fma(contrast, -domain_sum, range_sum) * num_pixels_of_blocks_inv;
//                    double a1 = fma(num_pixels, brightness, -sr_x_2);
//                    double a2 = fma(brightness, a1, range_sum_squared);
//                    double a3 = fma(brightness, sd_x_2, -rtd_270);
//                    double a4 = fma(contrast, domain_sum_squared, -rtd_270);
//                    double a5 = a3 + a4;
//                    double error = fma(a5, contrast, a2);
//                    error *= num_pixels_of_blocks_inv;
//                    __record_double_flops(18);
//
//                    if (error < best_error && contrast < 1.0 && contrast > -1.0) {
//                        best_error = error;
//                        best_domain_block_idx = idx_db;
//                        best_range_block_idx = curr_relative_rb_idx;
//                        best_contrast = contrast;
//                        best_brightness = brightness;
//                        best_angle = 270;
//                    }
//                }
            }

            if (best_error > error_threshold && current_quadtree_depth < MAX_QUADTREE_DEPTH &&
                range_blocks_size_next_iteration >= MIN_RANGE_BLOCK_SIZE &&
                range_blocks_size_next_iteration % 2 == 0) {

                quad3(range_blocks_idx_next_iteration + range_blocks_length_next_iteration,
                      curr_relative_rb_idx, range_blocks_size_current_iteration, image->size);

                range_blocks_length_next_iteration += 4;
            } else {
                struct transformation_t *best_transformation =
                    (struct transformation_t *)malloc(sizeof(struct transformation_t));

                int rb_rel_x =
                    BLOCK_CORD_REL_X(best_range_block_idx, range_blocks_size_current_iteration, image->size);
                int rb_rel_y =
                    BLOCK_CORD_REL_Y(best_range_block_idx, range_blocks_size_current_iteration, image->size);

                best_transformation->range_block =
                    make_block(rb_rel_x, rb_rel_y, range_blocks_size_current_iteration,
                               range_blocks_size_current_iteration);

                int db_rel_x =
                    BLOCK_CORD_REL_X(best_domain_block_idx, domain_block_size_current_iteration, image->size);
                int db_rel_y =
                    BLOCK_CORD_REL_Y(best_domain_block_idx, domain_block_size_current_iteration, image->size);

                best_transformation->domain_block =
                    make_block(db_rel_x, db_rel_y, domain_block_size_current_iteration,
                               domain_block_size_current_iteration);

                best_transformation->brightness = best_brightness;
                best_transformation->contrast = best_contrast;
                best_transformation->angle = best_angle;
                enqueue(transformations, best_transformation);
            }
        }
    }

    free(prep_domain_blocks);
    free(domain_block_sums);
    free(domain_block_sums_squared);
    free(range_blocks_idx_curr_iteration);
    free(range_blocks_idx_next_iteration);
    free(prepared_range_block);

    return transformations;
}

static void apply_transformation(struct image_t *image, const struct transformation_t *t) {
    assert(t->domain_block.width == t->domain_block.height);
    assert(t->range_block.width == t->range_block.height);

    double *scaled_domain_block = ALLOCATE(sizeof(double) * t->range_block.width * t->range_block.height);
    scale_block(scaled_domain_block, image->data, image->size, t->domain_block.rel_x, t->domain_block.rel_y,
                t->domain_block.height);
    double *rotated_domain_block = ALLOCATE(sizeof(double) * t->range_block.height * t->range_block.height);

    switch (t->angle) {
        case 0:
            rotate_raw_0(rotated_domain_block, scaled_domain_block, t->range_block.height);
            break;
        case 90:
            rotate_raw_90(rotated_domain_block, scaled_domain_block, t->range_block.height);
            break;
        case 180:
            rotate_raw_180(rotated_domain_block, scaled_domain_block, t->range_block.height);
            break;
        case 270:
            rotate_raw_270(rotated_domain_block, scaled_domain_block, t->range_block.height);
            break;
        default:
            assert(0);
            break;
    }

    for (int i = 0; i < t->range_block.height; ++i) {
        for (int j = 0; j < t->range_block.width; ++j) {
            double value = rotated_domain_block[i * t->range_block.height + j];
            int idx = get_index_in_image(&t->range_block, i, j, image);
            int new_pixel_value = value * t->contrast + t->brightness;
            if (new_pixel_value < 0) new_pixel_value = 0;
            if (new_pixel_value > 255) new_pixel_value = 255;
            image->data[idx] = new_pixel_value;
            __record_double_flops(2);
        }
    }

    free(scaled_domain_block);
    free(rotated_domain_block);
}

static void decompress(struct image_t *decompressed_image, const struct queue *transformations,
                       const int num_iterations) {
    for (int iter = 0; iter < num_iterations; ++iter) {
        const struct queue_node *current = transformations->front;
        while (current != transformations->back) {
            apply_transformation(decompressed_image, (const struct transformation_t *)current->data);
            current = current->next;
        }
    }
}

struct func_suite_t register_suite(void) {
    struct func_suite_t suite = {.compress_func = &compress, .decompress_func = &decompress};
    return suite;
}
