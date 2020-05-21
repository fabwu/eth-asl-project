#include <assert.h>
#include <float.h>
#include <immintrin.h>
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
#define MAX_QUADTREE_DEPTH 8
#define MIN_RANGE_BLOCK_SIZE 4

#define ALLOCATE(size) (aligned_alloc(32, size))
#define BLOCK_CORD_REL_X(block_id, block_size, image_size) \
    (block_id % (image_size / block_size)) * block_size
#define BLOCK_CORD_REL_Y(block_id, block_size, image_size) \
    ((int)(block_id / (image_size / block_size))) * block_size
#define BLOCK_CORD_X(block_id, block_size, image_size) \
    (block_id % (image_size / block_size))
#define BLOCK_CORD_Y(block_id, block_size, image_size) \
    ((int)(block_id / (image_size / block_size)))

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

void scale_block(double *out, const double *image, const int image_size,
                 const int block_rel_x, const int block_rel_y,
                 const int block_size) {
    assert(out != NULL);
    assert(image_size >= 4);
    assert(image_size % 4 == 0);

    size_t scaled_image_idx = 0;
    size_t original_image_idx = block_rel_y * image_size + block_rel_x;

    __m256d v_row1, v_row2, v_one_fourth;
    v_one_fourth = _mm256_set1_pd(0.25);
    for (int y = 0; y < block_size; y += 2) {
        for (int x = 0; x < block_size; x += 4) {
            v_row1 = _mm256_load_pd(image + original_image_idx);
            v_row2 = _mm256_load_pd(image + original_image_idx + image_size);
            v_row2 = _mm256_hadd_pd(v_row1, v_row2);
            v_row2 = _mm256_mul_pd(v_row2, v_one_fourth);

            out[scaled_image_idx] = v_row2[0] + v_row2[1];
            out[scaled_image_idx + 1] = v_row2[2] + v_row2[3];

            scaled_image_idx += 2;
            original_image_idx += 4;

            __record_double_flops(10);
        }
        original_image_idx += image_size;
    }
}

void prepare_domain_blocks_norotation(double *prepared_domain_blocks,
                                      double *sums, double *sums_squared,
                                      const struct image_t *image,
                                      const int domain_block_size,
                                      const int domain_blocks_length,
                                      const int range_block_size) {
    for (size_t i = 0; i < domain_blocks_length; ++i) {
        const int db_x = BLOCK_CORD_REL_X(i, domain_block_size, image->size);
        const int db_y = BLOCK_CORD_REL_Y(i, domain_block_size, image->size);
        scale_block(
            prepared_domain_blocks + i * range_block_size * range_block_size,
            image->data, image->size, db_x, db_y, domain_block_size);

        double sum = 0;
        double sum_squared = 0;
        for (size_t j = 0; j < range_block_size * range_block_size; j++) {
            double val = *(prepared_domain_blocks +
                           i * range_block_size * range_block_size + j);
            sum += val;
            sum_squared += val * val;
        }
        __record_double_flops(range_block_size * range_block_size * 3);

        sums[i] = sum;
        sums_squared[i] = sum_squared;
    }
}

void quad3(int *list, const int id, const int curr_block_size,
           const int image_size) {
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

void load_block(double *ret_out, double *ret_sum, double *ret_sum_squared,
                const int id, const int block_size,
                const struct image_t *image) {
    int block_rel_x = BLOCK_CORD_REL_X(id, block_size, image->size);
    int block_rel_y = BLOCK_CORD_REL_Y(id, block_size, image->size);

    double sum = 0;
    double sum_squared = 0;

    int idx = 0;
    int idx_in_image = block_rel_y * image->size + block_rel_x;
    for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
            double val = image->data[idx_in_image];
            sum += val;
            sum_squared = fma(val, val, sum_squared);
            ret_out[idx] = val;
            idx++;
            idx_in_image++;
        }
        idx_in_image += image->size - block_size;
    }

    __record_double_flops(block_size * block_size * 3);

    *ret_sum = sum;
    *ret_sum_squared = sum_squared;
}

double rtd_simd(const double *domain_block, const double *range_block,
                const int size) {
    __m256d v_rtd_sum_0 = _mm256_setzero_pd();
    __m256d v_rtd_sum_1 = _mm256_setzero_pd();
    __m256d v_rtd_sum_2 = _mm256_setzero_pd();
    __m256d v_rtd_sum_3 = _mm256_setzero_pd();

    for (int i = 0; i < size; i += 4) {
        for (int j = 0; j < size; j += 4) {
            const double *ri_start = range_block + j;
            const double *di_start = domain_block + j;

            __m256d v_ri_0 = _mm256_load_pd(ri_start + i * size);
            __m256d v_di_0 = _mm256_load_pd(di_start + i * size);
            v_rtd_sum_0 = _mm256_fmadd_pd(v_ri_0, v_di_0, v_rtd_sum_0);

            __m256d v_ri_1 = _mm256_load_pd(ri_start + (i + 1) * size);
            __m256d v_di_1 = _mm256_load_pd(di_start + (i + 1) * size);
            v_rtd_sum_1 = _mm256_fmadd_pd(v_ri_1, v_di_1, v_rtd_sum_1);

            __m256d v_ri_2 = _mm256_load_pd(ri_start + (i + 2) * size);
            __m256d v_di_2 = _mm256_load_pd(di_start + (i + 2) * size);
            v_rtd_sum_2 = _mm256_fmadd_pd(v_ri_2, v_di_2, v_rtd_sum_2);

            __m256d v_ri_3 = _mm256_load_pd(ri_start + (i + 3) * size);
            __m256d v_di_3 = _mm256_load_pd(di_start + (i + 3) * size);
            v_rtd_sum_3 = _mm256_fmadd_pd(v_ri_3, v_di_3, v_rtd_sum_3);
        }
    }

    // Now create the sum
    __m256d v_tmp1 = _mm256_add_pd(v_rtd_sum_0, v_rtd_sum_1);
    __m256d v_tmp2 = _mm256_add_pd(v_rtd_sum_2, v_rtd_sum_3);
    __m256d v_rtd = _mm256_add_pd(v_tmp1, v_tmp2);
    v_rtd = _mm256_hadd_pd(v_rtd, v_rtd);
    double rtd = v_rtd[0] + v_rtd[2];

    __record_double_flops(size * size * 2 + 2 + 16);
    return rtd;
}

struct queue *compress(const struct image_t *image, const int error_threshold) {
    const int initial_domain_block_size =
        image->size / (int)pow(2.0, (double)MIN_QUADTREE_DEPTH);
    const int initial_range_block_size = initial_domain_block_size / 2;

    const size_t initial_range_blocks_length =
        (image->size / initial_range_block_size) *
        (image->size / initial_range_block_size);

    // Forward declarations
    size_t domain_blocks_length;
    int domain_block_size_current_iteration;
    int num_pixels;
    double num_pixels_of_blocks_inv;

    // Queue for saving transformations
    struct queue *transformations =
        (struct queue *)ALLOCATE(sizeof(struct queue));
    *transformations = make_queue();

    int range_blocks_length_current_iteration = -1;
    int range_blocks_size_current_iteration = -1;

    int range_blocks_length_next_iteration = initial_range_blocks_length;
    int range_blocks_size_next_iteration = initial_range_block_size;

    const int upper_bound_domain_blocks = (int)pow(4.0, MAX_QUADTREE_DEPTH);
    const int upper_bound_range_blocks = (int)pow(4.0, MAX_QUADTREE_DEPTH + 1);

    double *rtd_rot = ALLOCATE(4 * upper_bound_domain_blocks * sizeof(double));

    double *prep_domain_blocks_0 =
        ALLOCATE(image->size * image->size * sizeof(double));
    double *prep_domain_blocks_90 =
        ALLOCATE(image->size * image->size * sizeof(double));
    double *prep_domain_blocks_180 =
        ALLOCATE(image->size * image->size * sizeof(double));
    double *prep_domain_blocks_270 =
        ALLOCATE(image->size * image->size * sizeof(double));

    double *domain_block_sums =
        ALLOCATE(sizeof(double) * upper_bound_domain_blocks);
    double *domain_block_sums_squared =
        ALLOCATE(sizeof(double) * upper_bound_domain_blocks);

    int *range_blocks_idx_curr_iteration =
        ALLOCATE(upper_bound_range_blocks * sizeof(int));
    int *range_blocks_idx_next_iteration =
        ALLOCATE(upper_bound_range_blocks * sizeof(int));
    for (int i = 0; i < initial_range_blocks_length; ++i) {
        range_blocks_idx_next_iteration[i] = i;
    }

    double *prepared_range_block =
        ALLOCATE(sizeof(double) * image->size * image->size / 4);

    for (int current_quadtree_depth = MIN_QUADTREE_DEPTH;
         current_quadtree_depth <= MAX_QUADTREE_DEPTH;
         ++current_quadtree_depth) {
        if (range_blocks_length_next_iteration == 0) break;

        range_blocks_length_current_iteration =
            range_blocks_length_next_iteration;
        memcpy(range_blocks_idx_curr_iteration, range_blocks_idx_next_iteration,
               sizeof(int) * range_blocks_length_next_iteration);
        range_blocks_size_current_iteration = range_blocks_size_next_iteration;
        num_pixels = range_blocks_size_current_iteration *
                     range_blocks_size_current_iteration;
        num_pixels_of_blocks_inv = 1.0 / num_pixels;

        range_blocks_length_next_iteration = 0;
        range_blocks_size_next_iteration /= 2;

        __record_double_flops(2);

        // Prepare the domain blocks
        domain_block_size_current_iteration =
            2 * range_blocks_size_current_iteration;
        assert(domain_block_size_current_iteration >= 2);
        domain_blocks_length =
            (image->size / domain_block_size_current_iteration) *
            (image->size / domain_block_size_current_iteration);

        prepare_domain_blocks_norotation(
            prep_domain_blocks_0, domain_block_sums, domain_block_sums_squared,
            image, domain_block_size_current_iteration, domain_blocks_length,
            range_blocks_size_current_iteration);

        for (int i = 0; i < domain_blocks_length; ++i) {
            const double *prep_db_rot_0 =
                prep_domain_blocks_0 + i * range_blocks_size_current_iteration *
                                           range_blocks_size_current_iteration;

            rotate_raw_90(prep_domain_blocks_90 +
                              i * range_blocks_size_current_iteration *
                                  range_blocks_size_current_iteration,
                          prep_db_rot_0, range_blocks_size_current_iteration);

            rotate_raw_180(prep_domain_blocks_180 +
                               i * range_blocks_size_current_iteration *
                                   range_blocks_size_current_iteration,
                           prep_db_rot_0, range_blocks_size_current_iteration);

            rotate_raw_270(prep_domain_blocks_270 +
                               i * range_blocks_size_current_iteration *
                                   range_blocks_size_current_iteration,
                           prep_db_rot_0, range_blocks_size_current_iteration);
        }
        __record_double_flops(4);

        // Process each range block
        for (size_t idx_rb = 0; idx_rb < range_blocks_length_current_iteration;
             ++idx_rb) {
            int curr_relative_rb_idx = range_blocks_idx_curr_iteration[idx_rb];
            double range_sum, range_sum_squared;
            load_block(prepared_range_block, &range_sum, &range_sum_squared,
                       curr_relative_rb_idx,
                       range_blocks_size_current_iteration, image);

            double best_error = DBL_MAX;

            int best_range_block_idx = -1;
            int best_domain_block_idx = -1;
            double best_contrast = -1;
            double best_brightness = -1;
            int best_angle = -1;

            const double sr_x_2 = 2 * range_sum;
            __record_double_flops(1);

            for (size_t idx_db = 0; idx_db < domain_blocks_length; ++idx_db) {
                double *prep_db_rot_0 =
                    prep_domain_blocks_0 +
                    idx_db * range_blocks_size_current_iteration *
                        range_blocks_size_current_iteration;
                rtd_rot[4 * idx_db] =
                    rtd_simd(prep_db_rot_0, prepared_range_block,
                             range_blocks_size_current_iteration);
                double *prep_db_rot_90 =
                    prep_domain_blocks_90 +
                    idx_db * range_blocks_size_current_iteration *
                        range_blocks_size_current_iteration;
                rtd_rot[4 * idx_db + 1] =
                    rtd_simd(prep_db_rot_90, prepared_range_block,
                             range_blocks_size_current_iteration);
                double *prep_db_rot_180 =
                    prep_domain_blocks_180 +
                    idx_db * range_blocks_size_current_iteration *
                        range_blocks_size_current_iteration;
                rtd_rot[4 * idx_db + 2] =
                    rtd_simd(prep_db_rot_180, prepared_range_block,
                             range_blocks_size_current_iteration);
                double *prep_db_rot_270 =
                    prep_domain_blocks_270 +
                    idx_db * range_blocks_size_current_iteration *
                        range_blocks_size_current_iteration;
                rtd_rot[4 * idx_db + 3] =
                    rtd_simd(prep_db_rot_270, prepared_range_block,
                             range_blocks_size_current_iteration);
            }

            for (size_t idx_db = 0; idx_db < domain_blocks_length; ++idx_db) {
                const double domain_sum = domain_block_sums[idx_db];
                const double domain_sum_squared =
                    domain_block_sums_squared[idx_db];

                const double ds_x_ds = domain_sum * domain_sum;
                const double num_pixels_x_dss = num_pixels * domain_sum_squared;
                const double denominator = num_pixels_x_dss - ds_x_ds;

                __m256d v_contrast, v_brightness, v_error, v_num_pixels,
                    v_minus_sd_x_sr, v_sd_x_2, v_denominator_inv, v_minus_rtd,
                    v_minus_domain_sum, v_domain_sum_squared, v_range_sum,
                    v_rtd, v_num_pixels_of_blocks_inv, v_a1, v_a2, v_a3, v_a4,
                    v_a5, v_minus_sr_x_2, v_range_sum_squared;

                v_rtd = _mm256_load_pd(rtd_rot + 4 * idx_db);

                if (denominator == 0) {
                    assert(0);
                } else {
                    v_denominator_inv = _mm256_set1_pd(1.0 / denominator);
                    __record_double_flops(1);

                    v_contrast =
                        _mm256_fmadd_pd(v_num_pixels, v_rtd, v_minus_sd_x_sr);
                    v_contrast = _mm256_mul_pd(v_contrast, v_denominator_inv);
                    __record_double_flops(12);
                }

                v_num_pixels = _mm256_set1_pd((double)num_pixels);
                v_minus_sd_x_sr = _mm256_set1_pd(-range_sum * domain_sum);
                __record_double_flops(1);

                // compute brightness
                v_minus_domain_sum = _mm256_set1_pd(-domain_sum);
                v_range_sum = _mm256_set1_pd(range_sum);
                v_num_pixels_of_blocks_inv =
                    _mm256_set1_pd(num_pixels_of_blocks_inv);
                v_brightness = _mm256_fmadd_pd(v_contrast, v_minus_domain_sum,
                                               v_range_sum);
                v_brightness =
                    _mm256_mul_pd(v_brightness, v_num_pixels_of_blocks_inv);
                __record_double_flops(12);

                // compute error
                v_minus_sr_x_2 = _mm256_set1_pd(-sr_x_2);
                v_minus_rtd = _mm256_sub_pd(_mm256_setzero_pd(), v_rtd);
                v_sd_x_2 = _mm256_set1_pd(2 * domain_sum);
                __record_double_flops(1);
                v_domain_sum_squared = _mm256_set1_pd(domain_sum_squared);
                v_range_sum_squared = _mm256_set1_pd(range_sum_squared);
                v_a1 =
                    _mm256_fmadd_pd(v_num_pixels, v_brightness, v_minus_sr_x_2);
                v_a2 = _mm256_fmadd_pd(v_brightness, v_a1, v_range_sum_squared);
                v_a3 = _mm256_fmadd_pd(v_brightness, v_sd_x_2, v_minus_rtd);
                v_a4 = _mm256_fmadd_pd(v_contrast, v_domain_sum_squared,
                                       v_minus_rtd);
                v_a5 = _mm256_add_pd(v_a3, v_a4);
                v_error = _mm256_fmadd_pd(v_a5, v_contrast, v_a2);
                v_error = _mm256_mul_pd(v_error, v_num_pixels_of_blocks_inv);
                __record_double_flops(48);

                int min_index = -1;
                double min_error = DBL_MAX;
                for (int i = 0; i < 4; i++) {
                    if (v_error[i] < min_error && v_contrast[i] < 1.0 &&
                        v_contrast[i] > -1.0) {
                        min_index = i;
                        min_error = v_error[i];
                    }
                }

                double contrast = v_contrast[min_index];
                double brightness = v_brightness[min_index];
                if (min_error < best_error) {
                    best_error = min_error;
                    best_domain_block_idx = idx_db;
                    best_range_block_idx = curr_relative_rb_idx;
                    best_contrast = contrast;
                    best_brightness = brightness;
                    best_angle = 90 * min_index;
                }
            }

            if (best_error > error_threshold &&
                current_quadtree_depth < MAX_QUADTREE_DEPTH &&
                range_blocks_size_next_iteration >= MIN_RANGE_BLOCK_SIZE &&
                range_blocks_size_next_iteration % 2 == 0) {
                assert(range_blocks_size_current_iteration >= 2);

                quad3(range_blocks_idx_next_iteration +
                          range_blocks_length_next_iteration,
                      curr_relative_rb_idx, range_blocks_size_current_iteration,
                      image->size);

                range_blocks_length_next_iteration += 4;
            } else {
                struct transformation_t *best_transformation =
                    (struct transformation_t *)ALLOCATE(
                        sizeof(struct transformation_t));

                int rb_rel_x = BLOCK_CORD_REL_X(
                    best_range_block_idx, range_blocks_size_current_iteration,
                    image->size);
                int rb_rel_y = BLOCK_CORD_REL_Y(
                    best_range_block_idx, range_blocks_size_current_iteration,
                    image->size);

                best_transformation->range_block = make_block(
                    rb_rel_x, rb_rel_y, range_blocks_size_current_iteration,
                    range_blocks_size_current_iteration);

                int db_rel_x = BLOCK_CORD_REL_X(
                    best_domain_block_idx, domain_block_size_current_iteration,
                    image->size);
                int db_rel_y = BLOCK_CORD_REL_Y(
                    best_domain_block_idx, domain_block_size_current_iteration,
                    image->size);
                best_transformation->domain_block = make_block(
                    db_rel_x, db_rel_y, domain_block_size_current_iteration,
                    domain_block_size_current_iteration);
                best_transformation->brightness = best_brightness;
                best_transformation->contrast = best_contrast;
                best_transformation->angle = best_angle;
                enqueue(transformations, best_transformation);
            }
        }
    }

    free(domain_block_sums);
    free(domain_block_sums_squared);
    free(prep_domain_blocks_90);
    free(prep_domain_blocks_180);
    free(prep_domain_blocks_270);
    free(prep_domain_blocks_0);
    free(rtd_rot);
    free(range_blocks_idx_curr_iteration);
    free(range_blocks_idx_next_iteration);
    free(prepared_range_block);

    return transformations;
}

void apply_transformation(struct image_t *image,
                          const struct transformation_t *t) {
    assert(t->domain_block.width == t->domain_block.height);
    assert(t->range_block.width == t->range_block.height);

    double *scaled_domain_block =
        ALLOCATE(sizeof(double) * t->range_block.width * t->range_block.height);
    scale_block(scaled_domain_block, image->data, image->size,
                t->domain_block.rel_x, t->domain_block.rel_y,
                t->domain_block.height);
    double *rotated_domain_block = ALLOCATE(
        sizeof(double) * t->range_block.height * t->range_block.height);

    switch (t->angle) {
        case 0:
            rotate_raw_0(rotated_domain_block, scaled_domain_block,
                         t->range_block.height);
            break;
        case 90:
            rotate_raw_90(rotated_domain_block, scaled_domain_block,
                          t->range_block.height);
            break;
        case 180:
            rotate_raw_180(rotated_domain_block, scaled_domain_block,
                           t->range_block.height);
            break;
        case 270:
            rotate_raw_270(rotated_domain_block, scaled_domain_block,
                           t->range_block.height);
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
