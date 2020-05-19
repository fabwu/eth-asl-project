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
#include <immintrin.h>
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

void precompute_sums(double *sums, double *sums_squared,
                     const struct block_t *blocks, const int blocks_length,
                     const struct image_t *image) {
    for (int i = 0; i < blocks_length; ++i) {
        const struct block_t *rb = blocks + i;
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
            range_blocks_next_iteration = (struct block_t *)malloc(
                4 * range_blocks_length_current_iteration *
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

            downsampled_domain_blocks = (struct image_t *)malloc(
                domain_blocks_length * sizeof(struct image_t));

            domain_block_sums =
                (double *)malloc(sizeof(double) * domain_blocks_length);
            domain_block_sums_squared =
                (double *)malloc(sizeof(double) * domain_blocks_length);

            prepare_domain_blocks_norotation(
                downsampled_domain_blocks, domain_block_sums,
                domain_block_sums_squared, image, domain_blocks,
                domain_blocks_length, range_blocks_size_current_iteration);

            __record_double_flops(4);
        }

        // Precomputations on domain / range blocks
        {
            range_block_sums = (double *)malloc(
                sizeof(double) * range_blocks_length_current_iteration);
            assert(range_block_sums != NULL);
            range_block_sums_squared = (double *)malloc(
                sizeof(double) * range_blocks_length_current_iteration);
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

            int best_range_block_rel_x = -1;
            int best_range_block_rel_y = -1;
            int best_domain_block_idx = -1;
            double best_contrast = -1;
            double best_brightness = -1;
            int best_angle = -1;

            const double range_sum = range_block_sums[idx_rb];
            const double range_sum_squared = range_block_sums_squared[idx_rb];
            const double sr_x_2 = 2 * range_sum;
            __record_double_flops(1);

            int rtd_start_rb =
                range_block->rel_y * image->size + range_block->rel_x;

            for (size_t idx_db = 0; idx_db < domain_blocks_length; ++idx_db) {
                assert(domain_blocks[idx_db].width == 2 * range_block->width);
                struct image_t *downsampled_db =
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

                /************************ BEGIN precompute rtd
                 * *************************/
                int rtd_idx_rb = rtd_start_rb;
                int dbs = downsampled_db->size;
                int dbs_dbs = dbs * dbs;

                __m256d v_rtd = _mm256_setzero_pd();
                for (int i = 0; i < dbs; i++) {
                    for (int j = 0; j < dbs; ++j) {
                        int idx_db_0 = i*dbs + j;
                        int idx_db_90 = dbs_dbs - j*dbs - dbs + i;
                        int idx_db_180 = dbs_dbs - i*dbs - j - 1;
                        int idx_db_270 = j*dbs + dbs - i - 1;

                        const __m256i offset_db = _mm256_set_epi64x(
                            idx_db_270, idx_db_180, idx_db_90, idx_db_0);
                        __m256d v_di = _mm256_i64gather_pd(downsampled_db->data,
                                                           offset_db, 8);

//                        __m256d v_di;
//                        v_di[0] = downsampled_db->data[idx_db_0];
//                        v_di[1] = downsampled_db->data[idx_db_90];
//                        v_di[2] = downsampled_db->data[idx_db_180];
//                        v_di[3] = downsampled_db->data[idx_db_270];

                        double ri = image->data[rtd_idx_rb];
                        const __m256d v_ri = _mm256_set1_pd(ri);

                        v_rtd = _mm256_fmadd_pd(v_ri, v_di, v_rtd);

                        rtd_idx_rb += 1;
                    }
                    rtd_idx_rb += image->size - dbs;
                }
                __record_double_flops(dbs * dbs * 8);

                /************************ END precompute rtd
                 * *************************/

                __m256d v_contrast, v_brightness, v_error, v_num_pixels,
                    v_minus_sd_x_sr, v_sd_x_2, v_denominator_inv, v_minus_rtd,
                    v_minus_domain_sum, v_domain_sum_squared, v_range_sum,
                    v_num_pixels_of_blocks_inv, v_a1, v_a2, v_a3, v_a4, v_a5,
                    v_minus_sr_x_2, v_range_sum_squared;

                v_num_pixels = _mm256_set1_pd((double)num_pixels);
                v_minus_sd_x_sr = _mm256_set1_pd(-sd_x_sr);
                /* TODO:l make this more efficient */
                v_denominator_inv = _mm256_set1_pd(denominator_inv);

                // compute contrast
                v_contrast =
                    _mm256_fmadd_pd(v_num_pixels, v_rtd, v_minus_sd_x_sr);
                v_contrast = _mm256_mul_pd(v_contrast, v_denominator_inv);
                __record_double_flops(12);

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
                v_minus_rtd = _mm256_sub_pd(_mm256_set1_pd(0.0), v_rtd);
                v_sd_x_2 = _mm256_set1_pd(sd_x_2);
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
                    if (v_error[i] < min_error && fabs(v_contrast[i]) <= 1) {
                        min_index = i;
                        min_error = v_error[i];
                    }
                }

                if (min_error < best_error) {
                    best_error = min_error;
                    best_domain_block_idx = idx_db;
                    best_range_block_rel_x = range_block->rel_x;
                    best_range_block_rel_y = range_block->rel_y;
                    best_contrast = v_contrast[min_index];
                    best_brightness = v_brightness[min_index];
                    best_angle = 90 * min_index;
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
                    (struct transformation_t *)malloc(
                        sizeof(struct transformation_t));
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
