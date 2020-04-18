#include <string>
#include "common.h"
#include <cassert>
#include <vector>
#include <limits>
#include <tuple>

using namespace std;

pair<double, double> solve_2x2_eq_system(double a11, double a12, double a21, double a22, double b1, double b2) {
    /// requires det(A) != 0
    double f1 = a22 == 0 ? 0.0 : a12 / a22;
    double f2 = a11 == 0 ? 0.0 : a21 / a11;

    double x = (b1 - f1 * b2) / (a11 - f1 * a21);
    double y = (b2 - a21 * x) / a22; //(b2 - f2 * b1) / (a22 - f2 * a11);

    return make_pair(x, y);
}


/**
 * Assumes image is a squared image
 * @param image
 * @param image_size
 * @param block_size
 * @return
 */
inline vector<block_t>
create_squared_blocks(const int image_size, const int block_size, int x_offset = 0, int y_offset = 0) {
    assert(image_size % block_size == 0);

    const auto num_blocks = (image_size / block_size) * (image_size / block_size);
    vector<block_t> blocks;
    blocks.reserve(num_blocks);
    for (int i = 0; i < image_size; i += block_size) {
        for (int j = 0; j < image_size; j += block_size) {
            blocks.emplace_back(j + x_offset, i + y_offset, block_size, block_size);
        }
    }

    return blocks;
}

double *compress_block(const image_t &image, block_t source_block, int width, int height) {
    assert(source_block.width >= width);
    assert(source_block.height >= height);
    assert(source_block.width == source_block.height); // just for simplicity
    assert(source_block.width % width == 0); // just for simplicity
    assert(source_block.height % height == 0); // just for simplicity
    assert(source_block.width % 2 == 0); // just for simplicity
    assert(source_block.height % 2 == 0); // just for simplicity
    assert(width == height); // just for simplicity

    const int n = source_block.width / width;
    const auto scaled_source_block = (double *) malloc(height * width * sizeof(double));
    const auto compression_blocks = create_squared_blocks(source_block.width, n, source_block.rel_x,
                                                          source_block.rel_y);
    int scaled_index = 0;
    for (auto b: compression_blocks) {
        double val = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                val += image.data[b.get_index_in_image(i, j, image.size)];
            }
        }
        val /= n * n;

        scaled_source_block[scaled_index] = val;
        scaled_index++;
    }

    return scaled_source_block;
}

pair<double, double> compute_brightness_and_contrast(const image_t &image, double *domain,
                                                     block_t target_block) {
    //assert(source_block.width == target_block.width);
    //assert(source_block.height == target_block.height);
    //assert(source_block.height == source_block.width);
    // todo asserts
    assert(target_block.height == target_block.width);

    // build equation system
    int n = target_block.width;
    double a11 = target_block.width * target_block.height;
    double a12 = 0.0;
    double a21 = 0.0;
    double a22 = 0.0;
    double b1 = 0.0;
    double b2 = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double di = domain[i * n + j];
            double ri = image.data[target_block.get_index_in_image(i, j, image.size)];
            a12 += di;
            a21 += di;
            a22 += di * di;
            b1 += ri;
            b2 += di * ri;
        }
    }

    return solve_2x2_eq_system(a11, a12, a21, a22, b1, b2);
}

inline double diff_block(const image_t &image, block_t source_block, block_t target_block) {
    // Need to compress source_block, such that size(source_block) == size(target_block) in order to compare difference!
    // That means that a square of n pixels need to be compressed to one pixel
    const int n = source_block.width / target_block.width;
    const auto scaled_source_block = compress_block(image, source_block, target_block.width,
                                                    target_block.height);

    double brightness, contrast;
    std::tie(brightness, contrast) = compute_brightness_and_contrast(image, scaled_source_block,
                                                                     target_block);

    double squared_error = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double target_val = image.data[target_block.get_index_in_image(i, j, image.size)];
            double source_val = scaled_source_block[i * n + j] * contrast + brightness;
            double diff = source_val - target_val;
            squared_error += diff * diff;
        }
    }

    free(scaled_source_block);
    return squared_error;
}

vector<transformation_t> compress(const image_t &image) {
    // Goal: Try to map blocks of size block_size_source to blocks of size block_size_target
    const int block_size_target = 8;
    const int block_size_source = 32;
    static_assert(block_size_source >= block_size_target);

    const auto target_blocks = create_squared_blocks(image.size, block_size_target);
    const auto source_blocks = create_squared_blocks(image.size, block_size_source);

    vector<pair<int, int>> block_mappings; // specifies best mapping from source to target as pair of (source_index, target_index)

    // Learn mappings from source to target blocks
    // That is, find a target block for every source block, such that their difference is minimal
    // todo: rotations, i.e. source block may be rotated (equivalently: target block may be rotated)
    for (int target_block_index = 0; target_block_index < target_blocks.size(); ++target_block_index) {
        const auto target_block = target_blocks[target_block_index];

        double best_error = numeric_limits<double>::max();
        int best_source_block_index = -1;

        for (int source_block_index = 0; source_block_index < source_blocks.size(); ++source_block_index) {
            const auto source_block = source_blocks[source_block_index];

            double target_block_error = diff_block(image, source_block, target_block);
            if (best_source_block_index == -1 || target_block_error < best_error) {
                best_source_block_index = source_block_index;
                best_error = target_block_error;
            }
        }

        block_mappings.emplace_back(target_block_index, best_source_block_index);
    }

    // Now we need to learn the transformations
    // That is, for each best mapping pair (i, j), we calculate transformation from block i to block j.
    // todo: rotation
    vector<transformation_t> transformations;
    for (auto mapping: block_mappings) {
        int source_block_index, target_block_index;
        std::tie(target_block_index, source_block_index) = mapping;
        auto source_block = source_blocks[source_block_index];
        auto target_block = target_blocks[target_block_index];
        double contrast, brightness;

        const auto scaled_source_block = compress_block(image, source_block, target_block.width,
                                                        target_block.height);
        std::tie(brightness, contrast) = compute_brightness_and_contrast(image, scaled_source_block,
                                                                         target_block);

        int target_block_x = target_block.rel_x;
        int target_block_y = target_block.rel_y;
        double scaling = (double) block_size_target /
                         (double) block_size_source; // we have squared matrices, so nothing special here

        transformation_t t(source_block, scaling, target_block_x, target_block_y, contrast, brightness);
        transformations.push_back(t);
    }

    return transformations;
}

void apply_transformation(image_t &image, const transformation_t &t) {
    assert(t.source_block.width == t.source_block.height);

    const int source_block_size = t.source_block.width;
    const int downsampled_block_size = ((double) source_block_size) * t.scaling;

    double *downsampled_source_block = compress_block(image, t.source_block, downsampled_block_size,
                                                      downsampled_block_size);
    for (int i = 0; i < downsampled_block_size; ++i) {
        for (int j = 0; j < downsampled_block_size; ++j) {
            image.data[(t.target_block_y + i) * image.size + t.target_block_x + j] = downsampled_source_block[
                                                                                             i *
                                                                                             downsampled_block_size +
                                                                                             j] * t.contrast +
                                                                                     t.brightness;
        }
    }

    free(downsampled_source_block);
}

void decompress(image_t &decompressed_image,
                const vector<transformation_t> &transformations,
                const int num_iterations) {
    for (int iter = 0; iter < num_iterations; ++iter) {
        for (const auto t: transformations) {
            apply_transformation(decompressed_image, t);
        }
    }
}

func_suite_t register_suite() {
    func_suite_t suite;
    suite.compress_func = &compress;
    suite.decompress_func = &decompress;
    return suite;
}
