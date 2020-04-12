#include <string>
#include <cassert>
#include <vector>
#include <limits>

#include "filehandler.cpp"
#include "rotate.h"

using namespace std;

pair<double, double> solve_2x2_eq_system(double a11, double a12, double a21, double a22, double b1, double b2) {
    /// requires det(A) != 0
    double f1 = a22 == 0 ? 0.0 : a12 / a22;
    double f2 = a11 == 0 ? 0.0 : a21 / a11;

    double x = (b1 - f1 * b2) / (a11 - f1 * a21);
    double y = (b2-a21*x)/a22; //(b2 - f2 * b1) / (a22 - f2 * a11);

    return make_pair(x, y);
}


class block_t {
public:
    const int rel_x, rel_y;
    const int width, height;

    block_t(int x, int y, int width, int height) : rel_x(x), rel_y(y), width(width), height(height) {}

    void print_block(const double *image, const int image_size) const {
        printf("Block rel_x: %d, rel_y: %d, width: %d, height: %d\n", rel_x, rel_y, width, height);
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                int index = get_index_in_image(i, j, image_size);
                printf("%.1f, ", image[index]);
            }
            printf("\n");
        }
    }

    int get_index_in_image(int y_rel_block, int x_rel_block, int full_image_size) const {
        return (rel_y + y_rel_block) * full_image_size + rel_x + x_rel_block;
    }

};

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

double *compress_block(const double *grayscale, const int image_size, block_t source_block, int width, int height) {
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
                val += grayscale[b.get_index_in_image(i, j, image_size)];
            }
        }
        val /= n * n;

        scaled_source_block[scaled_index] = val;
        scaled_index++;
    }

    return scaled_source_block;
}

class transformation_t {
public:
    const block_t source_block;
    const double scaling;
    const double contrast, brightness;
    const int target_block_x, target_block_y;

    transformation_t(block_t source_block, double scaling, int target_block_x, int target_block_y, double contrast,
                     double brightness) :
            source_block(source_block),
            scaling(scaling),
            target_block_x(target_block_x),
            target_block_y(target_block_y), contrast(contrast), brightness(brightness) {}

    void apply(double *image, int image_size) const {
        assert(source_block.width == source_block.height);

        const int source_block_size = source_block.width;
        const int downsampled_block_size = ((double) source_block_size) * scaling;

        double *downsampled_source_block = compress_block(image, image_size, source_block, downsampled_block_size,
                                                          downsampled_block_size);
        for (int i = 0; i < downsampled_block_size; ++i) {
            for (int j = 0; j < downsampled_block_size; ++j) {
                image[(target_block_y + i) * image_size + target_block_x + j] = downsampled_source_block[
                                                                                        i * downsampled_block_size +
                                                                                        j] * contrast + brightness;
            }
        }

        free(downsampled_source_block);
    }
};


pair<double, double> compute_brightness_and_contrast_naive(const double *image,
                                                           const int image_size, double *domain,
                                                           block_t target_block) {
    double contrast = 0.75;
    int n = target_block.width;

    // find average brightness for given contrast
    double brightness = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double di = domain[i * n + j];
            double ri = image[target_block.get_index_in_image(i, j, image_size)];
            brightness += ri - contrast*di;
        }
    }
    brightness /= n*n;

    return make_pair(brightness, contrast);
}
    // contrast = 0.75
    // brightness = (np.sum(D - contrast*S)) / D.size
    // return contrast, brightness

pair<double, double> compute_brightness_and_contrast(const double *image, const int image_size, double *domain,
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
            double ri = image[target_block.get_index_in_image(i, j, image_size)];
            a12 += di;
            a21 += di;
            a22 += di * di;
            b1 += ri;
            b2 += di * ri;
        }
    }

    return solve_2x2_eq_system(a11, a12, a21, a22, b1, b2);
}

inline double diff_block(const double *grayscale, const int image_size, block_t source_block, block_t target_block) {
    // Need to compress source_block, such that size(source_block) == size(target_block) in order to compare difference!
    // That means that a square of n pixels need to be compressed to one pixel
    const int n = source_block.width / target_block.width;
    const auto scaled_source_block = compress_block(grayscale, image_size, source_block, target_block.width,
                                                    target_block.height);

    double brightness, contrast;
    std::tie(brightness, contrast) = compute_brightness_and_contrast(grayscale, image_size, scaled_source_block,
                                                                     target_block);

    double squared_error = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double target_val = grayscale[target_block.get_index_in_image(i, j, image_size)];
            double source_val = scaled_source_block[i * n + j] * contrast + brightness;
            double diff = source_val - target_val;
            squared_error += diff * diff;
        }
    }

    free(scaled_source_block);
    return squared_error;
}

vector<transformation_t> compress(double *grayscale, int image_size) {
    // Goal: Try to map blocks of size block_size_source to blocks of size block_size_target
    const int block_size_target = 8;
    const int block_size_source = 16;
    static_assert(block_size_source >= block_size_target);

    const auto target_blocks = create_squared_blocks(image_size, block_size_target);
    const auto source_blocks = create_squared_blocks(image_size, block_size_source);

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

            double target_block_error = diff_block(grayscale, image_size, source_block, target_block);
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

        const auto scaled_source_block = compress_block(grayscale, image_size, source_block, target_block.width,
                                                        target_block.height);
        std::tie(brightness, contrast) = compute_brightness_and_contrast(grayscale, image_size, scaled_source_block,
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

double *decompress(int image_size,
                   vector<transformation_t> &transformations,
                   const int num_iterations) {
    // Generate random starting image
    double *decompressed_image = (double *) calloc(image_size * image_size, sizeof(double));

    for (int iter = 0; iter < num_iterations; ++iter) {
        for (const auto t: transformations) {
            t.apply(decompressed_image, image_size);
        }
    }

    return decompressed_image;
}

int main(int argc, char const *argv[]) {
    ios_base::sync_with_stdio(false);

    // TODO: remove
    // auto res = solve_2x2_eq_system(1.0, 3.0, 0.0, 2.0, 2.0, 3.0);
    // printf("x: %f, y: %f", res.first, res.second);

    // TODO: remove
    // for(int i=1;i<=16;++i){
    //     for(int j=1;j<=16;++j){
    //         printf("%d.%d;", i, j);
    //     }
    //     printf("\n");
    // }

    // argument processing
    if(argc != 2){
        printf("\nSpecify an image and nothing else\n");
        return 1;
    }
    string filename = argv[1];

    // read image
    int width;
    int height;
    double *grayscale_image = read_grayscale_file(filename, &height, &width);

    assert(width == height);
    auto transformations = compress(grayscale_image, width);
    double *converted_grayscale_image = decompress(width, transformations, 3);
    // write image again
    print_grayscale_file(converted_grayscale_image, height, width);
}
