#ifndef COMMON_H
#define COMMON_H

#include "filehandler.h"
#include "tsc_x86.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>


#define VERIFY_ALLOWED_SQUARED_ERROR_PER_PIXEL 6.0
#define VERIFY_DECOMPRESS_ITERATIONS 10
#define WARMUP_CYCLES_REQUIRED 1e8
#define BENCHMARK_REPETITIONS 50

struct image_t {
    double *data;
    int size;

    image_t(double *data, int size) : data(data), size(size) {}

    explicit image_t(int size) : size(size) {
        data = (double *) calloc(size * size, sizeof(double));
    }
};

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

struct transformation_t {
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
};

typedef std::vector<transformation_t>(*compress_func_type)(const image_t &image);

typedef void(*decompress_func_type)(image_t &image, const std::vector<transformation_t> &transformations,
                                    const int iterations);

struct func_suite_t {
    compress_func_type compress_func;
    decompress_func_type decompress_func;
};

func_suite_t register_suite();

inline double squared_error(const image_t &original, const image_t &converted) {
    assert(original.size == converted.size);

    double squared_error = 0.0;
    for (int i = 0; i < original.size; ++i) {
        double diff = original.data[i] - converted.data[i];
        squared_error += diff * diff;
    }

    return squared_error;
}

inline double verify_compress_decompress_error(const image_t &image, const func_suite_t &suite) {
    auto transformations = suite.compress_func(image);
    image_t decompressed_image(image.size);
    suite.decompress_func(decompressed_image, transformations, VERIFY_DECOMPRESS_ITERATIONS);

    return squared_error(image, decompressed_image);
}

/**
 * A function to be timed by some benchmark
 */
class benchmark_t {
public:
    virtual void perform() const = 0;
};

class benchmark_compress_t : public virtual benchmark_t {
private:
    const image_t image;
    const func_suite_t suite;
public:
    benchmark_compress_t(const image_t &image, const func_suite_t &suite) : image(image), suite(suite) {}

    void perform() const override {
        suite.compress_func(image);
    }
};

class benchmark_decompress_t : public virtual benchmark_t {
private:
    const image_t original_image;
    const func_suite_t suite;
    const std::vector<transformation_t> transformations;
    const size_t iterations;
public:
    benchmark_decompress_t(const image_t original_image,
                           std::vector<transformation_t> transformations,
                           const int decompression_iterations,
                           const func_suite_t &suite
    ) :
            original_image(original_image),
            suite(suite),
            iterations(decompression_iterations),
            transformations(std::move(transformations)) {}

    void perform() const override {
        image_t decompressed_image(original_image.size); // TODO: Always create this new?
        suite.decompress_func(decompressed_image, transformations, iterations);
    }

};

/**
 * Copied from Code Expert
 * Returns the number of repetitions
 * @param image
 * @param suite
 */
inline long warmup(const benchmark_t &benchmark) {
    long num_runs = 2;
    double multiplier = 1;
    myInt64 start, end;
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = (long) ((double) num_runs * multiplier);
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            benchmark.perform();
        }
        end = stop_tsc(start);

        auto cycles = (double) end;
        multiplier = (WARMUP_CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
    return num_runs;
}


inline double median(std::vector<double> &vec) {
    sort(vec.begin(), vec.end());
    auto size = vec.size();
    if (size % 2 == 0) {
        size_t i = size / 2 - 1;
        size_t j = i + 1;
        return (vec[i] + vec[j]) / 2;
    } else {
        return vec[size / 2];
    }
}

inline void benchmark_generic(const benchmark_t &benchmark) {
    std::cout << "\033[1m" << "WARMUP phase" << "\033[0m" << std::endl;
    std::cout << "\t" << "performing function for at least " << WARMUP_CYCLES_REQUIRED << " cycles" << std::endl;
    auto needed_runs = warmup(benchmark);
    std::cout << "\t" << needed_runs << " runs needed" << std::endl;


    std::cout << "\033[1m" << "BENCHMARK phase" << "\033[0m" << std::endl;
    std::vector<double> cycles;
    myInt64 start, end;
    for (size_t rep = 0; rep < BENCHMARK_REPETITIONS; ++rep) {

        start = start_tsc();
        for (size_t run = 0; run < needed_runs; ++run) {
            benchmark.perform();
        }
        end = stop_tsc(start);
        double cycles_run = ((double) end) / needed_runs;
        cycles.push_back(cycles_run);
    }

    auto median_cycles = median(cycles);
    std::cout << "\t" << "cycles (median): " << median_cycles << std::endl;
}

inline bool verify_suite(const func_suite_t &suite, const image_t &image) {
    std::cout << "\033[1m" << "VERIFICATION phase" << "\033[0m" << std::endl;
    double verification_error = verify_compress_decompress_error(image, suite);
    double verification_error_per_pixel = verification_error / image.size / image.size;
    bool verification_failed = verification_error_per_pixel > VERIFY_ALLOWED_SQUARED_ERROR_PER_PIXEL;
    std::cout << "\t" << "error² (Ø p.p): " << verification_error_per_pixel << " (allowed is up to "
              << VERIFY_ALLOWED_SQUARED_ERROR_PER_PIXEL << ")" << std::endl;
    if (verification_failed) {
        std::cout << "\t" << "\033[1;31m" << "✗ verification failed" << "\033[0m" << std::endl;
    } else {
        std::cout << "\t" << "\033[1;32m" << "✔ verification succeeded" << "\033[0m" << std::endl;
    }
    return !verification_failed;
}

inline void benchmark_compress(const std::string &image_path) {
    int width, height;
    double *original_image_data = read_grayscale_file(image_path, &height, &width);
    const image_t original_image(original_image_data, width);

    const auto suite = register_suite();
    if (!verify_suite(suite, original_image)) return;

    const benchmark_compress_t benchmark(original_image, suite);

    benchmark_generic(benchmark);
}

inline void benchmark_decompress(const std::string &image_path, const int decompression_iterations) {
    int width, height;
    double *original_image_data = read_grayscale_file(image_path, &height, &width);
    const image_t original_image(original_image_data, width);

    const auto suite = register_suite();
    if (!verify_suite(suite, original_image)) return;

    auto transformations = suite.compress_func(original_image);
    const benchmark_decompress_t benchmark(original_image, transformations, decompression_iterations, suite);

    benchmark_generic(benchmark);
}

#endif //COMMON_H
