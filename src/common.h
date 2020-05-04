#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <vector>
#include <cmath>

#include "filehandler.h"
#include "tsc_x86.h"
#include "performance.h"

#define VERIFY_MIN_PSNR 30.0
#define VERIFY_DECOMPRESS_ITERATIONS 10
#define WARMUP_CYCLES_REQUIRED 1e8
#define BENCHMARK_REPETITIONS 2

class params_t {
   public:
    std::string image_path;
    int block_size_range;
    int block_size_domain;
    int decompression_iterations;
    bool csv_output;
    std::string csv_output_path;

    params_t(std::string image_path, int block_size_range, int block_size_domain,
             int decompression_iterations, bool csv_output = false,
             std::string csv_output_path = std::string())
        : image_path(image_path),
          block_size_range(block_size_range),
          block_size_domain(block_size_domain),
          decompression_iterations(decompression_iterations),
          csv_output(csv_output),
          csv_output_path(csv_output_path) {}
};

class image_t {
public:
    double *data;
    int size;

    image_t(double *data, int size) : data(data), size(size) {}

    image_t(int size, bool randomize_data) : size(size) {
        data = (double *) malloc(size * size * sizeof(double));
        if (randomize_data) {
            for (int i = 0; i < size * size; ++i) data[i] = rand() % 256;
        }
    }

    ~image_t() {
        free(data);
    }

    double &operator[](int index) const {
        return data[index];
    }
};

class block_t {
public:
    int rel_x, rel_y;
    int width, height;

    block_t() {}

    block_t(int x, int y, int width, int height)
            : rel_x(x), rel_y(y), width(width), height(height) {}

    void print_block(const image_t &image, const int image_size) const {
        printf("Block rel_x: %d, rel_y: %d, width: %d, height: %d\n", rel_x, rel_y,
               width, height);
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                int index = get_index_in_image(i, j, image);
                printf("%.1f, ", image[index]);
            }
            printf("\n");
        }
    }

    int get_index_in_image(const int y_rel_block, const int x_rel_block,
                           const image_t &full_image) const {
        return (rel_y + y_rel_block) * full_image.size + rel_x + x_rel_block;
    }

    std::vector<block_t> quad() const {
        assert(width % 2 == 0);
        assert(height % 2 == 0);
        assert(width >= 2);
        assert(height >= 2);

        const int quad_width = width / 2;
        const int quad_height = height / 2;

        std::vector<block_t> quad_blocks;
        quad_blocks.push_back(block_t(rel_x, rel_y, quad_width, quad_height));
        quad_blocks.push_back(block_t(rel_x + quad_width, rel_y, quad_width, quad_height));
        quad_blocks.push_back(block_t(rel_x, rel_y + quad_height, quad_width, quad_height));
        quad_blocks.push_back(block_t(rel_x + quad_width, rel_y + quad_height, quad_width, quad_height));

        return quad_blocks;
    }

};

struct transformation_t {
    block_t domain_block, range_block;
    double contrast, brightness;
    int angle;
};

typedef std::vector<transformation_t>
(*compress_func_type)(const image_t &image,
                      const int block_size_range,
                      const int block_size_domain);

typedef void (*decompress_func_type)(
        image_t &image, const std::vector<transformation_t> &transformations,
        const int iterations);

struct func_suite_t {
    compress_func_type compress_func;
    decompress_func_type decompress_func;
};

func_suite_t register_suite();

inline double mean_squared_error(const image_t &original, const image_t &converted) {
    assert(original.size == converted.size);

    double squared_error = 0.0;
    for (int i = 0; i < original.size * original.size; ++i) {
        double diff = original[i] - converted[i];
        squared_error += diff * diff;
    }

    return squared_error / (original.size * original.size);
}

inline double psnr(const double mse) {
    // see https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio
    //
    // for 8 bits a typical value is between 30 and 50 dB (higher is better)
    double maxPixelValue = 255;
    return 20*std::log10(maxPixelValue) - 10*std::log10(mse);
}

inline double verify_compress_decompress_error(const image_t &image,
                                               const int block_size_range,
                                               const int block_size_domain,
                                               const func_suite_t &suite) {
    auto transformations = suite.compress_func(image,
                                               block_size_range,
                                               block_size_domain);
    image_t decompressed_image(image.size, true);
    suite.decompress_func(decompressed_image, transformations,
                          VERIFY_DECOMPRESS_ITERATIONS);

    double mse = mean_squared_error(image, decompressed_image);

    return psnr(mse);
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
    const image_t &image;
    const int block_size_range;
    const int block_size_domain;
    const func_suite_t suite;

public:
    benchmark_compress_t(const image_t &image,
                         const int block_size_range,
                         const int block_size_domain,
                         const func_suite_t &suite)
            : image(image),
              block_size_range(block_size_range),
              block_size_domain(block_size_domain),
              suite(suite) {}

  void perform() const override { suite.compress_func(image,
                                                      block_size_range,
                                                      block_size_domain); }
};

class benchmark_decompress_t : public virtual benchmark_t {
private:
    const image_t &original_image;
    const func_suite_t suite;
    const size_t iterations;
    const std::vector<transformation_t> transformations;

public:
    benchmark_decompress_t(const image_t &original_image,
                           std::vector<transformation_t> transformations,
                           const int decompression_iterations,
                           const func_suite_t &suite)
            : original_image(original_image),
              suite(suite),
              iterations(decompression_iterations),
              transformations(std::move(transformations)) {}

    void perform() const override {
        image_t decompressed_image(
                original_image.size, true);
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
        for (long i = 0; i < num_runs; i++) {
            benchmark.perform();
        }
        end = stop_tsc(start);

        auto cycles = (double) end;
        multiplier = (WARMUP_CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);
    return num_runs;
}

template <typename T>
inline double median(std::vector<T> &vec) {
    sort(vec.begin(), vec.end());
    auto size = vec.size();
    if (size % 2 == 0) {
        size_t i = size / 2 - 1;
        size_t j = i + 1;
        return (double)(vec[i] + vec[j]) / 2;
    } else {
        return (double)vec[size / 2];
    }
}

inline void benchmark_generic(const benchmark_t &benchmark, bool csv_output,
                              const std::string &csv_output_path) {
    std::cout << "\033[1m"
              << "WARMUP phase"
              << "\033[0m" << std::endl;
    std::cout << "\t"
              << "performing function for at least " << WARMUP_CYCLES_REQUIRED
              << " cycles" << std::endl;
    auto needed_runs = warmup(benchmark);
    std::cout << "\t" << needed_runs << " runs needed" << std::endl;

    std::cout << "\033[1m"
              << "BENCHMARK phase"
              << "\033[0m" << std::endl;
    std::vector<double> cycles;
    std::vector<long long> flops;
    myInt64 start, end;
    for (size_t rep = 0; rep < BENCHMARK_REPETITIONS; ++rep) {
        __reset_flop_counter();
        start = start_tsc();
        for (long run = 0; run < needed_runs; ++run) {
            benchmark.perform();
        }
        end = stop_tsc(start);
        double cycles_run = ((double)end) / needed_runs;
        cycles.push_back(cycles_run);
        flops.push_back(nbr_double_flops / needed_runs);
    }

    // CSV output has to be generated here before cycles gets sorted in median
    if (csv_output) {
        output_csv(cycles, flops, csv_output_path);
        std::cout << "\t"
                  << "Created csv file '" << csv_output_path << "'" << std::endl;
    }

    auto median_cycles = median(cycles);
    std::cout << "\t"
              << "cycles (median): " << median_cycles << std::endl;

#if ENABLE_PERF_COUNTER
    auto median_flops = median(flops);
    std::cout << "\t"
              << "flops: " << median_flops << std::endl;
    std::cout << "\t"
              << "perf [flops/cycle(median)]: "
              << (double)median_flops / median_cycles << std::endl;
#endif
}

inline bool verify_suite(const func_suite_t &suite,
                         const int block_size_range,
                         const int block_size_domain,
                         const image_t &image) {
    std::cout << "\033[1m"
              << "VERIFICATION phase"
              << "\033[0m" << std::endl;
    double psnr = verify_compress_decompress_error(image, block_size_range, block_size_domain, suite);
    bool verification_failed = psnr < VERIFY_MIN_PSNR;
    std::cout << "\t"
              << "PSNR (dB): " << psnr << " (at least " << VERIFY_MIN_PSNR
              << " dB / higher is better)" << std::endl;
    if (verification_failed) {
        std::cout << "\t"
                  << "\033[1;31m"
                  << "✗ verification failed"
                  << "\033[0m" << std::endl;
    } else {
        std::cout << "\t"
                  << "\033[1;32m"
                  << "✔ verification succeeded"
                  << "\033[0m" << std::endl;
    }
    return !verification_failed;
}

inline void benchmark_compress(const params_t &params) {
    int width, height;
    double *original_image_data =
        read_grayscale_file(params.image_path, &height, &width);
    const image_t original_image(original_image_data, width);

    const auto suite = register_suite();
    verify_suite(suite, params.block_size_range, params.block_size_domain,
                 original_image);

    const benchmark_compress_t benchmark(original_image,
                                         params.block_size_range,
                                         params.block_size_domain, suite);

    benchmark_generic(benchmark, params.csv_output, params.csv_output_path);
}

inline void benchmark_decompress(const params_t &params) {
    int width, height;
    double *original_image_data =
        read_grayscale_file(params.image_path, &height, &width);
    const image_t original_image(original_image_data, width);

    const auto suite = register_suite();
    if (!verify_suite(suite, params.block_size_range, params.block_size_domain,
                      original_image))
        return;

    auto transformations = suite.compress_func(
        original_image, params.block_size_range, params.block_size_domain);
    const benchmark_decompress_t benchmark(original_image, transformations,
                                           params.decompression_iterations,
                                           suite);

    benchmark_generic(benchmark, params.csv_output, params.csv_output_path);
}

inline void compress_decompress(const params_t &params) {
    int width, height;
    double *original_image_data =
        read_grayscale_file(params.image_path, &height, &width);
    const image_t image(original_image_data, width);

    const auto suite = register_suite();
    // if (!verify_suite(suite, image)) return;

    auto transformations = suite.compress_func(image, params.block_size_range,
                                               params.block_size_domain);
    image_t decompressed_image(width, true);
    suite.decompress_func(decompressed_image, transformations,
                          params.decompression_iterations);
    print_grayscale_file(decompressed_image.data, height, width);
}

#endif  // COMMON_H
