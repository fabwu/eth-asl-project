/**
 * \file Common C++ types and functions.
 */

#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "filehandler.hpp"
#include "tsc_x86.h"
extern "C" {
#include "performance.h"
#include "types.h"
}

#define VERIFY_MIN_PSNR 30.0
#define VERIFY_DECOMPRESS_ITERATIONS 10
#define WARMUP_CYCLES_REQUIRED 1e8

extern "C" struct func_suite_t register_suite();

class params_t {
   public:
    std::string image_path;
    int error_threshold;
    int decompression_iterations;
    int benchmark_repetitions;
    bool csv_output;
    std::string csv_output_path;

    params_t(std::string image_path, int error_threshold,
             int decompression_iterations, int benchmark_repetitions, bool csv_output = false,
             std::string csv_output_path = std::string())
        : image_path(image_path),
          error_threshold(error_threshold),
          decompression_iterations(decompression_iterations),
          benchmark_repetitions(benchmark_repetitions),
          csv_output(csv_output),
          csv_output_path(csv_output_path) {}
};

inline double mean_squared_error(const struct image_t &original,
                                 const struct image_t &converted) {
    assert(original.size == converted.size);

    double squared_error = 0.0;
    for (int i = 0; i < original.size * original.size; ++i) {
        double diff = original.data[i] - converted.data[i];
        squared_error += diff * diff;
    }

    return squared_error / (original.size * original.size);
}

inline double psnr(const double mse) {
    // see https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio
    //
    // for 8 bits a typical value is between 30 and 50 dB (higher is better)
    double maxPixelValue = 255;
    return 20 * std::log10(maxPixelValue) - 10 * std::log10(mse);
}

inline double verify_compress_decompress_error(const struct image_t &image,
                                               const int error_threshold,
                                               const func_suite_t &suite) {
    auto transformations = suite.compress_func(&image, error_threshold);
    struct image_t decompressed_image = make_image(image.size, true);
    suite.decompress_func(&decompressed_image, transformations,
                          VERIFY_DECOMPRESS_ITERATIONS);

    double mse = mean_squared_error(image, decompressed_image);
    free_image_data(&decompressed_image);

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
    const struct image_t &image;
    const int error_threshold;
    const func_suite_t suite;

   public:
    benchmark_compress_t(const struct image_t &image, const int error_threshold,
                         const func_suite_t &suite)
        : image(image), error_threshold(error_threshold), suite(suite) {}

    void perform() const override {
        auto transformations = suite.compress_func(&image, error_threshold);
        free_queue(transformations);
        free(transformations);
    }
};

class benchmark_decompress_t : public virtual benchmark_t {
   private:
    const struct image_t &original_image;
    const func_suite_t suite;
    const size_t iterations;
    const struct queue *transformations;

   public:
    benchmark_decompress_t(const struct image_t &original_image,
                           const struct queue *transformations,
                           const int decompression_iterations,
                           const func_suite_t &suite)
        : original_image(original_image),
          suite(suite),
          iterations(decompression_iterations),
          transformations(transformations) {}

    void perform() const override {
        struct image_t decompressed_image =
            make_image(original_image.size, true);
        suite.decompress_func(&decompressed_image, transformations, iterations);
        free_image_data(&decompressed_image);
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
        num_runs = (long)((double)num_runs * multiplier);
        start = start_tsc();
        for (long i = 0; i < num_runs; i++) {
            benchmark.perform();
        }
        end = stop_tsc(start);

        auto cycles = (double)end;
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

inline void benchmark_generic(const benchmark_t &benchmark,
                              int benchmark_repetitions, bool csv_output,
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
    for (int rep = 0; rep < benchmark_repetitions; ++rep) {
        __reset_flop_counter();
        start = start_tsc();
        for (long run = 0; run < needed_runs; ++run) {
            benchmark.perform();
        }
        end = stop_tsc(start);
        double cycles_run = ((double)end) / needed_runs;
        cycles.push_back(cycles_run);
        flops.push_back(__nbr_double_flops() / needed_runs);
    }

    // CSV output has to be generated here before cycles gets sorted in median
    if (csv_output) {
        output_csv(cycles, flops, csv_output_path);
        std::cout << "\t"
                  << "Created csv file '" << csv_output_path << "'"
                  << std::endl;
    }

    auto median_cycles = median(cycles);
    std::cout << "\t"
              << "cycles (median): " << median_cycles << std::endl;

#ifdef ENABLE_PERF_COUNTER
    auto median_flops = median(flops);
    std::cout << "\t"
              << "flops: " << median_flops << std::endl;
    std::cout << "\t"
              << "perf [flops/cycle(median)]: "
              << (double)median_flops / median_cycles << std::endl;
#endif
}

inline bool verify_suite(const func_suite_t &suite, const int error_threshold,
                         const image_t &image) {
    std::cout << "\033[1m"
              << "VERIFICATION phase"
              << "\033[0m" << std::endl;
    double psnr =
        verify_compress_decompress_error(image, error_threshold, suite);
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
    struct image_t original_image;
    original_image.data = original_image_data;
    original_image.size = width;

    const auto suite = register_suite();
    verify_suite(suite, params.error_threshold, original_image);

    const benchmark_compress_t benchmark(original_image, params.error_threshold,
                                         suite);


    benchmark_generic(benchmark, params.benchmark_repetitions,
                      params.csv_output, params.csv_output_path);
    free(original_image_data);
}

inline void benchmark_decompress(const params_t &params) {
    int width, height;
    double *original_image_data =
        read_grayscale_file(params.image_path, &height, &width);
    struct image_t original_image;
    original_image.data = original_image_data;
    original_image.size = width;

    const auto suite = register_suite();
    if (!verify_suite(suite, params.error_threshold, original_image)) return;

    auto transformations =
        suite.compress_func(&original_image, params.error_threshold);
    const benchmark_decompress_t benchmark(original_image, transformations,
                                           params.decompression_iterations,
                                           suite);

    benchmark_generic(benchmark, params.benchmark_repetitions,
                      params.csv_output, params.csv_output_path);
    free_queue(transformations);
    free(transformations);
    free(original_image_data);
}

inline void compress_decompress(const params_t &params) {
    int width, height;
    double *original_image_data =
        read_grayscale_file(params.image_path, &height, &width);
    struct image_t image;
    image.data = original_image_data;
    image.size = width;

    const auto suite = register_suite();
    // if (!verify_suite(suite, image)) return;

    auto transformations = suite.compress_func(&image, params.error_threshold);
    struct image_t decompressed_image = make_image(width, true);
    suite.decompress_func(&decompressed_image, transformations,
                          params.decompression_iterations);
    print_grayscale_file(decompressed_image.data, height, width);

    free_queue(transformations);
    free(transformations);
    free_image_data(&decompressed_image);
    free(original_image_data);
}

#endif  // COMMON_H
