#include <cassert>
#include <string>
#include "../lib/common.hpp"

#include <boost/program_options.hpp>
using namespace boost::program_options;

#include <iostream>
using namespace std;

int main(int argc, char const *argv[]) {
    ios_base::sync_with_stdio(false);

    options_description desc("Options");

    desc.add_options()
        ("help,h", "produce help message")
        ("benchmark,b", "benchmark")
        ("decompress,d", "decrompress image")
        ("compress,c", "crompress image")
        ("iterations,i",
         value<int>()->default_value(10),
         "number of decompression iterations")
        ("filename,f",
         value<string>(),
         "gray image file")
        ("error,e",
         value<int>()->default_value(100),
         "the error threshold, which may not be exceeded by any transformation")
        ("csv,o", value<string>(), "report output csv file")
        ;


    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);

    // show help
    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    };

    if (vm.count("filename") && vm.count("error")) {
        params_t params(
            vm["filename"].as<string>(),
            vm["error"].as<int>(), vm["iterations"].as<int>());

        if (vm.count("csv")) {
            params.csv_output = true;
            params.csv_output_path = vm["csv"].as<string>();
        }

        // benchmark compress
        if (vm.count("benchmark") && vm.count("compress")) {
            benchmark_compress(params);
            return 0;
        }

        // benchmark decompress
        if (vm.count("benchmark") && vm.count("decompress")) {
            benchmark_decompress(params);
            return 0;
        }

        // compress and decompress
        if (vm.count("compress") && vm.count("decompress")) {
            compress_decompress(params);
            return 0;
        }
    }

    // show help if not yet done
    cout << desc << "\n";

}
