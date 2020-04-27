#include <cassert>
#include <string>
#include "../common.h"

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
        ("sizes,s",
         value<vector<int> >()->multitoken()->default_value(vector<int>{8,16}, "8 16"),
         "size of range and domain blocks")
        ;


    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);

    // show help
    if(vm.count("help")) {
        cout << desc << "\n";
        return 0;
    };

    // benchmark compress
    if(vm.count("benchmark") &&
       vm.count("compress") &&
       vm.count("iterations")){
        benchmark_compress(vm["filename"].as<string>(),
                           vm["sizes"].as<vector<int> >()[0],
                           vm["sizes"].as<vector<int> >()[1]);
        return 0;
    }

    // benchmark decompress
    if(vm.count("benchmark") &&
       vm.count("decompress") &&
       vm.count("iterations")){
        benchmark_decompress(vm["filename"].as<string>(),
                             vm["sizes"].as<vector<int> >()[0],
                             vm["sizes"].as<vector<int> >()[1],
                             vm["iterations"].as<int>());
        return 0;
    }

    // compress and decompress
    if(vm.count("compress") &&
       vm.count("decompress") &&
       vm.count("filename") &&
       vm.count("iterations")){
        compress_decompress(vm["filename"].as<string>(),
                            vm["sizes"].as<vector<int> >()[0],
                            vm["sizes"].as<vector<int> >()[1],
                            vm["iterations"].as<int>());
        return 0;
    }

    // show help if not yet done
    cout << desc << "\n";

}
