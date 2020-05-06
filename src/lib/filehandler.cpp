#include "filehandler.hpp"
extern "C" {
#include "performance.h"
}

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cassert>

using namespace std;

/*
  Prints the GRAYSCALE_IMAGE as csv
*/
void print_grayscale_file(double *grayscale_image, int height, int width) {
    cout << setprecision(6) << std::fixed;
    cout << height << ";" << width << ";" << endl;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            cout << grayscale_image[i * width + j] << ";";
        }
        cout << endl;
    }
}

/*
  Reads file with FILENAME, sets height and width and returns the image.
*/
double *read_grayscale_file(const string &filename, int *height, int *width) {
    ifstream fin;
    fin.open(filename);

    if (!fin.is_open()) {
        throw std::runtime_error("Couldn't read grayscale file");
    }

    std::string line, word;
    int ctr = 0;

    // read image dimension
    fin >> line;
    stringstream s(line);
    if (getline(s, word, ';')) {
        const char *x = word.c_str();
        *height = stoi(x);
    }
    if (getline(s, word, ';')) {
        const char *x = word.c_str();
        *width = stoi(x);
    }

    // read image data
    auto grayscale_image =
            (double *) malloc((*width) * (*height) * sizeof(double));
    while (fin >> line) {
        stringstream s(line);
        while (getline(s, word, ';')) {
            const char *x = word.c_str();
            grayscale_image[ctr] = stod(x);
            ++ctr;
        }
    }
    return grayscale_image;
}

/* Write fic-file. */
void write_fic_file(string &filename, double *fic_image) {
    cout << "writing fic file (to be implemented)" << endl;
}

void output_csv(const vector<double> &cycles,
                       const vector<long long> &flops,
                       const string &csv_output_path) {
#ifdef ENABLE_PERF_COUNTER
    assert(cycles.size() == flops.size());
#endif
    ofstream fout;
    fout.open(csv_output_path);

    fout << "cycles";
#ifdef ENABLE_PERF_COUNTER
    fout << ";flops;flops/cycle";
#endif
    fout << endl;

    fout << std::fixed;
    for (size_t i = 0; i < cycles.size(); ++i) {
        fout << cycles[i];
#ifdef ENABLE_PERF_COUNTER
        fout << ';' << flops[i] << ';' << (double)flops[i] / cycles[i];
#endif
        fout << endl;
    }

    fout.close();
}
