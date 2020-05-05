#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include <string>
#include <vector>

using namespace std;

void print_grayscale_file(double *grayscale_image, int height, int width);

double *read_grayscale_file(const std::string &filename, int *height, int *width);

void write_fic_file(const std::string &filename, double *fic_image);

void output_csv(const vector<double> &cycles, const vector<long long> &flops, const string &csv_output_path);

#endif  // FILEHANDLER_H
