#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include <string>

void print_grayscale_file(double* grayscale_image, int height, int width);
double* read_grayscale_file(std::string filename, int* height, int* width);
void write_fic_file(std::string filename, double* fic_image);

#endif //FILEHANDLER_H
