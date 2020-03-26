#include <string>
#include "filehandler.cpp"

using namespace std;


int main(int argc, char const *argv[]) {
  ios_base::sync_with_stdio(false);

  // image settings
  string filename = "./earth.csv";

  // read image
  int width;
  int height;
  double* grayscale_image = read_grayscale_file(filename, &height, &width);

  // write image again
  print_grayscale_file(grayscale_image, height, width);

}
