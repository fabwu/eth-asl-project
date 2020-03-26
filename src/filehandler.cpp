#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;


/*
  Prints the GRAYSCALE_IMAGE as csv
*/
void print_grayscale_file(double* grayscale_image, int height, int width){
  cout << setprecision(6) << std::fixed;
  cout << height << ";" << width << ";" << endl;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      cout << grayscale_image[i*width+j] << ";";
    }
    cout << endl;
  }
}

/*
  Reads file with FILENAME, sets height and width and returns the image.
*/
double* read_grayscale_file(string filename, int* height, int* width){
  ifstream fin;
  fin.open(filename);
  std::string line, word;
  int ctr = 0;

  // read image dimension
  fin >> line;
  stringstream s(line);
  if(getline(s, word, ';')){
    const char* x = word.c_str();
    height[0] = stoi(x);
  }
  if(getline(s, word, ';')){
    const char* x = word.c_str();
    width[0] = stoi(x);
  }

  // read image data
  double* grayscale_image = (double*) malloc(width[0]*height[0]*sizeof(double));
  while(fin >> line){
    stringstream s(line);
    while(getline(s, word, ';')){
      const char* x = word.c_str();
      grayscale_image[ctr] = stod(x);
      ++ctr;
    }
  }
  return grayscale_image;
}



/* Write fic-file. */
void write_fic_file(string filename, double* fic_image){
  cout << "writing fic file (to be implemented)" << endl;
}
