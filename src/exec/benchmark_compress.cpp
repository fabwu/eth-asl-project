#include <string>

#include "../common.h"

using namespace std;

int main(int argc, char const *argv[]) {
  ios_base::sync_with_stdio(false);

  // argument processing
  if (argc != 2) {
    printf("\nSpecify an image and nothing else\n");
    return 1;
  }
  string filename = argv[1];

  benchmark_compress(filename);
}
