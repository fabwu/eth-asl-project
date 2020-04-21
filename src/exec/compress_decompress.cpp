#include <cassert>
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

    // read image
    int width;
    int height;
    double *image_data = read_grayscale_file(filename, &height, &width);
    image_t image(image_data, width);

    assert(width == height);

    func_suite_t suite = register_suite();

    auto transformations = suite.compress_func(image);
    image_t decompressed_image(width, true);
    suite.decompress_func(decompressed_image, transformations, 3);

    print_grayscale_file(decompressed_image.data, height, width);
}
