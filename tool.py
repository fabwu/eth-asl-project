import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import argparse, sys, os
import numpy as np

GRAY_DELIM = ";"

def __replace_file_extension(file, new_extension):
    return os.path.splitext(file)[0] + new_extension


def __rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2126, 0.7152, 0.0722])


def __read_gray(file):
    f = open(file, "r")
    dimensions = f.readline().split(GRAY_DELIM)
    height = int(dimensions[0])
    width = int(dimensions[1])

    img = np.zeros((height, width))
    for i in range(height):
        line_splitted = f.readline().split(GRAY_DELIM)
        line_splitted.remove(os.linesep)
        assert len(line_splitted) == width
        for j in range(width):
            img[i][j] = line_splitted[j]

    return img


def img_to_gray(file):
    rgb = mpimg.imread(file)

    gray = __rgb2gray(rgb)
    height = gray.shape[0]
    width = gray.shape[1]

    out_file = __replace_file_extension(file, ".gray")
    f = open(out_file, "w+")
    f.write("%d%s%d%s" % (height, GRAY_DELIM, width, os.linesep))
    for row in gray:
        for val in row:
            f.write("%f%s" % (val, GRAY_DELIM))
        f.write(os.linesep)
    f.close()


def gray_to_png(file):
    img = __read_gray(file)
    file_decoded = __replace_file_extension(file, ".png")
    mpimg.imsave(file_decoded, img, cmap="gray")


def show_gray(file):
    img = __read_gray(file)
    plt.imshow(img, cmap="gray")
    plt.show()


if __name__ == "__main__":
    METHOD_CREATE_GRAYSCALE_IMAGE = "togray"
    METHOD_TO_PNG = "topng"
    METHOD_SHOW_GRAY = "showgray"

    parser = argparse.ArgumentParser(description="FIC tool")
    parser.add_argument("command", help="The command to use", choices=[METHOD_CREATE_GRAYSCALE_IMAGE, METHOD_TO_PNG, METHOD_SHOW_GRAY])
    parser.add_argument("--file", type=str, required=True)
    args = parser.parse_args(sys.argv[1:])

    if args.command == METHOD_CREATE_GRAYSCALE_IMAGE:
        img_to_gray(args.file)
    elif args.command == METHOD_TO_PNG:
        gray_to_png(args.file)
    elif args.command == METHOD_SHOW_GRAY:
        show_gray(args.file)