import argparse
import configparser
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import sys

GRAY_DELIM = ";"


class Fic(object):
    @staticmethod
    def __replace_file_extension(file, new_extension):
        return os.path.splitext(file)[0] + new_extension

    @staticmethod
    def __rgb2gray(rgb):
        # uses PAL/NTSC conversion
        return np.dot(rgb[..., :3], [299, 587, 114]) / 1000

    @staticmethod
    def __parse_image(data):
        lines = data.splitlines()

        dimensions = lines.pop(0).split(GRAY_DELIM)
        height = int(dimensions[0])
        width = int(dimensions[1])

        img = np.zeros((height, width))
        for i, line in enumerate(lines):
            line_splitted = line.split(GRAY_DELIM)
            line_splitted.pop()
            assert len(line_splitted) == width
            for j in range(width):
                img[i][j] = line_splitted[j]

        return img

    @staticmethod
    def __read_gray(file):
        f = open(file, "r")
        return Fic.__parse_image(f.read())

    def togray(self):
        rgb = mpimg.imread(self.file)

        if rgb.ndim == 2:
            gray = rgb
            height = len(gray)
            width = len(gray[0])
        else:
            gray = self.__rgb2gray(rgb)
            height = gray.shape[0]
            width = gray.shape[1]

        out_file = self.__replace_file_extension(self.file, ".gray")
        f = open(out_file, "w+")
        f.write("%d%s%d%s" % (height, GRAY_DELIM, width, os.linesep))
        for row in gray:
            for val in row:
                f.write("%f%s" % (val, GRAY_DELIM))
            f.write(os.linesep)
        f.close()

    def topng(self):
        img = self.__read_gray(self.file)
        file_decoded = self.__replace_file_extension(self.file, ".png")
        mpimg.imsave(file_decoded, img, cmap="gray")

    def showgray(self):
        img = self.__read_gray(self.file)
        plt.figure(num=self.file)
        plt.imshow(img, cmap="gray")
        plt.show()

    def __call_fic(self, error_threshold, iterations):
        args = "build/fic -c -d -e {} -i {} -f {}".format(error_threshold, iterations, self.file)
        print("call: {}".format(args))
        p = subprocess.run(args.split(), capture_output=True, text=True)
        return self.__parse_image(p.stdout)

    @staticmethod
    def __read_num(label, default):
        while True:
            try:
                data = input("{} (default: {}): ".format(label, default))
                if not data:
                    return default
                return str(int(data))
            except ValueError:
                print("Oops!  That was no valid number.  Try again...")

    @staticmethod
    def __psnr(orig, new):
        mse = np.mean(np.square(new - orig))
        return 20 * np.log10(255) - 10 * np.log10(mse)

    def seq(self):
        config = configparser.ConfigParser()
        # default values
        config.read_dict({'FIC': {
            'error': 100,
            'iterations': '1 2 3 6 10 20',
        }})
        config.read('fic.ini')
        conf = config['FIC']

        conf['error'] = self.__read_num('error threshold', conf['error'])

        data = input('iteration sequence [e.g. 1 2 6 10] (default: {}): '.format(conf['iterations']))
        if data:
            conf['iterations'] = data
        iterations = conf['iterations'].split()
        assert len(iterations) > 0, 'Empty sequence found. Provide some numbers...'
        assert len(iterations) <= 6, "Can't show more than 6 plots"

        with open('fic.ini', 'w') as configfile:
            config.write(configfile)

        orig_image = self.__read_gray(self.file)

        fig = plt.figure()
        fig.suptitle("FIC error threshold {}".format(conf['error']))
        for idx, it in enumerate(iterations, start=1):
            img = self.__call_fic(conf.getint('error'), it)
            error = self.__psnr(img, orig_image)
            ax = fig.add_subplot(2, 3, idx)
            ax.imshow(img, cmap="gray", interpolation='nearest')
            ax.set_title('{} Iterations, {:.2f} psnr'.format(it, error))
            ax.set_yticks([])
            ax.set_xticks([])
        plt.show()

    def __init__(self):
        parser = argparse.ArgumentParser(description="FIC tool")
        parser.add_argument("command", help="The command to use", choices=[
            'togray',
            'topng',
            'showgray',
            'seq',
        ])

        parser.add_argument("file", type=str)
        args = parser.parse_args(sys.argv[1:3])

        self.file = args.file

        # call function with same name
        getattr(self, args.command)()


if __name__ == '__main__':
    Fic()
