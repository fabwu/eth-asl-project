/*
 * Fractal Image Compression. Copyright 2004 Alex Kennberg.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
//#include <windows.h>
using namespace std;

#include "Image.h"
#include "IFSTransform.h"
#include "Encoder.h"
#include "QuadTreeEncoder.h"
#include "Decoder.h"
#include "tsc_x86.h"
#include "performance.h"

int verb = 0;
bool useYCbCr = true;

void printUsage(char *exe);

void Convert(Encoder* enc, Image* source, int maxphases, int output)
{
	printf("Loading image...\n");
	source->Load();

	int width = source->GetWidth();
	int height = source->GetHeight();
	int imagesize = source->GetOriginalSize();

	printf("Encoding...\n");
   
    std::vector<std::pair<myInt64, myInt64>> cycles_flops; 
    myInt64 start, end;
    Transforms* transforms;
    
    // warm-up
    enc->Encode(source);
    
    for (size_t rep = 0; rep < 3; ++rep) {
        __reset_flop_counter();
        start = start_tsc();
        __record_double_flops(1);
	    transforms = enc->Encode(source);
        end = stop_tsc(start);
        cycles_flops.emplace_back(end, __nbr_double_flops());
    }
    
    printf("benchmark (cycles, flops, perf):\n");
    for(std::pair<myInt64, myInt64> p : cycles_flops) {
        printf("%lld;%lld;%f\n", p.first, p.second, p.second / (double) p.first);
    }
    printf("\n");

	int numTransforms = transforms->ch[0].size() +
		transforms->ch[1].size() + transforms->ch[2].size();

	printf("Number of transforms: %d\n", numTransforms);
	printf("Raw image bytes per transform: %d\n", imagesize/numTransforms);

	int transformSize = numTransforms * sizeof(IFSTransform);
	int bogocompressionratio = imagesize/(transformSize/sizeof(int));
	printf("Compression Ratio: %d:1\n", bogocompressionratio);

	printf("Decoding...\n");
	Decoder* dec = new Decoder(width, height);

	for (int phase = 1; phase <= maxphases; phase++)
	{
		dec->Decode(transforms);

		// Save all channels (note: channel 0 means all channels).
		if (output >= 2)
		{
			for (int ch = 0; ch <= transforms->channels; ch++)
			{
				string outName("output");
				outName += static_cast<char> ('0' + phase);
				outName += static_cast<char> ('0' + ch);
				outName += ".raw";

				// Generate output image per each phase
				Image* producer = dec->GetNewImage(outName, ch);
				producer->Save();
				delete producer;

				if (output == 2)
					break;
			}
		}
	}

	// Save the final image.
	if (output == 1)
	{
		Image* producer = dec->GetNewImage("output.raw", 0);
		producer->Save();
		delete producer;
	}

	delete dec;
	delete transforms;

	printf("Finished.\n");
}


int main(int argc, char **argv)
{
	QuadTreeEncoder* enc;
	Image *source;
	string fileName;
	int threshhold = 100;
	bool symmetry = false;
	int phases = 5;
	int output = 1;
	bool usage = true;

	// Load parameters
	for (int i = 1; i < argc && usage; i++)
	{
		string param(argv[i]);

		if (param == "-v" && i + 1 < argc)
			verb = atoi(argv[i + 1]);
		else if (param == "-t" && i + 1 < argc)
			threshhold = atoi(argv[i + 1]);
		else if (param == "-p" && i + 1 < argc)
			phases = atoi(argv[i + 1]);
		else if (param == "-o" && i + 1 < argc)
			output = atoi(argv[i + 1]);
		else if (param == "-f" && --i >= 0)
			symmetry = true;
		else if (param == "-r" && --i >= 0)
			useYCbCr = false;

		if (param.at(0) == '-')
		{
			i++;
		}
		else
		{
			fileName = param;
			usage = false;
		}
	}

	if (usage)
	{
		printUsage(argv[0]);
		return -1;
	}

	source = new Image(fileName);
	enc = new QuadTreeEncoder(threshhold, symmetry);

	Convert(enc, source, phases, output);

	delete enc;
	delete source;
}

void printUsage(char *exe)
{
	printf("Usage: %s [-v #] [-t #] [-p #] [-o #] [-f] [-r] filename\n"
		"\t-v 0    Verbous level (0-4)\n"
		"\t-t 100  Threshold (i.e. quality)\n"
		"\t-p 5    Number of decoding phases\n"
		"\t-o 1    1:Output final image,\n"
		"\t        2:Output at each phase,\n"
		"\t        3:Output at each phase & channel\n"
		"\t-f      Force symmetry operations during encoding\n"
		"\t-r      Enable RGB instead of YCbCr\n",
		exe
	);
}
