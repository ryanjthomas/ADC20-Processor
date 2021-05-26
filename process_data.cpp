#include <iostream>

#include <time.h>
#include <vector>
#include <numeric>
#include <string>
#include <cmath>
#include "ImageProcess.h"
#include "IdleProcess.h"
#include "ADC20File.h"

#include "tclap/CmdLine.h"

static const int DEFAULT_COLUMNS=4200;
static const int DEFAULT_ROWS=2100;

using namespace std;
int main (int argc, char **argv) {
	//Process our command line arguments first
	unsigned int window=25;
	int pre_cut=2;
	int post_cut=10;
	int rows=DEFAULT_ROWS;
	int columns=DEFAULT_COLUMNS;
	string filename;
	string outputDir;
	bool inverted=false;
	bool expose=true;
	bool writeFiles=true;
	bool idleFile=false;
	bool getOffset=true;
	int logLevel=0;
	
	try {
		TCLAP::CmdLine cmd("Image processor for a file taken with the new 20-bit ADC.", ' ', PROCESSOR_VERSION);
		TCLAP::ValueArg<unsigned int> windowArg("w","window", "Integration window/PSD window (for idle data)", false, window, "integer", cmd);
		TCLAP::ValueArg<int> preArg("p", "precut", "Bins before summing well to cut", false, pre_cut, "integer", cmd);
		TCLAP::ValueArg<int> postArg("P", "postcut", "Bins after summing well to cut", false, post_cut, "integer", cmd);
		TCLAP::ValueArg<int> rowArg("r", "rows", "Number of rows to process", false, rows, "integer", cmd);
		TCLAP::ValueArg<int> columnArg("c", "columns", "Number of columns in the image", false, columns, "integer", cmd);
		TCLAP::ValueArg<string> outputDirArg("o", "outputdir", "Directory to output fits/root file to", false, "none", "directory", cmd);
		TCLAP::ValueArg<int> logLevelArg("l", "loglevel", "Level of logging to report (0=errors only, 4=debug)",false, logLevel, "integer", cmd);
		TCLAP::SwitchArg getOffsetArg("g", "get-offset", "Switch to write the offsets to a CSV file (default yes)", cmd, getOffset);
		TCLAP::SwitchArg writeFilesArg("s", "save", "Switch to save fits/root files (default yes)", cmd, writeFiles);
		TCLAP::SwitchArg invertedArg("i", "inverted", "Whether the video is inverted or not (default no)", cmd, inverted/*default*/);
		TCLAP::SwitchArg exposeArg("e", "expose", "If the image was taken with an exposure or not (default yes)", cmd, expose);
		TCLAP::SwitchArg idleFileArg("I", "idle", "Whether to process the data as an idle file or not (default no)", cmd, idleFile);
		TCLAP::UnlabeledValueArg<string> filenameArg("filename", "Filename", true, "test.bin", "filename", cmd);

		cmd.parse(argc, argv);
		filename = filenameArg.getValue();
		inverted=invertedArg.getValue();
		window = windowArg.getValue();
		pre_cut=preArg.getValue();
		post_cut=postArg.getValue();
		expose=exposeArg.getValue();
		rows=rowArg.getValue();
		columns=columnArg.getValue();
		writeFiles=writeFilesArg.getValue();
		outputDir=outputDirArg.getValue();
		logLevel=logLevelArg.getValue();
		idleFile=idleFileArg.getValue();
		getOffset=getOffsetArg.getValue();

	} catch (TCLAP::ArgException &e ) {
		std::cerr << "Error: " << e.error() << " for arg " << e.argId() <<std::endl;
	}
 	/*	
	if (argc < 2 ) return 0;
	if (argc > 2 ) sscanf(argv[2], "%d", &window);
	if (argc ==5 ) {
		sscanf(argv[3], "%d", &pre_cut);
		sscanf(argv[4], "%d", &post_cut);
	}
	string filename=argv[1];*/
	ADC20File file(filename.c_str());
	if (!idleFile) {
		ImageProcess processor(rows, columns, &file);
		processor.setInverted(inverted);
		processor.setExpose(expose);
		processor.setCuts(pre_cut, post_cut);
		processor.setLogLevel(logLevel);
		processor.setComputeOffsets(getOffset);
		processor.init();
		processor.readImage(window);
		//Filename stuff
		size_t dir_end=filename.find_last_of("/");
		size_t extension_start=filename.find_last_of(".");
		//The base file name, w/o extension
		std::string basename = filename.substr(dir_end+1, extension_start-dir_end-1);
		//Our directory name
		std::string directory= filename.substr(0, dir_end+1);
		if (outputDir != "none" or outputDir.empty()) {
			directory=outputDir + std::string("/");
		} else {
			directory=directory+"processed/";
		}
		std::string fitsfilename = directory+basename+std::string(".fits");
		std::string rootfilename = directory+basename+std::string(".root");
		if (getOffset) {
			std::string offsetsfilename=directory+basename+std::string("-offsets.csv");
			processor.writeOffsets(offsetsfilename.c_str());
		};
		if (writeFiles) {
			processor.writeImage(fitsfilename.c_str());
		}
		//		processor.computeStatistics(writeFiles, rootfilename.c_str());
	} else {
		if (window < 23) window=std::pow(2,window);
		else window=std::pow(2,16);
		std::cout << "Processing idle data with window size of " << window << std::endl;
		IdleProcess processor(&file);
		processor.computePSD(window);
		processor.plotPSD();
	};
};
