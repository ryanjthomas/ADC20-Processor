#ifndef IMAGEPROCESS_H
#define IMAGEPROCESS_H
#include <stdio.h>
#include "utils.h"
#include "ADC20File.h"
#include <vector>

#define PROCESSOR_VERSION ".2"

class ImageProcess {
 public:
	ImageProcess(int rows, ADC20File *file);
	ImageProcess(int rows, int columns, ADC20File *file);
	~ImageProcess();

	off_t findStart();
	off_t getStart();
	bool gotoStart();
	
	int getColumns();	
	int findPixelLength();
	
	void readImage(int integration_window);
	bool writeImage(const char *filename);
	bool writeOffsets(const char *filename);
	void computeStatistics(bool writeRoot=true, const char *filename=NULL);
	void init();

	void setInverted(bool invert);
	void setRows(int r);
	void setColumns(int cols);
	void setBinning(int xbin, int ybin);
	void setExpose(bool expose);
	void setCleverStart(bool clever);
	void setCuts(int pre, int post);
	void setComputeOffsets(bool compute);

	void setLogLevel(int level);

 private:

	std::vector<std::vector<double> > image;
	std::vector<double> flat_image;
	std::vector<uint32> buffer;

	
	ADC20File *file;

	//Parameters for the CDS
	int rows, columns;
	int integration_samps;
	int pre_cut, post_cut;

	int pixel_length;
	int reset_height;
	int reset_width;

	int transition_window_size;

	bool inverted;
	bool exposure;
	bool cleverStart;
	int binning[2];
	int video_polarity;

	off_t imageStart;
	bool startSet;
	bool computeOffsets;

	bool isClocking(const std::vector<uint32> &data);
	bool isClocking(const std::vector<uint32> &data, int start, int end);
	
	int findTransition(const std::vector<uint32> &vec);
	int findTransition2(const std::vector<uint32> &vec);

	//Not working
	uint32 findNumColumns();

	std::vector<double> readRow(int integration_length, int pre_cut, int post_cut);
	double computeNoise(const std::vector<double> &pixels);
	std::vector<double> computeNoiseIntegrationTime(const std::vector<uint32> &vec, int pre_cut, int post_cut, int step=1);
	std::vector<double> rowNoiseIntegration(int row, int step=1);
	std::vector<double> rowAverageValue();
	std::vector<double> rowNoise();
	std::vector<off_t> rowStarts;
	std::vector<std::vector<off_t>> pixelStarts; //Offsets for the start of each pixel (for the last row read)
};	


#endif
