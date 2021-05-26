#include <numeric>
#include <algorithm>


#include <cmath>
#include <cstdlib>
#include "ImageProcess.h"
#include "ImageProcessDefs.h"
#include "vec_funcs.h"

//Root crap
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TROOT.h"


#include "stdio.h"
#include "fitsio.h"

#include <iostream>
#include <fstream>

#include "log.h"


//These are various helper functions for our ImageProcess class implementation
//TODO: put these into their own file
template <typename T>
std::vector<int> findPeaks(std::vector<T> &vec, const int peakWidth, const T peakThreshold) {
	std::vector<int> peaks;
	for (int i=peakWidth; i < vec.size()-peakWidth; i++ ) {
		if ((abs(vec[i]-vec[i+peakWidth])) < peakThreshold) {
			std::vector<T> peak_win(&vec[i-peakWidth], &vec[i+peakWidth]);
			if (vec[i]-vec[i+peakWidth] > 0) {
				peaks.push_back(getVectorArgMax<T>(peak_win)+(i-peakWidth));
			} else{
				peaks.push_back(getVectorArgMin<T>(peak_win)+(i-peakWidth));
			}
			//Assumes very well defined peaks
			i += peakWidth*2;
		};
	};
	return peaks;
};


template <typename T>
std::vector<int> findPeaks2(std::vector<T> &vec, const int peakWidth, const T peakThreshold, const char type='a') {
	std::vector<int> peaks;
	int winSize=2*peakWidth;
	typename std::vector<T>::iterator winStart=vec.begin();
	typename std::vector<T>::iterator winEnd=vec.begin()+winSize;
	int max, mean, min;
	bool findPositive = (type=='a' or type=='p');
	bool findNegative = (type=='a' or type=='n');
	
	for (int i=peakWidth; i < vec.size()-peakWidth; i++, winStart++, winEnd++){		
		max=*std::max_element(winStart, winEnd);
		min=*std::min_element(winStart, winEnd);
		mean=std::accumulate(winStart, winEnd, 0.0)/winSize;
		if (max > mean+peakThreshold) {
			if (findPositive) {peaks.push_back(std::distance(winStart, std::max_element(winStart, winEnd))+i);}
			i+=winSize;
			winStart+=winSize;
			winEnd+=winSize;
		} else if ( min < mean-peakThreshold ){
			if (findNegative) {peaks.push_back(std::distance(winStart, std::min_element(winStart, winEnd))+i);}
			i+=winSize;
			winStart+=winSize;
			winEnd+=winSize;
		};
	};
	return peaks;
};



//Finds the type of peaks given certain peak conditions
//True is a positive peak, false is a negative peak
template <typename T>
std::vector<bool> findPeaksType(std::vector<T> &vec, const int peakWidth, const T peakThreshold) {
	std::vector<bool> peaks;
	int winSize=2*peakWidth;
	typename std::vector<T>::iterator winStart=vec.begin();
	typename std::vector<T>::iterator winEnd=vec.begin()+winSize;
	int max, mean, min;
	
	for (int i=peakWidth; i < vec.size()-peakWidth; i++, winStart++, winEnd++){		
		max=*std::max_element(winStart, winEnd);
		min=*std::min_element(winStart, winEnd);
		mean=std::accumulate(winStart, winEnd, 0.0)/winSize;
		if (max > mean+peakThreshold) {
			peaks.push_back(true);
			i+=winSize;
			winStart+=winSize;
			winEnd+=winSize;
		} else if ( min < mean-peakThreshold) {
			peaks.push_back(false);
			i+=winSize;
			winStart+=winSize;
			winEnd+=winSize;
		};
	};
	return peaks;
};


//Smooths over a vector using a moving average (i.e. acts as a low-pass filter)
template <typename T>
std::vector<T> smoothVector(const std::vector<T> &vec, const int bins=3) {
	std::vector<T> smoothedVec(vec.size());
	int i=0;
	for (typename std::vector<T>::const_iterator it=vec.begin()+bins; it != vec.end()+bins; it++,i++) {
		smoothedVec[i]=std::accumulate(it-bins, it+bins, 0.0)/bins;
	};
	return smoothedVec;
};

//Crude code to check if we're in a "clocking" region of data, or a flat exposure/idle section
bool clockCheck(const std::vector<uint32> &data, int start, int end) {
	uint32 maximum = *std::max_element(data.begin()+start, data.begin()+end);
	uint32 minimum = *std::min_element(data.begin()+start, data.begin()+end);
	if ((maximum-minimum) > clock_threshold) {
		return true;
	} else {
		return false;
	};
};

ImageProcess::ImageProcess(int rows, ADC20File *file) : buffer(buffer_size), 
																												file(file), rows(rows), 
																												columns(default_columns), pre_cut(1), post_cut(10),
																												inverted(false), 
																												exposure(true), cleverStart(true),
																												startSet(false){
	
	binning[0]=1;
	binning[1]=1;
  #ifdef DEBUG
	logLevel=LogLevel_t(logDEBUG);
  #else
	logLevel=LogLevel_t(logERROR);
  #endif

};

ImageProcess::ImageProcess(int rows, int columns, ADC20File *file) : buffer(buffer_size), file(file), rows(rows),
																																		 columns(columns), pre_cut(1), post_cut(10),
																																		 inverted(false), 
																																		 exposure(true), cleverStart(true),
																																		 startSet(false){  
	binning[0]=1;
	binning[1]=1;
	#ifdef DEBUG
	logLevel=LogLevel_t(logDEBUG);
  #else
	logLevel=LogLevel_t(logERROR);
  #endif

};

ImageProcess::~ImageProcess(){};

//Initialization routine
void ImageProcess::init() {
	video_polarity=inverted? -1 : 1 ;
	imageStart=-1; 
	transition_window_size=base_transition_window_size*binning[1];
	findStart();
	findPixelLength();
};


//Interfaces to set various parameters
void ImageProcess::setColumns(int cols) {
	columns=cols;
};

void ImageProcess::setRows(int r) {
	rows=r;
};

void ImageProcess::setExpose(bool expose) {
	exposure=expose;
};

void ImageProcess::setCleverStart(bool clever) {
	cleverStart=clever;
};

void ImageProcess::setInverted(bool invert) {
	inverted=invert;
};

void ImageProcess::setComputeOffsets(bool compute) {
	computeOffsets=compute;
};

void ImageProcess::setBinning(int xbin, int ybin){
	binning[0]=xbin;
	binning[1]=ybin;
	transition_window_size=base_transition_window_size*ybin;
};

void ImageProcess::setCuts(int pre, int post) {
	post_cut=post;
	pre_cut=pre;
};

void ImageProcess::setLogLevel(int level) {
	logLevel=LogLevel_t(level);
};

int ImageProcess::getColumns() {
	return columns;
};

bool ImageProcess::gotoStart() {
	if (imageStart >0) {
		file->setPosition(imageStart);
		return true;
	} else {
		return false;
	};	
};

//Wrappers for clockCheck
bool ImageProcess::isClocking(const std::vector<uint32> &vec) {
	int start=0;
	int end=vec.size();
	return clockCheck(vec, start, end);
};

bool ImageProcess::isClocking(const std::vector<uint32> &vec, int start, int end) {
	return clockCheck(vec, start, end);
};


off_t ImageProcess::getStart() {
	return imageStart;
};	

/*
NOT WORKING
 */
uint32 ImageProcess::findNumColumns() {
	return 0;
	bool foundTransition=false;
	int i_transition=-1;
	uint32 numPixels=0;
	while (! foundTransition) {
		file->readData(buffer_size, buffer);
		//First look for a row transition
		i_transition=findTransition2(buffer);
		foundTransition= i_transition==-1 ? false : true;
		if (!foundTransition) {
			numPixels+=findPeaks2<uint32>(buffer, 20, 5000,'n').size();
		} else {
			std::vector<uint32> temp(&buffer[0], &buffer[i_transition]);
			numPixels+=findPeaks2<uint32>(temp, 20, 5000,'n').size();
		};
	}
	columns=numPixels;
	LOG(logDEBUG) << "There are " << numPixels << " columns" << std::endl;
	LOG(logDEBUG) << "Byte position of transition is: " << file->getLastPosition(i_transition) << std::endl;
	return numPixels;
};

off_t ImageProcess::findStart() {	
	bool foundExpose=false;
	bool foundStart=false;
	file->reset();
	//If we took the image with an exposure, look for the exposure where it stops clocking first	
	if (exposure) {
		while (not foundStart) {
			file->readData(buffer_size, buffer);
			for (int i=0; i < buffer_size/window_size-1; i++) {
				if (!foundExpose) {
					if (!isClocking(buffer,i*window_size, (i+1)*window_size)) {
						foundExpose=true;
						LOG(logDEBUG) << "Found start of exposure at: " << file->getLastPosition(0) << std::endl;
					}; 
				} else {
					if (isClocking(buffer, i*window_size, (i+1)*window_size)) {
						imageStart=file->getLastPosition(0);
						foundStart=true;
						break;
					};
				};
			};
		};
		//The above gives us a *rough* estimation of the start point, but we can do better
		//First look for a jump upwards
		for (int i=20; i < buffer_size; i++ ) {
			if (buffer[i] > buffer[i-20]+start_jump) {
				imageStart+=indexToBytes(i);
				break;
			};
		};
		//Next look for a jump between our little fast-clocking waveforms to a row-transition
		file->setPosition(imageStart);
		LOG(logINFO) << "Found jump at:" << imageStart << std::endl;
		file->readData(buffer_size, buffer);
		int transition=findTransition2(buffer);
		LOG(logDEBUG) << "Transiton index is: " << transition << std::endl;
		if (transition > 9000 or transition == -1) {
			transition = 7660;
			LOG(logWARNING) << "Warning: using hardcoded transition location" << std::endl;
		};
		imageStart+=indexToBytes(transition+transition_window_size);
		file->setPosition(imageStart);
		startSet=true;
	} else { // End code for starting with an exposure
		//Now do code for an "image" taking just by clocking the CCD
		bool foundTransition=false;
		while (!foundTransition) {
			int status = file->readData(buffer_size, buffer);
			int transition=findTransition2(buffer);
			if (transition !=-1) {
				imageStart=file->getLastPosition(transition);
				imageStart+=indexToBytes(transition_window_size);
				foundTransition=true;
				startSet=true;
				file->setPosition(imageStart);
			} else if (status == 0x3) {
				LOG(logERROR) << "Error, cannot find row transition" << std::endl;
				return -1;
			};
		};
	};
	LOG(logDEBUG) << "Image start bytes are: " << imageStart << std::endl;
	return imageStart;
};


//This code tries to find a row transition wavefunction
//It doesn't really work very well (matches with some pixel waveforms, for instance)	
int ImageProcess::findTransition(const std::vector<uint32> &vec) {
	//Looks for a transition in the data set, returns the integer of the start of the transition
	//Returns -1 if it can't find a transition
	int win_size=transition_window_size;
	std::vector<double> derivative;
	std::vector<uint32> smoothed;
	//	smoothed=smoothVector<uint32>(vec);
	derivative=finiteDifference<uint32>(vec);
	std::vector<int> peaks;
	peaks = findPeaks2(derivative, 7, 300.0);
	//	std::cout << "Numer of peaks is:" << peaks.size() << std::endl;
	//Now go through and find 2 peaks that are in the right place
	int jitter = 15; // Bins to allow in our derivative for jitter/derivative innaccuracy
	int a,b,c;
	int i_transition=-1;
	unsigned int i=1;
	for (; i < peaks.size(); i++ ) {
		if (i > peaks.size()-4) {return -1;} // Window not long enough
		if ( abs(peaks[i]-peaks[i-1] - win_size/6) < jitter ) {
			//NOw check that the next 3 also meet this condition
			a=abs(peaks[i+1]-peaks[i]);
			b=abs(peaks[i+2]-peaks[i+1]);
			c=abs(peaks[i+3]-peaks[i+2]);
			LOG(logDEBUG) << "Candidate peak at:" << peaks[i] << ":"<<a << ":" << b << ":" << c << std::endl;
			if ( (abs(peaks[i+1]-peaks[i] - win_size/6) < jitter ) and 
					 (abs(peaks[i+2]-peaks[i+1] - win_size/6) < jitter ) and
					 (abs(peaks[i+3]-peaks[i+2] - win_size/6) < jitter ) ){
				i_transition=peaks[i-1];
				break;

			};
		};
	};
	//Now do a little pattern matching: we want a +,-,+,-,+,-, then start tansition is done
	//First create a small window to look at
	/*
	std::vector<double> window(&vec[i_transition-win_size], &vec[i_transition+win_size]);
	double mean=std::accumulate(&vec[peaks[i-1]], &vec[peaks[i+1]],0.0)/(peaks[i+1]-peaks[i-1]);
	std::vector<int> gtmean(window.size());
	for (int j=0; j < gtmean.size(); j++) {
		if (window[j] >= mean) {
			gtmean[j]=1;
		} else {
			gtmean[j]=-1;
		};
	};
	std::vector<int> sum(win_size);
	//Now find the place where this sum is closest to zero, that should be the start of our transition
	for (int j=0; j < win_size; j++) {
		sum[j]=abs(std::accumulate(gtmean.begin()+j, gtmean.begin()+j+win_size, 0.0));
	};
	i_transition=i_transition-win_size+getVectorArgMin(sum);*/
	return i_transition;
};
		
int ImageProcess::findTransition2(const std::vector<uint32> &vec) {
	const double threshold=.05;
	int winSize=transition_window_size;
	std::vector<double> squareWave(winSize);
	//Set up our square wave
	for (int i=0; i < winSize; i++){
		if (i/(winSize/6)%2==0) {
			squareWave[i]=-1*video_polarity;
		} else{
			squareWave[i]=1*video_polarity;
		};
	};
	double normalization=vec.size()*winSize;
	//Now compute a normalized cross-correlation between these two matrices
	std::vector<double> correlation(vec.size());

	for (unsigned int i=0; i < vec.size(); i++ ){
		correlation[i]=0;
		if (i>winSize and i<vec.size()-winSize) {
			for (int j=0; j<winSize; j++ ) {
				correlation[i]+=vec[i-j]*squareWave[j];
			};
			
		};
		correlation[i]/=normalization;
	};
	//If we're below the threshold at our maximum, we don't have what we want anywhere;
	int maxCorr=getVectorArgMax<double>(correlation);
	if (correlation[maxCorr] < threshold) return -1;
	//otherwise let's go

	return maxCorr-177; // Hardcoded. Bad bad Ryan
};


//This goes through a section of the data and algorithmically determines the pixel length
//Due to jitter, we take the median result as the actual length
int ImageProcess::findPixelLength() {
	if (imageStart==-1) return -1;
	gotoStart();
	//Necessary because the 1st pixel is often much lower than the rest which messes with the threshold code
	file->setPosition(imageStart+indexToBytes(start_shift));
	file->readData(buffer_size, buffer);
	gotoStart();
	std::vector<int> peaks;
	std::vector<int> peakToPeaks;
	uint32 max=getVectorMax<uint32>(buffer);
	double mean=getVectorMean<uint32>(buffer);
	uint32 threshold = uint32((double(max)-mean)*peak_threshold);
	char reset_polarity = inverted ? 'n' : 'p';
	peaks = findPeaks2<uint32>(buffer, reset_pulse_width, threshold, reset_polarity);

	for (unsigned int i=1; i < peaks.size(); i++) {
		peakToPeaks.push_back(peaks[i]-peaks[i-1]);
	}; 
	int median=getVectorMedian(peakToPeaks);
	LOG(logDEBUG) << "Median length is: " << median << std::endl;
	pixel_length=median;

	//Now go through an characterize the reset pulses
	int resetIndex=0;
	//Start at at least the second peak
	for (	int i=1; i < peakToPeaks.size(); i++ ) {
		if (abs(peakToPeaks[i]-median)<pixel_jitter) {
			resetIndex=peaks[i];
			break;
		};
	};
	std::vector<int> heights;
	int height;
	std::vector<uint32>::iterator resetWindowStart=buffer.begin();
	std::vector<uint32>::iterator resetWindowEnd=buffer.begin();
	for ( int i=1; i < peaks.size()-1; i++ ) {
		if (buffer.size() < resetIndex+pixel_length*2) break;
		resetWindowStart=buffer.begin()+resetIndex+pixel_length-reset_pulse_width;
		resetWindowEnd=buffer.begin()+resetIndex+pixel_length+reset_pulse_width;
		if (inverted) {
			resetIndex=std::distance(buffer.begin(), std::min_element(resetWindowStart, resetWindowEnd));
		} else {
			resetIndex=std::distance(buffer.begin(), std::max_element(resetWindowStart, resetWindowEnd));
		};
		height=buffer[resetIndex]-buffer[resetIndex+reset_pulse_width];
		heights.push_back(height);
	};	
	return median;
};

//Actually performs the CDS on a single row
std::vector<double> ImageProcess::readRow(int integration_length, int pre_cut, int post_cut){
	std::vector<double> rowPixels(columns);
	std::vector<off_t> pixelOffsets(columns);
	//	int row_length_estimate=pixel_length*columns;
	
	int pixels_read=0;
	int pixels_per_buffer=buffer_size/pixel_length;
	double pedestal, signal;
	//char reset_polarity = inverted ? 'p' : 'n';
	if ((integration_length*2+pre_cut+post_cut) >= pixel_length) {
		LOG(logWARNING) << "Warning: integration time was too long for pixel length, reducing integration time..." << std::endl;
		integration_length=(pixel_length/2)-pre_cut-post_cut;
	};
	int i_sum_well=pixel_length/2;
	int i_ped_start=i_sum_well-integration_length-pre_cut;
	int i_sig_start=i_sum_well+post_cut;
	int i_ped_end=i_sum_well-pre_cut;
	int i_sig_end=i_sum_well+post_cut+integration_length;
	std::vector<uint32> pixel(pixel_length);
	off_t rowStart=file->getPosition();
	while (pixels_read<columns) {
		file->readData(buffer_size, buffer);
		std::vector<uint32>::iterator pixel_start=buffer.begin();
		for (int i=0; i < pixels_per_buffer-1; i++ ) {
			//Perform our CDS
			pedestal=std::accumulate(pixel_start+i_ped_start, pixel_start+i_ped_end, 0.0)/integration_length;
			signal=std::accumulate(pixel_start+i_sig_start, pixel_start+i_sig_end, 0.0)/integration_length;
			rowPixels[pixels_read]=signal-pedestal;
			if (computeOffsets) {
				pixelOffsets[pixels_read]=file->getLastPosition(std::distance(buffer.begin(), pixel_start));
			};
			if (rowPixels[pixels_read] < -1000) {
				//	std::cout << "Pixel at: " << file->getLastPosition(std::distance(buffer.begin(), pixel_start)) 
				//				<< "Pixel #: " << pixels_read << std::endl;
			};
			pixels_read++;
			if (pixels_read == columns) {
				pixel_start+=pixel_length;
				break;
			};
			if (!inverted) {
				pixel_start=std::max_element(pixel_start+pixel_length-pixel_jitter, pixel_start+pixel_length+pixel_jitter);
			} else {
				pixel_start=std::min_element(pixel_start+pixel_length-pixel_jitter, pixel_start+pixel_length+pixel_jitter);
			};
		};
		file->setPosition(file->getLastPosition(std::distance(buffer.begin(), pixel_start)));
	};
	//Skip over our transition window
	file->setPosition(file->getPosition()+indexToBytes(transition_window_size));
	if (computeOffsets) {
		pixelStarts.push_back(pixelOffsets);
	};
	//Compute the average for the row, if it's very large, spit it out
	//	double rowMean=std::accumulate(rowPixels.begin(), rowPixels.end(), 0.0)/columns;
	//if (rowMean > 1000) std::cout << "Row starting at " << rowStart << " bytes has mean " << rowMean ;
	//	std::cout << "Current file position is: " << file->getPosition() ;
	//std::cout << "Estimated length is: " << row_length_estimate ;
	return rowPixels;
};

//This function handles reading through the entire image, calling readRow to read row by row and forms the
//vector of the CDS result
void ImageProcess::readImage(int integration_window){
	LOG(logINFO) << "Processing image with parameters: \nInt Window: " << integration_window << "\nPre Cut:   " << pre_cut << "\nPost Cut:   " << post_cut << std::endl;
	for (int i=0; i < rows; i++ ) {
		//std::cout << "reading row: " << i ;
		if (i==1) {
			LOG(logDEBUG) << "Row 1 starts at: " << file->getPosition() << " bytes."  << std::endl;
		};
		rowStarts.push_back(file->getPosition());
		image.push_back(readRow(integration_window, pre_cut, post_cut));
	};
	flat_image=flatten<double>(image);
};

bool ImageProcess::writeImage(const char *filename) {
	fitsfile *fFile;
	int status=0;
	int bitpix=DOUBLE_IMG;
	long naxis=2;
	long naxes[2] = {columns, rows};	
	std::remove(filename);
	if (!fits_create_file(&fFile, filename, &status) ) {
		if (fits_create_img(fFile, bitpix, naxis, naxes, &status)) {LOG(logERROR) << "Error during image creation: " << status  << std::endl;};
		if (fits_write_img(fFile, TDOUBLE, 1,naxes[0]*naxes[1], &flat_image[0], &status)) {LOG(logERROR) << "Error during image writing: " << status << std::endl;};
		if (fits_close_file(fFile, &status)) {LOG(logERROR) << "Error during file closing: " << status << std::endl;}
	};
	LOG(logINFO) << "Writing to fits file: " << filename << std::endl;
	return true;
};

//Writes the byte offsets of all the pixels in the image to a file,
//so they can be used by external programs to test various things
bool ImageProcess::writeOffsets(const char *filename) {
	std::ofstream outputFile(filename);
	//First check if we actually computed our offsets
	if (pixelStarts.empty()) {
		LOG(logWARNING) << "Warning: got told to write offsets file without computing them" << std::endl;
		return false;
	};

	LOG(logINFO) << "Writing to root file: " << filename << std::endl;
	for (int i=0; i < rows; i++) {
		bool first = true; 
		std::vector<off_t> rowOffsets=pixelStarts[i];
		for (int j=0; j < columns;j++) { 
			bool last = i==columns-1;
			if (!first and !last) { outputFile << ","; } 
			first = false;
			outputFile << rowOffsets[j];
		};
		if (i != rows-1) {outputFile << std::endl;}
	};	
	return true;
};


std::vector<double> ImageProcess::computeNoiseIntegrationTime(const std::vector<uint32> &vec, int pre_cut, int post_cut, int step) {
	int num_pixels=vec.size()/pixel_length;
	int i_sum_well=pixel_length/2;
	int i_sig_start=i_sum_well+post_cut;
	int i_ped_end=i_sum_well-pre_cut;
	double pedestal, signal;
	std::vector<double> noise;
	std::vector<double> pixels;

	std::vector<uint32>::const_iterator first_pixel;
	if (!inverted) {
		first_pixel=std::max_element(vec.begin(), vec.begin()+pixel_length/2);
	} else {
		first_pixel=std::min_element(vec.begin(), vec.begin()+pixel_length/2);
	};
	for (int integration_length=1; integration_length < pixel_length/2; integration_length+=step) {		
		int i_ped_start=i_sum_well-integration_length-pre_cut;
		int i_sig_end=i_sum_well+post_cut+integration_length;
		std::vector<uint32>::const_iterator pixel_start=first_pixel;
		for (int i=0; i < num_pixels; i++ ) {
			pedestal=std::accumulate(pixel_start+i_ped_start, pixel_start+i_ped_end, 0.0)/integration_length;
			signal=std::accumulate(pixel_start+i_sig_start, pixel_start+i_sig_end, 0.0)/integration_length;
			pixels.push_back(signal-pedestal);
			if (!inverted) {
				pixel_start=std::max_element(pixel_start+pixel_length-pixel_jitter, pixel_start+pixel_length+pixel_jitter);
			} else {
				pixel_start=std::min_element(pixel_start+pixel_length-pixel_jitter, pixel_start+pixel_length+pixel_jitter);
			};
		};
		noise.push_back(computeNoise(pixels));
		pixels.clear();
	};
	return noise;
};

std::vector<double> ImageProcess::rowNoiseIntegration(int row, int step) {
	off_t rowStart=rowStarts[row];
	std::vector<double> noise;
	int maxWindow=pixel_length/2-2-pre_cut-post_cut;
	for (int window=1; window < maxWindow; window+=step){
		file->setPosition(rowStart);
		std::vector<double> pixels=readRow(window, pre_cut, post_cut);
		noise.push_back(computeNoise(pixels));
	};
	return noise;
};

double ImageProcess::computeNoise(const std::vector<double> &pixels) {
	double max= getVectorMax(pixels);
	double min= getVectorMin(pixels);
	int nbins=int((max-min)*5);
	
	TH1D* pixel_distribution_temp = new TH1D("pixel_distribution_temp", "pixel_distribution_temp", nbins, min, max);
	for (std::vector<double>::const_iterator it=pixels.begin(); it!= pixels.end(); it++ ) {
		pixel_distribution_temp->Fill(*it);
	};
	int maximum=pixel_distribution_temp->GetBinContent(pixel_distribution_temp->GetMaximumBin());
	double minrange=min;
	double maxrange=max;
	double median=getVectorMedian<double>(pixels);
	double meanGuess=median;
	double sigmaGuess;
	if (inverted) {
		sigmaGuess = median-minrange;
		maxrange = median+sigmaGuess;
	} else {
		sigmaGuess = maxrange-median;
		minrange = median-sigmaGuess;
	};

	TF1* g = new TF1("g", "gaus", minrange, maxrange);

	g->SetParameter(0, maximum);
	g->SetParameter(1, meanGuess);
	g->SetParameter(2, sigmaGuess);

	pixel_distribution_temp->Fit("g", "QLR", "",minrange, maxrange);

	double sum = std::accumulate(pixels.begin(), pixels.end(), 0.0);
	double mean = sum / pixels.size();
	mean = g->GetParameter(1);
	double sigma = g->GetParameter(2);

	delete g;
	delete pixel_distribution_temp;
	return sigma;
}; 


std::vector<double> ImageProcess::rowNoise() {
	std::vector<double> noise;
	for (int row=0; row < rows; row++ ) {
		noise.push_back(computeNoise(image[row]));
	};
	return noise;
};

std::vector<double> ImageProcess::rowAverageValue() {
	std::vector<double> noise;
	for (int row=0; row < rows; row++) {
		noise.push_back(getVectorMedian(image[row]));
	};
	return noise;
};

void ImageProcess::computeStatistics(bool writeRoot, const char *filename) {

	double max= getVectorMax(flat_image);
	double min= getVectorMin(flat_image);
	int nbins=int((max-min)*5);
	//	gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
	std::string rootFname=file->getFilename();
	if (filename) {
		rootFname=filename;
	} else {
		rootFname.replace(rootFname.end()-3, rootFname.end(), "root");
	};
	TFile *f1 = new TFile(rootFname.c_str(), "RECREATE");
	TH1D* pixel_distribution = new TH1D("pixel_distribution", "pixel_distribution", nbins, min, max);
	for (std::vector<double>::const_iterator it=flat_image.begin(); it!= flat_image.end(); it++ ) {
		pixel_distribution->Fill(*it);
	};
	int maximum=pixel_distribution->GetBinContent(pixel_distribution->GetMaximumBin());
	double minrange=min;
	double maxrange=max;
	double median=getVectorMedian<double>(flat_image);
	double meanGuess=median;
	double sigmaGuess;
	if (inverted) {
		sigmaGuess = median-minrange;
		maxrange = median+sigmaGuess;
	} else {
		sigmaGuess = maxrange-median;
		minrange = median-sigmaGuess;
	};

	TF1* g = new TF1("g", "gaus", minrange, maxrange);

	g->SetParameter(0, maximum);
	g->SetParameter(1, meanGuess);
	g->SetParameter(2, sigmaGuess);

	pixel_distribution->Fit("g", "QLR", "",minrange, maxrange);

	double sum = std::accumulate(flat_image.begin(), flat_image.end(), 0.0);
	double mean = sum / flat_image.size();
	mean = g->GetParameter(1);
	double sigma = g->GetParameter(2);
	std::cout << "Mean: " << mean <<std::endl;
	std::cout << "Std : " << sigma <<std::endl;
	//Now do noise vs. integration window
	LOG(logDEBUG) << "Computing noise vs integration window.";
	std::vector<double> noise_int=rowNoiseIntegration(rows-1);
	std::vector<double> x=arange<double>(1,noise_int.size()+1);
	LOG(logDEBUG) << "Length of noise v integration window array: " + noise_int.size() <<":" << x.size() << std::endl;

	std::vector<double> row_avg=rowAverageValue();
	TGraph* nint=new TGraph(noise_int.size(), &x[0],&noise_int[0]);
	nint->SetName("noise_int");
	nint->SetTitle("Noise vs Integration Window");

		
	std::vector<double> rNoise=rowNoise();
	x=arange<double>(0, rows);
	TGraph* rnoise = new TGraph(rNoise.size(), &x[0], &rNoise[0]);
	rnoise->SetName("row_noise");
	rnoise->SetTitle("Noise per row");
	TGraph* ravg=new TGraph(row_avg.size(), &x[0], &row_avg[0]);
	ravg->SetName("row_avg");
	ravg->SetTitle("Average value per row");

	if (writeRoot) {
		LOG(logINFO) << "Writing to root file: " << filename << std::endl;
		pixel_distribution->Write();
		nint->Write();
		rnoise->Write();
		ravg->Write();
	};
	if (pixel_distribution) delete pixel_distribution;
	if (nint) delete nint;
	if (rnoise) delete rnoise;
	if (ravg) delete ravg;
	if (g) delete g;
};
