#ifndef IDLEPROCESS_H
#define IDLEPROCESS_H

#include <vector>

#include "utils.h"
#include  "ADC20File.h"


const int default_PSD_winsize=8096;//*1024;

class IdleProcess {
 public: 
	IdleProcess(ADC20File *file);
	~IdleProcess();
	
	void process();
	
	void setMinIntWindow(int minimum=2);
  void setMaxIntWindow(int maximum=100);
	void plotPSD();
	std::vector<double> computePSD(int winsize=default_PSD_winsize);
	//	std::vector<double> computeNoiseIntegrationTime();

 private:

	int minWindow, maxWindow;

	std::vector<double> PSD;
	std::vector<uint32> buffer;

	ADC20File *file;
	

};


#endif
