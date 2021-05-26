#include <iostream>
#include "IdleProcess.h"
#include "fft.h"
#include "vec_funcs.h"
#include "gnuplot_i.hpp"

const int buffer_size=3;

IdleProcess::IdleProcess(ADC20File *file) : buffer(buffer_size), file(file) {
	setMinIntWindow();
	setMaxIntWindow();
};
IdleProcess::~IdleProcess(){};

void IdleProcess::setMinIntWindow(int minimum) {
	minWindow=minimum;
};

void IdleProcess::setMaxIntWindow(int maximum) {
	maxWindow=maximum;
};

std::vector<double> IdleProcess::computePSD(int winsize) {
	// CArray data(winsize);
  // PSD.resize(winsize/2);
	// CArray PSD_arr(Complex(0.0,0.0),winsize);
	// std::fill(PSD.begin(), PSD.end(), 0.0);
	// int loops=0;
	// bool atEnd=false;
	// uint32 mean=0;
	// //First establish a rough mean value
	// for (int i=0; i < 100; i++) {
	// 	file->readData(buffer_size, buffer);
	// 	mean+=buffer[0]+buffer[1]+buffer[2];
	// };
	// mean/=(100*3);
	// std::cout <<"Mean is: " << mean << std::endl;
	// while (!atEnd) {
	// 	data=Complex(0.0,0.0);
	// 	for (int i=0; i < winsize-3; i+=3) {
	// 		if (file->readData(buffer_size, buffer)!=0x3) {
	// 			data[i]=Complex(buffer[0],0.0)-Complex(mean,0.0);
	// 			data[i+1]=Complex(buffer[1],0.0)-Complex(mean,0.0);
	// 			data[i+2]=Complex(buffer[2],0.0)-Complex(mean,0.0);				 
	// 		} else {
	// 			atEnd=true;
	// 			break;
	// 		};
	// 	};
	// 	if (!atEnd) {
	// 		//if (loops==0) std::cout <<data[0]<<std::endl;
	// 		fft(data);
	// 		//	if (loops==0) std::cout << data[0] << std::endl;
	// 		data=std::abs(data);
	// 		if (loops==3) std::cout << "data[0]: "<< data[0] << std::endl;
	// 		data=std::pow(data, Complex(2.0,0.0));
	// 		if (loops==3) std::cout << "data[0]:" << data[0] << std::endl;
	// 		PSD_arr+=data;
	// 		//std::cout << PSD_arr[0] << std::endl;
	// 		loops++;
	// 	};
	// };

	// //Finally take the mean, then the root, of our PSD
	// std::cout << PSD_arr[0] << std::endl;
	// PSD_arr/=Complex(double(loops),0.0);
	// PSD_arr=std::sqrt(PSD_arr);
	// for (int i=0; i < winsize/2; i++ ) {
	// 	PSD[i]=PSD_arr[i].real();
	// };
	// return PSD;
};

void IdleProcess::plotPSD() {
	if (PSD.size()==0) return;
	std::vector<double>x=linspace<double>(0,800000,PSD.size());
	
	//	std::cout << PSD.size() << ":" << x.size() << std::endl;
	Gnuplot g1;

	g1.set_style("lines").set_ylogscale().set_xlogscale().plot_xy(x,PSD);
	
	std::cin.get();
};

