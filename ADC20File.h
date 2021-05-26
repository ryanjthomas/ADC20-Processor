#ifndef ADC20FILE_H
#define ADC20FILE_H

#include <stdio.h>
#include <string>
#include <vector>
#include "utils.h"

typedef enum {
	SUCCESS = 0x0,
	OPEN_ERROR = 0x01,
	READ_ERROR = 0x10
} ADC20FILE_STATUS;

class ADC20File {
 public:
	ADC20File(const char *filename);
	~ADC20File();
	int readData(int numSamples, std::vector<uint32> &pData);
	
	std::string getFilename();
	
	//Gets the position for the value at index from the *last* read
	off_t getLastPosition(int index);
	off_t getFileLength();
	off_t getPosition() {return curPos;}
	bool setPosition(off_t position);
	void reset();
	FILE* pFile;
	uint32 buffer32;
	uint64 buffer64;
	uint32 buffer_arr[2];


 private:
	std::string fname;
	int lastNumSamples;
	off_t fileLength;
	off_t lastPos;
	off_t curPos;
	bool readHeader();
};
#endif
