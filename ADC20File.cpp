#include "ADC20File.h"
#include <stdio.h>
#include <iostream>
#include <byteswap.h>
#include <string.h>


using namespace std;

/*
The masks for our MSB/LSB for each sample in each 2-word set of data read
 */
/*
//Old attempts using a single 64-bit read
//This didn't work, for unknown reasons
const uint64 MSB1 = 0x3ff;
const uint64 LSB1 = 0xffc00;
const uint64 MSB2 = 0x3ff00000;
const uint64 LSB2 = 0x3ff00000000;
const uint64 MSB3 = 0xffc0000000000;
const uint64 LSB3 = 0x3ff0000000000000;
*/

const uint32 MSB1 = 0x3ff;
const uint32 LSB1 = 0xffc00;
const uint32 MSB2 = 0x3ff00000;
const uint32 LSB2 = MSB1;
const uint32 MSB3 = LSB1;
const uint32 LSB3 = MSB2;


std::string ADC20File::getFilename(){return fname;};

ADC20File::ADC20File(const char* filename)
{
	pFile=NULL;
	//	cout << "Opening file " << filename << endl;
	pFile=fopen(filename, "rb");
	//Get the length of the file
	fseeko(pFile, 0L, SEEK_END);
	fileLength = ftell(pFile);
	fseeko(pFile, 0L, SEEK_SET);
	if (!readHeader() ) {
		cout << "Error reading header" << endl;
	};
	fname=string(filename);
};

ADC20File::~ADC20File() {
	if (pFile) {
		fclose(pFile);
	};
};

bool ADC20File::readHeader() {

	if (pFile) {
		//		cout << "Reading header" << endl;
		fread(&buffer32, sizeof(uint32), 1, pFile);
		//cout << "First word is " << buffer32 << endl;
		return true;
	} else {
		return false;
	};
};


int ADC20File::readData(int numSamples, vector<uint32> &pData) {
	//Reads numSamples samples.
	//Note that numSamples **MUST** be a multiple of 3
	if (numSamples%3 != 0) return 0x2;
	lastNumSamples=numSamples;
	int numReads = numSamples/3;
	size_t readWords;
	lastPos=ftello(pFile);
	for (int i=0; i < numReads; i++) {
		readWords = fread(buffer_arr, sizeof(uint32), 2, pFile);
		if (readWords != 2) {
			//			memset((void*)pData[i], 0, (numSamples-i*3)*sizeof(uint32));
			for (int j=i; j < numReads; j++) {
				pData[j*3]=0;
				pData[j*3+1]=0;
				pData[j*3+2]=0;
			};
			return 0x3;
		};
		//		if (i==0) printf("%016llx\n", buffer64);
		//buffer64=__bswap_64(buffer64);
		//if (i==0) printf("%016llx\n", buffer64);
		/*
		pData[i*3] = ((buffer64 & MSB1) << 10) + ((buffer64 & LSB1) >> 10 );
		pData[i*3+1] = ((buffer64 & MSB2) >> 10) + ((buffer64 & LSB2) >> 32 );
		pData[i*3+2] = ((buffer64 & MSB3) >> 32) + ((buffer64 & LSB3) >> 52 );
		*/
		buffer_arr[0]=__bswap_32(buffer_arr[0]);
		buffer_arr[1]=__bswap_32(buffer_arr[1]);
		pData[i*3] = ((buffer_arr[0] & MSB1) << 10) + ((buffer_arr[0] & LSB1) >> 10 );
		pData[i*3+1] = ((buffer_arr[0] & MSB2) >> 10) + ((buffer_arr[1] & LSB2) );
		pData[i*3+2] = ((buffer_arr[1] & MSB3) ) + ((buffer_arr[1] & LSB3) >> 20 );

	};
	curPos=ftell(pFile);
	return 0;
};

off_t ADC20File::getLastPosition(int index) {
	
	off_t bytes_from_start= ((off_t) index )/3*4*2;
	return lastPos+bytes_from_start;
}

off_t ADC20File::getFileLength() {
	return fileLength;
};

bool ADC20File::setPosition(off_t position) {
	//TODO: add error code so we can't set position past end of file
	fseeko(pFile, position, SEEK_SET);
	curPos=position;
	return true;
};
	
void ADC20File::reset() {
	fseeko(pFile, 0, SEEK_SET);
	readHeader();
};

