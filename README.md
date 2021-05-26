# ADC20-Processor

C++ code to process raw CCD video streams taken using a 20-bit ADC
developed for the DAMIC collaboration. The code contains signal
processing logic to attempt to find the start of the image (it assumes
the start of the video stream occured during an exposure, with no
video output). It will before digital CDS, and output a .fits file.

Code is provided as-is. It includes a lot of logic to find the start
of image taking patterns, and to detect vertical clocking, reset, and
summing well pulses, but it never worked very well, and was abandoned
in favor of on-FPGA signal processing, or adding pixel trigger signals
to the ADC stream to make locating pedestal/signal regions trivial.

The TCLAP library (http://tclap.sourceforge.net) is included to parse
command line options.

### Prerequisites

* C++03+ compiler (g++ is assumed)
* ROOT (for computing statistics of the output file, mainly the noise)
* cfitsio (for generating a .fits file)

