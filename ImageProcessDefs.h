#ifndef IMAGEPROCESSDEFS_H
#define IMAGEPROCESSDEFS_H

//Size of our initial reading buffer
const int buffer_size=33333;
//const int buffer_size=66666;
const int BUFFERSIZE=buffer_size;
const double peak_threshold = .7;
//const uint32 clock_threshold = 200;
const uint32 clock_threshold = 2000;
const int start_shift=100*86;
const int default_columns=4200;

//Window for determining various startup stuff (must be multiple of 3)
const int window_size=3333;
//Size required for our first pre-image "jump"
const int start_jump=300;
//Size of our pixel row transition, in samples (roughly)
const int base_transition_window_size=216;
const int pixel_jitter=5;
const int reset_pulse_width=10;

#define DEFAULT_LOGFILE "process_data.log"
//#define DEFAULT_LOGFILE "/dev/null"

#endif
