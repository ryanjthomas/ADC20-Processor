#ifndef MY_LOG_H
#define MY_LOG_H

enum LogLevel_t {
	logERROR,
	logWARNING,
	logINFO,
	logDEBUG
};

LogLevel_t logLevel;

#define LOG(level) \
	if (level > logLevel); \
	else std::cout


#endif
