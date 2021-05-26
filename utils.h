#ifndef UTILS_H
#include <sys/time.h>

typedef unsigned long long uint64;
typedef unsigned int uint32;
typedef int int32;
typedef unsigned short uint16;

extern uint16 bsw(uint16 value);

extern uint32 bsw(uint32 value);

extern void bit_print(uint16 value);

extern void bit_print(uint32 value);

extern void bit_print(uint64 value);

extern bool test_bit(uint32 value, uint32 bit);

double dt(timeval start, timeval end);
#endif
