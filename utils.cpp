#include <stdio.h>

#include "utils.h"

double dt(timeval start, timeval end)
{ return ((double)(end.tv_sec-start.tv_sec)*1e3 
	  + (double)(end.tv_usec-start.tv_usec)*1e-3); }

uint16 bsw(uint16 value)
{ return ( ((value&0xff00)>>8)+ ((value&0xff)<<8) );}

uint32 bsw(uint32 value)
{
  return ( ((value&0xff000000)>>24) + ((value&0x00ff0000)>>8) + 
           ((value&0x0000ff00)<<8)  + ((value&0x000000ff)<<24));
}

void bit_print(uint16 value)
{
  int i;
  printf("| ");
  for(i=15;i>=0;i--) {
    uint16 a = ((0x0001<<i)&value)>>i;
    printf("%d ",a);
    if(!(i%8)) printf("| ");
  }
  printf("\n");
}

void bit_print(uint32 value)
{
  int i;
  printf("| ");
  for(i=31;i>=0;i--) {
    uint32 a = ((0x00000001<<i)&value)>>i;
    printf("%d ",a);
    if(!(i%8)) printf("| ");
  }
  printf("\n");
}
bool test_bit(uint32 value, uint32 bit)
{
  return (bool( ((0x00000001<<bit)&value)>>bit ));
}
void bit_print(uint64 value)
{
  int i;
  uint64 b64 = 0x00000000000000001;
 
  printf("| ");
  for(i=63;i>=0;i--) {
    //uint64 a = ((0x0000000000000001<<i)&value)>>i;
    uint64 a = ((b64<<i)&value)>>i;
    printf("%llu ",a);
    if(!(i%8)) printf("| ");
  }
  printf("\n");
}

