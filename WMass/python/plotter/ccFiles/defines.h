#ifndef DEFINES_H
#define DEFINES_H

typedef enum {BToH=0, BToF, GToH} DataEra;
typedef enum {MC=0, Data} DataType;

bool isOddEvent(ULong64_t evt) {

  return (evt%2) ? 1 : 0;       

}

bool isEvenEvent(ULong64_t evt) {

  return (evt%2) ? 0 : 1;       

}


#endif
