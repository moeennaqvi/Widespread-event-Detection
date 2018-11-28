#ifndef _PREPROCESS_H
#define _PREPROCESS_H

#pragma warning(disable:4996)

typedef struct __parsedData {
	long long srcIp;
	long long dstIp;
}ParsedData;

void getDstSet(int* numberOfDst, long long** mainArray, std::set<long long>** mainSet, int option);
int dataCompare(const void* a, const void* b);
int getHashValue(long long inputData, int bitmapSize);
double LOG2(double x);
long long RandomNumber(int min, int max);

#endif // !_PREPROCESS_H
