#include <iostream>
#include <set>

#include "PREPROCESS.h"
#include "openssl\md5.h"

void getDstSet(int* numberOfDst, long long** mainArray, std::set<long long>** mainSet, int option) 
{

	FILE* readStream = 0x00L;
	readStream = fopen("src_dst_text.txt", "r");

	long long srcIp = 0, dstIp = 0;
	int packetCount = 0;

	//reading total packet count from file
	//printf("Reading data count from file...\n");
	while (!feof(readStream)) 
	{
		fscanf(readStream, "%lld %lld", &srcIp, &dstIp);
		++packetCount;
	}
	fclose(readStream);
	readStream = fopen("src_dst_text.txt", "r");

	std::set<long long> dstSet; //for temporary use
	std::set<long long>* finalSet = 0x00L; //finalMain set to return 
	ParsedData* dataSet = 0x00L; //for temp reason inside this function
	long long currentDstIp = 0; // current destination IP 
	long long* dstIpArray = 0x00L; //destination IP array
	int* countingArray = 0x00L; // no purpose
	dataSet = (ParsedData*)malloc(sizeof(ParsedData) * packetCount); //allocating memory to dataSet var
	

	packetCount = 0;
	//reading data from file and inserting into dataSet and dstSet
	//printf("Reading actual data from file...\n");
	while (!feof(readStream)) 
	{
		fscanf(readStream, "%lld %lld", &srcIp, &dstIp);
		long long finalDst = dstIp >> 5;
		dataSet[packetCount].srcIp = srcIp;
		dataSet[packetCount].dstIp = finalDst;
		/* for independent destination ip */
		dstSet.insert(finalDst);
		++packetCount;
	}
	//printf("Reading completed...\n");
	qsort(dataSet, packetCount, sizeof(ParsedData), dataCompare); //sorting data array
	//printf("Sorting data array...\n");
	//printf("Total set size is %d\n", dstSet.size());
	//printf("Total packet count is %d\n", packetCount);

	finalSet = new std::set<long long>[dstSet.size()]; //allocating memory size to final set

	dstIpArray = new long long[dstSet.size()]; //allocating memory to destination IP array

	int setFlag = 0;
	currentDstIp = dataSet[0].dstIp; //destination IP of first dataset
	dstIpArray[0] = dataSet[0].dstIp; //first Destination IP 

	//putting values in finalset and dstIpArray  
	//printf("Values started to put in output variables...\n");
	for (int i = 0; i < packetCount; ++i) 
	{
		if (currentDstIp != dataSet[i].dstIp) 
		{
			if (option == 1) 
			{
				if (finalSet[setFlag].size() <= 1000) 
				{
					++setFlag;
					dstIpArray[setFlag] = dataSet[i].dstIp;
				}
				else 
				{
					finalSet[setFlag].clear();
				}
			}
			else 
			{
				++setFlag;
				dstIpArray[setFlag] = dataSet[i].dstIp;
			}
			currentDstIp = dataSet[i].dstIp;
		}
		finalSet[setFlag].insert(dataSet[i].srcIp);
		//printf("Source IP: %lld\n", dataSet[i].srcIp);
	}

	*numberOfDst = (setFlag + 1);
	*mainArray = dstIpArray;
	*mainSet = finalSet;

	if (dataSet != 0x00L) 
	{
		free(dataSet);
		dataSet = 0x00L;
	}
	if (readStream != 0x00L) 
	{
		fclose(readStream);
		readStream = 0x00L;
	}
}

int dataCompare(const void* a, const void* b) 
{
	const ParsedData *m, *n;
	m = (const ParsedData*)a;
	n = (const ParsedData*)b;

	if (m->dstIp < n->dstIp) return -1;
	else if (m->dstIp == n->dstIp) return 0;
	else return 1;
}

int getHashValue(long long inputData, int bitmapSize) {
	unsigned char digest[MD5_DIGEST_LENGTH];
	MD5_CTX ctx;
	MD5_Init(&ctx);
	uint64_t frontResultValue = 0, endResultValue = 0;
	int lastValue = 0;
	char inputString[20];
	sprintf(inputString, "%lld", inputData);

	MD5_Update(&ctx, inputString, strlen(inputString));
	MD5_Final(digest, &ctx);

	frontResultValue = *((uint64_t*)digest);
	endResultValue = *((uint64_t*)(digest + 8));
	lastValue = (frontResultValue ^ endResultValue) % bitmapSize;
	return lastValue;
}
double LOG2(double x) {
	return log(x) / log(2.0);
}

long long RandomNumber(int min, int max)
{
	long long randomValue;
	long long RandMax = (long long) sqrt(max);
	do
	{
		long long r1 = rand() % (RandMax);
		long long r2 = rand() % (RandMax);
		randomValue = r1*RandMax + r2;

	} while (randomValue < min);
	return randomValue;
}

