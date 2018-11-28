#include <iostream>
#include <set>
#include <map>
#include <math.h>
#include <stdint.h>
#include <fstream>
#include <algorithm>
#include "BloomFilter.h"
#include "CombinableFilter.h"


void getHashedValue(long long inputData, unsigned char* hashedValue)
{
	unsigned char digest[MD5_DIGEST_LENGTH];
	MD5_CTX ctx;
	MD5_Init(&ctx);
	char inputString[20];
	sprintf_s(inputString, "%lld", inputData);
	MD5_Update(&ctx, inputString, strlen(inputString));
	MD5_Final(digest, &ctx); //computed hash
	hashedValue = digest;

}

long long getHashedIndex(unsigned int elem, long long bitmapSize)
{
	unsigned char hashedValue[MD5_DIGEST_LENGTH];
	getHashedValue(elem, hashedValue);
	uint64_t frontResultValue = 0, endResultValue = 0;
	frontResultValue = *((uint64_t*)hashedValue); //first part of hash
	//endResultValue = *((uint64_t*)(hashedValue + 8)); //last part of hash
	//printf("%" PRIu64 "\n", frontResultValue);
	long long index = frontResultValue % bitmapSize;
	return index;
}


//Method to generate bloom filter for each monitor
void GenerateFilters(int numberOfAgents, char scheme, long long* monitorsArray, set<long long>* monitorEventsSet, map<long long, BloomFilter*>** Filters)
{
	map<long long, BloomFilter*>* tempFilterPointer = new map<long long, BloomFilter*>();
	for (int i = 0; i <numberOfAgents; ++i)
	{
		set<long long>::iterator iter;

		//printf("\n %lld monitoring agent has %lld events", monitorsArray[i], monitorEventsSet[i].size());
		BloomFilter* bf = new BloomFilter(monitorEventsSet[i].size(), 99.9); //Created Bloom filter
		//values will be set according to scheme
		if (scheme == 'c') //this will be used for CFP scheme
		{
			double e = (100 - bf->GetAccuracyPercentage()) / 100.0;
			bf->SetVectorLength((int)pow(2, ceil(LOG2(bf->GetNumberofElements()*(-log(e) / (log(2)*log(2)))))));
			bf->SetNumberofHashes((int)ceil(-log(e) / log(2)));
		}
		else if(scheme == 'p') //this will be used for CFPS scheme
		{
			bf->SetVectorLength((int)pow(2, ceil(LOG2(bf->GetNumberofElements() / log(2)))));
			bf->SetNumberofHashes(1);
			bf->SetAccuracyPercentage(50);
		}
		for (iter = monitorEventsSet[i].begin(); iter != monitorEventsSet[i].end(); ++iter) //iterate through all events of each monitoring agent
		{			
			bf->insertElementByMD5(*iter);
		}
		if (tempFilterPointer->find(monitorsArray[i]) == tempFilterPointer->end())
		{
			tempFilterPointer->insert(std::make_pair(monitorsArray[i], bf));
		}		
	}
	*Filters = tempFilterPointer;
}

// This method will extend all varying size bloom filters to same size
void ExtendFilters(map<long long, BloomFilter*>* filter, int maxSize, map<long long, BloomFilter*>** ExtendedFilters)
{
	map<long long, BloomFilter*>* tempExtendedFilter = new map<long long, BloomFilter*>();
	map<long long, BloomFilter*>::iterator it_er;
	int i = 0;
	for (it_er = filter->begin(); it_er != filter->end(); ++it_er) //This will iterate through all the monitoring agents and their Bloom filters 
	{		
		int currentSize = it_er->second->GetVectorLength();
		vector<bool> bf_current(it_er->second->GetVector());
		vector<bool> bf_temp;
		bf_temp.reserve(maxSize);
		bf_temp.resize(maxSize);
		int iterations = 0;
		if(currentSize > 0)
			iterations = maxSize / currentSize;
		//Loop for copying small loop to larger
		for (int i = 0; i < iterations; i++)
		{
			for (int j = 0; j < currentSize; j++)
			{
				bf_temp[j + i*currentSize] = bf_current[j];
			}
		}
		BloomFilter* tempBF = new BloomFilter(it_er->second->GetNumberofElements(), it_er->second->GetAccuracyPercentage()); // Making a new bloom filter
		//setting all values for tempBF (it will work for both CFP and CFPS)
		tempBF->SetNumberofHashes(it_er->second->GetNumberofHashes());
		tempBF->SetVectorLength(maxSize);
		tempBF->SetVector(bf_temp);

		tempExtendedFilter->insert(std::make_pair(it_er->first, tempBF)); //insert new bloomfilter and its key
		
	}
	*ExtendedFilters = tempExtendedFilter;
}

void GetCombinedFilter(map<long long, BloomFilter*>* extendedFilters, int maxSize, vector<bool>* combinedFilter)
{
	vector<bool> temp_combined(maxSize, true); // initialize all the combined filter's values with 1.
	//outer loop on extended filter
	map<long long, BloomFilter*>::iterator it_er;
	for (it_er = extendedFilters->begin(); it_er != extendedFilters->end(); ++it_er)
	{
		vector<bool> t = it_er->second->GetVector();
		//inner loop to compare results
		for (int i = 0; i < maxSize; i++)
		{
			if (temp_combined[i] == false)
				continue;
			else if (t[i] == false)
				temp_combined[i] = false;
		}
	}
	*combinedFilter = temp_combined;
}

//This method will perform membership testing and return widespread event set
void PerformMemberShipTest(int numberOfAgents, long long* monitorsArray, set<long long>* monitorEventsSet, map<long long, BloomFilter*>* combinedFilter, set<long long>** widespreadEvents)
{
	set<long long>* temp_wide_events = 0x00L;
	temp_wide_events =  new set<long long>[numberOfAgents];

	map<long long, BloomFilter*>::iterator it_er;
	int i = 0;
	for (it_er = combinedFilter->begin(); it_er != combinedFilter->end() && i< numberOfAgents; ++it_er, ++i) //This will iterate through each monitoring  
	{
		set<long long>::iterator iter;
		for (iter = monitorEventsSet[i].begin(); iter != monitorEventsSet[i+504].end(); ++iter) //iterate through all events of each monitoring agent
		{
			if (it_er->second->isFoundElement(*iter))
			{
				temp_wide_events[i].insert(*iter);
			}
		}
	}

	/*
	for (int i = 0; i < numberOfAgents; i++) //outer loop for each monitoring agent
	{
		std::set<long long>::iterator iter;
		for (iter = monitorEventsSet[i].begin(); iter != monitorEventsSet[i].end(); ++iter) //iterate through all events of each monitoring agent
		{
			//Check if element exists in Combined bloom filter
			if (combinedBloomFilter->isFoundElement(*iter)) //if element is found in BF
			{
				temp_wide_events[numberOfAgents].insert(*iter);
			}
		}
	}
	*/
	*widespreadEvents = temp_wide_events;
}

//This method size of all bits of bloom filter and also the largest bloom filter size among all
long long CalculateBloomFilterSize(map<long long, BloomFilter*>* filter, int* maximumBF)
{
	long long tempSize = 0;
	int tempMax = 0;
	map<long long, BloomFilter*>::iterator it_er;
	for (it_er = filter->begin(); it_er != filter->end(); ++it_er) //This will count total number of real events that occur. 
	{
		if (it_er->second->GetVectorLength() > tempMax)
			tempMax = it_er->second->GetVectorLength();
		tempSize += it_er->second->GetVectorLength();
	}
	*maximumBF = tempMax;
	return tempSize;
}

void CFS(long long* data, int numberOfAgents, char scheme, long long* monitorsArray, set<long long>* monitorEventSet, set<long long>** widespreadEventSet)
{
	map<long long, BloomFilter*>* Filters;
	map<long long, BloomFilter*>* ExtendedFilters;
	vector<bool> combinedVector;
	set<long long>* WideEvents = 0x00L;
	//printf("\nAll agents are making their Bloom filters...\n");
	GenerateFilters(numberOfAgents, scheme, monitorsArray, monitorEventSet, &Filters);
	//printf("\nMonitors sending their Bloom filters to Coordinator..\n");
	
	//Calculate total data sent and find Maximum bloom filter size
	int maxBF = 0;
	long long dataSize = CalculateBloomFilterSize(Filters,&maxBF);
	//printf("\nTotal %lld bits are sent to coordinator..\n", dataSize);
	if (maxBF > 0)
	{
		if (scheme == 'p')
			printf("CSFP's ");
		else
			printf("CFS's ");
		printf("Largest BloomFilter size is %d..\n", maxBF);
	}
	
	//printf("\nCoordinator is now making Extended bloom filters..\n");
	
	//Make all Bloom filters with same size
	ExtendFilters(Filters, maxBF, &ExtendedFilters);	
	
	//Perform Bitwise AND to obtain Combined Filter
	combinedVector.reserve(maxBF);
	combinedVector.resize(maxBF);
	GetCombinedFilter(ExtendedFilters,maxBF,&combinedVector);
	//printf("\nCoordinator performed Bitwise AND to obtain Combined Filter..\n");
	
	//Calculate total size of all Bloom filters
	long long dataSize2 = combinedVector.size()*numberOfAgents;
	//printf("\nCoordinator sending Combined filter to all monitors. \n%lld bits are sent\n", dataSize2);
	//printf("Total communication: %lld\n", dataSize + dataSize2);
	
	//assign combined filter to all agents
	map<long long, BloomFilter*>::iterator it_er;
	for (it_er = ExtendedFilters->begin(); it_er != ExtendedFilters->end(); ++it_er) //It will iterate through each agent in Extended filter
	{
		it_er->second->SetVector(combinedVector); //change each vector to combined one
	}
	//Perform membership testing to obtain widespread events..
	/*
	BloomFilter* combinedBloomFilter = new BloomFilter(maxBF, 'e'); //creating a temp BF to pass for membership testing
	combinedBloomFilter->SetVector(combinedVector);
	*/
	PerformMemberShipTest(numberOfAgents, monitorsArray, monitorEventSet, ExtendedFilters, &WideEvents);
	//printf("\nMonitors performing membership testing to obtain widespread events..\n");
	//Get widespread events
	*widespreadEventSet = WideEvents;
	*data = dataSize + dataSize2;
	/*
	if (Filters != 0x00L)
	{
		free(Filters);
		//Filters = 0x00L;
	}
	if (ExtendedFilters != 0x00L)
	{
		free(ExtendedFilters);
		//ExtendedFilters = 0x00L;
	}
	*/
}

void CFPS(long long* data, int numberOfAgents, long long* monitorsArray, set<long long>* monitorEventSet, set<long long>** widespreadEventSet)
{
	long long outerTemp = 0;
	//printf("\nCFPS scheme: \n");
	set<long long>* WideEvents = 0x00L;
	double e = (100 - 99.9) / 100.0;
	int rounds = (int)ceil(-log(e) / log(2));
	for (int i = 0; i < rounds; i++)
	{
		long long tempData = 0;
		//printf("\nFor Iteration %d\n", i + 1);
		if (i == 0)
		{
			CFS(&tempData, numberOfAgents, 'p', monitorsArray, monitorEventSet, &WideEvents);
		}
		else 
		{
			CFS(&tempData, numberOfAgents, 'p', monitorsArray, WideEvents, &WideEvents);
		}
		outerTemp += tempData;
	}
	*widespreadEventSet = WideEvents;
	*data = outerTemp;
}

void CFS_Simulated(long long* data, int numberOfAgents, long long minEventSet, long long nCommonEvents, double Rint, double error, map<int, long long>* monitors, string scheme)
{
	long long dataLength = 0;
	map<int, long long>* bloomFilters = new map<int, long long>();
	map<int, long long>::iterator it_er;
	long long vectorLength = 0;
	for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
	{
		if (scheme == "real")
		{
			vectorLength = (int)pow(2, ceil(LOG2(it_er->second*(-log(error) / (log(2)*log(2))))));
		}
		if (scheme == "compact")
		{
			vectorLength = (int)(-(log(error) / (log(2)*log(2))*it_er->second));
		}		
		dataLength += vectorLength; //calculating data length
		bloomFilters->insert(std::make_pair(it_er->first, vectorLength));
	}
	//Bloom filters made 
	//now find max size filter
	long long maxLength = 0;
	for (it_er = bloomFilters->begin(); it_er != bloomFilters->end(); ++it_er)
	{
		if (it_er->second > maxLength)
			maxLength = it_er->second;
	}
	//Make extended filter and combined filter
	map<int, long long>* ExtendedbloomFilters = new map<int, long long>();
	for (it_er = bloomFilters->begin(); it_er != bloomFilters->end(); ++it_er)
	{
		dataLength += maxLength;
		ExtendedbloomFilters->insert(std::make_pair(it_er->first, maxLength));
	}
	*data = dataLength;
}

void CFSP_Simulated(long long* data, int numberOfAgents, long long minEventSet, long long nCommonEvents, double Rint, double error, map<int, long long>* monitors, string scheme)
{
	int rounds = (int)ceil(-log(error) / log(2.0));
	printf("Rounds= %d", rounds);
	long long dataLength = 0;
	map<int, long long>::iterator it_er;
	map<int, long long>* tempMonitor = new map<int, long long>();
	//inserting values into temp Monitor to further process
	for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
	{
		tempMonitor->insert(std::make_pair(it_er->first, it_er->second));
	}
	for (int i = 0; i < rounds; i++)
	{
		map<int, long long>* bloomFilters = new map<int, long long>();
		long long vectorLength = 0;
		for (it_er = tempMonitor->begin(); it_er != tempMonitor->end(); ++it_er)
		{
			if (scheme == "real")
			{
				vectorLength = (int)pow(2, ceil(LOG2(it_er->second / log(2)))); //calculate vector length according to size
			}
			if (scheme == "compact")
			{
				vectorLength = (int)it_er->second / log(2);
			}			
			dataLength += vectorLength; //calculating data length
			bloomFilters->insert(std::make_pair(it_er->first, vectorLength));
		}
		//Bloom filters made 
		//now find max size filter
		long long maxLength = 0;
		for (it_er = bloomFilters->begin(); it_er != bloomFilters->end(); ++it_er)
		{
			if (it_er->second > maxLength)
				maxLength = it_er->second;
		}
		//Make extended filter and combined filter
		map<int, long long>* ExtendedbloomFilters = new map<int, long long>();
		for (it_er = bloomFilters->begin(); it_er != bloomFilters->end(); ++it_er)
		{
			dataLength += maxLength;
			ExtendedbloomFilters->insert(std::make_pair(it_er->first, maxLength));
		}

		//Now the reduction in size of each monitor set. Since after each iteration, half of non-common events are eliminated. So, we reduce second of monitors.
		for (it_er = tempMonitor->begin(); it_er != tempMonitor->end(); ++it_er)
		{
			long long nonCommonEvents = it_er->second - nCommonEvents;
			//reduce it to half
			long long reduction = nonCommonEvents / 2;
			it_er->second -= reduction;
		}
	}
	*data = dataLength;
}
void RawData_Simulated(long long* data, int numberOfAgents, long long averageEventSize, long long nCommonEvents, map<int, long long>* monitors)
{
	long long dataLength = 0;
	map<int, long long>::iterator it_er;
	for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
	{
		dataLength += it_er->second*averageEventSize; //data that monitors send to coordinator
	}
	//Now coordinator sending back data to monitors..
	dataLength += numberOfAgents*nCommonEvents*averageEventSize;
	*data = dataLength;
}


void BDPThetaFalseEvents_Simulated(int* BDPDfalseEvents, int* BDPD90falseEvents, int* BDPD80falseEvents, int* BDPD70falseEvents, int numberOfAgents, long long minEventSet, long long nCommonEvents, long long nEvents90, long long nEvents80, long long nEvents70,  double Rint, int theta, long long bitmapSize, int digestSize, set<long long>* monitorEventSet)
{
	//step 1: Each monitor observing its events

	//step 2: each Monitor sending number of events to coordinator.
	

	//Step 3: Coordinator calculating appropriate bitmap size.. Already done.. (we know bitmap size before this).. 

	//step 4: coordinator sending back bitmap size to each monitor after determining it

	
	// step 5: every monitor makes a bitmap of its eventset according to the bitmap size told by coordinator
	char** bitmapArray = 0x00L; //arrays of each monitor
	bitmapArray = (char**)malloc(sizeof(char*) * numberOfAgents);
	for (int i = 0; i < numberOfAgents; ++i)
	{
		bitmapArray[i] = 0x00L;
		bitmapArray[i] = (char*)malloc(sizeof(char) * bitmapSize); //allocate memory to each bitmap
		memset(bitmapArray[i], '0', sizeof(char) * bitmapSize); //initialize each bitmap with '0'
	}
	
	printf("Bitmap size: %lld \n", bitmapSize);
	//Make hash table to store events against hash values
	vector<map<int, long long>*> hashTableArray;
	//hashTableArray.reserve(numberOfAgents);
	for (int i = 0; i < numberOfAgents; ++i)
	{
		std::set<long long>::iterator iter;
		map<int, long long>* tempHash = new map<int, long long>();
		for (iter = monitorEventSet[i].begin(); iter != monitorEventSet[i].end(); ++iter)
		{
			int hashValue = getHashValue(*iter, bitmapSize); //getting hash function's return.. 
			bitmapArray[i][hashValue] = 1; //Make a certain bit as 1 in 2D bitmap array (mapping bits to bit array)
			tempHash->insert(std::make_pair(hashValue, *iter));
		}
		hashTableArray.push_back(tempHash);
	}
	// Step 6: Monitors send their bitmaps to coordinator

	//step 7: Coordinator makes an integer array with value against each index. Those index with value more than theta..

	int* coordinator = 0x00L; //int array used to store coordinator array values after all bitmaps are received.. 
	coordinator = (int*)malloc(sizeof(int) * bitmapSize);
	memset(coordinator, 0, sizeof(int)*bitmapSize); //initialize coordinator with 0 values across all

	for (int i = 0; i < bitmapSize; ++i) //according to each bit of bitmap
	{
		for (int j = 0; j < numberOfAgents; ++j) //for each monitor
		{
			if (bitmapArray[j][i] == 1) //if a bit is 1. Subsequent bit is made 1. 
				++coordinator[i];
		}
	}

	//Now coordinator determine which indexes are above theta threshold..
	
	vector<int> indexes100;
	/*
	vector<int> indexes90;
	vector<int> indexes80;
	vector<int> indexes70;
	*/
	for (int i = 0; i < bitmapSize; ++i)
	{
		if (coordinator[i] >= theta)
		{
			indexes100.push_back(i);			
		}
		else
		{
			/*
			if (coordinator[i] >= (int)theta*0.9)
			{
				indexes90.push_back(i);				
			}
			else
			{
				if (coordinator[i] >= (int)theta*0.8)
				{
					indexes80.push_back(i);
				}
				else
				{
					if (coordinator[i] >= (int)theta*0.7)
					{
						indexes70.push_back(i);
					}
				}
			}
			*/
		}
	}


	//printf("W indexes: %d \n", indexes100.size());
	//printf("0.9 indexes: %d \n", indexes90.size());
	//printf("0.8 indexes: %d \n", indexes80.size());
	//printf("0.7 indexes: %d \n", indexes70.size());
	// Step 8: Coordinator sending back indexes of potential widespread events to monitors

	//Step 9: Each monitor determining its potential widespread events using the indexes provided by Coordinator.
	
	set<long long>* monitorsPotentialEvents100 = 0x00L; //finalMain set to return
	//set<long long>* monitorsPotentialEvents90 = 0x00L; //finalMain set to return
	//set<long long>* monitorsPotentialEvents80 = 0x00L; //finalMain set to return
	//set<long long>* monitorsPotentialEvents70 = 0x00L; //finalMain set to return
	monitorsPotentialEvents100 = new set<long long>[numberOfAgents]; //number of monitors
	//monitorsPotentialEvents90 = new set<long long>[numberOfAgents]; //number of monitors
	//monitorsPotentialEvents80 = new set<long long>[numberOfAgents]; //number of monitors
	//monitorsPotentialEvents70 = new set<long long>[numberOfAgents]; //number of monitors

	//Each monitor for Checking events for each event if its equal to indexes provided by coordinator

	for (int i = 0; i < numberOfAgents; i++)
	{
		std::set<long long>::iterator iter;
		for (iter = monitorEventSet[i].begin(); iter != monitorEventSet[i].end(); ++iter)
		{
			int hashValue = getHashValue(*iter, bitmapSize); //getting hash function's return..
			if (bitmapArray[i][hashValue] == 1)
			{
				for (std::vector<int>::iterator index100_it = indexes100.begin(); index100_it != indexes100.end(); ++index100_it)
				{
					if (hashValue == *index100_it)
					{
						monitorsPotentialEvents100[i].insert(*iter);
						break;
					}

				}
				/*
				for (std::vector<int>::iterator index90_it = indexes90.begin(); index90_it != indexes90.end(); ++index90_it)
				{
					if (hashValue == *index90_it)
					{
						monitorsPotentialEvents90[i].insert(*iter);
						break;
					}
				}
				for (std::vector<int>::iterator index80_it = indexes80.begin(); index80_it != indexes80.end(); ++index80_it)
				{
					if (hashValue == *index80_it)
					{
						monitorsPotentialEvents80[i].insert(*iter);
						break;
					}
				}
				for (std::vector<int>::iterator index70_it = indexes70.begin(); index70_it != indexes70.end(); ++index70_it)
				{
					if (hashValue == *index70_it)
					{
						monitorsPotentialEvents70[i].insert(*iter);
						break;
					}
				}
				*/
			}
		}
		printf("monitor done\n");
	}

	indexes100.clear();
	indexes100.shrink_to_fit();
	//indexes90.clear();
	//indexes80.clear();
	//indexes70.clear();
	//indexes90.shrink_to_fit();
	//indexes80.shrink_to_fit();
	//indexes70.shrink_to_fit();
	hashTableArray.clear();
	hashTableArray.shrink_to_fit();

	//Monitors generating an array of digests of the length provided

	set<int>* monitorsDigestEvents100 = 0x00L; //finalMain set to return
	monitorsDigestEvents100 = new set<int>[numberOfAgents]; //number of monitors
	/*
	set<int>* monitorsDigestEvents90 = 0x00L; //finalMain set to return
	set<int>* monitorsDigestEvents80 = 0x00L; //finalMain set to return
	set<int>* monitorsDigestEvents70 = 0x00L; //finalMain set to return
	
	monitorsDigestEvents90 = new set<int>[numberOfAgents]; //number of monitors
	monitorsDigestEvents80 = new set<int>[numberOfAgents]; //number of monitors
	monitorsDigestEvents70 = new set<int>[numberOfAgents]; //number of monitors
	*/
	std::set<long long>::iterator iter;
	int size = pow(2, digestSize);
	int digest = 0;
	for (int i = 0; i < numberOfAgents; i++)
	{
		
		for (iter = monitorsPotentialEvents100[i].begin(); iter != monitorsPotentialEvents100[i].end(); ++iter)
		{
			
			digest = getHashValue(*iter, size);
			monitorsDigestEvents100[i].insert(digest);
		}
		/*
		for (iter = monitorsPotentialEvents90[i].begin(); iter != monitorsPotentialEvents90[i].end(); ++iter)
		{
			digest = getHashValue(*iter, size);
			monitorsDigestEvents90[i].insert(digest);
		}
		
		for (iter = monitorsPotentialEvents80[i].begin(); iter != monitorsPotentialEvents80[i].end(); ++iter)
		{	
			digest = getHashValue(*iter, size);
			monitorsDigestEvents80[i].insert(digest);
		}
		
		for (iter = monitorsPotentialEvents70[i].begin(); iter != monitorsPotentialEvents70[i].end(); ++iter)
		{
			digest = getHashValue(*iter, size);
			monitorsDigestEvents70[i].insert(digest);
		}
		*/
	}
	printf("Digest done\n");
	//Step 10: Monitors sending the digested events to coordinator including false events
	

	//Step 11: Coordinator determines widespread event after comparing the digests to check which one have theta number of observance..

	map<int, int> digestCount100;
	//map<int, int> digestCount90;
	//map<int, int> digestCount80;
	//map<int, int> digestCount70;

	std::set<int>::iterator it_er;
	std::map<int, int>::iterator it;
	std::pair<std::map<int, int>::iterator, bool> ret;
	for (int i = 0; i < numberOfAgents; i++)
	{
		//printf("WSize: %d \n", monitorsDigestEvents100[i].size());
		for (it_er = monitorsDigestEvents100[i].begin(); it_er != monitorsDigestEvents100[i].end(); ++it_er)
		{
			it = digestCount100.find(*it_er);
			if (it == digestCount100.end())
			{				
				ret = digestCount100.insert(std::make_pair(*it_er, 1));
				if (ret.second == false)
				{
					printf("Failed to insert.. ");
				}
			}
			else
			{
				it->second += 1;
			}
		}
		//printf("0.9 Size: %d \n", monitorsDigestEvents90[i].size());
		/*
		for (it_er = monitorsDigestEvents90[i].begin(); it_er != monitorsDigestEvents90[i].end(); ++it_er)
		{
			it = digestCount90.find(*it_er);
			if (it == digestCount90.end())
			{
				digestCount90.insert(std::make_pair(*it_er, 1));
			}
			else
			{
				it->second += 1;
			}
		}
		//printf("0.8 Size: %d \n", monitorsDigestEvents80[i].size());
		for (it_er = monitorsDigestEvents80[i].begin(); it_er != monitorsDigestEvents80[i].end(); ++it_er)
		{
			it = digestCount80.find(*it_er);
			if (it == digestCount80.end())
			{
				digestCount80.insert(std::make_pair(*it_er, 1));
			}
			else
			{
				it->second += 1;
			}
		}
		//printf("0.7 Size: %d \n", monitorsDigestEvents70[i].size());
		for (it_er = monitorsDigestEvents70[i].begin(); it_er != monitorsDigestEvents70[i].end(); ++it_er)
		{
			it = digestCount70.find(*it_er);
			if (it == digestCount70.end())
			{
				digestCount70.insert(std::make_pair(*it_er, 1));
			}
			else
			{
				it->second += 1;
			}
		}
		*/
	}
	
	//Now checking for count
	int numberofEvents100 = 0;
	int numberofEvents90 = 0;
	int numberofEvents80 = 0;
	int numberofEvents70 = 0;



	for (it = digestCount100.begin(); it != digestCount100.end(); ++it)
	{
		if (it->second >= numberOfAgents)
		{
			numberofEvents100++;
		}
	}
	/*
	for (it = digestCount90.begin(); it != digestCount90.end(); ++it)
	{
		if (it->second >= (int)0.9*numberOfAgents)
		{
			numberofEvents90++;
		}
	}
	for (it = digestCount80.begin(); it != digestCount80.end(); ++it)
	{
		if (it->second >= (int)0.8*numberOfAgents)
		{
			numberofEvents80++;
		}
	}
	for (it = digestCount70.begin(); it != digestCount70.end(); ++it)
	{
		if (it->second >= (int)0.7*numberOfAgents)
		{
			numberofEvents70++;
		}
	}
	*/
	//Step 12: Coordinator send indexes of respective widespread events of Hi to respective monitors

	printf("Detected W events: %d\n", numberofEvents100);
	//printf("Detected 0.9 events: %d\n", (numberofEvents100 + numberofEvents90));
	//printf("Detected 0.8 events: %d\n", (numberofEvents100 + numberofEvents90 + numberofEvents80));
	//printf("Detected 0.7 events: %d\n", (numberofEvents100 + numberofEvents90 + numberofEvents80 + numberofEvents70));

	*BDPDfalseEvents = abs(numberofEvents100 - nCommonEvents);
	*BDPD90falseEvents = abs((numberofEvents100 + numberofEvents90) - (nCommonEvents+ nEvents90));
	*BDPD80falseEvents = abs((numberofEvents100 + numberofEvents90 + numberofEvents80) - (nCommonEvents + nEvents90 + nEvents80));
	*BDPD70falseEvents = abs((numberofEvents100 + numberofEvents90 + numberofEvents80 + numberofEvents70) - (nCommonEvents + nEvents90 + nEvents80 + nEvents70));

	
	for (int i = 0; i < numberOfAgents; i++)
	{
		monitorsPotentialEvents100[i].clear();
		monitorsDigestEvents100[i].clear();
		/*
		monitorsPotentialEvents90[i].clear();
		monitorsPotentialEvents80[i].clear();
		monitorsPotentialEvents70[i].clear();
		
		monitorsDigestEvents90[i].clear();
		monitorsDigestEvents80[i].clear();
		monitorsDigestEvents70[i].clear();
		*/
	}

	digestCount100.clear();
	
	//digestCount90.clear();
	//digestCount80.clear();
	//digestCount70.clear();
	

	if (bitmapArray != 0x00L)
	{
		free(bitmapArray[0]);
		free(bitmapArray);
		bitmapArray = 0x00L;
	}
	if (coordinator != 0x00L)
	{
		free(coordinator);
		coordinator = 0x00L;
	}

}

void BDPTheta_Simulated(long double* BDPdata, long double* BDPDdata, long double* BDP90data, long double* BDP80data, long double* BDP70data, long double* BDPD90data, long double* BDPD80data, long double* BDPD70data, int numberOfAgents, long long minEventSet, long long nCommonEvents, double Rint, int theta, long long bitmapSize, int digestSize, set<long long>* monitorEventSet)
{
	long long BDPdataLength = 0;
	long long BDPDdataLength = 0;
	long long BDP90dataLength = 0;
	long long BDP80dataLength = 0;
	long long BDP70dataLength = 0;
	long long BDPD90dataLength = 0;
	long long BDPD80dataLength = 0;
	long long BDPD70dataLength = 0;
	//step 1: Each monitor observing its events

	//step 2: each Monitor sending number of events to coordinator.
	double step2 = 0;
	for (int i = 0; i < numberOfAgents; i++)
	{
		step2 += LOG2(monitorEventSet[i].size()) + 1;
	}
	BDPdataLength += step2;
	BDPDdataLength += step2;
	BDP90dataLength += step2;
	BDP80dataLength += step2; 
	BDP70dataLength += step2; 
	BDPD90dataLength += step2; 
	BDPD80dataLength += step2; 
	BDPD70dataLength += step2;
	//printf("Step2: %lld ", BDPdataLength);
	//printf("Step 2: \n");
	//Step 3: Coordinator calculating appropriate bitmap size.. Already done.. (we know bitmap size before this).. 

	//step 4: coordinator sending back bitmap size to each monitor after determining it
	double step4 = 0;
	for (int i = 0; i < numberOfAgents; i++)
	{

		step4 += LOG2(bitmapSize) + 1;
	}
	BDPdataLength += step4;
	BDPDdataLength += step4;
	BDP90dataLength += step4;
	BDP80dataLength += step4;
	BDP70dataLength += step4;
	BDPD90dataLength += step4;
	BDPD80dataLength += step4;
	BDPD70dataLength += step4;

	//printf("Step4: %lld ", BDPdataLength);
	//printf("Step 4: \n");
	// step 5: every monitor makes a bitmap of its eventset according to the bitmap size told by coordinator
	char** bitmapArray = 0x00L; //arrays of each monitor
	bitmapArray = (char**)malloc(sizeof(char*) * numberOfAgents);
	for (int i = 0; i < numberOfAgents; ++i)
	{
		bitmapArray[i] = 0x00L;
		bitmapArray[i] = (char*)malloc(sizeof(char) * bitmapSize); //allocate memory to each bitmap
		memset(bitmapArray[i], '0', sizeof(char) * bitmapSize); //initialize each bitmap with '0'
	}
	//Make hash table to store events against hash values
	vector<map<int, long long>*> hashTableArray;
	//hashTableArray.reserve(numberOfAgents);
	for (int i = 0; i < numberOfAgents; ++i)
	{
		std::set<long long>::iterator iter;
		map<int, long long>* tempHash = new map<int, long long>();
		int j = 0;
		for (iter = monitorEventSet[i].begin(); iter != monitorEventSet[i].end(); ++iter)
		{
			int hashValue = getHashValue(*iter, bitmapSize); //getting hash function's return.. 
			bitmapArray[i][hashValue] = 1; //Make a certain bit as 1 in 2D bitmap array (mapping bits to bit array)
			tempHash->insert(std::make_pair(hashValue, *iter));
		}
		hashTableArray.push_back(tempHash);
		//tempHash->clear();
	}
	// Step 6: Monitors send their bitmaps to coordinator
	double step6 = 0;
	for (int i = 0; i < numberOfAgents; i++) //All schemes have it
	{

		step6 += (int)bitmapSize;
	}
	BDPdataLength += step6;
	BDPDdataLength += step6;
	BDP90dataLength += step6;
	BDP80dataLength += step6;
	BDP70dataLength += step6;
	BDPD90dataLength += step6;
	BDPD80dataLength += step6;
	BDPD70dataLength += step6;

	//printf("Step6: %lld ", BDPdataLength);
	//printf("Step 6: \n");
	//step 7: Coordinator makes an integer array with value against each index. Those index with value more than theta..

	int* coordinator = 0x00L; //int array used to store coordinator array values after all bitmaps are received.. 
	coordinator = (int*)malloc(sizeof(int) * bitmapSize);
	memset(coordinator, 0, sizeof(int)*bitmapSize); //initialize coordinator with 0 values across all

	for (int i = 0; i < bitmapSize; ++i) //according to each bit of bitmap
	{
		for (int j = 0; j < numberOfAgents; ++j) //for each monitor
		{
			if (bitmapArray[j][i] == 1) //if a bit is 1. Subsequent bit is made 1. 
				++coordinator[i];
		}
	}

	//Now coordinator determine which indexes are above theta threshold..
	int numberofIndex100 = 0;
	int numberofIndex90 = 0;
	int numberofIndex80 = 0;
	int numberofIndex70 = 0;
	vector<int> indexes100;
	vector<int> indexes90;
	vector<int> indexes80;
	vector<int> indexes70;
	for (int i = 0; i < bitmapSize; ++i)
	{
		if (coordinator[i] >= theta)
		{
			indexes100.push_back(i);
			numberofIndex100+= coordinator[i];
		}
		else
		{
			if (coordinator[i] >= theta*0.9)
			{
				indexes90.push_back(i);
				numberofIndex90+= coordinator[i];
			}
			else
			{
				if (coordinator[i] >= theta*0.8)
				{
					indexes80.push_back(i);
					numberofIndex80 += coordinator[i];
				}
				else
				{
					if (coordinator[i] >= theta*0.7)
					{
						indexes70.push_back(i);
						numberofIndex70 += coordinator[i];
					}
				}				
			}		
		}		
	}

	// Step 8: Coordinator sending back indexes of potential widespread events to monitors
	double step8 = 0;
	step8 += LOG2(numberofIndex100) + 1;
	BDPdataLength += LOG2(numberofIndex100) + 1;
	BDPDdataLength += LOG2(numberofIndex100) + 1;
	BDP90dataLength += LOG2(numberofIndex90 + numberofIndex100) + 1;
	BDP80dataLength += LOG2(numberofIndex80 + numberofIndex90 + numberofIndex100) + 1;
	BDP70dataLength += LOG2(numberofIndex70 + numberofIndex80 + numberofIndex90 + numberofIndex100) + 1;
	BDPD90dataLength += LOG2(numberofIndex90 + numberofIndex100) + 1;
	BDPD80dataLength += LOG2(numberofIndex80 + numberofIndex90 + numberofIndex100) + 1;
	BDPD70dataLength += LOG2(numberofIndex70 + numberofIndex80 + numberofIndex90 + numberofIndex100) + 1;
	
	//printf("Step8: %lld ", BDPdataLength);
	//printf("Step 8: \n");

	//Step 9: Each monitor determining its potential widespread events using the indexes provided by Coordinator.
	/*
	set<long long>* monitorsPotentialEvents100 = 0x00L; //finalMain set to return
	set<long long>* monitorsPotentialEvents90 = 0x00L; //finalMain set to return
	set<long long>* monitorsPotentialEvents80 = 0x00L; //finalMain set to return
	set<long long>* monitorsPotentialEvents70 = 0x00L; //finalMain set to return
	monitorsPotentialEvents100 = new set<long long>[numberOfAgents]; //number of monitors
	monitorsPotentialEvents90 = new set<long long>[numberOfAgents]; //number of monitors
	monitorsPotentialEvents80 = new set<long long>[numberOfAgents]; //number of monitors
	monitorsPotentialEvents70 = new set<long long>[numberOfAgents]; //number of monitors
	*/

	int numberofPotentialEvents100[52];
	int numberofPotentialEvents90[52];
	int numberofPotentialEvents80[52];
	int numberofPotentialEvents70[52];

	//assigning values
	for (int i = 0; i < 50; i++)
	{
		numberofPotentialEvents100[i] = 0;
		numberofPotentialEvents90[i] = 0;
		numberofPotentialEvents80[i] = 0;
		numberofPotentialEvents70[i] = 0;
	}

	//Each monitor for Checking events for each event if its equal to indexes provided by coordinator

	//printf("Monitor gonna check events..\n");
	//Monitor searches for events against index
	/*
	for (int i = 0; i < numberOfAgents; i++)
	{
		map<int, long long>* hashTable = hashTableArray[i];
		for (std::vector<int>::iterator index100_it = indexes100.begin(); index100_it != indexes100.end(); ++index100_it)
		{
			int index = *index100_it;
			if (bitmapArray[i][index] == 1)
			{
				std::map<int, long long>::iterator it;
				it = hashTable->find(index);
				if (it != hashTable->end())
				{
					long long potentialEvent = it->second;
					//printf("Potential event: %lld \n", potentialEvent);
					//monitorsPotentialEvents70[i].insert(potentialEvent);
					//monitorsPotentialEvents80[i].insert(potentialEvent);
					//monitorsPotentialEvents90[i].insert(potentialEvent);
					//monitorsPotentialEvents100[i].insert(potentialEvent);
					numberofPotentialEvents100[i]++;
					numberofPotentialEvents90[i]++;
					numberofPotentialEvents80[i]++;
					numberofPotentialEvents70[i]++;
				}
			}
			
		}
		for (std::vector<int>::iterator index90_it = indexes90.begin(); index90_it != indexes90.end(); ++index90_it)
		{
			int index = *index90_it;
			if (bitmapArray[i][index] == 1)
			{
				std::map<int, long long>::iterator it;
				it = hashTable->find(index);
				if (it != hashTable->end())
				{
					long long potentialEvent = it->second;
					//printf("Potential event: %lld \n", potentialEvent);
					//monitorsPotentialEvents70[i].insert(potentialEvent);
					//monitorsPotentialEvents80[i].insert(potentialEvent);
					//monitorsPotentialEvents90[i].insert(potentialEvent);
					numberofPotentialEvents90[i]++;
					numberofPotentialEvents80[i]++;
					numberofPotentialEvents70[i]++;
				}
			}
		}
		
		for (std::vector<int>::iterator index80_it = indexes80.begin(); index80_it != indexes80.end(); ++index80_it)
		{
			int index = *index80_it;
			if (bitmapArray[i][index] == 1)
			{
				std::map<int, long long>::iterator it;
				it = hashTable->find(index);
				if (it != hashTable->end())
				{
					long long potentialEvent = it->second;
					//printf("Potential event: %lld \n", potentialEvent);
					//monitorsPotentialEvents70[i].insert(potentialEvent);
					//monitorsPotentialEvents80[i].insert(potentialEvent);
					numberofPotentialEvents80[i]++;
					numberofPotentialEvents70[i]++;
				}
			}
		}
		for (std::vector<int>::iterator index70_it = indexes70.begin(); index70_it != indexes70.end(); ++index70_it)
		{
			int index = *index70_it;
			if (bitmapArray[i][index] == 1)
			{
				std::map<int, long long>::iterator it;
				it = hashTable->find(index);
				if (it != hashTable->end())
				{
					long long potentialEvent = it->second;
					//printf("Potential event: %lld \n", potentialEvent);
					//monitorsPotentialEvents70[i].insert(potentialEvent);
					numberofPotentialEvents70[i]++;
				}
			}
		}
		hashTable->clear();
	}
	
	indexes100.clear();
	indexes100.shrink_to_fit();
	indexes90.clear();
	indexes80.clear();
	indexes70.clear();
	indexes90.shrink_to_fit();
	indexes80.shrink_to_fit();
	indexes70.shrink_to_fit();	
	hashTableArray.clear();
	hashTableArray.shrink_to_fit();
	*/
	
	for (int i = 0; i < numberOfAgents; i++)
	{
		std::set<long long>::iterator iter;
		for (iter = monitorEventSet[i].begin(); iter != monitorEventSet[i].end(); ++iter)
		{
			int hashValue = getHashValue(*iter, bitmapSize); //getting hash function's return..
			if (bitmapArray[i][hashValue] == 1)
			{
				for (std::vector<int>::iterator index100_it = indexes100.begin(); index100_it != indexes100.end(); ++index100_it)
				{
					if (hashValue == *index100_it)
					{
						//monitorsPotentialEvents100[i].insert(*iter);
						numberofPotentialEvents100[i]++;
						numberofPotentialEvents90[i]++;
						numberofPotentialEvents80[i]++;
						numberofPotentialEvents70[i]++;
						break;
					}

				}
				//printf("100 complete \n");
				for (std::vector<int>::iterator index90_it = indexes90.begin(); index90_it != indexes90.end(); ++index90_it)
				{
					if (hashValue == *index90_it)
					{
						//monitorsPotentialEvents90[i].insert(*iter);
						numberofPotentialEvents90[i]++;
						numberofPotentialEvents80[i]++;
						numberofPotentialEvents70[i]++;
						break;
					}
				}
				//printf("90 complete \n");
				for (std::vector<int>::iterator index80_it = indexes80.begin(); index80_it != indexes80.end(); ++index80_it)
				{
					if (hashValue == *index80_it)
					{
						//monitorsPotentialEvents80[i].insert(*iter);
						numberofPotentialEvents80[i]++;
						numberofPotentialEvents70[i]++;
						break;
					}
				}
				//printf("80 complete \n");
				for (std::vector<int>::iterator index70_it = indexes70.begin(); index70_it != indexes70.end(); ++index70_it)
				{
					if (hashValue == *index70_it)
					{
						//monitorsPotentialEvents70[i].insert(*iter);
						numberofPotentialEvents70[i]++;
						break;
					}
				}
				//printf("70 complete \n");
			}
		}
		printf("monitor done\n");
	}
	
	//Monitors generating an array of digests of 20 bit length each
	//Step 10: Monitors sending the digested events to coordinator including false events
	double step10 = 0;
	for (int i = 0; i < numberOfAgents; i++)
	{
		//printf("100Indexes: %d, Events: %d, w: %lld \n", numberofIndex100, numberofPotentialEvents100[i], nCommonEvents);
		//printf("90Indexes: %d, Events: %d \n", (numberofIndex90+ numberofIndex100), numberofPotentialEvents90[i]);
		//printf("80Indexes: %d, Events: %d \n", (numberofIndex80 + numberofIndex90 + numberofIndex100), numberofPotentialEvents80[i]);
		//printf("70Indexes: %d, Events: %d \n", (numberofIndex70 + numberofIndex80 + numberofIndex90 + numberofIndex100), numberofPotentialEvents70[i]);
		step10 += numberofPotentialEvents100[i] * digestSize;
		BDPdataLength += numberofPotentialEvents100[i] * 64;
		BDPDdataLength += numberofPotentialEvents100[i] * digestSize;
		BDP90dataLength += numberofPotentialEvents90[i] * 64;
		BDP80dataLength += numberofPotentialEvents80[i] * 64;
		BDP70dataLength += numberofPotentialEvents70[i] * 64;
		BDPD90dataLength += numberofPotentialEvents90[i] * digestSize;
		BDPD80dataLength += numberofPotentialEvents80[i] * digestSize;
		BDPD70dataLength += numberofPotentialEvents70[i] * digestSize;
	}

	//printf("Step 10: %lld ", BDPdataLength);
	//dataLength += step10;
	//printf("Step 10: \n");
	//Step 11: Coordinator determines widespread event after comparing the digests to check which one have theta number of observance..

	//Step 12: Coordinator send indexes of respective widespread events of Hi to respective monitors
	double step12 = 0;
	for (int i = 0; i < numberOfAgents; i++)
	{
		step12 += LOG2(nCommonEvents) + 1;
	}
	//printf("Widespread: %d ", nCommonEvents);
	//printf("Step 12 value: %f ", step12);

	BDPdataLength += step12;
	BDPDdataLength += step12;
	BDP90dataLength += step12;
	BDP80dataLength += step12;
	BDP70dataLength += step12;
	BDPD90dataLength += step12;
	BDPD80dataLength += step12;
	BDPD70dataLength += step12;

	//printf("Step12: %lld ", BDPdataLength);
	//printf("Step 12: \n");

	*BDPdata = BDPdataLength;
	*BDPDdata = BDPDdataLength;
	*BDPD90data = BDPD90dataLength;
	*BDPD80data = BDPD80dataLength;
	*BDPD70data = BDPD70dataLength;
	*BDP90data = BDP90dataLength;
	*BDP80data = BDP80dataLength;
	*BDP70data = BDP70dataLength;
	
	
	if (bitmapArray != 0x00L)
	{
		free(bitmapArray[0]);
		free(bitmapArray);
		bitmapArray = 0x00L;
	}
	if (coordinator != 0x00L)
	{
		free(coordinator);
		coordinator = 0x00L;
	}
	
}

void BDP_Simulated(long long* data, int numberOfAgents, long long averageEventSize, long long minEventSet, long long nCommonEvents, double Rint, double error, map<int, long long>* monitors, string scheme, int bitArraySize, int digest)
{
	ofstream StepwiseCommunicationFile;
	StepwiseCommunicationFile.open("StepwiseCommunication.txt", std::ios_base::app);
	//ofstream FalsePositiveRatioFile;
	//FalsePositiveRatioFile.open("FalsePositiveRatio.txt", std::ios_base::app);
	long long dataLength = 0;
	double bitmapSize = 0;
	map<int, long long>::iterator it_er;
	double combinedProbability = 1;
	double averageEvents = 0;
	//step 1: Each monitor observing its events

	//step 2: each Monitor sending number of events to coordinator. 
	double step2 = 0;
	for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
	{
		step2 += LOG2(it_er->second) + 1;		
	}
	if (scheme == "renewed" || scheme == "extended" || scheme == "final")
	{
		StepwiseCommunicationFile << step2 << ",";
		//dataLength += sizeof(int) * 8;
		dataLength += step2;
	}

	//step 3: Coordinator calculating the appropriate bitmapSize
	for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
	{
		averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
	}
	//Check if bitArraySize value is passed by method. If so, replace bitmapSize with that

	if (bitArraySize != 0)
	{
		bitmapSize = bitArraySize;
	}
	else
	{
		//Use equation 4. 
		bitmapSize = (int)ceil(sqrt(nCommonEvents * (averageEvents - nCommonEvents) * averageEventSize));
		//use average events
		//bitmapSize = averageEvents;
	}
		

	//step 4: coordinator sending back bitmap size to each monitor after determining it
	double step4 = 0;
	for (int i = 0; i < numberOfAgents; i++)
	{

		step4 += LOG2(bitmapSize) + 1;
	}
	if (scheme == "renewed" || scheme == "extended" || scheme == "final")
	{
		StepwiseCommunicationFile << step4 << ",";
		//dataLength += sizeof(int) * 8;
		dataLength += step4;
	}	

	// step 5: every monitor makes a bitmap of its eventset according to the bitmap size told by coordinator

	// Step 6: Monitors send their bitmaps to coordinator
	double step6 = 0;
	for (int i = 0; i < numberOfAgents; i++) //All schemes have it
	{

		step6 += (int) bitmapSize;
	}
	StepwiseCommunicationFile << step6 << ",";
	dataLength += step6;

	

	//step 7: coordinator makes an integer array with value against each index. Those index with value more than theta..
	//here nCommonEvents are number of events monitored by all agents. So coordinator will send indexes of these events but with some false events attached

	

	// Step 8: Coordinator sending back indexes of potential widespread events to monitors
	double step8 = 0;
	double totalFalseEvents = 0;
	double eq5FalseEvents = 0;
	for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er) //All schemes have this..
	{
		//Determining number of false events
		double numberofFalseEvents = bitmapSize * (1 - exp(-(nCommonEvents / bitmapSize)))*(1 - exp(-((it_er->second - nCommonEvents) / bitmapSize)));
		eq5FalseEvents += nCommonEvents*((it_er->second - nCommonEvents) / bitmapSize);
		totalFalseEvents += numberofFalseEvents;
		/*
		ofstream falseEventsFile;
		double OurFalseEvents = bitmapSize * (1 - exp(-val)) * val2;
		printf("BitmapSize: %f, Actual false events: %f, Our false events: %f  \n", bitmapSize, numberofFalseEvents, OurFalseEvents);
		falseEventsFile.open("NumberofFalseEvents.txt", std::ios_base::app);
		falseEventsFile << bitmapSize << "," << numberofFalseEvents << "," << OurFalseEvents << "\n"; //writing values	
		falseEventsFile.close();
		*/

		//Each monitor will receive nCommonEvent indexes
		//dataLength += nCommonEvents * sizeof(int) * 8;

		step8 += LOG2(nCommonEvents + (int)numberofFalseEvents) + 1;
		//printf("Bitmap: %f, False: %f Step 8: %f \n", bitmapSize, numberofFalseEvents, step8);
	}
	StepwiseCommunicationFile << step8 << ",";
	dataLength += step8;

	if (scheme == "final")
	{
		//Step 9: Each monitor determining its potential widespread events using the indexes provided by Coordinator.
		//Monitors generating an array of digests of 20 bit length each 
		
		//Step 10: Monitors sending the digested events to coordinator including false events

		double digestSize = 0;
		
		//check if digest value or error is passed.
		if (digest != 0 || error != 0)
		{
			if (digest != 0)
			{
				digestSize = digest;
			}
			else
			{
				digestSize = ceil(LOG2(nCommonEvents / (bitmapSize*error))); //calculate digest according to error rate
			}			
		}
		else
		{
			digestSize = averageEventSize;
		}
					
		//printf("Digest size: %f, Error: %f \n", digestSize, error);
		double falsePositiveProbability = totalFalseEvents/pow(2,digestSize);
		double eq5falsePositives = eq5FalseEvents / pow(2, digestSize);
		//FalsePositiveRatioFile << digestSize << "," << falsePositiveProbability << "," << eq5falsePositiveProb << "\n";
		

		double step10 = 0;
		for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
		{
			//Determining number of false events
			double numberofFalseEvents = bitmapSize * (1 - exp(-(nCommonEvents / bitmapSize)))*(1 - exp(-((it_er->second - nCommonEvents) / bitmapSize)));
			step10 += (nCommonEvents + (int)numberofFalseEvents) * digestSize; //bit corresponding to digest size are transferred.
		}
		StepwiseCommunicationFile << step10 << ",";
		dataLength += step10;

		//Step 11: Coordinator determines widespread event after comparing the digests to check which one have theta number of observance.. 

		//Step 12: Coordinator send indexes of respective widespread events of Hi to respective monitors
		double step12 = 0;
		for (int i = 0; i < numberOfAgents; i++)
		{
			step12 += LOG2(nCommonEvents + (int)eq5falsePositives) + 1;
		}
		StepwiseCommunicationFile << step12 << ",";
		dataLength += step12;
	}
	else
	{
		//Each monitor sends the common event list to coordinator..
		for (int i = 0; i < numberOfAgents; i++)
		{
			dataLength += nCommonEvents * averageEventSize;
		}

		if (scheme == "extended")
		{
			//Coordinator determines common events and sends back to each monitor..

			for (int i = 0; i < numberOfAgents; i++)
			{
				dataLength += nCommonEvents * averageEventSize;
			}
		}
	}	
	StepwiseCommunicationFile << dataLength << "," << bitmapSize << "\n";
	StepwiseCommunicationFile.close();
	//FalsePositiveRatioFile.close();

	//finally monitor receives respective indices and find out the corresponding potential widespread events.. 
	*data = dataLength;

}


void SimulationforBitmapSize()
{
	ofstream BitmapSizeSettingsfile;
	ofstream DigestSizeCommunicationfile;
	int i = 0;
	long long numberofEvents;
	map<int, long long>* monitorsSizes = new map<int, long long>();
	int n = 10; //number of monitors
	while (i < n)
	{
		numberofEvents = RandomNumber(100000, 1000000);
		monitorsSizes->insert(std::make_pair(i, numberofEvents));
		i++;
	}

	//Now loop upon nodes to find minimum and maximum
	long long minEventSet = 1000000; //max possible value
	long long maxEventSet = 100000; //min possible value
	double averageEvents = 0;
	map<int, long long>::iterator it_er;
	for (it_er = monitorsSizes->begin(); it_er != monitorsSizes->end(); ++it_er)
	{
		averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
		if (minEventSet > it_er->second)
			minEventSet = it_er->second;
		if (maxEventSet < it_er->second)
			maxEventSet = it_er->second;
	}

	long long nCommonEvents = 0.1*minEventSet; //finding number of common events
	double Rint = (double)nCommonEvents / minEventSet; //finding intersection ratio
	
	BitmapSizeSettingsfile.open("BitmapSizeSettings.txt", std::ios_base::app);
	BitmapSizeSettingsfile << n << "," << minEventSet << "," << maxEventSet << "," << averageEvents << "," << nCommonEvents << "," << Rint << "\n"; //writing values	
	BitmapSizeSettingsfile.close();
	
	int bitmapSize = (int) ceil(sqrt(nCommonEvents * (averageEvents - nCommonEvents) * 64));
	int digestLength = 20;
	for (int i = 0; i < 15; i++)
	{
		//if (i != 0)
			//bitmapSize = (int) (averageEvents*((0.1*i)+0.1));
		bitmapSize = (int)averageEvents*(0.25)*(i + 1);
		long long BDPFinalLength = 0;		
		BDP_Simulated(&BDPFinalLength, n, 64, minEventSet, nCommonEvents, Rint, 0, monitorsSizes, "final", bitmapSize, digestLength);
		DigestSizeCommunicationfile.open("DigestSizeCommunication.txt", std::ios_base::app);
		DigestSizeCommunicationfile << digestLength << "," << BDPFinalLength << "\n"; //writing values	
		DigestSizeCommunicationfile.close();
	}
}

void SimulationFalseEventsBREDD()
{
	ofstream myfile;
	int i = 0;
	//long long numberofEvents;
	long long numberofEvents;
	double Rint = 1.0;
	map<int, long long>* monitorsSize = new map<int, long long>();
	int numberofMonitors = 10; //number of monitors
	while (i < numberofMonitors)
	{
		numberofEvents = RandomNumber(100000, 1000000);
		monitorsSize->insert(std::make_pair(i, numberofEvents));
		i++;
	}
	//ten monitors made
	//Now loop upon nodes to find minimum
	long long minEventSet = 1000000; //max possible value
	long long maxEventSet = 0;
	double averageEvents = 0;
	map<int, long long>::iterator it_er;
	for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er)
	{
		averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
		if (minEventSet > it_er->second)
			minEventSet = it_er->second;
		if (maxEventSet < it_er->second)
			maxEventSet = it_er->second;
	}

	long long nCommonEvents = (long long)Rint*minEventSet; //finding number of common events
	long long nEvents90 = 50000;
	long long nEvents80 = 60000;
	long long nEvents70 = 70000;
	//generate dataset
	set<long long>* monitors = 0x00L; //finalMain set to return
	monitors = new set<long long>[numberofMonitors]; //number of monitors
	//first put common values among all monitors
	for (int i = 0; i < nCommonEvents; i++)
	{
		long long commonEvent = RandomNumber(1, 100000000);
		for (int j = 0; j < numberofMonitors; j++)
		{
			monitors[j].insert(commonEvent);
		}
	}
	int nonCommonEvents = 0;
	//Put 0.9 events to 0.9 monitors
	int iterator90 = 0;
	while (iterator90 < nEvents90)
	{
		long long Event90 = RandomNumber(1, 100000000);
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
		{
			if (it_er->second != minEventSet)
			{
				monitors[it_er->first].insert(Event90);
				nonCommonEvents++;
				iterator90++;
			}
		}
	}

	//Put 0.8 events to 0.8 monitors
	int iterator80 = 0;
	while (iterator80 < nEvents80)
	{
		long long Event80 = RandomNumber(1, 100000000);
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
		{
			if (it_er->second != minEventSet && it_er->second != (int)(minEventSet + (800000 / numberofMonitors)))
			{
				monitors[it_er->first].insert(Event80);
				nonCommonEvents++;
				iterator80++;
			}
		}
	}

	//Put 0. events to 0.7 monitors
	int iterator70 = 0;
	while (iterator70 < nEvents70)
	{
		long long Event70 = RandomNumber(1, 100000000);
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
		{
			if (it_er->second != minEventSet && it_er->second != (int)(minEventSet + (800000 / numberofMonitors)) && it_er->second != (int)(minEventSet + 2 * (800000 / numberofMonitors)))
			{
				monitors[it_er->first].insert(Event70);
				nonCommonEvents++;
				iterator70++;
			}
		}
	}


	//Insert non-Common events
	for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
	{
		int remainingIndex = 0;
		if (it_er->second == minEventSet)
		{
			remainingIndex = it_er->second - nCommonEvents;
		}
		else
		{
			if (it_er->second == (minEventSet + (800000 / numberofMonitors)))
			{
				remainingIndex = it_er->second - nCommonEvents - nEvents90;
			}
			else
			{
				if (it_er->second == (minEventSet + 2 * (800000 / numberofMonitors)))
				{
					remainingIndex = it_er->second - nCommonEvents - nEvents90 - nEvents80;
				}
				else
				{
					remainingIndex = it_er->second - nCommonEvents - nEvents90 - nEvents80 - nEvents70;
				}
			}
		}
		for (int i = 0; i < remainingIndex; i++)
		{
			long long NoncommonEvent = RandomNumber(1, 100000000);
			monitors[it_er->first].insert(NoncommonEvent);
			nonCommonEvents++;
		}
	}


	//Find eq. 5 false +ve results

	printf("Widespread events: %lld\n", nCommonEvents);

	int BDPDfalseEvents = 0;
	int BDPD90falseEvents = 0;
	int BDPD80falseEvents = 0;
	int BDPD70falseEvents = 0;
	int digestSize = 10;
	std::set<long long>::iterator iter;
	//double eq5FalseEvents = 0;	
	while(digestSize<=30)
	{
		BDPThetaFalseEvents_Simulated(&BDPDfalseEvents, &BDPD90falseEvents, &BDPD80falseEvents, &BDPD70falseEvents, numberofMonitors, minEventSet, nCommonEvents, 0, 0, 0, Rint, numberofMonitors, averageEvents, digestSize, monitors);
		printf("False: %d, %d, %d, %d \n", BDPDfalseEvents, BDPD90falseEvents, BDPD80falseEvents, BDPD70falseEvents);
		

		printf("Digest: %d, FP: %llf \n", digestSize, BDPDfalseEvents);
		
		/*
		//Find eq. 5 false +ve results
		for (int i = 0; i < numberofMonitors; ++i)
		{
			for (iter = monitors[i].begin(); iter != monitors[i].end(); ++iter)
			{
				eq5FalseEvents += nCommonEvents*((*iter - nCommonEvents) / averageEvents);
			}
		}
		*/
		//double eq5falsePositives = eq5FalseEvents / pow(2, digestSize);
		myfile.open("NumberofFalseEvents.txt", std::ios_base::app);
		myfile << digestSize << "," << BDPDfalseEvents << "," << nonCommonEvents << "\n"; //writing values	
		myfile.close();
		printf("Done..\n");
		digestSize++;
	}
}

void SimulationFalseEvents()
{
	ofstream myfile;
	int digestSize = 10;


	int i = 0;
	//long long numberofEvents;
	long long EventSize = 15000;
	double Rint = 0.1;
	map<int, long long>* monitorsSize = new map<int, long long>();
	int numberofMonitors = 10; //number of monitors
	while (i < numberofMonitors)
	{
		EventSize += (80000 / numberofMonitors);
		monitorsSize->insert(std::make_pair(i, EventSize));
		i++;
	}
	//ten monitors made
	//Now loop upon nodes to find minimum
	long long minEventSet = 100000; //max possible value
	long long maxEventSet = 0;
	double averageEvents = 0;
	map<int, long long>::iterator it_er;
	for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er)
	{
		averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
		if (minEventSet > it_er->second)
			minEventSet = it_er->second;
		if (maxEventSet < it_er->second)
			maxEventSet = it_er->second;
	}

	long long nCommonEvents = Rint*minEventSet; //finding number of common events
	long long nEvents90 = 5000;
	long long nEvents80 = 6000;
	long long nEvents70 = 7000;

	long long nonCommon100 = 0;
	long long nonCommon90 = 0;
	long long nonCommon80 = 0;
	long long nonCommon70 = 0;

	//generate dataset
	set<long long>* monitors = 0x00L; //finalMain set to return
	monitors = new set<long long>[numberofMonitors]; //number of monitors
	//first put common values among all monitors
	for (int i = 0; i < nCommonEvents; i++)
	{
		long long commonEvent = RandomNumber(1, 10000000);
		for (int j = 0; j < numberofMonitors; j++)
		{
			monitors[j].insert(commonEvent);
		}
	}

	//Put 0.9 events to 0.9 monitors
	int iterator90 = 0;
	while (iterator90 < nEvents90)
	{
		long long Event90 = RandomNumber(1, 10000000);
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
		{
			if (it_er->second != minEventSet)
			{
				monitors[it_er->first].insert(Event90);
				nonCommon100++;
				iterator90++;
			}
		}
	}

	//Put 0.8 events to 0.8 monitors
	int iterator80 = 0;
	while (iterator80 < nEvents80)
	{
		long long Event80 = RandomNumber(1, 10000000);
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
		{
			if (it_er->second != minEventSet && it_er->second != (int)(minEventSet + (80000 / numberofMonitors)))
			{
				monitors[it_er->first].insert(Event80);
				nonCommon100++;
				nonCommon90++;
				iterator80++;
			}
		}
	}

	//Put 0. events to 0.7 monitors
	int iterator70 = 0;
	while (iterator70 < nEvents70)
	{
		long long Event70 = RandomNumber(1, 10000000);
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
		{
			if (it_er->second != minEventSet && it_er->second != (int)(minEventSet + (80000 / numberofMonitors)) && it_er->second != (int)(minEventSet + 2 * (80000 / numberofMonitors)))
			{
				monitors[it_er->first].insert(Event70);
				nonCommon100++;
				nonCommon90++;
				nonCommon80++;
				iterator70++;
			}
		}
	}


	//Insert non-Common events
	for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
	{
		int remainingIndex = 0;
		if (it_er->second == minEventSet)
		{
			remainingIndex = it_er->second - nCommonEvents;
		}
		else
		{
			if (it_er->second == (minEventSet + (800000 / numberofMonitors)))
			{
				remainingIndex = it_er->second - nCommonEvents - nEvents90;
			}
			else
			{
				if (it_er->second == (minEventSet + 2 * (800000 / numberofMonitors)))
				{
					remainingIndex = it_er->second - nCommonEvents - nEvents90 - nEvents80;
				}
				else
				{
					remainingIndex = it_er->second - nCommonEvents - nEvents90 - nEvents80 - nEvents70;
				}
			}
		}
		for (int i = 0; i < remainingIndex; i++)
		{
			long long NoncommonEvent = RandomNumber(1, 100000000);
			monitors[it_er->first].insert(NoncommonEvent);
			nonCommon100++;
			nonCommon90++;
			nonCommon80++;
			nonCommon70++;
		}
	}
	while (digestSize <= 30)
	{
		int BDPDfalseEvents = 0;
		int BDPD90falseEvents = 0;
		int BDPD80falseEvents = 0;
		int BDPD70falseEvents = 0;
		
		BDPThetaFalseEvents_Simulated(&BDPDfalseEvents, &BDPD90falseEvents, &BDPD80falseEvents, &BDPD70falseEvents, numberofMonitors, minEventSet, nCommonEvents, nEvents90, nEvents80, nEvents70, Rint, numberofMonitors, maxEventSet, digestSize, monitors);
		double fp100 = (double) BDPDfalseEvents / nonCommon100;
		double fp90 = (double)BDPD90falseEvents / nonCommon90;
		double fp80 = (double)BDPD80falseEvents / nonCommon80;
		double fp70 = (double)BDPD70falseEvents / nonCommon70;

		printf("False: %f, %f, %f, %f \n", fp100, fp90, fp80, fp70);
		myfile.open("NumberofFalseEvents.txt", std::ios_base::app);
		myfile << digestSize << "," << fp100 << "," << fp90 << "," << fp80 << "," << fp70 << "\n"; //writing values	
		myfile.close();
		digestSize++;
		printf("Done..\n");
	}
}

void SimulationCommunicationwrtMonitorTheta()
{
	ofstream myfile;
	double error = 0.001;
	long long averageEventSize = 64;
	for (int i = 2; i <= 52; i++)
	{
		double avgBDPFinal = 0;
		double avgBDPDigested = 0;
		double avgBDP90 = 0;
		double avgBDP80 = 0;
		double avgBDP70 = 0;
		double avgBDPDigested90 = 0;
		double avgBDPDigested80 = 0;
		double avgBDPDigested70 = 0;
		int numberofMonitors = i;
		for (int j = 0; j < 1; j++)
		{			
			map<int, long long>* monitorSize = new map<int, long long>();
			//long long numberofEvents;
			long long EventSize = 150000;
			for (int k = 0; k < numberofMonitors; k++)
			{
				EventSize += (800000 / numberofMonitors);
				//numberofEvents = RandomNumber(100000, 1000000);
				monitorSize->insert(std::make_pair(k, EventSize));
			}
			//Now loop upon nodes to find minimum and maximum
			long long minEventSet = 1000000; //max possible value
			long long maxEventSet = 100000; //min possible value
			double averageEvents = 0;
			map<int, long long>::iterator it_er;
			for (it_er = monitorSize->begin(); it_er != monitorSize->end(); ++it_er)
			{
				averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
				if (minEventSet > it_er->second)
					minEventSet = it_er->second;
				if (maxEventSet < it_er->second)
					maxEventSet = it_er->second;
			}
			long long nCommonEvents = 0.5*minEventSet; //finding number of common events
			double Rint = (double)nCommonEvents / minEventSet; //finding intersection ratio

															   //generate dataset
			set<long long>* monitors = 0x00L; //finalMain set to return
			monitors = new set<long long>[numberofMonitors]; //number of monitors
			//first put common values among all monitors
			int iterator = 0;
			while (iterator < nCommonEvents)
			{
				long long commonEvent = RandomNumber(1, 100000000);
				for (int cmi = 0; cmi < numberofMonitors; cmi++)
				{
					monitors[cmi].insert(commonEvent);
				}
				iterator++;
			}
			//Insert non-Common events
			for (it_er = monitorSize->begin(); it_er != monitorSize->end(); ++it_er) //for each monitor
			{
				int iterator = 0;
				while (iterator < (it_er->second - nCommonEvents))
				{
					long long NoncommonEvent = RandomNumber(1, 100000000);
					monitors[it_er->first].insert(NoncommonEvent);
					iterator++;
				}
			}	

			long double BDPFinalLength = 0;
			long double BDPDigestedLength = 0;
			long double BDP90Length = 0;
			long double BDP80Length = 0;
			long double BDP70Length = 0;
			long double BDPDigested90Length = 0;
			long double BDPDigested80Length = 0;
			long double BDPDigested70Length = 0;
			BDPTheta_Simulated(&BDPFinalLength, &BDPDigestedLength, &BDP90Length, &BDP80Length, &BDP70Length, &BDPDigested90Length, &BDPDigested80Length, &BDPDigested70Length, numberofMonitors, minEventSet, nCommonEvents, Rint, numberofMonitors, averageEvents, 20, monitors);
			avgBDPFinal = (avgBDPFinal*j + BDPFinalLength) / (j + 1);
			avgBDPDigested = (avgBDPDigested*j + BDPDigestedLength) / (j + 1);
			avgBDP90 = (avgBDP90*j + BDP90Length) / (j + 1);
			avgBDP80 = (avgBDP80*j + BDP80Length) / (j + 1);
			avgBDP70 = (avgBDP70*j + BDP70Length) / (j + 1);
			avgBDPDigested90 = (avgBDPDigested90*j + BDPDigested90Length) / (j + 1);
			avgBDPDigested80 = (avgBDPDigested80*j + BDPDigested80Length) / (j + 1);
			avgBDPDigested70 = (avgBDPDigested70*j + BDPDigested70Length) / (j + 1);
			printf("j: %d \n", j);			
			monitorSize->clear();			
			monitors->clear();			
		}
		//write
		myfile.open("CommunicationoverheadwrtMonitorsTheta.txt", std::ios_base::app);
		myfile << i << "," << avgBDPFinal << "," << avgBDPDigested << "," << avgBDP90 << "," << avgBDP80 << "," << avgBDP70 << "," << avgBDPDigested90 << "," << avgBDPDigested80 << "," << avgBDPDigested70 << "\n"; //writing values	
		myfile.close();
		printf("%d \n", i);		
	}
}


void SimulationCommunicationwrtRintTheta()
{
	ofstream myfile;
	double Rint = 0;
	long long averageEventSize = 64;
	while (Rint <= 1.0)
	{

		int i = 0;
		//long long numberofEvents;
		long long EventSize = 150000;
		map<int, long long>* monitorsSize = new map<int, long long>();
		int numberofMonitors = 10; //number of monitors
		while (i < numberofMonitors)
		{
			EventSize += (800000 / numberofMonitors);
			monitorsSize->insert(std::make_pair(i, EventSize));
			i++;
		}
		//ten monitors made
		//Now loop upon nodes to find minimum
		long long minEventSet = 1000000; //max possible value
		double averageEvents = 0;
		map<int, long long>::iterator it_er;
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er)
		{
			averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
			if (minEventSet > it_er->second)
				minEventSet = it_er->second;
		}

		long long nCommonEvents = Rint*minEventSet; //finding number of common events		

			//generate dataset
		set<long long>* monitors = 0x00L; //finalMain set to return
		monitors = new set<long long>[numberofMonitors]; //number of monitors
		//first put common values among all monitors
		for (int i = 0; i < nCommonEvents; i++)
		{
			long long commonEvent = RandomNumber(1, 100000000);
			for (int j = 0; j < numberofMonitors; j++)
			{
				monitors[j].insert(commonEvent);
			}
		}
		//Insert non-Common events
		for (it_er = monitorsSize->begin(); it_er != monitorsSize->end(); ++it_er) //for each monitor
		{
			for (int i = 0; i < (it_er->second - nCommonEvents); i++)
			{
				long long NoncommonEvent = RandomNumber(1, 100000000);
				monitors[it_er->first].insert(NoncommonEvent);
			}
		}
		long double BDPFinalLength = 0;
		long double BDPDigestedLength = 0;
		long double BDP90Length = 0;
		long double BDP80Length = 0;
		long double BDP70Length = 0;
		long double BDPDigested90Length = 0;
		long double BDPDigested80Length = 0;
		long double BDPDigested70Length = 0;
		BDPTheta_Simulated(&BDPFinalLength, &BDPDigestedLength, &BDP90Length, &BDP80Length, &BDP70Length, &BDPDigested90Length, &BDPDigested80Length, &BDPDigested70Length, numberofMonitors, minEventSet, nCommonEvents, Rint, numberofMonitors, averageEvents, 20, monitors);
		myfile.open("SimulatedAggregatedResultsTheta.txt", std::ios_base::app);
		myfile << Rint << "," << BDPFinalLength << "," << BDPDigestedLength << "," << BDP90Length << "," << BDP80Length << "," << BDP70Length << "," << BDPDigested90Length << "," << BDPDigested80Length << "," << BDPDigested70Length << "\n"; //writing values	
		myfile.close();
		Rint += 0.05;
		printf("Done..\n");
	}
}

void SimulationCommunicationwrtRint()
{
	ofstream myfile;
	double RintRange = 0;
	long long averageEventSize = 64;
	//for (int outer = 0; outer < 20; outer++)
	while (RintRange<1)
	{
		int iteration = 0;
		double avgRint = 0;
		double avgRaw = 0;
		/*
		double avgCFSCompact = 0;
		double avgCFSReal = 0;
		double avgCFSPCompact = 0;
		*/
		double avgCFSPReal = 0;
		/*
		double avgBDPOriginal = 0;
		double avgBDPRenewed = 0;
		double avgBDPExtended = 0;
		*/
		double avgBDPFinal = 0;
		//double avgBDPDigested = 0;
		//For simulation environment 
		//for (int iteration = 0; iteration < 1000; iteration++)
		while (iteration <1000)
		{
			int i = 0;
			long long numberofEvents;
			map<int, long long>* monitors = new map<int, long long>();
			int n = 10; //number of monitors
			while (i < n)
			{
				/*
				r1 = rand() % 1000;
				r2 = rand() % 1000;
				numberofEvents = r1 * 1000 + r2;
				if (numberofEvents >= 100000) //if value is less than 100,000 it will not insert
				{
				i++;
				monitors->insert(std::make_pair(i, numberofEvents));
				}
				*/
				numberofEvents = RandomNumber(100000, 1000000);
				monitors->insert(std::make_pair(i, numberofEvents));
				i++;
			}
			//ten monitors made
			//Now loop upon nodes to find minimum
			long long minEventSet = 1000000; //max possible value
			map<int, long long>::iterator it_er;
			for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
			{
				if (minEventSet > it_er->second)
					minEventSet = it_er->second;
			}

			long long nCommonEvents = RandomNumber(0, minEventSet); //finding number of common events	
			double Rint = (double)nCommonEvents / minEventSet; //finding intersection ratio			
			if (Rint >= RintRange && Rint <= RintRange + 0.05)
			{
				//printf("Iteration: %d , Rint: %G \n", iteration, Rint);
				double error = 0.01;
				long long RawDataLength = 0;
				//Calculate Rawdata
				RawData_Simulated(&RawDataLength, n, averageEventSize, nCommonEvents, monitors);

				/*
				long long BDPOriginalLength = 0;
				BDP_Simulated(&BDPOriginalLength, n, averageEventSize, minEventSet, nCommonEvents, Rint, monitors, "original");
				long long BDPRenewedLength = 0;
				BDP_Simulated(&BDPRenewedLength, n, averageEventSize, minEventSet, nCommonEvents, Rint, monitors, "renewed");
				long long BDPExtendedLength = 0;
				BDP_Simulated(&BDPExtendedLength, n, averageEventSize, minEventSet, nCommonEvents, Rint, monitors, "extended");
				*/
				long long BDPFinalLength = 0;
				BDP_Simulated(&BDPFinalLength, n, averageEventSize, minEventSet, nCommonEvents, Rint, 0, monitors, "final");
				//long long BDPDigestedLength = 0;
				//BDP_Simulated(&BDPDigestedLength, n, averageEventSize, minEventSet, nCommonEvents, Rint, 0, monitors, "final", 0, 20);

				/*
				long long CFSCompactLength = 0;
				//calculate CSF
				CFS_Simulated(&CFSCompactLength, n, minEventSet, nCommonEvents, Rint, error, monitors, "compact");
				long long CFSRealLength = 0;
				//calculate CSF
				CFS_Simulated(&CFSRealLength, n, minEventSet, nCommonEvents, Rint, error, monitors, "real");
				*/

				/*
				long long CFSPCompactLength = 0;
				//calculate CSFP
				CFSP_Simulated(&CFSPCompactLength, n, minEventSet, nCommonEvents, Rint, error, monitors, "compact");
				*/
				long long CFSPRealLength = 0;
				//calculate CSFP
				CFSP_Simulated(&CFSPRealLength, n, minEventSet, nCommonEvents, Rint, error, monitors, "real");

				//printf("File opened..\n");
				//Calculating runtime averages
				avgRint = (avgRint*iteration + Rint) / (iteration + 1);
				avgRaw = (avgRaw*iteration + RawDataLength) / (iteration + 1);
				//avgCFSCompact = (avgCFSCompact*iteration + CFSCompactLength) / (iteration + 1);
				//avgCFSReal = (avgCFSReal*iteration + CFSRealLength) / (iteration + 1);
				//avgCFSPCompact = (avgCFSPCompact*iteration + CFSPCompactLength) / (iteration + 1);
				avgCFSPReal = (avgCFSPReal*iteration + CFSPRealLength) / (iteration + 1);
				//avgBDPOriginal = (avgBDPOriginal*iteration + BDPOriginalLength) / (iteration + 1);
				//avgBDPRenewed = (avgBDPRenewed*iteration + BDPRenewedLength) / (iteration + 1);
				//avgBDPExtended = (avgBDPExtended*iteration + BDPExtendedLength) / (iteration + 1);
				avgBDPFinal = (avgBDPFinal*iteration + BDPFinalLength) / (iteration + 1);
				//avgBDPDigested = (avgBDPDigested*iteration + BDPDigestedLength) / (iteration + 1);
				//myfile2.open("SimulatedResults.txt", std::ios_base::app);
				//myfile2 << Rint << "," << RawDataLength << "," << CFSCompactLength << "," << CFSPCompactLength << "," << BDPOriginalLength << "," << minEventSet << "," << nCommonEvents << "\n"; //writing values	
				//myfile2.close();
				iteration++;
			}

		}
		RintRange += 0.05;
		myfile.open("SimulatedAggregatedResults.txt", std::ios_base::app);
		myfile << avgRint << "," << avgRaw << "," << avgCFSPReal << "," << avgBDPFinal << "\n"; //"," << avgBDPDigested << "\n"; //writing values	
		myfile.close();
	}
}

void SimulationCommunicationwrtMonitors()
{
	ofstream myfile;
	double error = 0.001;
	long long averageEventSize = 64;
	for (int i = 2; i <= 100; i++)
	{
		double avgRaw = 0;
		double avgCFSPReal = 0;
		double avgBDPFinal = 0;
		double avgBDPDigested = 0;

		for (int j = 0; j < 1000; j++)
		{
			map<int, long long>* monitors = new map<int, long long>();
			long long numberofEvents;
			for (int k = 0; k < i; k++)
			{
				numberofEvents = RandomNumber(100000, 1000000);
				monitors->insert(std::make_pair(k, numberofEvents));
			}
			//Now loop upon nodes to find minimum and maximum
			long long minEventSet = 1000000; //max possible value
			long long maxEventSet = 100000; //min possible value
			double averageEvents = 0;
			map<int, long long>::iterator it_er;
			for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
			{
				averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
				if (minEventSet > it_er->second)
					minEventSet = it_er->second;
				if (maxEventSet < it_er->second)
					maxEventSet = it_er->second;
			}
			long long nCommonEvents = 0.5*minEventSet; //finding number of common events
			double Rint = (double)nCommonEvents / minEventSet; //finding intersection ratio

			long long RawDataLength = 0;
			//Calculate Rawdata
			RawData_Simulated(&RawDataLength, i, averageEventSize, nCommonEvents, monitors);

			long long BDPFinalLength = 0;
			BDP_Simulated(&BDPFinalLength, i, averageEventSize, minEventSet, nCommonEvents, Rint, 0, monitors, "final");
			long long BDPDigestedLength = 0;
			BDP_Simulated(&BDPDigestedLength, i, averageEventSize, minEventSet, nCommonEvents, Rint, 0, monitors, "final", 0, 20);

			long long CFSPRealLength = 0;
			//calculate CSFP
			CFSP_Simulated(&CFSPRealLength, i, minEventSet, nCommonEvents, Rint, error, monitors, "real");


			avgRaw = (avgRaw*j + RawDataLength) / (j + 1);
			avgCFSPReal = (avgCFSPReal*j + CFSPRealLength) / (j + 1);
			avgBDPFinal = (avgBDPFinal*j + BDPFinalLength) / (j + 1);
			avgBDPDigested = (avgBDPDigested*j + BDPDigestedLength) / (j + 1);
		}
		//write
		myfile.open("CommunicationoverheadwrtMonitors.txt", std::ios_base::app);
		myfile << i << "," << avgRaw << "," << avgCFSPReal << "," << avgBDPFinal << "," << avgBDPDigested << "\n"; //writing values	
		myfile.close();
	}

}

void SimulationAccuracyComparison()
{
	ofstream AccuracyComparisonfile;
	int i = 0;
	long long numberofEvents;
	long long averageEventSize = 64;	
	double aRangeErrorAvg = 0;
	double bRangeErrorAvg = 0;
	double cRangeErrorAvg = 0;
	double dRangeErrorAvg = 0;
	double eRangeErrorAvg = 0;
	int aCount=0, bCount=0, cCount=0, dCount=0, eCount = 0;
	long long BDPDigestedLength;
	//long long CFSPRealLength;
	double communicationOverhead;
	int iteration = 0;
	while ((aCount < 1000 || bCount < 1000 || cCount < 1000 || dCount < 1000 || eCount < 100) && (iteration < 100))
	{
		map<int, long long>* monitors = new map<int, long long>();
		int n = 10; //number of monitors
		while (i < n)
		{
			numberofEvents = RandomNumber(100000, 1000000);
			monitors->insert(std::make_pair(i, numberofEvents));
			i++;
		}

		//Now loop upon nodes to find minimum and maximum
		long long minEventSet = 1000000; //max possible value
		long long maxEventSet = 100000; //min possible value
		double averageEvents = 0;
		map<int, long long>::iterator it_er;
		for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
		{
			averageEvents = (averageEvents*it_er->first + it_er->second) / (it_er->first + 1);
			if (minEventSet > it_er->second)
				minEventSet = it_er->second;
			if (maxEventSet < it_er->second)
				maxEventSet = it_er->second;
		}

		long long nCommonEvents = 0.5*minEventSet; //finding number of common events
		double Rint = (double)nCommonEvents / minEventSet; //finding intersection ratio
		long double error = 0.01;
		printf("Iteration: %d \n\n", iteration);
		do
		{
			error /= 10;
			BDPDigestedLength = 0;
			//CFSPRealLength = 0;
			//calculate CSFP
			//CFSP_Simulated(&CFSPRealLength, n, minEventSet, nCommonEvents, Rint, error, monitors, "real");
			BDP_Simulated(&BDPDigestedLength, n, averageEventSize, minEventSet, nCommonEvents, Rint, error, monitors, "final");
			communicationOverhead = BDPDigestedLength / (1024 * 1024);
			//printf("CFSP: %llf, Error: %llf \n", communicationOverhead, error);
			if (communicationOverhead >= 30 && communicationOverhead <= 40 && aCount<1000)
			{
				aRangeErrorAvg = (aRangeErrorAvg*aCount + error) / (aCount + 1);
				aCount++;
			}
			else if (communicationOverhead >= 40 && communicationOverhead <= 50 && bCount<1000)
			{
				bRangeErrorAvg = (bRangeErrorAvg*bCount + error) / (bCount + 1);
				bCount++;
			}
			else if (communicationOverhead >= 50 && communicationOverhead <= 60 && cCount<1000)
			{
				cRangeErrorAvg = (cRangeErrorAvg*cCount + error) / (cCount + 1);
				cCount++;
			}
			else if (communicationOverhead >= 60 && communicationOverhead <= 70 && dCount<1000)
			{
				dRangeErrorAvg = (dRangeErrorAvg*dCount + error) / (dCount + 1);
				dCount++;
			}
			else if (communicationOverhead >= 70 && communicationOverhead <= 80 && eCount<1000)
			{
				eRangeErrorAvg = (eRangeErrorAvg*eCount + error) / (eCount + 1);
				eCount++;
			}
			printf("aCount: %d, bCount: %d, cCount: %d, dCount: %d, eCount: %d \n", aCount, bCount, cCount, dCount, eCount);
		} while ((error > 0) && (aCount < 1000 || bCount < 1000 || cCount < 1000 || dCount < 1000 || eCount < 100));
		iteration++;
	}

	getchar();
	//write
	AccuracyComparisonfile.open("AccuracyComparisonFile.txt", std::ios_base::app);	
	printf("a: %llf, b: %llf, c:%llf, d=%llf, e=%llf \n", aRangeErrorAvg, bRangeErrorAvg, cRangeErrorAvg, dRangeErrorAvg, eRangeErrorAvg);
	AccuracyComparisonfile << aRangeErrorAvg << "\n" << bRangeErrorAvg << "\n" << cRangeErrorAvg << "\n" << dRangeErrorAvg << "\n" << eRangeErrorAvg; //writing values	
	AccuracyComparisonfile.close();
}

void SimulationforFalsePositivewrtCommunicationOverhead()
{
	ofstream FalsePositivewrtCommunicationfile;
	long long averageEventSize = 64;
	int numberofmonitors = 100; //number of monitors	
	double avgBDP = 0;
	double avgCFSP = 0;
	long long numberofEvents;
	double error = 1.0;
	for (int errorloop = 0; errorloop < 10; errorloop++)
	{
		int iteration = 0;
		error /= 10.0;
		while (iteration < 1000)
		{
			int i = 0;
			map<int, long long>* monitors = new map<int, long long>();
			while (i < numberofmonitors)
			{
				numberofEvents = RandomNumber(100000, 1000000);
				monitors->insert(std::make_pair(i, numberofEvents));
				i++;
			}
			//monitors made
			//Now loop upon nodes to find minimum
			long long minEventSet = 1000000; //max possible value
			map<int, long long>::iterator it_er;
			for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
			{
				if (minEventSet > it_er->second)
					minEventSet = it_er->second;
			}

			long long nCommonEvents = 0.5*minEventSet; //finding number of common events	
			double Rint = (double)nCommonEvents / minEventSet; //finding intersection ratio

			long long BDPFinalLength = 0;
			BDP_Simulated(&BDPFinalLength, numberofmonitors, averageEventSize, minEventSet, nCommonEvents, Rint, 0, monitors, "final");

			long long CFSPRealLength = 0;
			//calculate CSFP
			CFSP_Simulated(&CFSPRealLength, i, minEventSet, nCommonEvents, Rint, error, monitors, "real");

			avgBDP = (avgBDP*iteration + BDPFinalLength) / (iteration + 1);
			avgCFSP = (avgCFSP*iteration + CFSPRealLength) / (iteration + 1);
			iteration++;
		}
		FalsePositivewrtCommunicationfile.open("FalsePositivewrtCommunication.txt", std::ios_base::app);
		printf("Error = %llf, BRED = %llf, CFSP = %llf \n", error, avgBDP, avgCFSP);
		FalsePositivewrtCommunicationfile << error << "," << avgBDP << "," << avgCFSP << "\n"; //writing values	
		FalsePositivewrtCommunicationfile.close();
	}	
}









/*
double e = (100 - m_PERCENTAGE) / 100.0;
switch (scheme)
{
case 'c': //CFP (combinable filters is in place
m_vec_len = (int)pow(2, ceil(LOG2(m_ELEM_NUM*(-log(e)/(log(2)*log(2)))))); //formula for CFP
m_hash_num = (int)ceil(-log(e)/log(2));
break;
case 'p': //CFPS (progressive filtering) is in place
m_vec_len = (int)pow(2, ceil(LOG2(m_ELEM_NUM/log(2)))); //formula for CFPS
m_hash_num = 1;
m_PERCENTAGE = 50;
break;
*/


/*
//Determining the probability of all the bits to be one together..
for (it_er = monitors->begin(); it_er != monitors->end(); ++it_er)
{

double currentProbability = (1 - pow(1 - (1.0 / bitmapSize), it_er->second));
//printf("Bitmap size: %d, Events: %lld, Current Probability: %f \n", bitmapSize, it_er->second, currentProbability);
combinedProbability *= currentProbability;
}

long long numberofIndexes = (long long)bitmapSize*combinedProbability; //total number of indexes that coordinator will send to each monitor

printf("Common E: %lld, Combined prob: %f, Coordinator Events: %lld \n", nCommonEvents, combinedProbability, numberofIndexes);
myfile.open("CombinedandProbable.txt", std::ios_base::app);
myfile << Rint << "," << nCommonEvents << "," << numberofIndexes << "\n"; //writing values
myfile.close();
*/