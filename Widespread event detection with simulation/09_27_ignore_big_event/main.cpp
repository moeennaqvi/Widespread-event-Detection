#include <iostream>
#include <set>
#include <map>
#include <math.h>
#include <stdint.h>
#include <fstream>
#include "PREPROCESS.h"
#include "BLOOMFILTER.h"
#include "CombinableFilter.h"
#include <cstdlib>
#pragma warning(disable:4996)
//#include "openssl\md5.h"


int main(int argc, char* argv[]) 
{
	/*---------------------------------------------------------------------*/
	/* numberOfNodes : */
	/* thetaRatio : ratio of the node we want to check */
	ofstream myfile;
	ofstream myfile2;
	int numberOfNodes = atoi(argv[1]);
	int bitmapSize = atoi(argv[2]);
	double thetaRatio = atof(argv[3]);
	int option = atoi(argv[4]);
	int averageEventSize = atoi(argv[5]);
	int alpha = atoi(argv[6]);
	int acceptRatio = (int)((double)numberOfNodes * (double)thetaRatio);
	double totalMemorySize = 0;
	int numberOfDst = 0;
	long long* dstArray = 0x00L;
	std::set<long long>* finalSet = 0x00L;
	numberOfNodes = 2;
	//for experiment
	bitmapSize = 1024;
	

	//SimulationFalseEventsBREDD();
	SimulationforFalsePositivewrtCommunicationOverhead();
	//SimulationforBitmapSize();
	//SimulationCommunicationwrtRint();
	//SimulationCommunicationwrtMonitors();
	//SimulationAccuracyComparison();
	//SimulationCommunicationwrtRintTheta();
	//SimulationFalseEvents();
	//SimulationCommunicationwrtMonitorTheta();
	getDstSet(&numberOfDst, &dstArray, &finalSet, option);
	printf("%lld \n", numberOfDst);
	printf("ended and file close\n");
	//for calculating number of monitors for each events
	/*
	map<long long, int>* eventsCount = new map<long long, int>();
	for (int i = 0; i < numberOfDst; i++)
	{
		set<long long>::iterator iter;
		for (iter = finalSet[i].begin(); iter != finalSet[i].end(); ++iter)
		{
			if (eventsCount->find(*iter) == eventsCount->end()) 
			{
				eventsCount->insert(std::make_pair(*iter, 1));
			}
			else
			{
				eventsCount->find(*iter)->second+=1;
			}
		}
	}
	printf("Data into map..\n");
	map<long long, int>::iterator it_er;
	printf("File opened..\n");
	myfile2.open("monitorsPerEventtoptenthousand.txt", std::ios_base::app);
	for (it_er = eventsCount->begin(); it_er != eventsCount->end(); ++it_er)
	{
		if(it_er->second > 30)
		myfile2 << it_er->first << "," << it_er->second << "\n"; //writing values
	}
	myfile2.close();
	printf("ended and file close\n");
	*/

	//for calculating number of events for each monitor
	/*	
	printf("File opened..\n");
	myfile2.open("eventpernode.txt", std::ios_base::app);
	for (int i = 0; i < numberOfDst; i++)
	{		
		myfile2 << i + 1 << "," << finalSet[i].size() << "\n"; //writing values		
	}
	myfile2.close();
	printf("ended and file close\n");
	*/
	
	//for calculating monitor event scatter plot
	/*	
	printf("File opened..\n");
	myfile2.open("eventmonitorscatter.txt", std::ios_base::app);
	for (int i = 0; i < numberOfDst; i++)
	{
		set<long long>::iterator iter;
		for (iter = finalSet[i].begin(); iter != finalSet[i].end(); ++iter)
		{

		myfile2 << i+1 << "," << *iter << "\n"; //writing values

		}		
	}
	myfile2.close();
	printf("ended and file close\n");
	*/
	
	for (numberOfNodes; numberOfNodes <= 1024; numberOfNodes++)
	{
		acceptRatio = (int)((double)numberOfNodes * (double)thetaRatio); // for temporary
		std::cout << "m : " << bitmapSize << " , " << "k : " << numberOfNodes << " , " << "Theta : " << thetaRatio << " , " << "alpha : " << alpha << std::endl;
		//printf(".\nChecking if the paramters fullfil requirement...\n");
		if (numberOfDst < numberOfNodes) //if agents we got from data are less than pre-defined nodes.  
		{
			fprintf(stderr, "Need more node\n");
			return 0;
		}
		if (argc < 5) //if number of arguments passed by main function is less than 5
		{
			fprintf(stderr, "Usage: [numberOfNodes] [bitmapSize] [thetaRatio] [option] [averageEventSize]\n");
			return 0;
		}

		long long totalCFPdata = 0;
		//Calling CFS Chen's scheme
		set<long long>* widespreadEventSetCFS = 0x00L;
		CFS(&totalCFPdata, numberOfNodes, 'c', dstArray, finalSet, &widespreadEventSetCFS);
		//printf("\nCFS: %lld\n", totalCFPdata);

		//calling CFPS chen's scheme
		long long totalCFPSdata = 0;
		set<long long>* widespreadEventSetCFPS = 0x00L;
		CFPS(&totalCFPSdata, numberOfNodes, dstArray, finalSet, &widespreadEventSetCFPS);
		//printf("\nCFPS: %lld\n", totalCFPSdata);





		std::map<long long, int> check; //contains sourceIP (events) on first value and its count on second
		int realEventCount = 0;
		long long totalEventCount = 0;
		//printf(".\nExtracting all distinct events and thier counts from raw dataset..\n");
		int maxEventByAgent = 0;
		for (int i = 0; i < numberOfNodes; ++i) //outer loop iterates one complete Monitoring agent (destination IP) 
		{
			int temp_count = 0;
			std::set<long long>::iterator iter;
			totalEventCount += finalSet[i].size();

			for (iter = finalSet[i].begin(); iter != finalSet[i].end(); ++iter) //inner loop for each event in a particular agent (Source IP). Creates a new value in Hash table for each new event and increases its count if its found already
			{
				temp_count++;
				if (check.find(*iter) == check.end())
				{
					check[*iter] = 1;
				}
				else
				{
					check[*iter] += 1;
				}
			}
			if (temp_count > maxEventByAgent)
				maxEventByAgent = temp_count;
		}

		//printf("\nMax event by an agent is: %d\n", maxEventByAgent);
		//printf("\nTotal events are: %lld\n", totalEventCount); //total number of events
		std::map<long long, int>::iterator it_er;
		for (it_er = check.begin(); it_er != check.end(); ++it_er) //This will count total number of real events that occur. 
		{
			if (it_er->second >= acceptRatio) //if a specific event (sourceIP) has more occurances than threshold, its an event. 
			{
				++realEventCount;
			}
		}
		//printf("\nEvents that crosses %d threshold are %d...\n", acceptRatio, realEventCount);
		char** bitmapArray = 0x00L; //arrays of each monitor
		int* coordinator = 0x00L; //int array used to store coordinator array values after all bitmaps are received..  

		//printf("\nInitializing bitmap and coordinator..\n");
		bitmapArray = (char**)malloc(sizeof(char*) * numberOfNodes);
		coordinator = (int*)malloc(sizeof(int) * bitmapSize);
		memset(coordinator, 0, sizeof(int)*bitmapSize); //initialize coordinator with 0 values across all 

		//printf("initialize each bitmap with 0..\n");
		for (int i = 0; i < numberOfNodes; ++i)
		{
			bitmapArray[i] = 0x00L;
			bitmapArray[i] = (char*)malloc(sizeof(char) * bitmapSize); //allocate memory to each bitmap
			memset(bitmapArray[i], '0', sizeof(char) * bitmapSize); //initialize each bitmap with '0'
		}
		/*
		step 1. encoding all events in the bitmap
		*/
		//printf("\nPutting values in bitmap after taking hash of each value..\n");
		for (int i = 0; i < numberOfNodes; ++i)
		{
			std::set<long long>::iterator iter;
			for (iter = finalSet[i].begin(); iter != finalSet[i].end(); ++iter)
			{
				int hashValue = getHashValue(*iter, bitmapSize); //getting hash function's return.. 
				bitmapArray[i][hashValue] = 1; //Make a certain bit as 1 in 2D bitmap array (mapping bits to bit array)   
			}
		}
		/* step 3. count all bits in the bitmap for BitwiseAdd */
		//printf("\nIncrementing coordinator according to multiple occurances in bitmap..\n");
		for (int i = 0; i < bitmapSize; ++i) //outer loop for each bitmap
		{
			for (int j = 0; j < numberOfNodes; ++j) //inner loop iterates for each item of each bitmap and check if it's one to store number of appearences in Coordinator
			{
				if (bitmapArray[j][i] == 1) ++coordinator[i];
			}
		}
		/* m * k */
		totalMemorySize = bitmapSize * numberOfNodes;
		double temp = totalMemorySize;
		//printf("\nstep 2 Memory size: \t%lf\n", totalMemorySize);
		//printf("%lf ", totalMemorySize);
		/* step4. check all bits in coordinator*/
		//printf("\nChecking all bytes in coordinator and calculating memory...\n");
		int count = 0;
		for (int i = 0; i < bitmapSize; ++i)
		{
			if (acceptRatio <= coordinator[i])
			{
				++count;
				for (int j = 0; j < numberOfNodes; ++j)
				{
					if (bitmapArray[j][i] == 1)
						totalMemorySize += LOG2(bitmapSize);
				}
			}
		}
		//std::cout << "total count : " << count << std::endl;
		//printf("\nstep 4 Memory size: \t%lf\n", (totalMemorySize - temp));
		//printf("%lf ", (totalMemorySize - temp));
		temp = totalMemorySize;
		//printf("\nTotal communication upto step 4 is: %lf..\n", totalMemorySize);
		long long BDPTotal = totalMemorySize;

		/* get events in coordinator array */
		std::set<int> overNode;
		std::map<long long, int> realEventMap;
		std::set<long long> tmpSet;
		int cArray[1025];
		memset(cArray, 0, sizeof(int) * 1025);

		for (int i = 0; i < bitmapSize; ++i) // for each bit
		{
			if (acceptRatio <= coordinator[i]) //those index that have widespread event
			{
				for (int j = 0; j < numberOfNodes; ++j) //for every monitor
				{
					if ((cArray[j] == 0) && (bitmapArray[j][i] == 1)) //if bit is one and cArray has zero
					{
						int numberOfEvents = 0;
						std::set<long long>::iterator iter;

						for (iter = finalSet[j].begin(); iter != finalSet[j].end(); ++iter) //iterating through finalset of each node
						{
							int hashValue = getHashValue(*iter, bitmapSize);

							if (hashValue == i)
							{
								++numberOfEvents;
							}
						}
						if (numberOfEvents < alpha)
						{
							totalMemorySize += averageEventSize * numberOfEvents;
							std::set<long long>::iterator it_Er;
							for (it_Er = finalSet[j].begin(); it_Er != finalSet[j].end(); ++it_Er)
							{
								int hashValue = getHashValue(*it_Er, bitmapSize);
								if (hashValue == i)
								{
									if (realEventMap.find(*it_Er) == realEventMap.end()) realEventMap[*it_Er] = 1;
									else realEventMap[*it_Er] += 1;
								}
							}
						}
						else
						{
							overNode.insert(j);
							if (cArray[j] == 0) cArray[j] = 1;
						}
					}
				}
			}
		}

		//printf("step 6 : \t%lf\n", (totalMemorySize - temp));
		//printf("%lf ", (totalMemorySize-temp));
		temp = totalMemorySize;
		int overCount = 0;
		std::map<long long, int>::iterator iter;
		std::set<long long> estimatedEventSet;

		/* check recieved events (in each monitoring node) */
		for (iter = realEventMap.begin(); iter != realEventMap.end(); ++iter)
		{
			if (iter->second >= (acceptRatio - overNode.size()))
			{
				++overCount;
				estimatedEventSet.insert(iter->first);
			}
		}

		totalMemorySize += (averageEventSize * overCount * overNode.size()) + overNode.size();
		std::set<int>::iterator overIter;
		for (overIter = overNode.begin(); overIter != overNode.end(); ++overIter)
		{

			std::set<long long>::iterator srcIter;
			for (srcIter = finalSet[*overIter].begin(); srcIter != finalSet[*overIter].end(); ++srcIter)
			{
				std::set<long long>::iterator esIter;
				for (esIter = estimatedEventSet.begin(); esIter != estimatedEventSet.end(); ++esIter)
				{
					if (*srcIter == *esIter) realEventMap[*esIter] += 1;
				}
			}
		}
		int estimatedCount = 0;
		for (iter = realEventMap.begin(); iter != realEventMap.end(); ++iter)
		{
			if (iter->second >= acceptRatio)
			{
				++estimatedCount;
			}
		}
		//printf("step 8 : \t%lf\n", (totalMemorySize - temp));
		//printf("%lf ", (totalMemorySize-temp));
		//printf("Total memory size : %lf\n", totalMemorySize);
		//if (realEventCount != estimatedCount) printf("%d %d False\n", realEventCount, estimatedCount);
		//else printf("ture\n");
		//printf("Total overhead : %lf\n", totalMemorySize);
		printf("\nBDP: %lld\t CFS: %lld\t CFPS: %lld\n\n", BDPTotal, totalCFPdata, totalCFPSdata);
		myfile.open("output.txt", std::ios_base::app);
		myfile << numberOfNodes<<"\t"<<BDPTotal << "\t" << totalCFPdata << "\t" << totalCFPSdata << "\n"; //writing values
		myfile.close();
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
		/*
		if (widespreadEventSetCFS != 0x00L)
		{
			free(widespreadEventSetCFS);
			widespreadEventSetCFS = 0x00L;
		}
		if (widespreadEventSetCFPS != 0x00L)
		{
			free(widespreadEventSetCFPS);
			widespreadEventSetCFPS = 0x00L;
		}
		*/
	}
		if (dstArray != 0x00L)
		{
			delete[] dstArray;
			dstArray = 0x00L;
		}
		if (finalSet != 0x00L)
		{
			delete[] finalSet;
			finalSet = 0x00L;
		}
		
	
	
	return 0;
}