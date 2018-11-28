#include <vector>
#include <cmath>
#include <cstdlib>
#include <set>
#include <map>
#include <inttypes.h>

#include <openssl/md5.h>
#include "PREPROCESS.h"
using namespace std;

class BloomFilter{
private:
		int m_ELEM_NUM;     // 
		double m_PERCENTAGE;// false positive percentage
		int m_hash_num;     // the number of hash function
		int m_vec_len;      // the length of the bloomfilter bit array
		vector<bool> m_vec; // the bloomfilter bit array
public:
	BloomFilter(int elem_num, double percentage) :m_ELEM_NUM(elem_num), m_PERCENTAGE(percentage)
	{
		init();
	}

	void init()
	{
		m_vec_len = (int)abs(m_ELEM_NUM*log(m_PERCENTAGE)) / (log(2)*log(2)) + 1;
		m_hash_num = (int)(log(2)*((double)m_vec_len / m_ELEM_NUM));
		m_vec.reserve(m_vec_len);
		m_vec.resize(m_vec_len);
		for (size_t i = 0; i < m_vec.size(); ++i)
		{
			m_vec[i] = false;
		}
	}
		
		/*
		 *@func: reset bloom filer to initialize status.
		 */
		void resetBloomFilter()
		{
			for (size_t i = 0; i < m_vec.size(); ++i)
			{
				m_vec[i] = false;
			}
		}
		
		/*
		 * @func: insert element which type is int into the bloom filter
		 * @params: elem {unsigned int} the elemnt inserted into bloom filter
		 */
		void insertElementByMD5(unsigned int elem)
		{			
			unsigned char hashedValue[MD5_DIGEST_LENGTH];
			getHashedValue(elem, hashedValue);
			uint64_t frontResultValue = 0, endResultValue = 0;
			frontResultValue = *((uint64_t*)hashedValue); //first part of hash
			endResultValue = *((uint64_t*)(hashedValue + 8)); //last part of hash
			
			//printf("\nfront value is %" PRIu64 "\n", frontResultValue);
			//printf("\nEnd value is %" PRIu64 "\n", endResultValue);			

			for (int i = 0; i < m_hash_num; i++)
			{
				int index = 0; //initialize index
				if (i == 0)
				{
					index = frontResultValue % m_vec_len;
				}
				else if (i == 1)
				{
					index = endResultValue % m_vec_len;
				}
				else
				{
					index = (frontResultValue + i*endResultValue) % m_vec_len;
				}

				//printf("\nIndex value for %dth hash is:%d\n ", i, index);
				m_vec[index] = true;
			}
			//printf("Value inserted in BF successfully");
		}
		/*
		 *@func: check the element which type is int locate in the bloom filter
		 *@param: elem {int}, the element to be found 
		 *@ret: true, the element is in the bloom filter
		 *      false, the element is not in the bloom filter
		 */
		
		bool isFoundElement(int elem)
		{
			unsigned char hashedValue[MD5_DIGEST_LENGTH];
			getHashedValue(elem, hashedValue);
			uint64_t frontResultValue = 0, endResultValue = 0;
			frontResultValue = *((uint64_t*)hashedValue); //first part of hash
			endResultValue = *((uint64_t*)(hashedValue + 8)); //last part of hash
			for (int i = 0; i < m_hash_num; i++)
			{
				int index = 0; //initialize index
				if (i == 0)
				{
					index = frontResultValue % m_vec_len;
				}
				else if (i == 1)
				{
					index = endResultValue % m_vec_len;
				}
				else
				{
					index = (frontResultValue + i*endResultValue) % m_vec_len;
				}
				if (false == m_vec[index])
					return false;
			}
			
			return true;
		}
		
		void getHashedValue(long long inputData, unsigned char* hashedValue) 
		{
			uint64_t frontResultValue = 0, endResultValue = 0;
			int lastValue = 0;
			unsigned char digest[MD5_DIGEST_LENGTH];
			MD5_CTX ctx;
			MD5_Init(&ctx);
			char inputString[20];
			sprintf_s(inputString, "%lld", inputData);
			MD5_Update(&ctx, inputString, strlen(inputString));
			MD5_Final(digest, &ctx); //computed hash
			hashedValue = digest;
			
		}

		//Getter events
		
		int GetNumberofElements()
		{
			return m_ELEM_NUM;
		}

		double GetAccuracyPercentage()
		{
			return m_PERCENTAGE;
		}

		int GetVectorLength()
		{
			return m_vec_len;
		}

		vector<bool> GetVector()
		{
			return m_vec;
		}

		int GetNumberofHashes()
		{
			return m_hash_num;
		}

		//setter events

		void SetNumberofElements(int n)
		{
			m_ELEM_NUM = n;
		}

		void SetAccuracyPercentage(double p)
		{
			m_PERCENTAGE = p;
		}

		void SetVectorLength(int len)
		{
			if (m_vec.size() != len)
			{
				m_vec.reserve(len);
				m_vec.resize(len);
			}
			m_vec_len = len;
		}

		void SetVector(vector<bool> bf)
		{
			m_vec_len = bf.size(); //change vector length
			m_vec.reserve(bf.size());
			m_vec.resize(bf.size());
			m_vec = bf;
		}

		void SetNumberofHashes(int k)
		{
			m_hash_num = k;
		}



};
