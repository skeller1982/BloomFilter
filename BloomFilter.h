// Bloom filter
// implementation according to https://en.wikipedia.org/wiki/Bloom_filter

/*
The MIT License(MIT)

Copyright 2017 Steffen Keller

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files(the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/* usage:
provide the type of your data as a first template parameter T (you can put in any type you want)
provide the size of the Hash returned by your Hashfunction as a second parameter HashSize (32 or 64 bit integral)

for the constructor provide:
- std::function<HashSize(T, uint32_t)> hashFunc		- a hashing function that takes your data T plus a uint32_t salt and returning HashSize
- uint32_t maxExpectedEntries						- the maximum number of entries you expect (deviations from this value will affect your false positive rate)
- double falsePositivePropability					- ]0,1] the acceptable propability of a false positive

remarks on how to choose parameters:
- if you put in more values than expected your false-positive rate will increase, but everything will work fine, if you put in less, your false positive rate will be lower
- falsePositivePropability will only affect the size and performance logarithmically
- maxExpectedEntries will affect the size linearly

You should choose a fast hashing like MurmurHash3 for optimal performance. a bad hash distribution will have a negative effect on the false positive rate of this filter.

example: BloomFilter<std::string, uint64_t>(MyHashStringFunction, 10000, 0.01);
*/

#pragma once

#include <cmath>
#include <random>
#include <functional>
#include <limits>

namespace
{
uint8_t mask[8] = { 0x01, 0x02, 0x04 , 0x08, 0x10 , 0x20 , 0x40 , 0x80 };

static uint32_t GetOptHashes(const double& falsePositivePropability)
{
	double optHashes = -1.0 * std::log(falsePositivePropability) / std::log(2);
	int32_t numHashes = static_cast<int32_t>(std::rint(optHashes));
	return std::max(numHashes, 1);
}

static size_t GetOptSize(uint32_t maxExpectedEntries, const double& falsePositivePropability)
{
	double optSize = (-1.0 * maxExpectedEntries * std::log(falsePositivePropability)) / (std::log(2) * std::log(2));
	size_t setSize = static_cast<size_t>(std::rint(optSize));
	return std::max(setSize, static_cast<size_t>(64u));
}
}

namespace bloom
{

template <typename T, typename HashSize>
class BloomFilter
{
public:

	using HashFunc = std::function<HashSize(T, uint32_t)>;

	void InsertValue(const T& val)
	{
		for (const auto& salt: m_salts)
		{
			HashSize hash = m_hashfunc(val, salt);
			uint64_t bit = hash % m_setSize;
			set(bit);
		}
	}

	bool ContainsValue(const T& val) const
	{
		for (const auto& salt : m_salts)
		{
			HashSize hash = m_hashfunc(val, salt);
			uint64_t bit = hash % m_setSize;
			
			if (!test(bit))
			{
				return false;
			}
		}
		return true;
	}

	BloomFilter(HashFunc hashFunc, uint32_t maxExpectedEntries, double falsePositivePropability)
		: m_hashfunc(hashFunc)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<uint32_t> dist;

		// make sure propability is within boundaries
		falsePositivePropability = std::min(1.0, falsePositivePropability);
		falsePositivePropability = std::max(std::numeric_limits<double>::denorm_min(), falsePositivePropability);

		// calculate optimal parameters
		m_setSize = GetOptSize(maxExpectedEntries, falsePositivePropability);
		uint32_t numHashes  = GetOptHashes(falsePositivePropability);

		//generate salts according to the number of different hashes needed
		for (uint32_t i = 0; i < numHashes; ++i)
		{
			m_salts.push_back(dist(gen));
		}

		m_bitVector.resize((m_setSize / 8) + 1);
	}

private:

	void set(const uint64_t& bit)
	{
		size_t vecpos = static_cast<size_t>(bit / 8);
		uint8_t bitPos = bit % 8;
		m_bitVector[vecpos] |= mask[bitPos];
	}

	bool test(const uint64_t& bit) const
	{
		size_t vecpos = static_cast<size_t>(bit / 8);
		uint8_t bitPos = bit % 8;
		return m_bitVector[vecpos] & mask[bitPos];
	}

	HashFunc m_hashfunc;
	size_t m_setSize;
	std::vector<uint32_t> m_salts;
	std::vector<uint8_t> m_bitVector;
};

} //namespace bloom
