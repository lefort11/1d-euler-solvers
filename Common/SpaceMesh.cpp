#include "SpaceMesh.h"
#include "assert.h"

const SpaceMesh operator-(SpaceMesh const& left, SpaceMesh const& right)
{
	assert(left.size() == right.size());

	SpaceMesh result(left.size());

	for(unsigned long i = 0; i < left.size(); ++i)
	{
		result[i].density = left[i].density - right[i].density;
		result[i].velocity = left[i].velocity - right[i].velocity;
		result[i].pressure = left[i].pressure - right[i].pressure;

	}

	return result;
}