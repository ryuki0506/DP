#include <iostream>
#include <cmath>
#include "culc.hpp"
#include "field.hpp"
#include "output.hpp"

using namespace std;

double limited_average(double *data, int data_len)
{
	double sum = 0;
	int non0_count = 0;
	int inf_count = 0;

	for (int i = 0; i < data_len; i++)
	{
		if (data[i] != log(0))
		{
			if (data[i] != 0)
			{
				sum += data[i];
				non0_count++;
			}
		}
		else
		{
			inf_count++;
		}
	}
	if (non0_count != 0)
	{
		return sum / non0_count;
	}
	else if (inf_count == data_len)
	{
		return log(0);
	}
	else
	{
		return 0;
	}
};

double average(double *data, int data_len)
{
	double sum = 0;
	for (size_t i = 0; i < data_len; i++)
	{
		sum += data[i];
	}
	return sum / data_len;
}
