#include "field.hpp"
#include "culc.hpp"
#include "output.hpp"

#include <iostream>
#include <fstream>
using namespace std;

const int lenmax = 100;//ポリマーの長さの最大値
const double pmax = 1;//サイトがopenな確率の最大
const double pmin = 0.5;//サイトがopenな確率の最小
const int steps = 100; //>1 pの
const int shots = 10;

int main()
{
	double Dp = (pmax - pmin) / (steps - 1);
	ofstream ofs("../result/result.txt");

	for (int len = 0; len <= lenmax; len++)
	{
		if (len < lenmax)
		{
			ofs << len << ", ";
		}
		else
		{
			ofs << len << endl;
		}
	}

	for (int step = 0; step < steps; step++)
	{
		double p = pmin + step * Dp;
		ofs << p << ", ";
		for (int len = 1; len <= lenmax; len++)
		{
			double shots_data[shots];
			for (int shot = 0; shot < shots; shot++)
			{
				Field *field;
				field = new Field(len);
				field->set_potential(p);

				shots_data[shot] = field->get_growth_rate();
				delete field;
			}
			double re = limited_average(shots_data, shots);

			if (len < lenmax)
			{
				ofs << re << ", ";
			}
			else
			{
				ofs << re << endl;
			}
		}
	}
	ofs.close();

	return 0;
}