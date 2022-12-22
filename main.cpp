#include "field.hpp"
#include "culc.hpp"
#include "output.hpp"

#include <iostream>
#include <fstream>
using namespace std;

const int lenmax = 200;//ポリマーの長さの最大値
const double pmax = 1;//サイトがopenな確率の最大
const double pmin = 0;//サイトがopenな確率の最小
const int steps = 100; //>1 pの刻み数
const int shots = 100;//試行回数

const int noise_mode=2;//計算するノイズの種類
/*
noize_mode==1 :Bernulli分布
noize_mode==2 :geometric分布
noize_mode==3 :exponential分布
noize_mode==4 :log-gamma分布
*/
const int calc_mode=2;
/*
calc_mode==1 :growth rate
calc_mode==2 :entropy
*/
const bool Isfixed=false;
const bool parcolation=true;//parcolationとして計算するか？
const bool show_in_terminal=false;//ターミナルに表示するか？

int main()
{
	double Dp = (pmax - pmin) / (steps - 1);
	ofstream ofs("../result/result.txt");

	ofs<<noise_mode<<endl;
	ofs<<calc_mode<<endl;
	ofs<<Isfixed<<endl;
	ofs<<parcolation<<endl;

	ofs << 0 << ", ";
	for (int len = lenmax; len <= lenmax; len++)
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
		for (int len = lenmax; len <= lenmax; len++)
		{
			double shots_data[shots];
			for (int shot = 0; shot < shots; shot++)
			{
				Field *field;
				field = new Field(len);
				field->set_potential(p,noise_mode);
				field->set_partition_function();

				if (calc_mode==1)
				{
					show_field(field->get_partition_function(),len,lenmax,show_in_terminal);
					shots_data[shot] = field->get_growth_rate(parcolation,Isfixed);
				}else if (calc_mode==2)
				{
					show_field(field->get_num_of_least_energy_pathes(),len,lenmax,show_in_terminal);
					shots_data[shot] = field->get_entropy(parcolation,Isfixed);
				}
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