#include "field.hpp"
#include "culc.hpp"
#include "output.hpp"

#include <iostream>
#include <fstream>
using namespace std;

const int lenmax = 3;//ポリマーの長さの最大値
const double pmax = 1;//サイトがopenな確率の最大
const double pmin = 0;//サイトがopenな確率の最小
const int steps = 2; //>1 pの刻み数
const int shots = 1;//試行回数

const int noise_mode=1;//計算するノイズの種類
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
const bool parcolation=false;//parcolationとして計算するか？
const bool show_in_terminal=true;//ターミナルに表示するか？

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
			double re = limited_average(simulation(len,lenmax,p,shots,noise_mode,calc_mode,show_in_terminal,Isfixed,parcolation), shots);

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