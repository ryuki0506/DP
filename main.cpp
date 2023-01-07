#include "field.hpp"
#include "culc.hpp"
#include "output.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int lenmax = 500; // ポリマーの長さの最大値
const double pmax = 1;	// サイトがopenな確率の最大
const double pmin = 0.6;	// サイトがopenな確率の最小
const int steps = 100;	//>1 pの刻み数
const int shots = 100;	// 試行回数

const int Emax = 21;

const int noise_mode = 1; // 計算するノイズの種類
/*
noize_mode==1 :Bernulli分布
noize_mode==2 :geometric分布
noize_mode==3 :exponential分布
noize_mode==4 :log-gamma分布
*/
const int calc_mode = 2;
/*
calc_mode==1 :growth rate
calc_mode==2 :entropy
*/
const bool Isfixed = false;
const bool parcolation = false;		 // parcolationとして計算するか？
const bool show_in_terminal = false; // ターミナルに表示するか？

int main()
{
	double Dp = (pmax - pmin) / (steps - 1);
	ofstream ofs("../result/result.txt");

	ofs << noise_mode << endl;
	ofs << calc_mode << endl;
	ofs << Isfixed << endl;
	ofs << parcolation << endl;

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
/*
						double SofE_all_path[Emax];
						SofE(SofE_all_path,Emax,len,p,shots,noise_mode,Isfixed);
						output_SofE(SofE_all_path,Emax,show_in_terminal);

						for (size_t E = 0; E < Emax; E++)
						{
							if (E<Emax-1)
							{
								ofs<<log(SofE_all_path[E])/len<<",";
							}else{
								ofs<<log(SofE_all_path[E])/len<<endl;
							}
						}
*/

/* 			
						double _Sofp=Sofp(len,lenmax,p,shots,noise_mode,calc_mode,show_in_terminal,Isfixed,parcolation);
			//			double re = limited_average(_simulation, shots);

						if (len < lenmax)
						{
							ofs << log(_Sofp)/len << ", ";
						}
						else
						{
							ofs << log(_Sofp)/len << endl;
						}
 */			

			double _SofE_min[Emax];
			SofE_min(_SofE_min, Emax, len, p, shots, noise_mode, Isfixed,parcolation);

			for (size_t E = 0; E < Emax; E++)
			{
				if (E < Emax - 1)
				{
					ofs << log(_SofE_min[E]) / len << ",";
				}
				else
				{
					ofs << log(_SofE_min[E]) / len << endl;
				}
			}


		}
	}

	ofs.close();

	return 0;
}