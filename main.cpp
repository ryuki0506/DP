#include "field.hpp"
#include "culc.hpp"
#include "output.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int len = 500;	 // ポリマーの長さの最大値
const double pmax = 1;	 // サイトがopenな確率の最大
const double pmin = 0.6; // サイトがopenな確率の最小
const int steps = 100;	 //>1 pの刻み数
const int shots = 100;	 // 試行回数

const int Emax = 21;
const double temp =0;

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
const int output_mode = 2;
/*
output_mode==1 :A vs p
output_mode==2 :S vs E
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
	ofs << temp << endl;
	ofs << len << endl;

	for (int step = 0; step < steps; step++)
	{
		double p = pmin + step * Dp;
		ofs << p << ", ";

		double FPT = 0;
		double W = 0;
		double WofE[Emax];
		for (size_t E = 0; E < Emax; E++)
		{
			WofE[E] = 0;
		}

		for (size_t shot = 0; shot < shots; shot++)
		{
			Field *field;
			field = new Field(len, Emax);
			field->set_potential(p, noise_mode);
			field->time_evolution(temp);

			if (temp == 0)
			{
				double FPT_shot = field->calc_FPT(parcolation, Isfixed);
				double W_shot = field->calc_W_Emin(parcolation, Isfixed);

				double _Emin = field->calc_Emin(parcolation, Isfixed);
				int Emin = int(_Emin);

				FPT += FPT_shot;
				W += W_shot;
				if (Emin < Emax)
				{
					//WofE[Emin] += W_shot;
					WofE[Emin] += 1;
				}
			}
			else if (temp > 0)
			{
				FPT+=field->get_Z();
				for (size_t E = 0; E < Emax; E++)
				{
					WofE[E]+=field->get_WofE()[E];
				}				
			}

			delete field;
		}

		FPT /= shots;
		W /= shots;
		for (size_t E = 0; E < Emax; E++)
		{
			WofE[E] /= shots;
		}

		if (output_mode == 1)
		{
			if (calc_mode == 1)
			{
				ofs << FPT / len << endl;
			}
			else if (calc_mode == 2)
			{
				ofs << log(W) / len << endl;
			}
		}
		else if (output_mode == 2)
		{
			for (size_t E = 0; E < Emax; E++)
			{
				if (E < Emax - 1)
				{
					ofs << WofE[E] << ",";
				}
				else
				{
					ofs << WofE[E] << endl;
				}
			}
		}
	}

	ofs.close();

	return 0;
}
