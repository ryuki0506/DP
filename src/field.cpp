#include <iostream>
#include <cmath>
#include <random>
#include "field.hpp"
using namespace std;

Field::Field(int len, int E) : field_size(len), Emax(E)
{
	Z = 0;
	field = new double[field_size * field_size];
	FPT = new double[field_size * field_size];
	W_Emin = new double[field_size * field_size];
	WofE = new double[Emax];

	for (int i = 0; i < field_size; i++) // 分配関数を初期化
	{
		for (int j = 0; j < field_size; j++)
		{
			FPT[field_size * i + j] = 0;
			W_Emin[field_size * i + j] = 0;
		}
	}

	for (size_t E = 0; E < Emax; E++)
	{
		WofE[E] = 0;
	}
}

Field::~Field()
{
	delete[] field;
	delete[] FPT;
	delete[] W_Emin;
	delete[] WofE;
}

void Field::set_size(int size)
{
	field_size = size;
}

int Field::get_size()
{
	return field_size;
}

void Field::set_potential(double p, int mode)
{
	random_device seed_gen;
	default_random_engine engine(seed_gen());

	if (mode == 1)
	{
		bernoulli_distribution dist(1 - p);
		for (int i = 0; i < field_size; i++)
		{
			for (int j = 0; j < field_size; j++)
			{
				field[field_size * i + j] = dist(engine);
			}
		}
	}
	else if (mode == 2)
	{
		geometric_distribution<int> dist(p);
		for (int i = 0; i < field_size; i++)
		{
			for (int j = 0; j < field_size; j++)
			{
				field[field_size * i + j] = dist(engine);
			}
		}
	}
	else if (mode == 3)
	{
		exponential_distribution<double> dist(p);
		for (int i = 0; i < field_size; i++)
		{
			for (int j = 0; j < field_size; j++)
			{
				field[field_size * i + j] = -dist(engine);
			}
		}
	}
	else if (mode == 4)
	{
		gamma_distribution<double> dist(p, 1.0);
		for (int i = 0; i < field_size; i++)
		{
			for (int j = 0; j < field_size; j++)
			{
				field[field_size * i + j] = exp(-log(dist(engine)));
			}
		}
	}
}

double *Field::get_potential()
{
	return field;
}

void Field::time_evolution(double temp)
{
	if (temp == 0)
	{
		for (int i = 0; i < field_size; i++) // 分配関数を漸化式に従って計算
		{
			for (int j = 0; j <= i; j++)
			{
				if (i == 0 and j == 0)
				{
					FPT[field_size * (i - j) + j] = field[field_size * (i - j) + j];
					W_Emin[field_size * (i - j) + j] = 1.0;
				}
				else if (j == 0)
				{
					FPT[field_size * (i - j) + j] = field[field_size * (i - j) + j] + FPT[field_size * (i - j - 1) + j];
					W_Emin[field_size * (i - j) + j] = W_Emin[field_size * (i - j - 1) + j];
				}
				else if (j == i)
				{
					FPT[field_size * (i - j) + j] = field[field_size * (i - j) + j] + FPT[field_size * (i - j) + j - 1];
					W_Emin[field_size * (i - j) + j] = W_Emin[field_size * (i - j) + j - 1];
				}
				else
				{
					FPT[field_size * (i - j) + j] = field[field_size * (i - j) + j] + min(FPT[field_size * (i - j - 1) + j], FPT[field_size * (i - j) + j - 1]);
					if (FPT[field_size * (i - j - 1) + j] > FPT[field_size * (i - j) + j - 1])
					{
						W_Emin[field_size * (i - j) + j] = W_Emin[field_size * (i - j) + j - 1];
					}
					else if (FPT[field_size * (i - j - 1) + j] < FPT[field_size * (i - j) + j - 1])
					{
						W_Emin[field_size * (i - j) + j] = W_Emin[field_size * (i - j - 1) + j];
					}
					else
					{
						W_Emin[field_size * (i - j) + j] = W_Emin[field_size * (i - j) + j - 1] + W_Emin[field_size * (i - j - 1) + j];
					}
				}
			}
		}
	}
	else if (temp > 0)
	{
		set_WofE();
		set_Z(temp);
	}
}

double *Field::calc_WofE(int pos, int depth)
{
	double SofE[Emax];
	for (size_t i = 0; i < Emax; i++)
	{
		SofE[i] = 0;
	}

	if (depth == 0)
	{
		double pot = field[field_size * (depth - pos) + pos];
		int E = (int)pot;

		SofE[E] = 1;
	}
	else if (depth > 0)
	{
		if (pos == 0)
		{
			double pot_u = field[field_size * (depth - 1 - pos) + pos];
			int E_u = (int)pot_u;

			double *_SofE_u = calc_WofE(pos, depth - 1);

			for (size_t j = 0; j < Emax; j++)
			{
				if (E_u + j < Emax)
				{
					SofE[E_u + j] += _SofE_u[j];
				}
			}
		}
		else if (pos == depth)
		{
			double pot_l = field[field_size * (depth - 1 - (pos - 1)) + pos - 1];
			int E_l = (int)pot_l;

			double *_SofE_l = calc_WofE(pos - 1, depth - 1);

			for (size_t j = 0; j < Emax; j++)
			{
				if (E_l + j < Emax)
				{
					SofE[E_l + j] += _SofE_l[j];
				}
			}
		}
		else if (pos > 0 && pos < depth)
		{
			double pot_u = field[field_size * (depth - 1 - pos) + pos];
			int E_u = (int)pot_u;
			double pot_l = field[field_size * (depth - 1 - (pos - 1)) + pos - 1];
			int E_l = (int)pot_l;

			double *_SofE_u = calc_WofE(pos, depth - 1);
			double *_SofE_l = calc_WofE(pos - 1, depth - 1);

			for (size_t j = 0; j < Emax; j++)
			{
				if (E_u + j < Emax)
				{
					SofE[E_u + j] += _SofE_u[j];
				}
				if (E_l + j < Emax)
				{
					SofE[E_l + j] += _SofE_l[j];
				}
			}
		}
	}
	return SofE;
}

double *Field::get_FPT()
{
	return FPT;
}

double *Field::get_W_Emin()
{
	return W_Emin;
}

void Field::set_WofE()
{
	for (size_t i = 0; i < field_size; i++)
	{
		for (size_t E = 0; E < Emax; E++)
		{
			WofE[E] += calc_WofE(i, field_size - 1)[E];
		}
	}
}

double *Field::get_WofE()
{
	return WofE;
}

void Field::set_Z(double temp)
{
	for (size_t E = 0; E < Emax; E++)
	{
		Z += WofE[E] * exp(-E / temp);
	}
}

double Field::get_Z()
{
	return Z;
}

double Field::calc_pysical_quantity(int calc_mode, bool parcolation, bool Isfixed)
{
	double sum = 0;
	double min_FPT = calc_Emin(parcolation, Isfixed);

	if (!Isfixed)
	{
		for (int j = 0; j < field_size; j++) // 最小のもののみ計算
		{
			if (min_FPT == FPT[field_size * (field_size - j - 1) + j])
			{
				if (calc_mode == 1)
				{
					sum += W_Emin[field_size * (field_size - j - 1) + j] * FPT[field_size * (field_size - j - 1) + j];
				}
				else if (calc_mode == 2)
				{
					sum += W_Emin[field_size * (field_size - j - 1) + j];
				}
			}
		}
		return sum;
	}
	else
	{
		for (int j = field_size / 2 - 1; j <= field_size / 2; j++) // 最大のもののみ計算
		{
			if (min_FPT -FPT[field_size * (field_size - j - 1) + j]<0.1)
			{
				if (calc_mode == 1)
				{
					sum += W_Emin[field_size * (field_size - j - 1) + j] * FPT[field_size * (field_size - j - 1) + j];
				}
				else if (calc_mode == 2)
				{
					sum += W_Emin[field_size * (field_size - j - 1) + j];
				}
			}
		}
		return sum;
	}
}

double Field::calc_Emin(bool parcolation, bool Isfixed)
{
	double min_FPT;
	if (parcolation == true) // 分配関数の最大値の初期化(parcolationの時はノイズの最小値)
	{
		min_FPT = 0;
	}
	else
	{
		if (!Isfixed)
		{
			min_FPT = FPT[field_size * (field_size - 1)];
			for (int j = 0; j < field_size; j++) // 長さfield_sizeの各分配関数の最大値を計算
			{
				if (min_FPT > FPT[field_size * (field_size - j - 1) + j])
				{
					min_FPT = FPT[field_size * (field_size - j - 1) + j];
				}
			}
		}
		else
		{
			min_FPT = FPT[field_size * (field_size / 2) + field_size / 2-1];
			for (int j = field_size / 2 - 1; j <= field_size / 2; j++) // 長さfield_sizeの各分配関数の最大値を計算
			{
				if (min_FPT > FPT[field_size * (field_size - j - 1) + j])
				{
					min_FPT = FPT[field_size * (field_size - j - 1) + j];
				}
			}
		}
	}
	return min_FPT;
}

double Field::calc_FPT(bool parcolation, bool Isfixed)
{
	double FPT = calc_pysical_quantity(1, parcolation, Isfixed);
	return FPT;
}

double Field::calc_W_Emin(bool parcolation, bool Isfixed)
{
	double W = calc_pysical_quantity(2, parcolation, Isfixed);
	return W;
}

double Field::calc_entropy(bool parcolation, bool Isfixed)
{
	double W = calc_pysical_quantity(2, parcolation, Isfixed);
	return log(W);
}
