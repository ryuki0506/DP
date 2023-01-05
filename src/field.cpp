#include <iostream>
#include <cmath>
#include <random>
#include "field.hpp"
using namespace std;

Field::Field(int len) : field_size(len)
{
	field = new double[field_size * field_size];
	LPT = new double[field_size * field_size];
	num_of_least_energy_pathes = new double[field_size * field_size];

	for (int i = 0; i < field_size; i++) //分配関数を初期化
	{
		for (int j = 0; j < field_size; j++)
		{
			LPT[field_size * i + j] = 0;
			num_of_least_energy_pathes[field_size * i + j] = 0;
		}
	}
}

Field::~Field()
{
	delete[] field;
	delete[] LPT;
	delete[] num_of_least_energy_pathes;
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
				field[field_size * i + j] = dist(engine);
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
				field[field_size * i + j] = dist(engine);
			}
		}
	}
}

double *Field::get_potential()
{
	return field;
}

void Field::set_LPT()
{
	for (int i = 0; i < field_size; i++) //分配関数を漸化式に従って計算
	{
		for (int j = 0; j <= i; j++)
		{
			if (i == 0 and j == 0)
			{
				LPT[field_size * (i - j) + j] = field[field_size * (i - j) + j];
				num_of_least_energy_pathes[field_size * (i - j) + j] = 1.0;
			}
			else if (j == 0)
			{
				LPT[field_size * (i - j) + j] = field[field_size * (i - j) + j] + LPT[field_size * (i - j - 1) + j];
				num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j - 1) + j];
			}
			else if (j == i)
			{
				LPT[field_size * (i - j) + j] = field[field_size * (i - j) + j] + LPT[field_size * (i - j) + j - 1];
				num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j) + j - 1];
			}
			else
			{
				LPT[field_size * (i - j) + j] = field[field_size * (i - j) + j] + max(LPT[field_size * (i - j - 1) + j], LPT[field_size * (i - j) + j - 1]);
				if (LPT[field_size * (i - j - 1) + j] < LPT[field_size * (i - j) + j - 1])
				{
					num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j) + j - 1];
				}
				else if (LPT[field_size * (i - j - 1) + j] > LPT[field_size * (i - j) + j - 1])
				{
					num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j - 1) + j];
				}
				else
				{
					num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j) + j - 1] + num_of_least_energy_pathes[field_size * (i - j - 1) + j];
				}
			}
		}
	}
}

double *Field::get_LPT()
{
	/* 
	for (int i = 0; i < field_size; i++)
	{
		for (int j = 0; j < field_size; j++)
		{
			LPT[field_size * i + j] *= num_of_least_energy_pathes[field_size * i + j];
		}
	}
 */
	return LPT;
}

double *Field::get_num_of_least_energy_pathes()
{
	return num_of_least_energy_pathes;
}

double Field::calc_pysical_quantity(int calc_mode, bool parcolation, bool Isfixed)
{
	if (!Isfixed)
	{
		double sum = 0;
		double max_LPT;
		if (parcolation == true) //分配関数の最大値の初期化(parcolationの時はノイズの最小値)
		{
			max_LPT = 1.0;
		}
		else
		{
			max_LPT = LPT[field_size * (field_size - 1)];
		}

		for (int j = 0; j < field_size; j++) //長さfield_sizeの各分配関数の最大値を計算
		{
			if (max_LPT < LPT[field_size * (field_size - j - 1) + j])
			{
				max_LPT = LPT[field_size * (field_size - j - 1) + j];
			}
		}

		for (int j = 0; j < field_size; j++) //最大のもののみ計算
		{
			if (max_LPT == LPT[field_size * (field_size - j - 1) + j])
			{
				if (calc_mode == 1)
				{
					sum += num_of_least_energy_pathes[field_size * (field_size - j - 1) + j] * LPT[field_size * (field_size - j - 1) + j];
				}
				else if (calc_mode == 2)
				{
					sum += num_of_least_energy_pathes[field_size * (field_size - j - 1) + j];
				}
			}
		}
		return sum;
	}
	else
	{
		double max_LPT;
		if (parcolation == true) //分配関数の最大値の初期化(parcolationの時はノイズの最小値)
		{
			max_LPT = 1.0;
		}
		else
		{
			max_LPT = LPT[field_size * (field_size/2 - 1)+field_size/2];
		}

		if (max_LPT < LPT[field_size * (field_size/2 - 1)+field_size/2])
			{
				if (calc_mode == 1)
				{
					return num_of_least_energy_pathes[field_size * (field_size/2 - 1)+field_size/2] * LPT[field_size * (field_size/2 - 1)+field_size/2];
				}
				else if (calc_mode == 2)
				{
					return num_of_least_energy_pathes[field_size * (field_size/2 - 1)+field_size/2];
				}
			}else{
				if (calc_mode == 1)
				{
					return num_of_least_energy_pathes[field_size * (field_size/2)+field_size/2-1] * LPT[field_size * (field_size/2)+field_size/2-1];
				}
				else if (calc_mode == 2)
				{
					return num_of_least_energy_pathes[field_size * (field_size/2)+field_size/2-1];
				}
			}
	}
}

double Field::get_growth_rate(bool parcolation, bool Isfixed)
{
	double G = calc_pysical_quantity(1, parcolation, Isfixed);
	return G / field_size;
}

double Field::get_entropy(bool parcolation, bool Isfixed)
{
	double W = calc_pysical_quantity(2, parcolation, Isfixed);
	return log(W) / field_size;
}
