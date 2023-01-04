#include <iostream>
#include <cmath>
#include <random>
#include "field.hpp"
using namespace std;

Field::Field(int len,int Energy) : field_size(len)
{
	field = new double[field_size * field_size];
	partition_function = new double[field_size * field_size];
	num_of_least_energy_pathes = new double[field_size * field_size];

	for (int i = 0; i < field_size; i++) // 分配関数を初期化
	{
		for (int j = 0; j < field_size; j++)
		{
			partition_function[field_size * i + j] = 0;
			num_of_least_energy_pathes[field_size * i + j] = 0;
		}
	}
	
}

Field::~Field()
{
	delete[] field;
	delete[] partition_function;
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
				field[field_size * i + j] = exp(-dist(engine));
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
				field[field_size * i + j] = exp(-dist(engine));
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

void Field::set_partition_function()
{
	for (int i = 0; i < field_size; i++) // 分配関数を漸化式に従って計算
	{
		for (int j = 0; j <= i; j++)
		{
			if (i == 0 and j == 0)
			{
				partition_function[field_size * (i - j) + j] = field[field_size * (i - j) + j];
				num_of_least_energy_pathes[field_size * (i - j) + j] = 1.0;
			}
			else if (j == 0)
			{
				partition_function[field_size * (i - j) + j] = field[field_size * (i - j) + j] * partition_function[field_size * (i - j - 1) + j];
				num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j - 1) + j];
			}
			else if (j == i)
			{
				partition_function[field_size * (i - j) + j] = field[field_size * (i - j) + j] * partition_function[field_size * (i - j) + j - 1];
				num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j) + j - 1];
			}
			else
			{
				partition_function[field_size * (i - j) + j] = field[field_size * (i - j) + j] * max(partition_function[field_size * (i - j - 1) + j], partition_function[field_size * (i - j) + j - 1]);
				if (partition_function[field_size * (i - j - 1) + j] < partition_function[field_size * (i - j) + j - 1])
				{
					num_of_least_energy_pathes[field_size * (i - j) + j] = num_of_least_energy_pathes[field_size * (i - j) + j - 1];
				}
				else if (partition_function[field_size * (i - j - 1) + j] > partition_function[field_size * (i - j) + j - 1])
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

double *Field::SofE_all_path(int Emax,int pos, int depth)
{
	double SofE[Emax];
	for (size_t i = 0; i < Emax; i++)
	{
		SofE[i]=0;
	}
	
	if (depth == 0)
	{
		//cout<<"pos="<<pos<<","<<"depth="<<depth<<endl;
		double pot = field[field_size * (depth - pos) + pos];
		int E = (int)pot;

		SofE[E] = 1;
	}
	else if(depth>0)
	{	
		if (pos == 0)
		{
			//cout<<"pos="<<pos<<","<<"depth="<<depth<<endl;
			double pot_u = field[field_size * (depth-1 - pos) + pos];
			int E_u = (int)pot_u;

			double *_SofE_u=SofE_all_path(Emax,pos, depth-1);

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
			//cout<<"pos="<<pos<<","<<"depth="<<depth<<endl;
			double pot_l = field[field_size * (depth-1 - (pos-1)) + pos-1];
			int E_l = (int)pot_l;

			double *_SofE_l=SofE_all_path(Emax,pos-1, depth-1);

			for (size_t j = 0; j < Emax; j++)
			{
				if (E_l + j < Emax)
				{
					SofE[E_l + j] += _SofE_l[j];
				}
			}
		}
		else if(pos>0 && pos<depth)
		{
			//cout<<"pos="<<pos<<","<<"depth="<<depth<<endl;
			double pot_u = field[field_size * (depth-1 - pos) + pos];
			int E_u = (int)pot_u;
			double pot_l = field[field_size * (depth-1 - (pos-1)) + pos-1];
			int E_l = (int)pot_l;
			
			double *_SofE_u=SofE_all_path(Emax,pos, depth-1);
			double *_SofE_l=SofE_all_path(Emax,pos-1, depth-1);

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

double *Field::get_partition_function()
{
	return partition_function;
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
		double max_partition_func;
		if (parcolation == true) // 分配関数の最大値の初期化(parcolationの時はノイズの最小値)
		{
			max_partition_func = 1.0;
		}
		else
		{
			max_partition_func = partition_function[field_size * (field_size - 1)];
		}

		for (int j = 0; j < field_size; j++) // 長さfield_sizeの各分配関数の最大値を計算
		{
			if (max_partition_func < partition_function[field_size * (field_size - j - 1) + j])
			{
				max_partition_func = partition_function[field_size * (field_size - j - 1) + j];
			}
		}

		for (int j = 0; j < field_size; j++) // 最大のもののみ計算
		{
			if (max_partition_func == partition_function[field_size * (field_size - j - 1) + j])
			{
				if (calc_mode == 1)
				{
					sum += num_of_least_energy_pathes[field_size * (field_size - j - 1) + j] * partition_function[field_size * (field_size - j - 1) + j];
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
		double sum = 0;
		double max_partition_func;
		if (parcolation == true) // 分配関数の最大値の初期化(parcolationの時はノイズの最小値)
		{
			max_partition_func = 1.0;
		}
		else
		{
			max_partition_func = partition_function[field_size * (field_size / 2 - 1) + field_size / 2];
		}

		for (int j = field_size / 2 - 1; j <= field_size / 2; j++) // 長さfield_sizeの各分配関数の最大値を計算
		{
			if (max_partition_func < partition_function[field_size * (field_size - j - 1) + j])
			{
				max_partition_func = partition_function[field_size * (field_size - j - 1) + j];
			}
		}

		for (int j = field_size / 2 - 1; j <= field_size / 2; j++) // 最大のもののみ計算
		{
			if (max_partition_func == partition_function[field_size * (field_size - j - 1) + j])
			{
				if (calc_mode == 1)
				{
					sum += num_of_least_energy_pathes[field_size * (field_size - j - 1) + j] * partition_function[field_size * (field_size - j - 1) + j];
				}
				else if (calc_mode == 2)
				{
					sum += num_of_least_energy_pathes[field_size * (field_size - j - 1) + j];
				}
			}
		}
		return sum;
	}
}

double Field::get_growth_rate(bool parcolation, bool Isfixed)
{
	double Z = calc_pysical_quantity(1, parcolation, Isfixed);
	return log(Z) / field_size;
}

double Field::get_entropy(bool parcolation, bool Isfixed)
{
	double W = calc_pysical_quantity(2, parcolation, Isfixed);
	return log(W) / field_size;
}
