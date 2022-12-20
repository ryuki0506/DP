#include <iostream>
#include <cmath>
#include <random>
#include "field.hpp"
using namespace std;

Field::Field(int len) : field_size(len)
{
	field = new double[field_size * field_size];
	partition_function = new double[field_size * field_size];
	num_of_least_energy_pathes = new double[field_size * field_size];

	for (int i = 0; i < field_size; i++) //分配関数を初期化
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
				field[field_size * i + j] = exp(-dist(engine));
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
	}else if (mode == 3){
		exponential_distribution<double> dist(p);
		for (int i = 0; i < field_size; i++)
		{
			for (int j = 0; j < field_size; j++)
			{
				field[field_size * i + j] = exp(-dist(engine));
			}
		}
	}else if (mode == 4){
		gamma_distribution<double> dist(p,1.0);
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
	for (int i = 0; i < field_size; i++) //分配関数を漸化式に従って計算
	{
		for (int j = 0; j <= i; j++)
		{
			if (i == 0 and j == 0)
			{
				partition_function[field_size * (i - j) + j] = field[field_size * (i - j) + j];
				num_of_least_energy_pathes[field_size * (i - j) + j] = 1;
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

double *Field::get_partition_function()
{
	for (int i = 0; i < field_size; i++)
	{
		for (int j = 0; j < field_size; j++)
		{
			partition_function[field_size * i + j] *= num_of_least_energy_pathes[field_size * i + j];
		}
	}

	return partition_function;
}

double *Field::get_num_of_least_energy_pathes(){
	return num_of_least_energy_pathes;
}

double Field::calc_pysical_quantity(int calc_mode,bool parcolation)
{
	/* for (int i = 0; i < field_size; i++)//分配関数は縮退度がかかっているので、取り除く
	{
		for (int j = 0; j < field_size; j++)
		{
			if (num_of_least_energy_pathes[field_size * i + j]!=0)
			{
				partition_function[field_size * i + j] /= num_of_least_energy_pathes[field_size * i + j];				
			}else{
				partition_function[field_size * i + j]=0;
			}
		}
	} */

	double sum = 0;
	double max_partition_func;
	if (parcolation == true)//分配関数の最大値の初期化(parcolationの時はノイズの最小値)
	{
		max_partition_func = 1;
	}
	else
	{
		max_partition_func = partition_function[field_size * (field_size - 1)];
	}

	for (int j = 0; j < field_size; j++)//長さfield_sizeの各分配関数の最大値を計算
	{
		if (max_partition_func < partition_function[field_size * (field_size - j - 1) + j])
		{
			max_partition_func = partition_function[field_size * (field_size - j - 1) + j];
		}
	}

	for (int j = 0; j < field_size; j++)//最大のもののみ計算
	{
		if (max_partition_func == partition_function[field_size * (field_size - j - 1) + j])
		{
			if (calc_mode==1)
			{
				sum += num_of_least_energy_pathes[field_size * (field_size - j - 1) + j] * partition_function[field_size * (field_size - j - 1) + j];
			}else if(calc_mode==2){
				sum += num_of_least_energy_pathes[field_size * (field_size - j - 1) + j];
			}
		}
	}
/* 
	for (int i = 0; i < field_size; i++)//最初に取り除いた縮退度を戻す
	{
		for (int j = 0; j < field_size; j++)
		{
			partition_function[field_size * i + j] *= num_of_least_energy_pathes[field_size * i + j];
		}
	}
 */
	return sum;
}

double Field::get_growth_rate(bool parcolation)
{
	double Z = calc_pysical_quantity(1,parcolation);
	if (Z != 0)
	{
		return log(Z) / field_size;
	}
	else
	{
		return 0;
	}
}

double Field::get_entropy(bool parcolation){
	double W=calc_pysical_quantity(2,parcolation);
	return log(W)/field_size;
}
