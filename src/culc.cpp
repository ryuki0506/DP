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

double *simulation(int len, int lenmax, double p, int shots, int noise_mode, int calc_mode, bool show_in_terminal, bool Isfixed, bool parcolation)
{
	int Emax=0;
	double shots_data[shots];
	for (int shot = 0; shot < shots; shot++)
	{
		Field *field;
		field = new Field(len,Emax);
		field->set_potential(p, noise_mode);
		field->set_partition_function();

		if (calc_mode == 1)
		{
			show_field(field->get_partition_function(), len, lenmax, show_in_terminal);
			shots_data[shot] = field->get_growth_rate(parcolation, Isfixed);
		}
		else if (calc_mode == 2)
		{
			show_field(field->get_num_of_least_energy_pathes(), len, lenmax, show_in_terminal);
			shots_data[shot] = field->get_entropy(parcolation, Isfixed);
		}
		delete field;
	}
	return shots_data;
}

void SofE(double *SofE_all_path,int Emax, int len, double p, int shots, int noise_mode, bool Isfixed)
{
	for (size_t E = 0; E < Emax; E++)
	{
		SofE_all_path[E]=0;
	}
	
	Field *field;
	field = new Field(len,Emax);
	field->set_potential(p, noise_mode);

	for (size_t i = 0; i < len; i++)
	{
		for (size_t E = 0; E < Emax; E++)
		{
			SofE_all_path[E]+= field->calc_SofE_all_path(i, len-1)[E];
			field->clear_SofE();
		}
	}

	delete field;
}
