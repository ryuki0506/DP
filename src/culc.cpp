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

double average(double *data, int data_len)
{
	double sum = 0;
	for (size_t i = 0; i < data_len; i++)
	{
		sum += data[i];
	}
	return sum / data_len;
}

double Sofp(int len, int lenmax, double p, int shots, int noise_mode, int calc_mode, bool show_in_terminal, bool Isfixed, bool parcolation)
{
	double data=0;
	for (int shot = 0; shot < shots; shot++)
	{
		Field *field;
		field = new Field(len);
		field->set_potential(p, noise_mode);
		field->set_FPT();

		if (calc_mode == 1)
		{
			show_field(field->get_FPT(), len, lenmax, show_in_terminal);
			data += field->get_FPT(parcolation, Isfixed);
		}
		else if (calc_mode == 2)
		{
			show_field(field->get_num_of_least_energy_pathes(), len, lenmax, show_in_terminal);
			data += field->get_entropy(parcolation, Isfixed);
		}
		delete field;
	}

	return data/shots;
}

void SofE(double *SofE_all_path, int Emax, int len, double p, int shots, int noise_mode, bool Isfixed)
{
	for (size_t E = 0; E < Emax; E++)
	{
		SofE_all_path[E] = 0;
	}

	for (size_t shot = 0; shot < shots; shot++)
	{
		Field *field;
		field = new Field(len);
		field->set_potential(p, noise_mode);

		for (size_t i = 0; i < len; i++)
		{
			double *SofE = field->SofE_all_path(Emax, i, len - 1);
			for (size_t E = 0; E < Emax; E++)
			{
				// cout<< field->SofE_all_path(Emax,i, len-1)[E]<<endl;
				SofE_all_path[E] += SofE[E];
			}
		}

		delete field;
	}

	for (size_t E = 0; E < Emax; E++)
	{
		SofE_all_path[E] /= shots;
	}
}

void SofE_min(double *SofE_min,int Emax,int len, double p, int shots, int noise_mode, bool Isfixed,bool parcolation){
	for (size_t E = 0; E < Emax; E++)
	{
		SofE_min[E] = 0;
	}

	for (size_t shot = 0; shot < shots; shot++)
	{
		Field *field;
		field = new Field(len);
		field->set_potential(p, noise_mode);
		field->set_FPT();

		double _Emin=field->calc_Emin(parcolation);
		int Emin=int(_Emin);
		if (Emin<Emax)
		{
			SofE_min[Emin]+=field->get_entropy(parcolation, Isfixed);
		}		

		delete field;
	}

	for (size_t E = 0; E < Emax; E++)
	{
		SofE_min[E] /= shots;
	}
}
