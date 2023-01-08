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
	double data = 0;
	for (int shot = 0; shot < shots; shot++)
	{
		Field *field;
		field = new Field(len);
		field->set_potential(p, noise_mode);
		field->time_evolution();

		if (calc_mode == 1)
		{
			show_field(field->get_FPT(), len, lenmax, show_in_terminal);
			data += field->calc_FPT(parcolation, Isfixed);
		}
		else if (calc_mode == 2)
		{
			show_field(field->get_W_Emin(), len, lenmax, show_in_terminal);
			data += field->calc_W_Emin(parcolation, Isfixed);
		}
		delete field;
	}

	return data / shots;
}

void SofE(double *WofE_all_path, int Emax, int len, double p, int shots, int noise_mode, bool Isfixed)
{
	for (size_t E = 0; E < Emax; E++)
	{
		WofE_all_path[E] = 0;
	}

	for (size_t shot = 0; shot < shots; shot++)
	{
		Field *field;
		field = new Field(len);
		field->set_potential(p, noise_mode);

		for (size_t E = 0; E < Emax; E++)
		{
			WofE_all_path[E] += field->WofE_all_path(Emax)[E];
		}

		delete field;
	}

	for (size_t E = 0; E < Emax; E++)
	{
		WofE_all_path[E] /= shots;
	}
}

void WofE_min(double *WofE_min, int Emax, int len, double p, int shots, int noise_mode, bool Isfixed, bool parcolation)
{
	for (size_t E = 0; E < Emax; E++)
	{
		WofE_min[E] = 0;
	}

	for (size_t shot = 0; shot < shots; shot++)
	{
		Field *field;
		field = new Field(len);
		field->set_potential(p, noise_mode);
		field->time_evolution();

		double _Emin = field->calc_Emin(parcolation, Isfixed);
		int Emin = int(_Emin);
		if (Emin < Emax)
		{
			WofE_min[Emin] += field->calc_W_Emin(parcolation, Isfixed);
		}

		delete field;
	}

	for (size_t E = 0; E < Emax; E++)
	{
		WofE_min[E] /= shots;
	}
}
