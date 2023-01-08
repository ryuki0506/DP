#ifndef __CULC_H__
#define __CULC_H__

#include <vector>
using namespace std;

double limited_average(double *data, int data_len);
// limit_average(double data[],int shots);

double average(double *data, int data_len);

double Sofp(int len, int lenmax, double p, int shots, int noise_mode, int calc_mode, bool show_in_terminal, bool Isfixed, bool parcolation);

void WofE(double *SofE_all_path,int Emax,int len, double p, int shots, int noise_mode, bool Isfixed);

void WofE_min(double *SofE_min,int Emax,int len, double p, int shots, int noise_mode, bool Isfixed,bool parcolation);


#endif