#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <vector>
using namespace std;

void show_field(double *data,int size,int max_size,bool Isshow);
void output_settings(int noise_mode,int calc_mode);
void output_SofE(double *data,int size,bool Isshow);
#endif