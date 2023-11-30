/**
 * @file LineMandelCalculator.cc
 * @author Martin Zmitko <xzmitk01@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 2023-11-05
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

#include <stdlib.h>


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)aligned_alloc(64, height * width * sizeof(int));
	line_real = (float *)aligned_alloc(64, width * sizeof(float));
	line_imag = (float *)aligned_alloc(64, width * sizeof(float));
	line_real_start = (float *)aligned_alloc(64, width * sizeof(float));
	for (int i = 0; i < width; i++)
	{
		line_real_start[i] = x_start + i * dx;
	}
	for (int i = 0; i < width * height; i++)
	{
		data[i] = limit;
	}
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	free(line_real);
	free(line_imag);
	free(line_real_start);
	data = NULL;
	line_real = NULL;
	line_imag = NULL;
	line_real_start = NULL;
}


int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *pline_real = line_real, *pline_imag = line_imag, *pline_real_start = line_real_start;
	for (int i = 0; i < height / 2; i++)
	{
		float zImag = y_start + i * dy; // current imaginary value
		int data_addr = i * width;
		std::memcpy(pline_real, pline_real_start, width * sizeof(float));
		for (int j = 0; j < width; j++)
		{
			pline_imag[j] = zImag;
		}
		for (int n = 0; n < limit; n++)
		{
			int flag = width;
			#pragma omp simd aligned(pdata, pline_real, pline_imag, pline_real_start:64) reduction(-:flag) simdlen(64)
			for (int j = 0; j < width; ++j)
			{
				float r2 = pline_real[j] * pline_real[j];
				float i2 = pline_imag[j] * pline_imag[j];

				if (r2 + i2 > 4.0f && pdata[data_addr + j] == limit) {
					flag--;
					pdata[data_addr + j] = n;
				}

				pline_imag[j] = 2.0f * pline_real[j] * pline_imag[j] + zImag;
				pline_real[j] = r2 - i2 + pline_real_start[j];
			}
			if (flag == 0)
				break;
		}
	}
	pdata += width * (height / 2);
	for (int i = 0; i < height / 2; i++) {
		std::memcpy(pdata, data + (height / 2 - i - 1) * width, width * sizeof(int));
		pdata += width;
	}
	return data;
}
