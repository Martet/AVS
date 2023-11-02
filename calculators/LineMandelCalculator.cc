/**
 * @file LineMandelCalculator.cc
 * @author Martin Zmitko <xzmitk01@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 2023-11-03
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
	data = (int *)malloc(height * width * sizeof(int));
	line_real = (float *)aligned_alloc(64, width * sizeof(float));
	line_imag = (float *)aligned_alloc(64, width * sizeof(float));
	line_real_start = (float *)aligned_alloc(64, width * sizeof(float));
	line = (int *)aligned_alloc(64, width * sizeof(int));
	for (int i = 0; i < width; i++)
	{
		line_real_start[i] = x_start + i * dx;
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
	for (int i = 0; i < height / 2; i++)
	{
		float zImag = y_start + i * dy; // current imaginary value
		std::memcpy(line_real, line_real_start, width * sizeof(float));
		std::memset(line_imag, 0, width * sizeof(float));
		for (int j = 0; j < width; j++)
		{
			line[j] = limit;
		}
		for (int n = 0; n < limit; n++)
		{
			int flag = width;
			#pragma omp simd aligned(line_real, line_imag, line_real_start: 64) reduction(-:flag)
			for (int j = 0; j < width; ++j)
			{
				float r2 = line_real[j] * line_real[j];
				float i2 = line_imag[j] * line_imag[j];

				if (r2 + i2 > 4.0f && line[j] == limit) {
					flag--;
					line[j] = n;
				}
				if (line[j] == limit) {
					line_imag[j] = 2.0f * line_real[j] * line_imag[j] + zImag;
					line_real[j] = r2 - i2 + line_real_start[j];
				}
			}
			if (flag == 0)
				break;
		}
		std::memcpy(pdata, line, width * sizeof(int));
		pdata += width;
	}
	for (int i = 0; i < height / 2; i++) {
		std::memcpy(pdata, data + (height / 2 - i - 1) * width, width * sizeof(int));
		pdata += width;
	}
	return data;
}
