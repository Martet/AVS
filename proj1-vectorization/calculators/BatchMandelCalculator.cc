/**
 * @file BatchMandelCalculator.cc
 * @author Martin Zmitko <xzmitk01@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 2023-11-05
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)aligned_alloc(64, height * width * sizeof(int));
	block_real = (float *)aligned_alloc(64, block_size * sizeof(float));
	block_imag = (float *)aligned_alloc(64, block_size * sizeof(float));
	block_real_start = (float *)aligned_alloc(64, width * sizeof(float));
	for (int i = 0; i < width; i++)
	{
		block_real_start[i] = x_start + i * dx;
	}
	for (int i = 0; i < width * height; i++)
	{
		data[i] = limit;
	}
}

BatchMandelCalculator::~BatchMandelCalculator() {
	free(data);
	free(block_real);
	free(block_imag);
	free(block_real_start);
	data = NULL;
	block_real = NULL;
	block_imag = NULL;
	block_real_start = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *pblock_real = block_real, *pblock_imag = block_imag, *pblock_real_start = block_real_start;
	for (int i = 0; i < height / 2; i++)
	{
		float zImag = y_start + i * dy; // current imaginary value
		for (int j_block = 0; j_block < width / block_size; j_block++)
		{
			int data_addr = j_block * block_size + i * width;
			int start_addr = j_block * block_size;
			std::memcpy(pblock_real, pblock_real_start + start_addr, block_size * sizeof(float));
			for (int j = 0; j < block_size; j++)
			{
				pblock_imag[j] = zImag;
			}
			for (int n = 0; n < limit; n++)
			{
				int flag = block_size;
				#pragma omp simd aligned(pdata, pblock_real, pblock_imag, pblock_real_start:64) reduction(-:flag) simdlen(64) safelen(64)
				for (int j = 0; j < block_size; ++j)
				{

					float r2 = pblock_real[j] * pblock_real[j];
					float i2 = pblock_imag[j] * pblock_imag[j];

					if (r2 + i2 > 4.0f && pdata[j + data_addr] == limit) {
						flag--;
						pdata[j + data_addr] = n;
					}

					pblock_imag[j] = 2.0f * pblock_real[j] * pblock_imag[j] + zImag;
					pblock_real[j] = r2 - i2 + pblock_real_start[j + start_addr];
				}
				if (flag == 0)
					break;
			}
		}
	}
	pdata = data + (height / 2) * width;
	for (int i = 0; i < height / 2; i++)
	{
		std::memcpy(pdata, data + (height / 2 - i - 1) * width, width * sizeof(int));
		pdata += width;
	}
	return data;
}
