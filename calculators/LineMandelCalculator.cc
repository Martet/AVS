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

#include <stdlib.h>


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	// @TODO allocate & prefill memory
}

LineMandelCalculator::~LineMandelCalculator() {
	// @TODO cleanup the memory
}


int * LineMandelCalculator::calculateMandelbrot () {
	// @TODO implement the calculator & return array of integers
	return NULL;
}
