//
//  mvchelpers.h
//  ip2015
//
//  Created by ANGELA ZHOU on 5/3/15.
//  Copyright (c) 2015 Z Sweedyk. All rights reserved.
//

#ifndef __ip2015__mvchelpers__
#define __ip2015__mvchelpers__

#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>
#include "main.h"
#include "vmath.h"

double mvc_calculate_w (long double alpha, long double beta, Point p, Point x);
double mvc_calculate_angle(Point x, Point p1, Point p2);
vector<double> mvc_calculate_lambda(Point x, vector<Point> list, double& sum);
void mvc_calculate_lambda_array(Point x, vector<Point> list, vector<double>& matrix);
#endif /* defined(__ip2015__mvchelpers__) */
