//
//  helpfunctions.h
//  ip2015
//
//  Created by ANGELA ZHOU on 5/3/15.
//  Copyright (c) 2015 Z Sweedyk. All rights reserved.
//

#ifndef __ip2015__helpfunctions__
#define __ip2015__helpfunctions__

#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>
#include "main.h"
#include "vmath.h"

vector<Point> find_boundary(int width, int height);
void find_vertices_square(int width, int height, Point *top_left, Point *top_right, Point *bot_left, Point *bot_right);
#endif /* defined(__ip2015__helpfunctions__) */
