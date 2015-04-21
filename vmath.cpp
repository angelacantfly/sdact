//
//  vmath.cpp
//  ip2015
//
//  Created by ANGELA ZHOU on 4/21/15.
//  Copyright (c) 2015 Z Sweedyk. All rights reserved.
//

#include "vmath.h"

double dotProduct(vector<double> v1, vector<double> v2)
{
    double product = 0;
    unsigned long size = v1.size();
    
    assert(v1.size() == v2.size());
    for (int i = 0; i < size; ++i)
        product += v1[i] * v2[i];
    return product;
}

double length(vector<double> v)
{
    double len = 0;
    for (int i = 0; i < v.size(); ++i)
        len += pow(v[i],2);
    len = sqrt(len);
    
    return len;
}