
//
//  mvchelpers.cpp
//  ip2015
//
//  Created by ANGELA ZHOU on 5/3/15.
//  Copyright (c) 2015 Z Sweedyk. All rights reserved.
//

#include "mvchelpers.h"

double mvc_calculate_w (long double alpha, long double beta, Point p, Point x){
    // tan(alpha/2) + tan(beta/2)
    double top = tan(alpha/2) + tan(beta/2);
    // ||p - x||, distance between two points
    double bot = sqrt(pow(p.x - x.x, 2) + pow(p.y- x.y, 2));
    
    return top/bot;
}

double mvc_calculate_angle(Point x, Point p1, Point p2)
{
    if (x.x == p1.x && x.y == p1.y) {
        return 0;
    }
    if (x.x == p2.x && x.y == p2.y) {
        return 0;
    }
    vector<double> p1_x;
    p1_x.push_back(x.x - p1.x);
    p1_x.push_back(x.y - p1.y);
    vector<double> x_p2;
    x_p2.push_back(p2.x - x.x);
    x_p2.push_back(p2.y - x.y);
    double cosTheta = dotProduct(p1_x, x_p2)/(length(p1_x)* length(x_p2));
    if (x.x ==125) {
         cout << " acosTheta" << acos(cosTheta) << endl;
    }
    return acos(cosTheta);
}

vector<double> mvc_calculate_lambda(Point x, vector<Point> list, double& sum)
{
    vector<double> result;
    unsigned long size = list.size();
    result.reserve(size);
    double alpha, beta;
    sum = 0;    // set sum to 0 to count;
    
    alpha = mvc_calculate_angle(x, list[size-1], list[0]);
    beta = mvc_calculate_angle(x, list[0], list[1]);
    result[0] = mvc_calculate_w(alpha, beta, list[0], x);
    sum += result[0];
    
    alpha = mvc_calculate_angle(x, list[size-2], list[size-1]);
    beta = mvc_calculate_angle(x, list[size-1], list[0]);
    result[size-1] = mvc_calculate_w(alpha, beta, list[size-1], x);
    sum+= result[size-1];
    
    
    for (int i = 1 ; i < size-2; i++) {
        alpha = mvc_calculate_angle(x, list[i-1], list[i]);
        beta = mvc_calculate_angle(x, list[i], list[i+1]);
        result[i] = mvc_calculate_w(alpha, beta, list[i], x);
        sum += result[i];
    }
    
    vector<double> lambda;
    lambda.reserve(size);
    for (int k = 0; k < size -1; k++) {
        lambda[k] = result[k] / sum;
    }
    
    return lambda;
}

void mvc_calculate_lambda_array(Point x, vector<Point> list, vector<double>& matrix)
{
    unsigned long size = list.size();
    double alpha, beta;
    double sum = 0;
    
    alpha = mvc_calculate_angle(x, list[size-1], list[0]);
    beta = mvc_calculate_angle(x, list[0], list[1]);
    matrix[0] = mvc_calculate_w(alpha, beta, list[0], x);
    sum += matrix[0];
    
    alpha = mvc_calculate_angle(x, list[size-2], list[size-1]);
    beta = mvc_calculate_angle(x, list[size-1], list[0]);
    matrix[size-1] = mvc_calculate_w(alpha, beta, list[size-1], x);
    sum+= matrix[size-1];
    
    for (int i = 1 ; i < size-2; i++) {
        alpha = mvc_calculate_angle(x, list[i-1], list[i]);
        beta = mvc_calculate_angle(x, list[i], list[i+1]);
        matrix[i] = mvc_calculate_w(alpha, beta, list[i], x);
        sum += matrix[i];
    }
    
    for (int k = 0; k < size -1; k++) {
        matrix[k] = matrix[k] / sum;
    }
}
