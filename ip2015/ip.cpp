#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>
#include "vmath.h"
#include <array>

#define RED	0
#define GREEN	1
#define BLUE	2



/*
* use a mask image for a per-pixel alpha value to perform
* interpolation with a second image
*/
Image* ip_composite (Image* src1, Image* src2, 
					 Image* mask)
{
    int width = src1->getWidth();
    int height = src1->getHeight();
    
    Image* result = new Image(width, height);
    double value;
    double maskValue;
    
    for (int w = 0; w < width; ++w)
        for (int h = 0; h < height; ++h)
            for (int c= 0; c < 3; ++c) {
                maskValue = mask->getPixel(w, h, c);
                value = (1 -maskValue)* src1->getPixel(w, h, c) +  maskValue * src2->getPixel(w, h, c);
                result->setPixel(w, h, c, value);
            }
    return result;
}
/*
* create a new image with values equal to the psychosomatic intensities
* of the source image
*/
Image* ip_grey (Image* src)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    double grey;
    double RED_CONSTANT = .2126;
    double GREEN_CONSTANT = .7152;
    double BLUE_CONSTANT = .0722;
    
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            grey = RED_CONSTANT * src->getPixel(w, h, RED) +
                    GREEN_CONSTANT * src->getPixel(w, h, GREEN) +
                    BLUE_CONSTANT * src->getPixel(w, h, BLUE);
            newImage->setPixel(w, h, GREEN, grey);
            newImage->setPixel(w, h, BLUE, grey);
            newImage->setPixel(w, h, RED, grey);
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}

double mvc_calculate_w (double alpha, double beta, Point p, Point x){
    // tan(alpha/2) + tan(beta/2)
    double top = tan(alpha/2) + tan(beta/2);
    // ||p - x||, distance between two points
    double bot = sqrt(pow(p.x - x.x, 2) + pow(p.y- x.y, 2));
    
    return top/bot;
}

double mvc_calculate_angle(Point x, Point p1, Point p2)
{
    vector<double> p1_x = vector<double>(x.x - p1.x, x.y - p1.y);
    vector<double> x_p2 = vector<double>(p2.x - x.x, p2.y - x.y);
    double cosTheta = dotProduct(p1_x, x_p2)/(length(p1_x)* length(x_p2));
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

void mvc_cloning(Image* source, Image* patch)
{
    
}
