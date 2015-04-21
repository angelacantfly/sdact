#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>
#include "vmath.h"

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

void mvc_calculate_lambda(Point x, Point list[], double* result)
{
    
    
}

