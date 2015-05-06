#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>
#include "vmath.h"
#include <array>
#include "helpfunctions.h"
#include "mvchelpers.h"

#define RED	0
#define GREEN	1
#define BLUE	2


double correctChannel(double value)
{
    if (value > 1) return 1;
    if (value < 0) return 0;
    return value;
}

/*
 * convolve an image with a kernel
 */
Image* ip_convolve (Image* src, int size, double* kernel )
{
    int width = src->getWidth();
    int height = src->getHeight();
    assert(size%2 ==1); // makes sure that the size is odd number
    
    double currentRed;
    double currentGreen;
    double currentBlue;
    int countNeightbor = 0;
    
    int currentX;
    int currentY;
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            currentRed = 0;
            currentGreen = 0;
            currentBlue = 0;
            countNeightbor = 0;
            
            for (int nx = -(size-1)/2; nx < (size-1)/2 + 1; ++nx )
                for (int ny = -(size-1)/2; ny < (size-1)/2 + 1; ++ny)
                {
                    currentX = w + nx;
                    currentY = h + ny;
                    if (currentX < 0 || currentX >= width || currentY < 0 || currentY >= height) {
                        continue;
                    }
                    currentRed += src->getPixel(currentX, currentY, RED) * kernel[countNeightbor];
                    currentGreen += src->getPixel(currentX, currentY, GREEN) * kernel[countNeightbor];
                    currentBlue += src->getPixel(currentX, currentY, BLUE) * kernel[countNeightbor];
                    countNeightbor++;
                }
            
            newImage->setPixel(w, h, RED, correctChannel(currentRed) );
            newImage->setPixel(w, h, GREEN, correctChannel(currentGreen));
            newImage->setPixel(w, h, BLUE, correctChannel(currentBlue));
            
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}


/*
 * convolve with an edge detection kernel
 */
Image* ip_edge_detect (Image* src)
{
    //Normal Edge Detect
    double* kernel = new double[9];
    
    for(int i = 0; i < 9; i++){
        kernel[i] = -1;
    }
    
    kernel[4] = 8;
    
    return ip_convolve(src, 3, kernel);
    
}

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

Image* ip_select_center_square(Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* newImage =  new Image(*src);
    
    vector<Point> bounds = find_boundary(width, height);
    Pixel pink = Pixel(1,0,0.87);
    int x,y;
    for (int i = 0; i < bounds.size(); i++) {
        x = bounds[i].x;
        y = bounds[i].y;
        newImage->setPixel(x, y, pink);
    }
    return newImage;
}


void mvc_cloning(Image* source, Image* patch)
{
    int width = source->getWidth();
    int height = source->getHeight();
    Image* newImage =  new Image(*source);
    
    // Source patch image
    vector<Point> boundary_source = find_boundary(width, height);
    vector<Point> all_pixel_source;
    unsigned long size_all_pixel_source = all_pixel_source.size();
    
    vector<vector<double>> all_pixel_result; // FIXME: good idea?
    
    //Preprocessing stage
    for (int i = 0; i < size_all_pixel_source; i++) {
        double sum = 0;
        all_pixel_result[i] = mvc_calculate_lambda(all_pixel_source[i], boundary_source, sum);
    }
    
    // Target Image
    // Compute the differences along the boundary
    vector<Point> boundary_target;
    unsigned long size_bt = boundary_target.size();
    for (int k = 0; k < size_bt; ++k) {
        // Compute the differences along the boundary
        
        // Evaluate the mean-value interpolant at x
    }
}


Image* mvc_cloning_square(Image* source, Image* patch)
{
    int width = source->getWidth();
    int height = source->getHeight();
    Image* newImage =  new Image(*source);
    vector<Point> boundary_source = find_boundary(width, height);
    vector<Point> list = boundary_source;
//    cout << "size: "<< boundary_source.size() << endl;
//    for (int i = 0; i < boundary_source.size(); i++)
//    {
//        cerr << boundary_source[i].x << ","<< boundary_source[i].y <<  endl ;
//    }

    Point top_left, top_right, bot_left, bot_right;
    find_vertices_square(width, height, &top_left, &top_right, &bot_left, &bot_right);
    top_left.x +=1;
    top_left.y +=1;
    top_right.x -=1;
    top_right.y +=1;
    bot_left.x +=1;
    bot_left.y -=1;
    bot_right.x -=1;
    bot_right.y -=1;
    double size_of_square = top_right.x - top_left.x;
    double size_of_inner_pixels = size_of_square -2;
    cout << "size of inner pixels :" << size_of_inner_pixels<< endl;
    cout << "top left x: " << top_left.x << ", y :" << top_left.y << endl;
    cout << "bot_left.y : " << bot_left.y << endl;
    cout << "top_right.x : " << top_right.x << endl;
    
    vector<vector<vector<double>>> lambda_result;
    lambda_result.resize(size_of_inner_pixels);
    for (int i = 0; i <= size_of_inner_pixels; i++) {
        lambda_result[i].resize(size_of_inner_pixels);
        for (int j = 0; j <= size_of_inner_pixels ; j++) {
            lambda_result[i][j].resize(boundary_source.size());
        }
    }
    Point current;
    double row, col;
    unsigned long size = boundary_source.size();
    for (int i = top_left.x+1 ; i < top_right.x; i++) {
        row = i - (top_left.x+1);
        for (int j = top_left.y+1; j < bot_left.y; j++) {
            col = j - (top_left.y+1);
//            cout << "row : " << row << " , col :" << col << endl;
            current.x = i;
            current.y = j;
            long double alpha, beta;
            long double sum = 0;

           
            alpha = mvc_calculate_angle(current, boundary_source[size-1], boundary_source[0]);
            beta = mvc_calculate_angle(current, boundary_source[0], boundary_source[1]);
            lambda_result[row][col][0] = mvc_calculate_w(alpha, beta, boundary_source[0], current);
            sum += lambda_result[col][col][0];
            
            alpha = mvc_calculate_angle(current, boundary_source[size-2], boundary_source[size-1]);
            beta = mvc_calculate_angle(current, boundary_source[size-1], boundary_source[0]);
            lambda_result[row][col][size-1] = mvc_calculate_w(alpha, beta, boundary_source[size-1], current);
            sum+= lambda_result[row][col][size-1];
//            cout << i << ","<< j << ","<< size-1 << ": "<<lambda_result[row][col[size-1]<< endl;
//            cout << "SUM :" << sum << endl;

//            if (lambda_result[row][col[size-1] != lambda_result[row][col[size-1]) {
//                cout << i << ", " << j << ": "<< lambda_result[row][col[0] << endl;
//            }
            
            for (int k = 1 ; k < size-2; k++) {
                alpha = mvc_calculate_angle(current, list[k-1], list[k]);
                beta = mvc_calculate_angle(current, list[k], list[k+1]);
                lambda_result[row][col][k] = mvc_calculate_w(alpha, beta, list[k], current);
                sum += lambda_result[row][col][k];
                

                if (lambda_result[row][col][k] != lambda_result[row][col][k]) {
                    
                    cout << "alpha :" << alpha <<endl;
                    cout << "beta :" << beta << endl;
                    cout << "SUM :" << sum << endl;
                    cout << i << ","<< j << ","<< k << ": "<<lambda_result[row][col][k]<< endl;
                    cout << "NANANAN!" << endl;
                }
            }
            
            for (int k = 0; k < size -1; k++) {
                lambda_result[row][col][k] = lambda_result[row][col][k] / sum;
                if (lambda_result[row][col][k] != lambda_result[row][col][k]) {
                    cout << i << ", "<< j << ", " << k  << " - "<<  lambda_result[row][col][k] << sum << endl;
                }
            }
            
        }
        cout << row << "/" << size_of_inner_pixels << endl;
    }
    cout << "Mvc calculated"<< endl;
    
    return newImage;
}
