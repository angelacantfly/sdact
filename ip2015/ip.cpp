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

void computeLambda()
{
    cout << "HERE" << endl;
    
    // for fixed square  in duck
    Point upperLeft, upperRight, lowerLeft, lowerRight;
    upperLeft.x=0;
    upperLeft.y=0;
    upperRight.x=10;
    upperRight.y=0;
    lowerLeft.x=0;
    lowerLeft.y=10;
    lowerRight.x=10;
    lowerRight.y = 10;
    
    int boundaryLength=11;
    int perimeter=4*(boundaryLength-1);
    Point boundary[(boundaryLength-1)*4];
    for (int i=0;i<boundaryLength;i++) {
        // left
        boundary[i].x=0;
        boundary[i].y=i;
        
        // bottom
        boundary[i+boundaryLength].x=i;
        boundary[i+boundaryLength].y=boundaryLength;
        
        // right
        boundary[i+boundaryLength*2].x=boundaryLength;
        boundary[i+boundaryLength*2].y = boundaryLength-i;
        
        // top
        boundary[i+boundaryLength*3].x=boundaryLength-i;
        boundary[i+boundaryLength*3].y=0;
    }
    for (int i=0; i<perimeter; i++)
        cout << boundary[i].x << " " << boundary[i].y << endl;
//
//    Point interior[(boundaryLength-2)*(boundaryLength-2)];
//    for (int i=0; i<boundaryLength-1; i++) {
//        for (int j=0; j<boundaryLength-1; j++) {
//            interior[i*(boundaryLength-1)+j].x = i;
//            interior[i*(boundaryLength-1)+j].y = j;
//        }
//    }
//    
//    double lambda[(boundaryLength-2)*(boundaryLength-2)][perimeter];
//    
//    for (int i=0; i<(boundaryLength-2)*(boundaryLength-2); i++) {
//        double sum = 0;
//        for (int j=0; j<perimeter; j++) {
//            // pi-1 = p0, pi = p1, pi+1=p2
//            Point p0 = boundary[j];
//            Point p1 = boundary[(j+1)%perimeter];
//            Point p2 = boundary[(j+2)%perimeter];
//            
//            Point inner = interior[i];
//            
//            vector<double> ip0, ip1, ip2;
//            ip0.push_back(inner.x - p0.x);
//            ip0.push_back(inner.y - p0.y);
//            ip1.push_back(inner.x - p1.x);
//            ip1.push_back(inner.y - p1.y);
//            ip2.push_back(inner.x - p2.x);
//            ip2.push_back(inner.y - p2.y);
//            
//            if (i==80 && j==34) {
//                cout << "stop" << endl;
//            }
//            
//            
//            double alpha = acos(dotProduct(ip0, ip1)/(length(ip0)* length(ip1)));
//            double beta = acos(dotProduct(ip1, ip2)/(length(ip1)* length(ip2)));
//            lambda[i][j] = (tan(alpha/2) + tan(beta/2)) / (length(ip1));
//            sum += lambda[i][(j+1)%perimeter];
//            
//            
//        }
//        
//        for (int j=0; j<perimeter; j++) {
//            lambda[i][j] = lambda[i][j] / sum;
//            cout << i << " " << j << " " << lambda[i][j] << endl;
//        }
//        
//    }
//
    
}


Image* mvc_cloning_square(Image* source, Image* patch)
{
    cout << "HERE";
    int width = source->getWidth();
    int height = source->getHeight();
    Image* newImage =  new Image(*source);
    
    // for fixed square  in duck
    Point upperLeft, upperRight, lowerLeft, lowerRight;
    upperLeft.x=0;
    upperLeft.y=0;
    upperRight.x=10;
    upperRight.y=0;
    lowerLeft.x=0;
    lowerLeft.y=10;
    lowerRight.x=10;
    lowerRight.y = 10;
    
    int boundaryLength=11;
    int perimeter=4*(boundaryLength-1);
    Point boundary[(boundaryLength-1)*4];
    for (int i=0;i<boundaryLength;i++) {
        // left
        boundary[i].x=0;
        boundary[i].y=i;
        
        // bottom
        boundary[i+boundaryLength].x=i;
        boundary[i+boundaryLength].y=boundaryLength;
        
        // right
        boundary[i+boundaryLength*2].x=boundaryLength;
        boundary[i+boundaryLength*2].y = boundaryLength-i;
        
        // top
        boundary[i+boundaryLength*3].x=boundaryLength-i;
        boundary[i+boundaryLength*3].y=0;
    }
    for (int i=0; i<perimeter; i++)
        cout << boundary[i].x << " " << boundary[i].y << endl;
    
    Point interior[(boundaryLength-2)*(boundaryLength-2)];
    for (int i=0; i<boundaryLength-1; i++) {
        for (int j=0; j<boundaryLength-1; j++) {
            interior[i*(boundaryLength-1)+j].x = i;
            interior[i*(boundaryLength-1)+j].y = j;
        }
    }
    
    double lambda[(boundaryLength-2)*(boundaryLength-2)][perimeter];
    
    for (int i=0; i<(boundaryLength-2)*(boundaryLength-2); i++) {
        double sum = 0;
        for (int j=0; j<perimeter; j++) {
            // pi-1 = p0, pi = p1, pi+1=p2
            Point p0 = boundary[j];
            Point p1 = boundary[(j+1)%perimeter];
            Point p2 = boundary[(j+2)%perimeter];
            
            Point inner = interior[i];
            
            vector<double> ip0, ip1, ip2;
            ip0.push_back(inner.x - p0.x);
            ip0.push_back(inner.y - p0.y);
            ip1.push_back(inner.x - p1.x);
            ip1.push_back(inner.y - p1.y);
            ip2.push_back(inner.x - p2.x);
            ip2.push_back(inner.y - p2.y);
    
            if (i==80 && j==34) {
                cout << "stop" << endl;
            }
            
            
            double alpha = acos(dotProduct(ip0, ip1)/(length(ip0)* length(ip1)));
            double beta = acos(dotProduct(ip1, ip2)/(length(ip1)* length(ip2)));
            lambda[i][j] = (tan(alpha/2) + tan(beta/2)) / (length(ip1));
            sum += lambda[i][(j+1)%perimeter];
            
            
        }
        
        for (int j=0; j<perimeter; j++) {
            lambda[i][j] = lambda[i][j] / sum;
            cout << i << " " << j << " " << lambda[i][j] << endl;
        }
        
    }

       return newImage;
}
