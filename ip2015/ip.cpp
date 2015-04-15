#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>

#define RED	0
#define GREEN	1
#define BLUE	2




/*
* convolve with a box filter
*/
Image* ip_blur_box (Image* src, int size)
{
    int kernelSize = size*size;
    double* kernel = new double[kernelSize];
    
    for(int i = 0; i < kernelSize; i++){
        kernel[i] = ((double)1)/kernelSize;
        cerr << kernel[i] << endl;
    }
    
    return ip_convolve(src, size, kernel);
    
}


/*
* convolve with a gaussian filter
*/
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    int kernelSize = size*size;
    double* kernel = new double[kernelSize];
    int count = 0;
    double sum = 0;
    
    for (int nx = -(size-1)/2; nx < (size-1)/2 + 1; ++nx )
        for (int ny = -(size-1)/2; ny < (size-1)/2 + 1; ++ny)
        {
            double coord = pow(nx, 2)+ pow(ny, 2);
            kernel[count] = exp(-coord/(2*pow(sigma, 2)));
            sum += kernel[count];
            cerr<< kernel[count] << endl;
            ++count;
        }
    
    for (int i = 0; i < kernelSize; ++i) {
        kernel[i] = kernel[i]/sum;
        cerr<< kernel[count] << endl;
    }
    
    return ip_convolve(src, size, kernel);
}


/*
* convolve with a triangle filter
*/
Image* ip_blur_triangle (Image* src, int size)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* interpolate with a black image
*/
Image* ip_brighten (Image* src, double alpha)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    Image* blackImage = new Image(width,height);
    
    if (alpha > 0 && alpha <= 1)
        return ip_interpolate(src, blackImage, alpha);
    
    if (alpha > 1)
        return ip_interpolate(src, blackImage, alpha);

    return NULL;
    
}


/*
* shift colors
*/
Image* ip_color_shift(Image* src)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            newImage->setPixel(w, h, GREEN, src->getPixel(w, h, RED));
            newImage->setPixel(w, h, BLUE, src->getPixel(w, h, GREEN));
            newImage->setPixel(w, h, RED, src->getPixel(w, h, BLUE));
        }
    }
    
	cerr << "Done!" << endl;
	return newImage;
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
* interpolate with the average intensity of the src image
*/
Image* ip_contrast (Image* src, double alpha)
{
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* grey = new Image(width, height);
    

    for (int w = 0; w <width; ++w)
        for (int h = 0; h<height; ++h )
        {
            grey->setPixel(w, h, RED, 0.5);
            grey->setPixel(w, h, GREEN, 0.5);
            grey->setPixel(w, h, BLUE, 0.5);
        }

    
    return ip_interpolate(src, grey, alpha);

}

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
*  create cropped version of image
*/
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    // Create image of proper size
    Image* newImage =  new Image(x1- x0, y1-y0);
    
    cerr<< "size is" << x1-x0 << endl;
    
    for(int x = x0; x < x1; x++){
        for(int y = y0; y < y1; y++){
           
            newImage->setPixel(x-x0, y-y0, RED, src->getPixel(x, y, RED));
            newImage->setPixel(x-x0, y-y0, GREEN, src->getPixel(x, y, GREEN));
            newImage->setPixel(x-x0, y-y0, BLUE, src->getPixel(x, y, BLUE));
            
        }
    }
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
* extract channel of input image
*/
Image* ip_extract (Image* src, int channel)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            for (int c = 0; c < 3; ++c)
                if (c == channel) newImage->setPixel(w, h, c, src->getPixel(w, h, c));
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}


/*
* create your own fun warp
*/
Image* ip_fun_warp (Image* src, int samplingMode)
{
    cerr << "Sampling method is " << samplingMode << endl;
    int width = src->getWidth();
    int height = src->getHeight();
    Image* result = new Image(width, height);
    
    Pixel currentPixel;
    double x;
    double y;
    double radius;
    
    for (int w = 0; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            
            // Convert coordinate to be on a scale [0,1]
            x = 2 * ((double)w ) / (width - 1);
            y = 2 * ((double) h) / (height -1);
            radius = sqrt(pow(x - 0.5 , 2) + pow(y - 0.5,2));

            double a = atan2(x - 0.5, y - 0.5);
            double rn = pow(radius, 2)/2;
            
            y = rn * cos(a) + 0.5;
            x = rn * sin(a) + 0.5;
            
            x = x * (width - 1)/2;
            y = y * (height -1)/2;
            
            if (samplingMode == I_BILINEAR)
                currentPixel = ip_resample_bilinear(src, x, y);
            else if (samplingMode == I_GAUSSIAN)
                currentPixel = ip_resample_gaussian(src, x, y, 3, 1.0);
            else
                currentPixel = ip_resample_nearest(src, x, y);
            
            result->setPixel(w, h, currentPixel);
        }
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


/*
*  shift image by dx and dy (modulo width & height)
*/

Image* ip_image_shift (Image* src, double dx, double dy)
{
//    (i+dx)%width, (j+dy)%height.
    
    int width = src->getWidth();
    int height = src->getHeight();
    
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            for (int c = 0; c < 3; ++c) {
                
                int x = ( w + (int) dx ) % width;
                int y = ( h + (int) dy ) % height;
                newImage->setPixel(w, h,c,src->getPixel(x,y,c) );
            }
            
            
        }}
	
    return newImage;
}
/*
* interpolate an image with another image
*/
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    //FIXME: negative alpha value freaks out and doesn't end
    // Only seems to be a problem with brighten function
    // works for saturate and contrast
    
    // get width and height
    int width = src1->getWidth();
    int height = src1->getHeight();
    Pixel src1Pixel;
    Pixel src2Pixel;
    double currentChannel;
    
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            src1Pixel = src1->getPixel(w, h);
            src2Pixel = src2->getPixel(w, h);
            for (int c = 0; c < 3; ++c)
            {
                currentChannel = alpha * src1Pixel.getColor(c) + (1.0-alpha)* src2Pixel.getColor(c);
                newImage->setPixel(w, h, c, correctChannel(currentChannel));

            }
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}
/*
* invert input image
*/
Image* ip_invert (Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    

    Image* grey = new Image(width,height);
    
    for (int w =0; w<width; ++w)
        for (int h=0; h<height; ++h)
        {
            grey->setPixel(w, h, RED, 0.5);
            grey->setPixel(w, h, GREEN, 0.5);
            grey->setPixel(w, h, BLUE, 0.5);
        }
    
    return ip_interpolate(src, grey, -1);
}


/*
* define your own filter
* you need to request any input paraters here, not in control.cpp
*/

Image* ip_misc(Image* src, double gamma)
{
    
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            for (int c = 0; c < 3; ++c) {
                newImage->setPixel(w, h, c, pow(src->getPixel(w, h, c),1.0/gamma));
            }
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;


}

Image* ip_misc_gamma(Image* src, double gamma)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            for (int c = 0; c < 3; ++c) {
                newImage->setPixel(w, h, c, pow(src->getPixel(w, h, c),1.0/gamma));
            }
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}

Image* ip_misc_tv(Image* src)
{
    
    double* kernel = new double[9];
    kernel[0] = -1;
    kernel[1] = 0;
    kernel[2] = 1;
    kernel[3] = -2;
    kernel[4] = 0;
    kernel[5] = 2;
    kernel[6] = -1;
    kernel[7] = 0;
    kernel[8] = 1;
    
    Image* intermediate = ip_convolve(src, 3, kernel);
    
    double* kernel2 = new double[9];
    kernel[0] = -1;
    kernel[1] = -2;
    kernel[2] = -1;
    kernel[3] = 0;
    kernel[4] = 0;
    kernel[5] = 0;
    kernel[6] = 1;
    kernel[7] = 2;
    kernel[8] = 1;
    
    return ip_convolve(intermediate, 3, kernel2);
}

Image* ip_misc_sobel(Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    
    double* kernel = new double[9];
    kernel[0] = -1;
    kernel[1] = 0;
    kernel[2] = 1;
    kernel[3] = -2;
    kernel[4] = 0;
    kernel[5] = 2;
    kernel[6] = -1;
    kernel[7] = 0;
    kernel[8] = 1;
    
    Image* g_x = ip_convolve(src, 3, kernel);
    
    double* kernel2 = new double[9];
    kernel[0] = -1;
    kernel[1] = -2;
    kernel[2] = -1;
    kernel[3] = 0;
    kernel[4] = 0;
    kernel[5] = 0;
    kernel[6] = 1;
    kernel[7] = 2;
    kernel[8] = 1;
    
    Image* g_y = ip_convolve(src, 3, kernel2);
    
    Image* result = new Image(width, height);
    double resultPixel = 0;
    
    for (int w = 0; w < width; ++w)
        for (int h = 0; h < height; ++h)
            for (int c = 0; c < 3; ++c)
            {
                resultPixel = sqrt(pow(g_x->getPixel(w, h, c), 2) + pow(g_y->getPixel(w, h, c), 2));
                result->setPixel(w, h, c, correctChannel(resultPixel));
            }
    return result;
}
/*
* round each pixel to the nearest value in the new number of bits
*/
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
    
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image* newImage = new Image(width, height, bitsPerChannel);
    Pixel currentPixel;
    
    for (int w =0; w<width; ++w)
        for (int h=0; h<height; ++h)
        {
            newImage->setPixel(w, h, src->getPixel(w, h, currentPixel));
        }
    
    return newImage;
}

/*
 *
 *
 */
double quantize_fs_helper(double value, double* threshold){
    
    int i = 0;
    double newValue;
    
    while (value > threshold[i])        i++;
    
    if((abs(value - threshold[i-1])) <= (abs(value - threshold[i])))
        newValue = threshold[i-1];
    else
        newValue = threshold[i];
    
    return newValue;
}



/*
* dither each pixel to the nearest value in the new number of bits
* using a static 2x2 matrix
*/

Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{

//    int numPixelIntensity = pow(2, bitsPerChannel);
//    double intensityArray[numPixelIntensity];
//    intensityArray[0] = 0;
//    intensityArray[numPixelIntensity-1]= 1;
//    for (int i = 1; i < numPixelIntensity - 1; ++i){
//        intensityArray[i] = ((double) i )/(numPixelIntensity-1);
//    }
//
//    for (int i = 0; i < numPixelIntensity; ++i)
//        cerr << "intensityArray at "<< i << " " << intensityArray[i] << endl;
//    
//    int numBlocks = (pow(2, bitsPerChannel) -2)*4 + 5;
//    double* blockOutputLevels = new double[numBlocks];
//    blockOutputLevels[0]= 0;
//    blockOutputLevels[numBlocks-1] = 1.00;
//    for (int i = 1; i < numBlocks -1; ++i)
//        blockOutputLevels[i] = ((double) i)/(numBlocks -1);
//    for (int i = 0; i < numBlocks; ++ i)
//        cerr << "blockOutputLevels at "<< i << " " << blockOutputLevels[i] << endl;
//
//    
//    double* blockThreshold = new double[numBlocks + 1];
//    blockThreshold[0] = 0;
//    
//    for(int i = 1; i < numBlocks  ; i++)
//        blockThreshold[i] = (blockOutputLevels[i-1] + blockOutputLevels[i])/2;
//    for (int i = 0; i < numBlocks ; i++)
//        cerr << "blockThreshold at "<< i << " " << blockThreshold[i] << endl;
//
//    
//    int width = src->getWidth();
//    int height = src->getHeight();
//    int channels = 3;
//
//    int currentChannel = 0;
//    
//    
////    Image* newImage = new Image(width, height, bitsPerChannel);
////    Pixel currentPixel;
////    for (int w = 0; w < width; w = w + 2) {
////        for (int h = 0; h < height; h = h + 2)
////            for (int c = 0; c < 3; ++c) {
////                currentChannel = src->getPixel(w, h, c);
////            }
////    }
//    
////    return newImage;
//    cerr << "This function is not implemented." << endl;
    return NULL;
}


/*
* dither each pixel to the nearest value in the new number of bits
* using error diffusion
*/
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
    int width = src->getWidth();
    int height = src->getHeight();
    int channels = 3;
    
    double matrix [width][height][channels];

    for (int w = 0; w <width; ++w)
        for (int h = 0; h<height; ++h)
            for (int c = 0; c < 3; ++c) {
                matrix[w][h][c] = 0;
            }

    int outputLevels = pow(2, bitsPerChannel);
    double threshold[outputLevels+1];
    
    double outputArray[outputLevels];
    for(int j = 0; j < outputLevels; j++){
        outputArray[j] = j/(outputLevels-1);


    }
    
    threshold[0] = 0;
    threshold[outputLevels] = 1;
    
    for(int i = 1; i < outputLevels; i++){
        threshold[i] = (outputArray[i-1] + outputArray[i])/2;


    }
    
    
    Image* newImage = new Image(width, height, bitsPerChannel);
    Pixel currentPixel;
    
    for (int w =0; w<width; ++w)
        for (int h=0; h<height; ++h)
        {
            for(int c = 0; c < 3; c++){

                // Calculating the value for Pixel(w,h)
                double value = src->getPixel(w, h, c);
                double getError = matrix[w][h][c];
                if (w >= width or h>= height)
                    getError = 0;
                double quantized = quantize_fs_helper(value, threshold) + getError;
                newImage->setPixel(w, h, c, correctChannel(quantized));
                
                // Calculate errors for Pixel(w,h) and its neighbors
                double difference = src->getPixel(w, h, c) - newImage->getPixel(w, h, c);
                
                // Right
                if (w < width - 1 and h < height)
                    matrix[w+1][h][c] += difference * (7/16);
                // Down
                if (w < width and h < height -1)
                    matrix[w][h+1][c] += difference * (5/16);
                // Down left
                if (0 < w < width and h < height -1)
                    matrix[w-1][h+1][c] += difference * (3/16);
                // Down right
                if (w < width -1 and h < height -1)
                    matrix[w+1][h+1][c] += difference * (1/16);

            }
            
            newImage->setPixel(w, h, src->getPixel(w, h, currentPixel));
        }
    
    return newImage;
}


/*
* nearest neighbor sample
*/
Pixel ip_resample_nearest(Image* src, double x, double y) {
    double roundx = floor(x);
    double roundy = floor(y);
    if (x - roundx >= 0.5)
        roundx += 1;
    if (y - roundy >= 0.5)
        roundy +=1;
    
    Pixel resultPixel;
    if (roundx >= src->getWidth() or roundy >= src->getHeight()or roundx< 0 or roundy <0)
        return Pixel(0,0,0);
	return src->getPixel(roundx, roundy, resultPixel);
}

void ip_resample_nearest2( double x, double y, double& xx, double& yy)
{
    double roundx = floor(x);
    double roundy = floor(y);
    if (x - roundx >= 0.5)
        roundx += 1;
    if (y - roundy >= 0.5)
        roundy +=1;
    
    xx = roundx;
    yy = roundy;
}
/*
* bilinear sample
*/

Pixel ip_resample_bilinear(Image* src, double x, double y) {
    Pixel resultPixel;
    double roundx = floor(x);
    double roundy = floor(y);
    double leftTop = 0;
    double rightTop = 0;
    double leftBot = 0;
    double rightBot = 0;
    double value1 = 0;
    double value2 = 0;
    double finalValue = 0;
    
    int width = src->getWidth();
    int height = src->getHeight();
    
    for (int c = 0; c <3; ++c) {
        // Bound checking
        if (roundx < width and roundy < height)
        {
            leftTop = src->getPixel(roundx, roundy, c);

            if (roundx + 1 < width)
            {
                rightTop = src->getPixel(roundx+1, roundy, c);
                if (roundy + 1 < height)
                    rightBot = src->getPixel(roundx+1, roundy+1, c);
                else
                    rightBot = 0;
            }
            else
            {
                rightTop = 0;
                rightBot = 0;
                leftTop =0;
            }
            
            if (roundy + 1 < height)
                leftBot = src->getPixel(roundx, roundy+1, c);
            else
                leftBot = 0;
        }

        
        
        value1 = leftTop * (roundx + 1 - x) + rightTop * (x - roundx);
        value2 = leftBot * (roundx + 1 - x) + rightBot * (x - roundx);
        
        finalValue = value1 * (roundy + 1 - y) + value2 * (y - roundy);
        resultPixel.setColor(c, finalValue);
    }
    
    if (roundx >= src->getWidth() or roundy >= src->getHeight())
        return Pixel(0,0,0);
    
	return resultPixel;
}

/*
* gausian sample
*/
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma)
{
    Pixel resultPixel;
    int i;
    int j;
    int roundx = floor(x);
    int roundy = floor(y);
    
    for (int c = 0; c < 3; ++c) {
        double value = 0;
        double sum = 0;
        for (int dx = -((double) size)/2; dx < ((double) size)/2; ++dx) {
            for (int dy = -((double) size)/2; dy < ((double) size)/2; ++dy) {
                i = roundx + dx;
                j = roundy + dy;
                sum += exp(-(pow(x-i, 2) + pow(y-j, 2))/(2 * pow(sigma, 2)));
                value += src->getPixel_(i, j, c) * exp(-(pow(x-i, 2) + pow(y-j, 2))/(2 * pow(sigma, 2)));
                
            }
        }
        resultPixel.setColor(c, value/sum);
    }
    if (roundx >= src->getWidth() or roundy >= src->getHeight())
        return Pixel(0,0,0);
    return resultPixel;
}

/*
* rotate image using one of three sampling techniques
*/
Image* ip_rotate (Image* src, double theta, int x, int y, int samplingMode, 
				  int gaussianFilterSize, double gaussianSigma)
{
    
    int width = src->getWidth();
    int height = src->getHeight();
    double xx;
    double yy;
    theta = theta / 180 * M_PI;
    cerr << "Theta is " << theta << endl;
    
    Image* newImage =  new Image(width, height);

    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            xx = cos(theta) * ( w - x) + sin(theta)* (y - h)  + x;
            yy = -((- sin(theta) * ( w - x) + cos(theta)* (y - h) ) - y) ;
            
            Pixel currentPixel;
            
            if ( xx > 0 and xx < width and yy > 0 and yy < height) {
                if (samplingMode == I_BILINEAR)
                    currentPixel = ip_resample_bilinear(src, xx, yy);
                else if (samplingMode == I_GAUSSIAN)
                    currentPixel = ip_resample_gaussian(src, xx, yy, gaussianFilterSize, gaussianSigma);
                else
                {
                    currentPixel = ip_resample_nearest(src, xx, yy);
                }
            }
            else{
                currentPixel = Pixel(0, 0, 0);
            }
            
            newImage->setPixel(w, h, currentPixel);
        }
    }

    return newImage;
    
}


/*
* change saturation
*/
Image* ip_saturate (Image* src, double alpha)
{
    Image* greyImage = ip_grey(src);

    return ip_interpolate(src, greyImage, alpha);

}


/*
* scale image using one of three sampling techniques
*/
Image* ip_scale (Image* src, double xFac, double yFac, int samplingMode, 
				 int gaussianFilterSize, double gaussianSigma)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* threshold image
*/
Image* ip_threshold (Image* src, double cutoff)
{
    // get width and height
    int width = src->getWidth();
    int height = src->getHeight();
    
    
    Image* newImage =  new Image(width, height);
    
    for (int w = 0 ; w < width; ++w) {
        for (int h = 0; h < height; ++h) {
            for (int c = 0; c < 3; ++ c)
                if (src->getPixel(w, h, c) > cutoff) {
                    newImage->setPixel(w, h, c, 1);
                }
        }
    }
    
    cerr << "Done!" << endl;
    return newImage;
}


