#include "main.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "control.h"
#include "vmath.h"

/*
 * IMPORTANT - DO NOT CHANGE THIS FILE - IMPORTANT
 */


int  window_width  = 300;
int  window_height = 300;

Image* currentImage  = NULL;
Image* originalImage = NULL;

bool quietMode = false;
bool textMode  = false;


int main (int argc, char** argv)
{
    // initialize parameters
    char* toLoad = init(argc, argv);
    
    if (textMode)
    {
        if (toLoad)
            image_load(toLoad);
        textMenuLoop();
    }
    else
    {
        // set up the window
        glutInit(&argc, &argv[0]);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
        glutInitWindowPosition(100,100);
        glutInitWindowSize(window_width, window_height);
        glutCreateWindow("Mean Value Cloning");
        
        // register call back functions
        glutDisplayFunc(display);
        glutReshapeFunc(unreshape);
        
        glClearColor(0.0,0.0,0.0,0.0);
        glDisable(GL_DEPTH_TEST);
        
        // setup main menu
        make_menu();
        
        // register keyboard callback function
        glutKeyboardFunc(keyboard_func);
        
        if (toLoad)
            image_load(toLoad);
        
        // wait for something to happen
        try {
            glutMainLoop();
        }
        catch (const std::length_error& le) {
            std::cerr << "Length error: " << le.what() << '\n';
        }
    }
    return 0;
}


char* init (int argc, char** argv)
{
    // init random number generator
    //srand48(time(0));
    
    char* toLoad = NULL;
    
    // parse the command line options
    bool noMoreArgs  = false;
    bool noMoreFlags = false;
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            if (noMoreArgs)
                usage();
            
            if (!noMoreFlags && argv[i][0] == '-')
            {
                switch (argv[i][1])
                {
                    case 't':
                        textMode = true;
                        break;
                        
                    case 'q':
                        quietMode = true;
                        break;
                        
                    case '-':
                        if (argv[i][2] == '\0')
                            noMoreFlags = true;
                        else
                            usage();
                        break;
                        
                    default:
                        usage();
                }
            }
            else
            {
                noMoreArgs = true;
                //        image_load(argv[i]);
                toLoad = argv[i];
            }
        }
    }
    
    return toLoad;
}


void usage ()
{
    cerr << "usage: ./ip [ -t ] [ -q ] [ -- ] [ file ]" << endl;
    exit(-1);
}


void display ()
{
    if (textMode)
        return;
    
    // check if there have been any openGL problems
    GLenum errCode = glGetError();
    if (errCode != GL_NO_ERROR)
    {
        const GLubyte* errString = gluErrorString(errCode);
        cerr << "OpenGL error: " << errString << endl;
    }
    
    // clear the frame buffer
    glClear(GL_COLOR_BUFFER_BIT);
    
    // draw the image
    if (currentImage)
        currentImage->glDrawPixelsWrapper();
    
    // swap buffers
    glutSwapBuffers();
}


void unreshape (int width, int height)
{
    // don't allow user to manuall resize the window
    reshape(window_width, window_height);
}


void reshape (int width, int height)
{
    // set window height and width
    window_width  = max(width,  64);
    window_height = max(height, 64); 
    
    if (textMode)
        return;
    
    // change the actual window's size
    glutReshapeWindow(window_width, window_height);
    
    // the lower left corner of the viewport is 0,0
    // the upper right corner is width, height
    glViewport(0, 0, (GLint) window_width, (GLint) window_height);  
    
    // setup orthographic projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, window_width, 0.0, window_height);
    
    // default mode should be modelview
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

Image* computeLambda(Image* background, Image* patch)
{
    Image* result = new Image(*background);
    // for fixed square  in duck
    Point upperLeft, upperRight, lowerLeft, lowerRight;
    upperLeft.x=0;
    upperLeft.y=0;
    upperRight.x=3;
    upperRight.y=0;
    lowerLeft.x=0;
    lowerLeft.y=3;
    lowerRight.x=3;
    lowerRight.y = 3;
    
    const int boundaryLength=61;
    const int perimeter=4*(boundaryLength-1);
//    int boundary[perimeter][2];
    vector<vector<int>> boundary;
    boundary.resize(perimeter);
    for (int i = 0; i < perimeter; i++) {
        boundary[i].resize(2);
    }
//    int interiorPoints[(boundaryLength-2)][(boundaryLength-2)][2];
    vector<vector<vector<int>>> interiorPoints;
    interiorPoints.resize(boundaryLength-2);
    for (int a = 0 ; a < boundaryLength-2; a++) {
        interiorPoints[a].resize(boundaryLength-2);
        for (int b = 0 ; b < boundaryLength-2; b++) {
            interiorPoints[a][b].resize(2);
        }
    }
    
//    for (int c=0; c<perimeter; c++) {
//        cout <<  c << "("  << boundary[c][0] << "," << boundary[c][1] << ") " << endl;
//    }
    
    // finding boundary points
    for (int a = 0; a < boundaryLength-1; a++) {
        boundary[a][0] = boundaryLength - 1 - a; // top
        boundary[a][1] = 0;
        
        boundary[(boundaryLength-1) + a][0] = 0; // left
        boundary[(boundaryLength-1) + a][1] = a;
        
        boundary[2 * (boundaryLength - 1) + a][0] = a; // bottom
        boundary[2 * (boundaryLength - 1) + a][1] = boundaryLength -1;
        
        boundary[3 * (boundaryLength - 1) + a][0] = boundaryLength - 1; // right
        boundary[3 * (boundaryLength - 1) + a][1] = boundaryLength - 1 -a;
    }
    for (int c=0; c<perimeter; c++) {
        cout <<  c << "("  << boundary[c][0] << "," << boundary[c][1] << ") " << endl;
    }

    // setting interior pts
    for (int i=0; i<boundaryLength-2; i++) {
        for (int j=0; j<boundaryLength-2; j++) {
            interiorPoints[i][j][0] = i+1;
            interiorPoints[i][j][1] = j+1;
        }
    }
    
    // calculate lambda
    double lambda[(boundaryLength-2)][(boundaryLength-2)][perimeter];
    for (int i=0; i<(boundaryLength-2); i++) {
        for (int j=0;j<(boundaryLength-2);j++) {
            
            long double sum = 0;
            
            for (int c=0; c<perimeter; c++) {
                // pi-1 = p0, pi = p1, pi+1=p2
                Point p0, p1, p2;
                p1.x = boundary[c][0];
                p1.y = boundary[c][1];
                p0.x = boundary[(c-1 + perimeter)%(perimeter)][0];
                p0.y = boundary[(c-1 + perimeter)%(perimeter)][1];
                p2.x = boundary[(c+1)%(perimeter)][0];
                p2.y = boundary[(c+1)%(perimeter)][1];

                Point inner;
                inner.x= interiorPoints[i][j][0];
                inner.y= interiorPoints[i][j][1];
                    
                vector<double> ip0, ip1, ip2;
                ip0.push_back(inner.x - p0.x);
                ip0.push_back(inner.y - p0.y);
                ip1.push_back(inner.x - p1.x);
                ip1.push_back(inner.y - p1.y);
                ip2.push_back(inner.x - p2.x);
                ip2.push_back(inner.y - p2.y);
                    
                double alpha = acos(dotProduct(ip0, ip1)/(length(ip0)* length(ip1)));
                double beta = acos(dotProduct(ip1, ip2)/(length(ip1)* length(ip2)));
                if (alpha != alpha || tan(alpha/2) != tan(alpha/2))
                    cout << "alpha is nan" << endl;
                if (beta != beta || tan(beta/2) != tan(beta/2))
                    cout << "beta is nan" << endl;
                if (length(ip1) ==0)
                    cout << "length of ip1 is 0" << endl;
                lambda[i][j][(c+1)%perimeter] = (tan(alpha/2) + tan(beta/2)) / (length(ip1));
                sum += lambda[i][j][(c+1)%perimeter];

                }
//            cout << "Sum : " << sum << endl;
            for (int c=0; c<perimeter; c++) {
                lambda[i][j][c] = lambda[i][j][c] / sum;
//                cout << "(" <<interiorPoints[i][j][0] << "," << interiorPoints[i][j][1] << ") ("  << boundary[c][0] << "," << boundary[c][1] << ") " << lambda[i][j][c] << endl;
            }
        }
    }
    

    
    //Compute the differences along the boundary
    double diff_boundary[perimeter][3];
    for (int c = 0; c < perimeter; c++) {
        double x,y;
        x = boundary[c][0];
        y = boundary[c][1];
        Pixel Px_patch, Px_bg;
        Px_patch = patch->getPixel(x, y);
        Px_bg = background->getPixel(x, y);
        for (int channel = 0; channel < 3; channel++) {
            diff_boundary[c][channel] = Px_bg.getColor(channel) - Px_patch.getColor(channel);
        }
    }
    
    //Evaluate the mean-value interpolant at x
    for (int i = 0; i < boundaryLength - 2; i++)
        for (int j = 0; j < boundaryLength - 2; j++)
        {
            int x, y;
            x = interiorPoints[i][j][0];
            y = interiorPoints[i][j][1];
            double rx_red = 0;
            double rx_green = 0;
            double rx_blue = 0;
            
            for (int c = 0 ; c < perimeter; c++) {
                rx_red += lambda[i][j][c] * diff_boundary[c][0];
                rx_green += lambda[i][j][c] * diff_boundary[c][1];
                rx_blue += lambda[i][j][c] * diff_boundary[c][2];
            }
            
            result->setPixel(x, y, RED, correctColor(patch->getPixel(x, y, RED) + rx_red));
            result->setPixel(x, y, GREEN, correctColor(patch->getPixel(x, y, GREEN) + rx_green));
            result->setPixel(x, y, BLUE, correctColor(patch->getPixel(x, y, BLUE) + rx_blue));
        }
    return result;
}



double correctColor(double value)
{
    if (value > 1)
        return 1;
    else if (value < 0)
        return 0;
    else
        return value;
}

