#include "control.h"
#include "ip.h"
#include "main.h"
#include <stdlib.h>


/*
* IMPORTANT - DO NOT CHANGE THIS FILE - IMPORTANT
*/


enum {
	M_QUIT = 0,
	M_HELP = 1,

	M_FILE_OPEN =2,
    M_FILE_OPEN_DEFAULT = 33,
	M_FILE_SAVE =3,
	M_FILE_INFO =4,
	M_FILE_REVERT =5,

    
    M_MVC_SQUARE = 6,

} MENU_ITEMS;


int make_menu ()
{
	int file = glutCreateMenu(menu_func);
	glutAddMenuEntry( "Open...",		M_FILE_OPEN);
    glutAddMenuEntry("Open the background ...", M_FILE_OPEN_DEFAULT);
	glutAddMenuEntry( "Save...",		M_FILE_SAVE);
	glutAddMenuEntry( "Get Image Info",		M_FILE_INFO);
	glutAddMenuEntry( "Revert",		M_FILE_REVERT);


	int process = glutCreateMenu(menu_func);
    glutAddMenuEntry("MVC with square", M_MVC_SQUARE);

	int main = glutCreateMenu(menu_func);
	glutAddSubMenu(   "File",		file);
	glutAddSubMenu(   "Process",		process);
	glutAddMenuEntry( "Help",		M_HELP);
	glutAddMenuEntry( "Quit",		M_QUIT);

	glutAttachMenu(GLUT_RIGHT_BUTTON);

	return main;
}


static inline void checkStream (const istream& in)
{
	if (in.fail())
	{
		cerr << "Fatal error: stream failed!" << endl;
		exit(-1);
	}
}


void menu_func (int value)
{
	// variables used in the switch statement
	char filename[MAX_LINE];

	switch (value)
	{
	case M_QUIT:  // enum #0
		exit(0);
		break;



	case M_HELP:  // enum #1
		menu_help();
		break;



	case M_FILE_OPEN:   // enum #2
		if (!quietMode)
			cerr << "Open file (string - no spaces) : ";
		cin  >> filename;
		checkStream(cin);
		image_load(filename);
		break;
            
    case M_FILE_OPEN_DEFAULT:
        image_load("/Users/owlroro/Desktop/sdact/sand.bmp");
        break;
            

	case M_FILE_SAVE:   // enum #3
		if (!quietMode)
			cerr << "Save as (string - no spaces) : ";
		cin  >> filename;
		checkStream(cin);
		image_save(filename);
		break;


	case M_FILE_INFO:  // enum #4
		image_print_info();
		break;


	case M_FILE_REVERT:  // enum #5
		image_revert();
		break;

		default:
		process_func(value);
	}
	return;
}

void process_func (int value)
{
    Image* resultImage = NULL;
    int x_offset = getInt("x : ");
    int y_offset = getInt("y : ");
    
    if (x_offset > currentImage->getWidth() || y_offset > currentImage->getHeight()) {
        cout << "x or y is too large. Using default values instead."<<endl;
        x_offset = 200;
        y_offset = 120;
    }
    
//    char filename[MAX_NAME];
//    if (!quietMode)
//        cerr << "Enter filename of patch image (string - no spaces) : ";
//    cin  >> filename;
    Image* secondImage = new Image();
//    secondImage->read(filename);
    secondImage->read("/Users/owlroro/Desktop/sdact/polarbear2.bmp");
    resultImage = computeLambda(currentImage, secondImage, x_offset, y_offset);
    delete secondImage;
    

    if (resultImage != NULL)
    {
        delete currentImage;
        currentImage = resultImage;
        
        if (currentImage->getWidth()  != window_width    ||
            currentImage->getHeight() != window_height)
            reshape(currentImage->getWidth(), currentImage->getHeight());
        
        if (!quietMode)
            cerr << "done!" << endl;
        
        if (!textMode)
            glutPostRedisplay();
    }
}

void keyboard_func (unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'H':
	case 'h':
		menu_help();
		break;
		;;

	case 'Q':
	case 'q':
		exit(0);
		break;
		;;
	}
}


void menu_help ()
{
	cerr << endl
		<< "hmc cs155 image processor" << endl
		<< "please see the ip manual for usage and algorithm information" << endl
		<< "http://www.cs.hmc.edu/courses/2002/fall/cs155/proj1/doc/ip_manual.html"
		<< endl << endl;
}


#define MENUOP(num, tag)	cerr << " " << num << ") " << tag << endl;



void textMenuLoop ()
{
	char command[MAX_LINE];


	while (true)
	{
		if (!quietMode)
			cerr << endl
			<< "selection > " << flush;
		cin  >> command;

		switch (command[0])
		{
		case '\n':
		case '\0':
			//printMenu();
			break;

		case 'Q':
		case 'q':
			menu_func(M_QUIT);
			break;

		case 'H':
		case 'h':
			menu_func(M_HELP);
			break;

		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			menu_func(atoi(command));
			break;

		default:
			//printMenu();
			break;
		}
	}
}


void image_load (const char* filename)
{
	if (currentImage)
		delete currentImage;
	if (originalImage)
		delete originalImage;
	currentImage  = NULL;
	originalImage = NULL;

	originalImage = new Image();
	originalImage->read(filename);

	if (originalImage->good())
	{  
		currentImage = new Image(*originalImage);
		reshape(currentImage->getWidth(), currentImage->getHeight());
	}
	else
	{
		delete originalImage;  
		originalImage = NULL;
		cerr << "Couldn't load image " << filename << "!" << endl;
		return;
	}

	if (!textMode)
		glutPostRedisplay();

	if (!quietMode)
		cerr << "done!" << endl;
}  


void image_save (const char* filename)
{
	if (currentImage)
	{
		if (currentImage->write(filename) == 0)
		{
			//delete originalImage;
			//originalImage = new Image(*currentImage);
		}
	}  
	else if (originalImage)
	{
		originalImage->write(filename);
	}
	else
	{
		cerr << "No image!" << endl;
		return;
	}

	if (!quietMode)
		cerr << "done!" << endl;
}


void image_print_info ()
{  
	if (currentImage) {
		cerr << "width:    " << currentImage->getWidth() << endl
			<< "height:   " << currentImage->getHeight() << endl
			<< "bits:     " << currentImage->getBits() << endl;
	}
	cerr << "done!" << endl;
}


void image_revert ()
{
	if (currentImage)
		delete currentImage;

	if (originalImage)
	{
		currentImage = new Image(*originalImage);

		if (window_width  != currentImage->getWidth() ||
			window_height != currentImage->getHeight())
			reshape(currentImage->getWidth(), currentImage->getHeight());
	}
	else
	{
		cerr << "No image!" << endl;
		return;
	}

	if (!textMode)
		glutPostRedisplay();

	if (!quietMode)
		cerr << "done!" << endl;
}  

int getFilterSize()
{
	int filtersize;
	if (!quietMode)
		cerr << "Enter filter size (positive, odd integer) : ";
	cin  >> filtersize;
	if (filtersize % 2 !=1 || filtersize<=0)
	{
		cerr << "Sorry, the filter size must be a positive, odd integer." << endl;
		filtersize=0;
	}
	checkStream(cin);
	return filtersize;
}
double getDouble(const char* message)
{
	double value;
	if (!quietMode)
		cerr << "Enter " << message << "(double): ";
	cin  >> value;
	checkStream(cin);
	return value;
}

double getPositiveDouble(const char* message)
{
	double value;
	if (!quietMode)
		cerr << "Enter positive " << message << "(double): ";
	cin  >> value;
	checkStream(cin);
	return value;
}

int getInt(const char* message)
{

	int value;
	if (!quietMode)
		cerr << "Enter " << message << " (integer): ";
	cin  >> value;
	checkStream(cin);
	return value;
}
