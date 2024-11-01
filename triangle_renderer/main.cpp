///////////////////////////////////////////////////
//
//	Hamish Carr
//	January, 2018
//
//	------------------------
//	main.cpp
//	------------------------
//	
///////////////////////////////////////////////////

#include <QtWidgets/QApplication>
#include "GeometricWidget.h"
#include <stdio.h>

int main(int argc, char **argv)
	{ // main()
	// initialize QT
	QApplication app(argc, argv);

	// the geometric surface
	GeometricSurfaceFaceDS surface;

	// check the args to make sure there's an input file
	if (argc == 2)
		{ // two parameters - read a file
		if (!surface.ReadFileTriangleSoup(argv[1]))
			{ // surface read failed 
			printf("Read failed for file %s\n", argv[1]);
			exit(0);
			} // surface read failed
		else
			{ // surface read succeeded
			//	create a window
			GeometricWidget aWindow(&surface, NULL);
	
			// 	set the initial size
			aWindow.resize(600, 600);

			// show the window
			aWindow.show();

			// set QT running
			return app.exec();
			} // surface read succeeded			
		} // two parameters - read a file
	else 
		{ // too many parameters 
		printf("Usage: %s filename\n", argv[0]); 
		exit (0);
		} // too many parameters 

	// paranoid return value
	exit(0);
	} // main()
