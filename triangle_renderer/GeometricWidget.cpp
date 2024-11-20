///////////////////////////////////////////////////
//
//	Hamish Carr
//	January, 2018
//
//	------------------------
//	GeometricWidget.h
//	------------------------
//	
//	The main widget that shows the geometry
//	
///////////////////////////////////////////////////

#include <math.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "GeometricWidget.h"
static GLfloat light_position[] = {0.0, 0.0, 1.0, 0.0};							

// constructor
GeometricWidget::GeometricWidget(GeometricSurfaceFaceDS *newSurface, QWidget *parent)
	: _GEOMETRIC_WIDGET_PARENT_CLASS(parent)
	{ // constructor
	// store pointer to the model
	surface = newSurface;
	
	// initialise arcballs to 80% of the widget's size
	Ball_Init(&lightBall);		Ball_Place(&lightBall, qOne, 0.80);
	Ball_Init(&objectBall);		Ball_Place(&objectBall, qOne, 0.80);

	// initialise translation values
	translate_x = translate_y = 0.0;
	
	// and set the button to an arbitrary value
	whichButton = -1;
	
	} // constructor

// destructor
GeometricWidget::~GeometricWidget()
	{ // destructor
	// nothing yet
	} // destructor																	

// called when OpenGL context is set up
void GeometricWidget::initializeGL()
	{ // GeometricWidget::initializeGL()
	// enable Z-buffering
	glEnable(GL_DEPTH_TEST);
	
	// set lighting parameters
	glShadeModel(GL_FLAT);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

        glEnable(GL_CULL_FACE);
	
	// background is pink
	glClearColor(1.0, 0.7, 0.7, 1.0);
	} // GeometricWidget::initializeGL()

// called every time the widget is resized
void GeometricWidget::resizeGL(int w, int h)
	{ // GeometricWidget::resizeGL()
	// reset the viewport
	glViewport(0, 0, w, h);
	
	// set projection matrix to be glOrtho based on zoom & window size
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	// retrieve the scale factor
	float size = surface->boundingSphereSize;
	
	// compute the aspect ratio of the widget
	float aspectRatio = (float) w / (float) h;
	
	// depending on aspect ratio, set to accomodate a sphere of radius = diagonal without clipping
    float zoom = 2.0;
	if (aspectRatio <= 1.0)
		glOrtho(-aspectRatio * size / zoom, aspectRatio * size / zoom, -size / zoom, size / zoom, -size / zoom, size / zoom);
	else
		glOrtho(-size, size, -size/aspectRatio, size/aspectRatio, -size, size);

	} // GeometricWidget::resizeGL()
	
// called every time the widget needs painting
void GeometricWidget::paintGL()
	{ // GeometricWidget::paintGL()
	// clear the buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set lighting on
	glEnable(GL_LIGHTING);

	// set model view matrix based on stored translation, rotation &c.
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// set light position first
	// retrieve rotation from arcball & apply
	GLfloat mNow[16];
	Ball_Value(&lightBall, mNow);
	glMultMatrixf(mNow);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
			 
	// apply translation for interface control
	glLoadIdentity();
	glTranslatef(translate_x, translate_y, 0.0);

	// apply rotation matrix from arcball
	Ball_Value(&objectBall, mNow);
	glMultMatrixf(mNow);

	// now draw the surface
	surface->Render();
	} // GeometricWidget::paintGL()

// mouse-handling
void GeometricWidget::mousePressEvent(QMouseEvent *event)
	{ // GeometricWidget::mousePressEvent()
	// store the button for future reference
	whichButton = event->button();

	// find the minimum of height & width	
	float size = (width() > height()) ? height() : width();

	// convert to the ArcBall's vector type
	HVect vNow;

	// scale both coordinates from that
	vNow.x = (2.0 * event->x() - size) / size;
	vNow.y = (size - 2.0 * event->y() ) / size;

	// now either translate or rotate object or light
	switch(whichButton)
		{ // button switch
		case Qt::RightButton:
			// save the last x, y
			last_x = vNow.x; last_y = vNow.y;
			// and update
			_GL_WIDGET_UPDATE_CALL();
			break;
		case Qt::MiddleButton:
			// pass the point to the arcball code	
			Ball_Mouse(&lightBall, vNow);
			// start dragging
			Ball_BeginDrag(&lightBall);
			// update the widget
			_GL_WIDGET_UPDATE_CALL();
			break;
		case Qt::LeftButton:
			// pass the point to the arcball code	
			Ball_Mouse(&objectBall, vNow);
			// start dragging
			Ball_BeginDrag(&objectBall);
			// update the widget
			_GL_WIDGET_UPDATE_CALL();
			break;
		} // button switch
	} // GeometricWidget::mousePressEvent()
	
void GeometricWidget::mouseMoveEvent(QMouseEvent *event)
	{ // GeometricWidget::mouseMoveEvent()
	// find the minimum of height & width	
	float size = (width() > height()) ? height() : width();

	// convert to the ArcBall's vector type
	HVect vNow;

	// scale both coordinates from that
	vNow.x = (2.0 * event->x() - size) / size;
	vNow.y = (size - 2.0 * event->y() ) / size;

	// now either translate or rotate object or light
	switch(whichButton)
		{ // button switch
		case Qt::RightButton:
			// subtract the translation
			translate_x += vNow.x - last_x;
			translate_y += vNow.y - last_y;
			last_x = vNow.x; 
			last_y = vNow.y;
			// update the widget
			_GL_WIDGET_UPDATE_CALL();
			break;
		case Qt::MiddleButton:
			// pass it to the arcball code	
			Ball_Mouse(&lightBall, vNow);
			// start dragging
			Ball_Update(&lightBall);
			// update the widget
			_GL_WIDGET_UPDATE_CALL();
			break;
		case Qt::LeftButton:
			// pass it to the arcball code	
			Ball_Mouse(&objectBall, vNow);
			// start dragging
			Ball_Update(&objectBall);
			// update the widget
			_GL_WIDGET_UPDATE_CALL();
			break;
		} // button switch
	} // GeometricWidget::mouseMoveEvent()
	
void GeometricWidget::mouseReleaseEvent(QMouseEvent *event)
	{ // GeometricWidget::mouseReleaseEvent()
	// now either translate or rotate object or light
	switch(whichButton)
		{ // button switch
		case Qt::RightButton:
			// just update
			_GL_WIDGET_UPDATE_CALL();
			break;
		case Qt::MiddleButton:
			// end the drag
			Ball_EndDrag(&lightBall);
			// update the widget
			_GL_WIDGET_UPDATE_CALL();
			break;
		case Qt::LeftButton:
			// end the drag
			Ball_EndDrag(&objectBall);
			// update the widget
			_GL_WIDGET_UPDATE_CALL();
			break;
		} // button switch
	} // GeometricWidget::mouseReleaseEvent()
	
