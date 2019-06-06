#ifndef PARADEFINE_H
#define PARADEFINE_H
#include<opencv2\calib3d\calib3d_c.h>
#include<opencv\cv.h>
#include<opencv2\core\core_c.h>
#include<opencv\highgui.h>
#include<opencv2\highgui\highgui.hpp>
#include<opencv2\highgui\highgui_c.h>
typedef struct camera {
	int width, height;
	int distortion_dim;
	// intrinsic
	float fx, fy, cx, cy;
	float distortions[8];
	// extrinsic
	float r[3], t[3];
} camera_t;

typedef struct RotaTran {
	float rx, ry, rz;
} D2G;
/*template <typename T> struct Matrix {
	int rows, cols;
	int channels;
	T *data;
};*/
typedef struct {
	int rows;
	int cols;
	float ** data;
} Matrix;

typedef struct{
	int rows;
	int cols;
	// int channnel;
	float **data;
} DepthImage;

typedef struct {
	int rows;
	int cols;
	int channel = 3;
	int ** data;
}RGBImage;

typedef struct {
	int rows;
	int cols;
	float ** data;
} CloudImage;

typedef struct {
	int length;
	char** data;
}FilenameList;

typedef char* FilePath;

#endif 
