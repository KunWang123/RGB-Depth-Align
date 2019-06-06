#include<iostream>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
// #include<iostream>
 #include<opencv\cv.h>
#include<opencv\highgui.h>
#include<opencv2\calib3d\calib3d_c.h>
#include<opencv2\highgui\highgui.hpp>
#include "preprocess.h"
#include"ParaDefine.h"

//typedef char* FilePath;
using namespace std;
using namespace cv;

int main() {
	FilePath json_path = "../rgbd_preprocessing_data/sedna_bs20_date190418_intrinsic.json";
	FilePath img_path = "../rgbd_preprocessing_data/";	
	FilePath Out_folder = "./out";
/* // sample of save name of image
	FilenameList name_list;
	name_list.length = 3;
	name_list.data = (char**)malloc(sizeof(char*) * name_list.length);
	for (int i = 0; i < name_list.length; i++)
	{
		name_list.data[i] = Out_folder;
	}

	printf("%s\n", name_list.data[0]);
	printf("%s\n", name_list.data[1]);
	printf("%s\n", name_list.data[2]);*/

	
	//printf("%s\n", json_path);
	process_all(json_path, img_path, Out_folder);
	system("Pause");
	return 0;
}
