#pragma once
#ifndef FILESIO_H
#define FILESIO_H

#include<iostream>
#include<stdio.h>
#include<io.h>
#include<string.h>
#include"ParaDefine.h"
#include"preprocess.h"
#include<opencv2\highgui\highgui.hpp>
using namespace cv;

// load json files // need be improved
void loadjson(camera_t *cam_depth, camera_t *cam_rgb, D2G *Rots, D2G *Trans, FilePath json_path) {
	FILE *fp = fopen(json_path, "r");
	if (fp == NULL) {
		printf("Error\n");
		return;
	}
	char string[150], symble1[12], se[12], symble2[12];
	float var = 0.0;
	cam_depth->distortions[6] = var; cam_depth->distortions[7] = var;
	cam_rgb->distortions[6] = var; cam_rgb->distortions[7] = var;
	//while (1)
	//{
	fgets(string, 150, fp);
	fgets(string, 150, fp);
	fgets(string, 150, fp);
	fgets(string, 150, fp);

	int res = fscanf(fp, "%s", string);
	if (res == EOF) {
		//break;
	}
	if (strcmp(string, "\"DEPTH\"") == 0) {
		// printf("start DEPTH\n");
		fgets(string, 150, fp);
		//fscanf(fp, "%s %s", symble1, se);
		fscanf(fp, "%s", string);
		//printf("%s\n", string);
		// printf("%s\n", string);
		while (strcmp(string, "},") != 0) {
			//fscanf(fp, "%s %f", symble1, &var);
			//fgets(symble2, 12, fp);
			if (strcmp(string, "\"Cx\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->cx = var;
			}
			if (strcmp(string, "\"Cy\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->cy = var;
			}
			if (strcmp(string, "\"Fx\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->fx = var;
			}
			if (strcmp(string, "\"Fx\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->fy = var;
			}
			if (strcmp(string, "\"ImgHeight\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->height = var;
			}
			if (strcmp(string, "\"ImgWidth\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->width = var;
			}
			if (strcmp(string, "\"K1\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[0] = var;
			}
			if (strcmp(string, "\"K2\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[1] = var;
			}
			if (strcmp(string, "\"K3\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[4] = var;
			}
			if (strcmp(string, "\"K4\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[5] = var;
			}
			if (strcmp(string, "\"K5\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[6] = var;
			}
			if (strcmp(string, "\"K6\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[7] = var;
			}
			if (strcmp(string, "\"P1\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[2] = var;
			}
			if (strcmp(string, "\"P2\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_depth->distortions[3] = var;
			}
			fgets(string, 150, fp);
			fscanf(fp, "%s", string);
		}
		fscanf(fp, "%s", string);
	}// end if 	
	if (strcmp(string, "\"Depth2RGB_Rotate\"") == 0) {
		// printf("start Depth2RGB_Rotate\n");
		fscanf(fp, "%s %s", symble1, se);
		fscanf(fp, "%s", string);
		while (strcmp(string, "},") != 0) {
			fscanf(fp, "%s %f", symble1, &var);
			fgets(symble2, 12, fp);
			if (strcmp(string, "\"RX\"") == 0)
				Rots->rx = var;
			if (strcmp(string, "\"RY\"") == 0)
				Rots->ry = var;
			if (strcmp(string, "\"RZ\"") == 0)
				Rots->rz = var;
			fscanf(fp, "%s", string);
		}
		fscanf(fp, "%s", string);
	}//end if 

	if (strcmp(string, "\"Depth2RGB_Trans\"") == 0) {
		// printf("start Depth2RGB_Trans\n");
		fscanf(fp, "%s %s", symble1, se);
		fscanf(fp, "%s", string);
		while (strcmp(string, "},") != 0) {
			fscanf(fp, "%s %f", symble1, &var);
			fgets(symble2, 12, fp);
			if (strcmp(string, "\"TX\"") == 0)
				Trans->rx = var;
			if (strcmp(string, "\"TY\"") == 0)
				Trans->ry = var;
			if (strcmp(string, "\"TZ\"") == 0)
				Trans->rz = var;
			fscanf(fp, "%s", string);
		}
		fscanf(fp, "%s", string);
	}//end if
	// printf("%s\n", string);
	fgets(string, 150, fp);
	fgets(string, 150, fp);
	// printf("%s\n", string);
	fscanf(fp, "%s", string);

	if (strcmp(string, "\"RGB\"") == 0) {
		// printf("start RGB\n");
		fgets(string, 150, fp);
		//fscanf(fp, "%s %s", symble1, se);
		fscanf(fp, "%s", string);
		while (strcmp(string, "},") != 0) {
			// fscanf(fp, "%s %f", symble1, &var);
			// fgets(symble2, 12, fp);
			if (strcmp(string, "\"Cx\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->cx = var;
			}
			if (strcmp(string, "\"Cy\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->cy = var;
			}
			if (strcmp(string, "\"Fx\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->fx = var;
			}
			if (strcmp(string, "\"Fx\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->fy = var;
			}
			if (strcmp(string, "\"ImgHeight\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->height = var;
			}
			if (strcmp(string, "\"ImgWidth\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->width = var;
			}
			if (strcmp(string, "\"K1\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[0] = var;
			}
			if (strcmp(string, "\"K2\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[1] = var;
			}
			if (strcmp(string, "\"K3\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[4] = var;
			}
			if (strcmp(string, "\"K4\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[5] = var;
			}
			if (strcmp(string, "\"K5\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[6] = var;
			}
			if (strcmp(string, "\"K6\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[7] = var;
			}
			if (strcmp(string, "\"P1\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[2] = var;
			}
			if (strcmp(string, "\"P2\"") == 0) {
				fscanf(fp, "%s %f", symble1, &var);
				cam_rgb->distortions[3] = var;
			}
			fgets(string, 150, fp);
			fscanf(fp, "%s", string);
		}
		fscanf(fp, "%s", string);
	}// end if 
	// printf("END");
	//}
	fclose(fp);
	return;
}

// get all files
int loadfiles(FilePath filepath, FilenameList namergb_list, FilenameList namedep_list) {
	// init all name list in this function
	
	int num = 0;//the num of all images
	return num;
}

void load_rgbimg(FilePath file_path, RGBImage* rgb) {
	Mat img = imread(file_path, CV_LOAD_IMAGE_UNCHANGED);
	int height = img.rows;
	int width = img.cols;
	*rgb = alloc_RGB(height, width);
	for (int i = 0; i < height; i++)
	{
		unsigned char *data = img.data + i*width*img.channels();
		for (int j = 0; j < width; j++)
		{
			rgb->data[i][j * 3] = *(data + j*img.channels());
			rgb->data[i][j * 3 + 1] = *(data + j*img.channels() + 1);
			rgb->data[i][j * 3 + 2] = *(data + j*img.channels() + 2);
			// cout << rgb->data[i][j * 3] << " " << rgb->data[i][j * 3 + 1] << " " << rgb->data[i][j * 3 + 2] << endl;
		}
	}
	return;
}

void load_depimg(FilePath file_path, DepthImage* depth) { //
	
	// load_png_image(); 
	CvMat *img = cvLoadImageM(file_path, CV_LOAD_IMAGE_ANYDEPTH);
	int height = img->rows; int width = img->cols;
	CvMat *img_float = cvCreateMat(height, width, CV_32FC1);
	cvConvert(img, img_float);
	*depth = alloc_depth(height, width);
	float*p = (float*)(img_float->data.fl);
	// FILE *fp2 = fopen("out_depth.txt", "w");
	for (int y = 0; y < height; y++) {
		int d_start = y * width;
		for (int x = 0; x < width; x++) {
			depth->data[y][x] = p[d_start + x] / 10.0;
			// fprintf(fp2, "%7.2f ", depth->data[y][x]);
		}
		// fprintf(fp2, "\n");
	}
	// fclose(fp2);
	cvReleaseMat(&img);
	return;
}

char *myStrncpy(char *dest, const char *src, size_t n) {
	int size = sizeof(char)*(n + 1);
	char *tmp = (char*)malloc(size);  // 开辟大小为n+1的临时内存tmp
	if (tmp) {
		memset(tmp, '\0', size);  // 将内存初始化为0
		memcpy(tmp, src, size - 1);  // 将src的前n个字节拷贝到tmp
		memcpy(dest, tmp, size);  // 将临时空间tmp的内容拷贝到dest
		free(tmp);  // 释放内存
		return dest;
	}
	else {
		return NULL;
	}
}

void save_one( DepthImage new_depth, FilePath out_path, char* name) {
	char new_name[81];
	int length = strlen(name);
	myStrncpy(new_name, name, length - 4);
	//strcat(new_name, name);
	strcat(new_name, ".raw");
	FILE *fp = fopen(new_name, "w");
	FILE *fp1 = fopen("new_depth.txt", "w");
	float *a = (float*)malloc(sizeof(float)*new_depth.cols);
	for (int i = 0; i < new_depth.rows; i++){
		for (int j = 0; j < new_depth.cols; j++){
			a[j] = new_depth.data[i][j];
			fprintf(fp1, "%7.2f ", a[j]);
		}
		fwrite(a, sizeof(float)*new_depth.cols, 1, fp);
		fprintf(fp1, "\n");
	}
	fclose(fp);
	fclose(fp1);

	/*cv::Mat img = cv::Mat(new_depth.rows, new_depth.cols, CV_32FC1, new_depth.data);
	cv::Mat img2;
	img.convertTo(img2, CV_8UC1);
	cv::namedWindow("MyWindow", CV_WINDOW_AUTOSIZE);
	cv::imshow("MyWindow", img2);
	cv::waitKey(0);
	cv::destroyWindow("MyWindow");*/

	return;
}

void save_cld(Matrix Cld_3D, Matrix rgb_copy, FilePath out_path, char* name) {
	
	FILE *fp1 = fopen("xyzrgb_3D.obj", "w");
	int end = Cld_3D.cols - 1;
	for (int i = 0; i < Cld_3D.cols; i++) {
		fprintf(fp1, "v ");
		fprintf(fp1, "%f %f %f ", Cld_3D.data[0][i], Cld_3D.data[1][i], Cld_3D.data[2][i]);

		int b = rgb_copy.data[0][i];
		int g = rgb_copy.data[1][i];
		int r = rgb_copy.data[2][i];
		fprintf(fp1, "%d %d %d", r, g, b);
		fprintf(fp1, "\n");
	}
	fclose(fp1);
	return;
}

void save_rgbpixel(Matrix rgb_copy, int rows, int cols, FilePath out_path, char*name) {
	Matrix save_matrix1 = alloc_matrix(rows, cols);
	Matrix save_matrix2 = alloc_matrix(rows, cols);
	Matrix save_matrix3 = alloc_matrix(rows, cols);
	for (int i = 0; i < cols; i++){
		for (int j = 0; j < rows; j++){
			save_matrix1.data[j][i] = rgb_copy.data[0][i*rows + j];//bgr
			save_matrix2.data[j][i] = rgb_copy.data[1][i*rows + j];
			save_matrix3.data[j][i] = rgb_copy.data[2][i*rows + j];
		}

	}

	Mat mat(rows, cols, CV_8UC3);
	//mat.create(CV_8UC3)
	for (int i = 0; i < mat.rows; i++)
	{
		unsigned char *data = mat.data + i * 640 * mat.channels();
		for (int j = 0; j < mat.cols; j++)
		{
			*(data + j*mat.channels()) = save_matrix1.data[i][j];//r
			*(data + j*mat.channels() + 1) = save_matrix2.data[i][j];//g
			*(data + j*mat.channels() + 2) = save_matrix3.data[i][j];//b
		}
	}
	imwrite("1.jpg", mat);
	free_matrix(save_matrix1); free_matrix(save_matrix2); free_matrix(save_matrix3);
	return;
}
#endif // !FILESIO_H

