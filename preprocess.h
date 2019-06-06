#pragma once
#ifndef PREPROCESS_H
#define PREPROCESS_H

#include<stdio.h>
#include<stdbool.h>
#include<io.h>

#include"ParaDefine.h"
#include "MatrixOperation.h"
#include "Filesio.h"


void compute_cs_matrix(D2G Rots, float*sin_Rots, float *cos_Rots) {
	sin_Rots[0] = sin(Rots.rx); cos_Rots[0] = cos(Rots.rx);
	sin_Rots[1] = sin(Rots.ry); cos_Rots[1] = cos(Rots.ry);
	sin_Rots[2] = sin(Rots.rz); cos_Rots[2] = cos(Rots.rz);
	return;
}

void compute_Rxyz_matrix(Matrix Rx, Matrix Ry, Matrix Rz, float *sin_Rots, float *cos_Rots) {
	Rx.data[0][0] = 1; Rx.data[0][1] = 0; Rx.data[0][2] = 0;
	Rx.data[1][0] = 0; Rx.data[1][1] = cos_Rots[0]; Rx.data[1][2] = -sin_Rots[0];
	Rx.data[2][0] = 0; Rx.data[2][1] = sin_Rots[0]; Rx.data[2][2] = cos_Rots[0];
	// R{y} = [c(i), 0, s(i); 0, 1, 0;-s(i), 0, c(i)];
	Ry.data[0][0] = cos_Rots[1]; Ry.data[0][1] = 0; Ry.data[0][2] = sin_Rots[1];
	Ry.data[1][0] = 0; Ry.data[1][1] = 1; Ry.data[1][2] = 0;
	Ry.data[2][0] = -sin_Rots[1]; Ry.data[2][1] = 0; Ry.data[2][2] = cos_Rots[1];
	//R{z} = [c(i),-s(i), 0; s(i), c(i), 0; 0, 0, 1]
	Rz.data[0][0] = cos_Rots[2]; Rz.data[0][1] = -sin_Rots[2]; Rz.data[0][2] = 0;
	Rz.data[1][0] = sin_Rots[2]; Rz.data[1][1] = cos_Rots[2]; Rz.data[1][2] = 0;
	Rz.data[2][0] = 0; Rz.data[2][1] = 0; Rz.data[2][2] = 1;
	return;
}

void parametersintoMatrix(camera_t cam_depth, camera_t cam_rgb, Matrix *f_depth,
	Matrix *f_rgb, D2G Rots, D2G Trans, Matrix *D2RGB, char* type) {
	*f_depth = alloc_matrix(3, 3);
	*f_rgb = alloc_matrix(3, 3);
	*D2RGB = alloc_matrix(3, 4);

	f_depth->data[0][0] = cam_depth.fx;
	f_depth->data[0][2] = cam_depth.cx;
	f_depth->data[1][1] = cam_depth.fy;
	f_depth->data[1][2] = cam_depth.cy;
	f_depth->data[2][2] = 1;

	f_rgb->data[0][0] = cam_rgb.fx;
	f_rgb->data[0][2] = cam_rgb.cx;
	f_rgb->data[1][1] = cam_rgb.fy;
	f_rgb->data[1][2] = cam_rgb.cy;
	f_rgb->data[2][2] = 1;
	
	D2RGB->data[0][3] = Trans.rx;// Tx
	D2RGB->data[1][3] = Trans.ry;
	D2RGB->data[2][3] = Trans.rz;
	/*CvMat src,  dst;
	float src_r[] = { Rots.rx, Rots.ry, Rots.rz };
	float dsr_R[9] = {0,0,0,0,0,0,0,0,0};
	cvInitMatHeader(&src, 1, 3, CV_32FC1, src_r, CV_AUTOSTEP);
	cvInitMatHeader(&dst, 3, 3, CV_32FC1, dsr_R, CV_AUTOSTEP);
	cvRodrigues2(&src, &dst, 0);*/
	/*Matrix dst = alloc_matrix(3, 3);
	Matrix dst2 = alloc_matrix(3, 3);
	Matrix Rxyz_var = alloc_matrix(3, 3);
	set_identity_matrix(dst);
	for (int iter = 0; iter < 3; iter++) {
		char axis = type[iter];
		if (strcmp(&axis,"x")==0){

		}
		else if (strcmp(&axis, "y") == 0){

		} 
		else if(strcmp(&axis, "z") == 0)
		{

		}
	}*/
	float sin_Rots[3] = { 0 };
	float cos_Rots[3] = { 0 };
	Matrix Rx = alloc_matrix(3, 3);
	Matrix Ry = alloc_matrix(3, 3);
	Matrix Rz = alloc_matrix(3, 3);
	Matrix Rotation_Matrix1 = alloc_matrix(3, 3);
	Matrix Rotation_Matrix = alloc_matrix(3, 3);
	// compute the sin_Rots and cos_Rots
	compute_cs_matrix(Rots, sin_Rots, cos_Rots);
	// computer Rx Ry Rz
	compute_Rxyz_matrix(Rx, Ry, Rz, sin_Rots, cos_Rots);
	// can enmurate ("zxy","xyz","xzy"...)
	if (strcmp(type, "zxy")== 0) {
		multiply_matrix(Rz, Rx, Rotation_Matrix1);
		multiply_matrix(Rotation_Matrix1, Ry, Rotation_Matrix);
	}
	// copy_matrix(Rz, Rotation_Matrix);
	
	for (int D2G_i = 0; D2G_i < 3; D2G_i++)
		for (int D2G_j = 0; D2G_j < 3; D2G_j++)	
			D2RGB->data[D2G_i][D2G_j] = Rotation_Matrix.data[D2G_i][D2G_j];
	
	/*D2RGB->data[0][0] = cos_Rots[1] * cos_Rots[2] - sin_Rots[0] * sin_Rots[1] * sin_Rots[2];
	D2RGB->data[0][1] = -cos_Rots[0] * sin_Rots[2];
	D2RGB->data[0][2] = cos_Rots[2]*sin_Rots[1] + cos_Rots[1]* sin_Rots[0] * sin_Rots[2];

	D2RGB->data[1][0] = cos_Rots[2] * sin_Rots[0] * sin_Rots[1] + cos_Rots[1] * sin_Rots[2];
	D2RGB->data[1][1] = cos_Rots[0] * cos_Rots[2];
	D2RGB->data[1][2] = -cos_Rots[1] * cos_Rots[2] * sin_Rots[0] + sin_Rots[1] * sin_Rots[2];

	D2RGB->data[2][0] = -cos_Rots[0] * sin_Rots[1];
	D2RGB->data[2][1] = sin_Rots[0];
	D2RGB->data[2][2] = cos_Rots[0] * cos_Rots[1];*/

	//D2RGB; //define
	// printf("memory free\n");
	free_matrix(Rx); free_matrix(Ry); free_matrix(Rz);
	free_matrix(Rotation_Matrix1); free_matrix(Rotation_Matrix);
	return;
}

/* phase plane to camera plane */
// [uv1] location for depth 
void UV_matrix(int height, int wideth, Matrix p1) {
	for (int i = 0; i < wideth; i++){
		for (int j = 0; j < height; j++) {
			p1.data[0][height*i + j] = i;//u 480*1+480*2...+480*640
			p1.data[1][height*i + j] = j;//1..480+1..480+...+1..480
			p1.data[2][height*i + j] = 1;
		}
	}
	return;
}
// f Matrix inversion f_inv
void invert_matrix(Matrix in, Matrix *out) {
	*out = alloc_matrix(in.cols, in.rows);
	destructive_invert_matrix(in, *out);
	return;
}

// init variables
void init_allvar( Matrix *p2, Matrix* p3, Matrix* p4, int height, int weidth) {
	
	*p2 = alloc_matrix(3, height*weidth);
	*p3 = alloc_matrix(3, height*weidth);
	*p4 = alloc_matrix(4, height*weidth);
	// *p3_3 = alloc_matrix(3, height*weidth);
	return;
}

// f intenstic var p2 = f_invk * p1
void mult_f_uv1(Matrix UV_depth, Matrix f_depth,
				Matrix p2, int rows, int cols) {
	
	Matrix f_depinv, f_depthcopy;
	f_depthcopy = alloc_matrix(3, 3);
	copy_matrix(f_depth, f_depthcopy);
	
	UV_matrix(rows, cols, UV_depth);
	invert_matrix(f_depthcopy, &f_depinv);
	
	multiply_matrix(f_depinv, UV_depth, p2);

	free_matrix(f_depinv); free_matrix(f_depthcopy);
	return;
}
// 2D Point [xy1] to [XYZ]
void Point2Dto3D(DepthImage depth, Matrix p2) {
	Matrix depth_1D;
	depth_1D = alloc_matrix(1, depth.cols*depth.rows);
	resize_Matrix(depth, depth.rows, depth.cols, depth_1D);
	scale_mul_matrix(p2, depth_1D, true);
	free_matrix(depth_1D);
	return;
}

bool check_in_range(int u, int v, int height, int width) {
	if (u < width && u >= 0 && v < height && v >= 0)
		return true;
	return false;	
}

void new_depthmap(Matrix UV_depth, Matrix uv_newdep, DepthImage depth, DepthImage *new_depth) {

	*new_depth = alloc_depth(depth.rows, depth.cols);// need to be writen 

	for (int i = 0; i < uv_newdep.cols; i++)
	{	
		int u = uv_newdep.data[0][i];
		int v = uv_newdep.data[1][i];
		if (check_in_range(u, v, depth.rows, depth.cols)) {
			int row = UV_depth.data[1][i];
			int col = UV_depth.data[0][i];
			new_depth->data[v][u] = depth.data[row][col];
		}
	}
	return;
}

void coeff_of_P(Matrix p_x, Matrix p_y, Matrix* C, Matrix* r2, Matrix* r4, 
			Matrix* r6, Matrix* x2, Matrix* y2, Matrix* xy, float * dist) {
	*x2 = dot_multipy(p_x, p_x);
	*y2 = dot_multipy(p_y, p_y);
	*xy = dot_multipy(p_x, p_y);
	*r2 = add_matrix(*x2, *y2);
	*r4 = dot_multipy(*r2, *r2);
	*r6 = dot_multipy(*r2, *r4);
	
	*C = alloc_matrix(1, p_x.cols);
	for (int i = 0; i < p_x.cols; i++) {
		C->data[0][i] = (1 + dist[0] * r2->data[0][i] + dist[1] * r4->data[0][i] + dist[4] * r6->data[0][i]) /
			(1 + dist[5] * r2->data[0][i] + dist[6] * r4->data[0][i] + dist[7] * r6->data[0][i]);
	}
	return;
}

void isdistortion(Matrix *p_x, Matrix *p_y, Matrix src, float* dist) {
	Matrix C, r2, r4, r6, x2, y2, xy;
	// init vals
	*p_x = alloc_matrix(1, src.cols);
	*p_y = alloc_matrix(1, src.cols);
	for (int j = 0; j < src.cols; j++) {
		p_x->data[0][j] = src.data[0][j];
		p_y->data[0][j] = src.data[1][j];
	}

	coeff_of_P(*p_x, *p_y, &C, &r2, &r4, &r6, &x2, &y2, &xy, dist);
	for (int i = 0; i < src.cols; i++) {
		p_x->data[0][i] = p_x->data[0][i] * C.data[0][i] + 2 * dist[2] * xy.data[0][i] + 
						dist[3] * (r2.data[0][i] + 2 * x2.data[0][i]);
		p_y->data[0][i] = p_y->data[0][i] * C.data[0][i] + 2 * dist[3] * xy.data[0][i] +
			dist[2] * (r2.data[0][i] + 2 * y2.data[0][i]);
	}
	return;
}
void invdistortion(Matrix *p_x, Matrix *p_y, Matrix src, float* dist) {
	Matrix C, r2, r4, r6, x2, y2, xy;
	// init x y
	*p_x = alloc_matrix(1, src.cols);
	*p_y = alloc_matrix(1, src.cols);
	for (int j = 0; j < src.cols; j++) {
		p_x->data[0][j] = src.data[0][j];
		p_y->data[0][j] = src.data[1][j];
	}
	// x'' y'' iter to x' y'
	for (int iter = 0; iter < 5; iter++) {
		coeff_of_P(*p_x, *p_y, &C, &r2, &r4, &r6, &x2, &y2, &xy, dist);
		for (int i = 0; i < src.cols; i++){
			p_x->data[0][i] = p_x->data[0][i] - 2 * dist[2] * xy.data[0][i] -
				dist[3] * (r2.data[0][i] + 2 * x2.data[0][i]);
			p_x->data[0][i] = p_x->data[0][i] / C.data[0][i];
			p_y->data[0][i] = p_y->data[0][i] - 2 * dist[3] * xy.data[0][i] -
				dist[2] * (r2.data[0][i] + 2 * y2.data[0][i]);
			p_y->data[0][i] = p_y->data[0][i] / C.data[0][i];
		}
	}

	return;
}

void distortion(Matrix src, float dist[], Matrix dst, bool flag) {
	// dist = [k1 k2 p1 p2 k3 k4 k5 k6]
	Matrix p_x, p_y;
	// distortion
	if (flag)
		isdistortion(&p_x, &p_y, src, dist);
	else//invdistortion
		invdistortion(&p_x, &p_y, src, dist);	

	// after distortion 
	for (int i = 0; i < dst.cols; i++) {
		dst.data[0][i] = p_x.data[0][i];
		dst.data[1][i] = p_y.data[0][i];
		dst.data[2][i] = 1.0;
	}

	return;
}

void make_lastDto1(Matrix src) {
	int last_line = src.rows - 1;
	for (int i = 0; i < src.cols; i++)
	{
		src.data[last_line][i] = 1.0;
	}
	return;
}

void depth_align(Matrix f_rgb, Matrix f_depth, Matrix D2RGB, DepthImage depth, float *dist_dep, 
	float *dist_rgb, Matrix UV_depth, Matrix Cld_3D, Matrix uv_newdep, Matrix uv_newrgb) {
	/*
	* for each pixel in dpeth
	* p1 = (u, v, 1)
	* p2 = cam_invk * p1
	*
	* p2_2 = distortion(p2)
	*
	* p2_2 = d * p2_2;
	* p4 = [p2_2;1]
	* p3 = cam_R * p4 + cam_T
	* p' = (p3.x/p3.z, p3.y/p3.z, 1)
	*
	* p_2' = distotion(p')
	*
	* (u', v') = cam_k * p_2'
	* check (u', v') is in range
	* depth.at(u', v') = d;
	* cld.at(u', v') = p3
	*/
	Matrix p2, p3, p4;
	init_allvar(&p2, &p3, &p4, depth.rows, depth.cols);
	
	mult_f_uv1(UV_depth, f_depth, p2, depth.rows, depth.cols);// f_depth changed 

	distortion(p2, dist_dep, Cld_3D, false);

	Point2Dto3D(depth, Cld_3D);
	// save p2_2 3_dim [xyz]
	/*FILE *fp_p2_2 = fopen("p2_2xyz.txt","w");
	for (int i = 0; i < p2_2.cols; i++)
	{
		fprintf(fp_p2_2,"v ");
		for (int j = 0; j < p2_2.rows; j++)
		{
			fprintf(fp_p2_2, "%f ", p2_2.data[j][i]);
		}
		fprintf(fp_p2_2, "\n");
	}
	fclose(fp_p2_2);*/

	// p4 -> [x y z 1]
	add_oneD(Cld_3D, p4);
	multiply_matrix(D2RGB, p4, p3);// save p3
	/*FILE *fp_p3 = fopen("p3xyz.txt","w");
	for (int i = 0; i < p3.cols; i++)
	{
		fprintf(fp_p3,"v ");
		for (int j = 0; j < p3.rows; j++)
		{
			fprintf(fp_p3, "%f ", p3.data[j][i]);
		}
		fprintf(fp_p3, "\n");
	}
	fclose(fp_p3);*/

	// copy_matrix(p3, p3_3);
	scale_mul_matrix(p3, p3, false);//if z==0 x/z = 0 or 1
	make_lastDto1(p3);

	// distortion(p3_3, dist_dep, p3_2, true); // È¥depth»û±ä

	multiply_matrix(f_depth, p3, uv_newdep);// resuse p2 variable
	multiply_matrix(f_rgb, p3, uv_newrgb); // to rgb save [u'',v'']
	
	//FILE*fp1 = fopen("p1_2.txt","w");// check the p1 original [U V]
	//FILE*fp2 = fopen("p2_2.txt", "w"); // check the p2 destination [U' V']
	//for (int i = 0; i < p2.cols; i++)
	//{
	//	fprintf(fp1, "%7.4f ", p1.data[0][i]);
	//	fprintf(fp1, "%7.4f ", p1.data[1][i]);
	//	fprintf(fp1, "\n");
	//	fprintf(fp2, "%7.4f ", p2.data[0][i]);
	//	fprintf(fp2, "%7.4f ", p2.data[1][i]);
	//	fprintf(fp2, "\n");
	//}
	//fclose(fp1);
	//fclose(fp2);

	free_matrix(p2);
	free_matrix(p3); free_matrix(p4);
}

void new_rgbmap(Matrix uv_newrgb, RGBImage rgb, Matrix *rgb_copy) {
	*rgb_copy =alloc_matrix(uv_newrgb.rows, uv_newrgb.cols);
	int endu = rgb.cols - 1;
	int endv = rgb.rows - 1;
	for (int i = 0; i < uv_newrgb.cols; i++)
	{
		int u2 = uv_newrgb.data[0][i];
		int v2 = uv_newrgb.data[1][i];
		if (check_in_range(u2, v2, rgb.rows, rgb.cols)) {
			//rgb_copy->data[0][i] = rgb.data[endv - v2][(endu-u2) * 3];//bgr
			//rgb_copy->data[1][i] = rgb.data[endv - v2][(endu - u2) * 3 + 1];
			//rgb_copy->data[2][i] = rgb.data[endv - v2][(endu - u2) * 3 + 2];
			rgb_copy->data[0][i] = rgb.data[v2][u2 * 3];//bgr
			rgb_copy->data[1][i] = rgb.data[v2][u2 * 3 + 1];
			rgb_copy->data[2][i] = rgb.data[v2][u2 * 3 + 2];
		}
	}
	return;
}

void process_one(Matrix UV_depth, Matrix Cld_3D, Matrix uv_newdep, Matrix uv_newrgb,
	FilePath depth_path, FilePath rgb_path) {
	
	DepthImage depth, new_depth;
	RGBImage rgb, cld; Matrix rgb_copy;
	FilePath out_path = "./out";
	char name[45] = "1.png";
	load_depimg(depth_path, &depth);//
	load_rgbimg(rgb_path, &rgb);//
	new_depthmap(UV_depth, uv_newdep, depth, &new_depth);
	save_one(new_depth, out_path, name);
	new_rgbmap(uv_newrgb, rgb, &rgb_copy);
	save_cld(Cld_3D, rgb_copy, out_path, name);
	save_rgbpixel(rgb_copy, depth.rows, depth.cols, out_path, name);

	// save 480*640rgb
	// save_newrgb(uv_newdepth, uv_newrgb, out_path, name)

	return;
}

void init_map_val(Matrix *UV_depth, Matrix *Cld_3D, Matrix *uv_newdep, Matrix *uv_newrgb,
	int height, int weidth) {
	*UV_depth = alloc_matrix(3, height*weidth);
	*Cld_3D = alloc_matrix(3, height*weidth);
	*uv_newdep = alloc_matrix(3, height*weidth);
	*uv_newrgb = alloc_matrix(3, height*weidth);
	return;
}

void process_all(FilePath json_path, FilePath img_path, FilePath out_path){
	camera_t cam_depth, cam_rgb;
	D2G Rots, Trans;
	Matrix f_depth, f_rgb, D2RGB;
	Matrix UV_depth, Cld_3D, uv_newdep, uv_newrgb;
	FilenameList namergb_list, namedep_list;

	loadjson(&cam_depth, &cam_rgb, &Rots, &Trans, json_path);
	parametersintoMatrix(cam_depth, cam_rgb, &f_depth, &f_rgb, Rots, Trans, &D2RGB, "zxy");
	float *dist_dep = cam_depth.distortions;
	float *dist_rgb = cam_rgb.distortions;

	//D2RGB.data[0][3] = 0.012991820462048054;// Tx
	//D2RGB.data[1][3] = 300.472934722900391;
	//D2RGB.data[2][3] = -0.31010499596595764;

/*
	int num = 0;
	num = loadfiles(img_path, namergb_list, namedep_list);

    DepthImage depth, new_depth;
	RGBImage rgb;
	CloudImage cld;

	for(int i = 0; i < num; i++ ){
		load_rgbimg(namergb_list.data[i], rgb);
		load_depimg(namedep_list.data[i], &depth);
		process_one(f_rgb, f_depth, D2RGB, rgb, depth, cld, &new_depth);
		char * name = namedep_list.data[i];

		save_one(rgb, cld, new_depth, out_path, name);
	}
*/
	DepthImage depth, new_depth;
	RGBImage rgb;
	//RGBImage rgb;
	//CloudImage cld;
	// rgb = alloc_RGB();
	char *name = "1.png";
	// load one image to generate the map
	load_depimg(name, &depth);
	init_map_val(&UV_depth, &Cld_3D, &uv_newdep, &uv_newrgb, depth.rows, depth.cols);
	depth_align(f_rgb, f_depth, D2RGB, depth, dist_dep, dist_rgb, 
		UV_depth, Cld_3D, uv_newdep, uv_newrgb);
	
	/*printf("D2RGB\n");
	print_matrix(D2RGB);*/

	// load the rest images using the map
	//for (int i = 0; i < image_length; i++)
	//{
	//	// process_each_one(P1, P2, P3, P4);
	//	// equal to process one
	//	/*
	//	*UV_depth UV of depth
	//	*Cld_3D xyz of the 3D
	//	*uv_newdep u'v' of the new depth
	//	*uv_newrgb u''v'' of the new rgb
	//	*/
	//}
	
	// load each image using the parameters to generate ...
	FilePath depth_path = "1.png";
	FilePath rgb_path = "11.raw";
	process_one(UV_depth, Cld_3D, uv_newdep, uv_newrgb, depth_path, rgb_path);
	

// free_all()

}
#endif
