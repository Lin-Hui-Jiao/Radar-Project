//
// Created by LiYang on 2021/6/2
//
#pragma once
#include "interface.hpp"
using namespace std;

class caltools
{
public:
	float rad(float d);																										//角度转弧度
	void de2sec(float de, float *sec);																						//将度转化为度分秒
	float normalize_angle(float angle);																						// resolve clockwise angle (degree) to [-180,180)
	float cal_phased_Eangle(double angle, float n, float d, float f);														//计算相控阵雷达的对应表
	float cal_sew_Eangle(double angle, float d, float f);																	//计算缝阵雷达的对应表
	float cal_bamu_Etheta(double theta);
	float cal_bamu_Ephi(double phi);
	void make_angle_table(RadarType type, float n, float d, float f, float *table,bool flag);								//建立特定型号雷达角度与Eangle之间的对应表
	float cal_distance(float lat1, float lng1, float lat2, float lng2);														//根据两点经纬度计算出两点距离
	float cal_angle(float lat1, float lng1, float lat2, float lng2);														//根据两点经纬度计算出两点的角度
	// void lon2theta(float *Rlon, float *Rlat, float Ralt, const float *xlon, const float *xlat, float xalt, float *thephiR); //将某点x与雷达的经纬度绝对位置转化为以雷达为原点的相对位置
	void lon2xyz(float *Rlon, float *Rlat, float Ralt, const float *xlon, const float *xlat, float xalt, float *thephiR);   //将某点p与雷达的经纬度绝对位置转化为以雷达为原点的相对位置xyz
	void trans_axes(float *xyz,float theta,float phi);																		//加上雷达旋转角后的xyz值
	void xyz2tpr(float *xyz,float *thephiR);																					//直角坐标系转球标系
	float choose_rDB(float f, float phi);																					//匹配对应仰角的衰减值
	float cal_alpha(float reduce_db);																						//计算α值
	float cal_La(float alpha, float remote);																				//计算La值
	float cal_Gthetaphi(float G, float etheta, float ephi);																	//计算增益
	float cal_Qt(float _Pt, float remote, float La, float Gthetaphi);														//计算该点的功率密度
};