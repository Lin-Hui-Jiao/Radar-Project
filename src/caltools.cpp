//
// Created by LiYang on 2021/6/2
//
#include "hiradar/caltools.h"
#include<iostream>
#include<cmath>
#include<map>

#define PRECISION 1000 //角度精度
#define PI 3.1415926
#define EARTH_RADIUS 6371.193 //地球半径 KM

//这个地方时之后要改正的重点，因为原来的雷达能量的计算公式是针对大范围的，但是现在是在城市条件下，所以需要进行改正
/**
 * 将角度（Degree）转换为弧度（Radian）。
 */
float caltools::rad(float d)
{
	return d * PI / 180.0;
}

/**
 * 将度数（Degree）转换为度分秒（Degree-Minute-Second）格式。
 */
void caltools::de2sec(float de,float* sec)
{
	int D = de;
	int F = (de - D) * 60;
	int E = ((de - D) * 60 - F) * 60;
	sec[0] = D;
	sec[1] = F;
	sec[2] = E;
}
/**
 * 规范化角度。将任意角度（如370°或-90°）转换到 [-180, 180) 这个标准区间内，方便统一处理。
 */
float caltools::normalize_angle(float angle)
{
	auto angle_360 = fmod(angle, 360);
	return angle_360 > 180 ? 360 - angle : angle;
}

/**
 * 根据不同的雷达天线类型（相控阵、缝阵、八木），分别计算其方向性系数（E-angle）。这个系数是一个0到1之间的值，
 * 描述了当信号偏离天线正中心时，能量会衰减多少。
 * n是阵元的个数，d是阵元间距，f是工作频率 (Operating Frequency)，单位是兆赫兹(MHz)。
 */
float caltools::cal_phased_Eangle(double angle,float n,float d, float f)
{
	float Eangle;
	float lamda;
	lamda = static_cast<float>(300) / f;
	float aaa = sin(angle / 180 * PI);
	float bbb = sin(n * (PI * d / lamda) * aaa);
	float ccc = n * sin((PI * d / lamda) * aaa);
	Eangle = bbb / ccc;
	return Eangle;
}

float caltools::cal_sew_Eangle(double angle,float d,float f)
{
	float Eangle;
	float lamda;
	lamda = static_cast<float>(300) / f;
	float aaa = sin(angle / 180 * PI);
	float bbb = sin((PI * d / lamda) * aaa);
	float ccc = (PI * d / lamda) * aaa;
	Eangle = bbb / ccc;
	return Eangle;
}

float caltools::cal_bamu_Etheta(double theta)
{
	float Etheta;
	float A=1.3894;
	float thetab=7;
	float aaa=sin(A*theta/thetab);
	float bbb=A*theta/thetab;
	Etheta = aaa / bbb;
	return Etheta;
}

float caltools::cal_bamu_Ephi(double phi)
{
	float Ephi;
	float A=1.3894;
	float phib=20;
	float aaa=sin(A*phi/phib);
	float bbb=A*phi/phib;
	Ephi = aaa / bbb;
	return Ephi;
}

/**
 * 建立预计算查找表，这是一个核心的性能优化。它会在程序启动时，预先计算好所有角度（以极高精度）对应的方向性系数值并存入内存。
 * 在后续海量计算中，程序只需查表，而无需重复进行耗时的三角函数计算。
 */
void caltools::make_angle_table(RadarType type, float n, float d, float f, float *table,bool flag)
{
	switch (type)
	{
	case SEW:
		for (int i = 1; i < 360 * PRECISION; ++i)
		{
			if(flag==true)
			{
				table[i] = cal_sew_Eangle(i / static_cast<double>(PRECISION),d,f);
			}
			else//如果theta和ephi的计算方式不同可以换函数
			{
				table[i] = cal_sew_Eangle(i / static_cast<double>(PRECISION),d,f);
			}
		}
		break;
	case PHASED:
		for (int i = 1; i < 360 * PRECISION; ++i)
		{
			if(flag==true)
			{
				table[i] = cal_phased_Eangle(i / static_cast<double>(PRECISION),n,d,f);
			}
			else//如果theta和ephi的计算方式不同可以换函数
			{
				table[i] = cal_phased_Eangle(i / static_cast<double>(PRECISION),n,d,f);
			}
		}
		break;
	case YAGI:
		for (int i = 1; i < 360 * PRECISION; ++i)
		{
			if(flag==true)
			{
				table[i] = cal_bamu_Etheta(i / static_cast<double>(PRECISION));
			}
			else
			{
				table[i] = cal_bamu_Ephi(i / static_cast<double>(PRECISION));
			}
		}
		break;
	default:
		break;
	}
	//由洛必达法则及共方向时大气损耗最小，反方向时大气损耗最大知:
	table[0]=1;
	table[180*PRECISION]=0;
}

/**
 * 根据两个点的经纬度，计算它们在地球表面的大圆距离。
 * 
 * 
 */
float caltools::cal_distance(float lat1, float lng1, float lat2, float lng2)//lat1第一个点纬度,lng1第一个点经度,lat2第二个点纬度,lng2第二个点经度
{
	float a;
	float b;
	float radLat1 = rad(lat1);
	float radLat2 = rad(lat2);
	a = radLat1 - radLat2;
	b = rad(lng1) - rad(lng2);
	float s = 2 * asin(sqrt(pow(sin(a / 2), 2) + cos(radLat1) * cos(radLat2) * pow(sin(b / 2), 2)));
	s = s * EARTH_RADIUS;
	s = s * 1000;
	return s;
}

//根据两点经纬度计算出两点的角度差，与雷达本身方位角求差或求和即得方位角
//以正北方向为正方向，顺时针转过的角度大小
float caltools::cal_angle(float lat1, float lng1, float lat2, float lng2)
{ 
	float x = lat1 - lat2;
	float y = lng1 - lng2;
	int flag = 0;
	float angle = 0;
	if (y == 0 && x > 0) { angle = 0; flag = 1; }
	if (y == 0 && x < 0) { angle = 180; flag = 1; }
	if (x == 0 && y > 0) { angle = 90; flag = 1; }
	if (x == 0 && y < 0) { angle = 270; flag = 1; }
	if (flag == 0)
	{
		float dislat = cal_distance(lat1, lng2, lat2, lng2);
		float dislng = cal_distance(lat2, lng1, lat2, lng2);
		if (x > 0 && y > 0) angle = atan2(dislng, dislat) / PI * 180;
		if (x < 0 && y > 0) angle = atan2(dislat, dislng) / PI * 180 + 90;
		if (x < 0 && y < 0) angle = atan2(dislng, dislat) / PI * 180 + 180;
		if (x > 0 && y < 0) angle = atan2(dislat, dislng) / PI * 180 + 270;
	}
	return angle;
}

//将某点x与雷达的经纬度绝对位置转化为以雷达为原点的相对位置
//Input：雷达经纬度及高程坐标与空间某点x的经纬高
//		经纬度单位为度分秒，为一个数组，高程坐标为m
//Output：一个数组，thephiR，即点x与雷达的相对位置关系theta、phi、R
//		theta、phi以度为单位，R以米为单位
//假设在赤道附近，纬度差转距离误差忽略
// void caltools::lon2theta(float *Rlon, float *Rlat, float Ralt, const float *xlon, const float *xlat, float xalt, float *thephiR)
// {
// 	//static float tpr[3] = {};//保存x与雷达的相对位置关系
// 	////注：C++ 不支持在函数外返回局部变量的地址，除非定义局部变量为 static 变量。
// 	float Rlondu = Rlon[0] + Rlon[1] / 60 + Rlon[2] / 3600;  //将度秒分统一转换为度，方便计算
// 	float Rlatdu = Rlat[0] + Rlat[1] / 60 + Rlat[2] / 3600;
// 	float xlondu = xlon[0] + xlon[1] / 60 + xlon[2] / 3600;
// 	float xlatdu = xlat[0] + xlat[1] / 60 + xlat[2] / 3600;

// 	float dxy = cal_distance(Rlatdu, Rlondu, xlatdu, xlondu);
// 	//cout << "dxy" << dxy <<endl;
// 	float dz = xalt - Ralt;
// 	float xR = sqrt(pow(dxy, 2) + pow(dz, 2));
// 	float xtheta = cal_angle(xlatdu, xlondu, Rlatdu, Rlondu);//以角度为单位
// 	float xphi = dxy == 0 ? (dz >= 0 ? 90 : -90) : atan(dz / dxy) / PI * 180;
// 	thephiR[0] = xtheta;//方位角
// 	thephiR[1] = xphi;//俯仰角
// 	thephiR[2] = xR;//距离，以米为单位
// }

//计算p点在以正北方向为x轴，正西方向为y轴，正上方为z轴的xyz坐标值（方便旋转）
//假设在赤道附近，纬度差转距离误差忽略
void caltools::lon2xyz(float *Rlon, float *Rlat, float Ralt, const float *xlon, const float *xlat, float xalt, float *xyz)
{
	float Rlondu = Rlon[0] + Rlon[1] / 60 + Rlon[2] / 3600;  //将度秒分统一转换为度，方便计算
	float Rlatdu = Rlat[0] + Rlat[1] / 60 + Rlat[2] / 3600;
	float xlondu = xlon[0] + xlon[1] / 60 + xlon[2] / 3600;
	float xlatdu = xlat[0] + xlat[1] / 60 + xlat[2] / 3600;
	//以雷达为原点的坐标系
	float x=(Rlatdu<xlatdu)?cal_distance(Rlatdu, Rlondu, xlatdu, Rlondu):(-cal_distance(Rlatdu, Rlondu, xlatdu, Rlondu));//经度不变,纬度差转成距离
	float y=(Rlondu>xlondu)?cal_distance(Rlatdu, Rlondu, Rlatdu, xlondu):(-cal_distance(Rlatdu, Rlondu, Rlatdu, xlondu));//纬度不变，经度差转成距离
	//p在R北方则x为正，否则为负
	//p在R西方则y为正，否则为负
	float z=xalt-Ralt;//上方为正，下方为负
	xyz[0]=x;
	xyz[1]=y;
	xyz[2]=z;
}

//加上雷达旋转角后的xyz值
//theta为雷达相对于正北方向的偏移角,绕z轴逆时针转动角度，即雷达指向逆时针偏移正北方向的角度（如北偏西30度为+30，北偏东30度为-30）
//phi为雷达相对于平视方向的偏移角，绕y轴逆时针转动角度，即雷达指向向下偏移平视方向的角度（如向下30度为+30，向上30度为-30）
void caltools::trans_axes(float *xyz,float theta,float phi)
{
	float x=xyz[0];
	float y=xyz[1];
	float z=xyz[2];

	// 将度转换为弧度
	float theta_rad = theta * PI / 180.0;
	float phi_rad = phi * PI / 180.0;
	//绕z轴转动导致的坐标改变
	xyz[0]=x*cos(theta_rad)-y*sin(theta_rad);
	xyz[1]=x*sin(theta_rad)+y*cos(theta_rad);
	xyz[2]=xyz[2];

	x=xyz[0];
	y=xyz[1];
	z=xyz[2];
	//绕y轴转动导致的坐标改变
	xyz[0]=z*sin(phi_rad)+x*cos(phi_rad);
	xyz[1]=xyz[1];
	xyz[2]=z*cos(phi_rad)-x*sin(phi_rad);
}
/**
 * 第三步：将旋转后的XYZ坐标，转换为雷达物理公式真正需要的球坐标格式，
 * 即相对于天线朝向的方位角（theta）、俯仰角（phi）和距离（R）。
 */
void caltools::xyz2tpr(float *xyz,float *thephiR)
{
	//注意:这里的theta phi与球坐标系还是不一样的
	float r=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
	//根据象限确定theta角,theta角范围是(-180,180]
	if((xyz[0]>=0)&&(xyz[1]>0))//第一象限
	{
		thephiR[0]=(xyz[0]==0)?90:atan(xyz[1]/xyz[0])*180/PI;
	}
	else if((xyz[0]>=0)&&(xyz[1]<0))//第四象限
	{
		thephiR[0]=(xyz[0]==0)?-90:atan(xyz[1]/xyz[0])*180/PI;
	}
	else if((xyz[0]<=0)&&(xyz[1]>0))//第二象限
	{
		thephiR[0]=(xyz[0]==0)?90:((atan(xyz[1]/xyz[0])*180/PI)+180);
	}
	else if((xyz[0]<=0)&&(xyz[1]<0))//第三象限
	{
		thephiR[0]=(xyz[0]==0)?-90:((atan(xyz[1]/xyz[0])*180/PI)-180);
	}

	if((xyz[0]>=0)&&(xyz[1]==0)){thephiR[0]=0;}
	if((xyz[0]<=0)&&(xyz[1]==0)){thephiR[0]=180;}

	thephiR[1]=(r==0)?0: asin(xyz[2]/r)*180/PI;//phi的范围是[-90,90]
	thephiR[2]=r;
}



float caltools::choose_rDB(float f,float phi)
{
	//建立不同频率下的损耗表，如有新加频率，再增加表的内容
	if(f < 1000 || f > 20000)return 0;//忽略小于1000Mhz下的损耗 //这个地方还是要修改的，因为频率很高的时候，不知道损耗是怎样的，先假设为0
	float rDBat1000[]={3.4,2.7,2.2,1.55,0.8,0.4,0.55,0.0,0.0};
	float rDBat3000[]={4.4,3.4,2.8,1.9,0.95,0.5,0.08,0.0,0.0};
	float rDBat5000[]={5.0,3.8,3.0,2.1,1.0,0.54,0.19,0.11,0.0};
	float rDBat10000[]={6.8,5.0,4.0,2.6,1.4,0.7,0.24,0.14,0.12};
	std::map<float,float*> choose_rDB_table;
	choose_rDB_table[1000]=rDBat1000;
	choose_rDB_table[3000]=rDBat3000;
	choose_rDB_table[5000]=rDBat5000;
	choose_rDB_table[10000]=rDBat10000;

	float rDB;
	if (phi < 0.5)
	{
		rDB = choose_rDB_table[f][0] - (phi / 0.5) * (choose_rDB_table[f][0] - choose_rDB_table[f][1]);
	}
	else if ((0.5 <= phi) && (phi < 1))
	{
		rDB = choose_rDB_table[f][1] - ((phi - 0.5) / (1 - 0.5)) * (choose_rDB_table[f][1] - choose_rDB_table[f][2]);
	}
	else if ((1 <= phi) && (phi < 2))
	{
		rDB = choose_rDB_table[f][2] - ((phi - 1.0) / (2.0 - 1.0)) * (choose_rDB_table[f][2] - choose_rDB_table[f][3]);
	}
	else if ((2 <= phi) && (phi < 5))
	{
		rDB = choose_rDB_table[f][3] - ((phi - 2.0) / (5.0 - 2.0)) * (choose_rDB_table[f][3] - choose_rDB_table[f][4]);
	}
	else if ((5 <= phi) && (phi < 10))
	{
		rDB = choose_rDB_table[f][4] - ((phi - 5.0) / (10.0 - 5.0)) * (choose_rDB_table[f][4] - choose_rDB_table[f][5]);
	}
	else if ((10 <= phi) && (phi < 30))
	{
		rDB = choose_rDB_table[f][5] - ((phi - 10.0) / (30.0 - 10.0)) * (choose_rDB_table[f][5] - choose_rDB_table[f][6]);
	}
	else if ((30 <= phi) && (phi < 60))
	{
		rDB = choose_rDB_table[f][6] - ((phi - 30.0) / (60.0 - 30.0)) * (choose_rDB_table[f][6] - choose_rDB_table[f][7]);
	}
	else if ((60 <= phi) && (phi < 90))
	{
		rDB = choose_rDB_table[f][7] - ((phi - 60.0) / (90.0 - 60.0)) * (choose_rDB_table[f][7] - choose_rDB_table[f][8]);
	}
	else
	{
		rDB = 0;
	}
	return rDB;
}

float caltools::cal_alpha(float reduce_db)
{
	// float alpha;
	// alpha = reduce_db / 2 / 350 / 1.852;
	// return alpha;
	return reduce_db;
}
//这是什么玩意
float caltools::cal_La(float alpha,float remote)
{
	float La_DB;
	float La;
	// La_DB = 2 * alpha * remote;
	La_DB = alpha * remote / 1000.0;
	La = pow(10, (La_DB/10));
	return La;
}

float caltools::cal_Gthetaphi(float G,float etheta, float ephi)
{
	float Gthetaphi;
	float G_linear = pow(10, G/10.0);  // 将dB转换为线性值
	Gthetaphi = G_linear * abs(etheta) * abs(ephi);
	return Gthetaphi;
}

float caltools::cal_Qt(float _Pt,float remote,float La, float Gthetaphi)
{
	float Qt;
	float Pt_watts = _Pt * 1000.0;  // kW → W
    // float remote_m = remote * 1000.0;  // km → m
	Qt = (Pt_watts * Gthetaphi) / (4 * PI * remote * remote * La);
	return Qt;
}
