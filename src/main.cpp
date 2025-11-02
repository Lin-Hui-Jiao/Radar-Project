//
// Created by LiYang on 2021/5/23
//

#include<string>
#include <iostream>
#include "hiradar/caltools.h"
#include "hiradar/grid.hpp"
#include "hiradar/radar_pool.hpp"
#include <vector>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include "hiradar/radar.hpp"
#include "hiradar/writer.hpp"

using namespace std;

#define PRECISION 1000 //角度精度
#define PI 3.1415926
#define EARTH_RADIUS 6371.193 //地球半径 KM

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	caltools ct; //工具类对象

	//根据参数建立雷达
	string apg_66 = "AN/APG-66"; //雷达名字
	Position center = {
		{113, 10, 10}, //雷达经度
		{9, 10, 10},   //雷达纬度
		5000,
	};
	float range = 50000; //雷达作用距离
	auto radar = ANAPG66(center, 0, 0);
	unsigned short level = 16;

	//建立以雷达为中心的网格（共100万个）
	auto point_list = CreateGrid(center, 50000, level);
	auto points = point_list.first;
	auto points_n = point_list.second;
	GeoSotChPowerDensityWriter writer("127.0.0.1", 19000, 42, level);
	radar->BindWriter(&writer);

	radar->PowerDensity(points, points_n);

	radar->CapablePowerDensity(points,points_n);
	free(points);
	MPI_Finalize();
	return 0;
}

