#pragma once

#include"DataStructure.h"
#include"ClassInterpolate.h"
#include"include.h"

namespace ClassInter {

	vector<SourceData> GetData(string filename);	//获取几何轨道数据
	vector<Point> GetKnownPoint(vector<SourceData>& data,int interval);	//获取已知点
	int FindIndex(int fp,	//寻找星历t的内插位置
		int ep,
		Time t, 
		vector<Point> data);
	void GetOrbitPosition(vector<Point>& orbit_a,	//获取动力学轨道数据
		vector<Point>&orbit_b,
		string filename);
	Time Normalization(vector<Point>& data,		//归一化处理
		int begin, 
		int staff, 
		Time t);
}

	/*void CmpError(Interpolate& inter,	
		int staff,
		vector<SourceData>& data, 
		const vector<Point>& knownp,
		const vector<Point>& ob,
		ostream& os,
		ostream& os2, 
		ostream& os3);*/

