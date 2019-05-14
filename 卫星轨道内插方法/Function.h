#pragma once
#include"DataStructure.h"
#include"include.h"
#include"ClassInterpolate.h"

namespace ClassInter
{
	int FindIndex(int fp, int ep, Time t, vector<Point> data);
	vector<SourceData> GetData(string filename);
	vector<Point> GetKnownPoint(vector<SourceData>& data, int interval);
	void GetOrbitPosition(vector<Point>& orbit_a, vector<Point>&orbit_b, string filename);
	Time Normalization(vector<Point>& data, int begin, int staff, Time t);

}
