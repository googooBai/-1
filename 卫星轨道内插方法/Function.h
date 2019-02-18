#pragma once

#include"DataStructure.h"
#include"ClassInterpolate.h"
#include"include.h"




vector<SourceData> GetData(string filename);
vector<Point> GetKnownPoint(vector<SourceData> data);
int FindIndex(int fp, int ep,Time t, vector<Point> data);
void CmpError(Interpolate& inter, int staff, vector<SourceData> data, vector<Point> knownp,ostream& os);
Point Spline(vector<Point> data, Time t, int staff);//已知点 内插时间点 阶数