#pragma once
#include<vector>
#include<string>
#include"DataStructure.h"

using std::string;
using std::vector;


vector<SourceData> GetData(string filename);
vector<Point> SetKnownPoint(vector<SourceData> data);
Point LagrangeInterpolating(vector<SourceData> data, Time t);