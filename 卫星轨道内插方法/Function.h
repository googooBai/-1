#pragma once
#include<vector>
#include<string>
#include"DataStructure.h"

using std::string;
using std::vector;


vector<SourceData> GetData(string filename);
vector<Point> GetKnownPoint(vector<SourceData> data);
Point LagrangeInterpolating(vector<Point> data, Time t);
Point ChebyshevInterpolating(vector<Point> data, Time t);