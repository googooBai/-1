#pragma once
#include<vector>
#include<string>
#include"DataStructure.h"

using std::string;
using std::vector;


vector<SourceData> GetData(string filename);
Point LagrangeInterpolating(vector<SourceData> data, Time t);
Point ChebyshevInterpolating(vector<SourceData> data, Time t);