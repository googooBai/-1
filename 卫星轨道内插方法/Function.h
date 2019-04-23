#pragma once

#include"DataStructure.h"
#include"ClassInterpolate.h"
#include"include.h"

namespace ClassInter {

	vector<SourceData> GetData(string filename);
	vector<Point> GetKnownPoint(vector<SourceData>& data,int interval);
	int FindIndex(int fp,
		int ep,
		Time t, 
		vector<Point> data);
	void CmpError(Interpolate& inter,
		int staff,
		vector<SourceData>& data, 
		const vector<Point>& knownp,
		const vector<Point>& ob,
		ostream& os,
		ostream& os2, 
		ostream& os3);
	Time Normalization(vector<Point>& data,
		int index, 
		int staff, 
		Time t);
	void GetOrbitPosition(vector<Point>& orbit_a,
		vector<Point>&orbit_b,
		string filename);
}


