#include<iostream>
#include<fstream>
#include"DataStructure.h"
#include"Function.h"
#include<vector>

using namespace std;

int main()
{
	string filename{ "D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-02.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> knownp = GetKnownPoint(data);
	double start, end;
	start = knownp.front().time;
	end = knownp.back().time;
	vector<Point> interpolating;
	/*
	for (auto i = data.begin(); i !=data.end(); i++)
	{
		interpolating.push_back(LagrangeInterpolating(data, i->time));
	}
	*/
	while (start <= end)
	{
		//interpolating.push_back(LagrangeInterpolating(knownp, start));
		interpolating.push_back(ChebyshevInterpolating(knownp, start));
		start += min/6; //5分钟
	}
	ofstream fout;
	fout.open("output.txt");
	if (!fout.is_open())
		exit(EXIT_FAILURE);
	
	fout.setf(ios::fixed);
	for (auto i = interpolating.begin(); i != interpolating.end(); i++)
	{
		fout.precision(15);
		fout << i->time;
		fout.precision(4);
		fout<< " " << i->x << " " << i->y << " " << i->z << endl;
	}
	fout.close();
	system("pause");
	return 0;
}