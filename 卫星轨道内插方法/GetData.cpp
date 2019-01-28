#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include"DataStructure.h"

using namespace std;

vector<SourceData> GetData(string filename)
{
	ifstream fid;
	fid.open(filename);
	if (!fid.is_open())
	{
		exit(EXIT_FAILURE);
	}
	
	string strline;
	getline(fid, strline);
	getline(fid, strline);
	vector<SourceData> data;
	while (!fid.eof())
	{
		SourceData temp;
		fid >> temp.time;
		fid >> temp.x;
		fid >> temp.y;
		fid >> temp.z;
		fid >> temp.xx;
		fid >> temp.yy;
		fid >> temp.zz;
		fid >> temp.xy;
		fid >> temp.xz;
		fid >> temp.yz;
		data.push_back(temp);		
	}
	data.pop_back();
	fid.close();
	return data;
}
