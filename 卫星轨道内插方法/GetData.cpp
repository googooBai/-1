#include"Function.h"

namespace ClassInter {

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
			fid >> temp.time_g;
			fid >> temp.x;
			fid >> temp.y;
			fid >> temp.z;
			fid >> temp.xx;
			fid >> temp.yy;
			fid >> temp.zz;
			fid >> temp.xy;
			fid >> temp.xz;
			fid >> temp.yz;
			temp.time_d = temp.time_g*86400.0-4453444800.0;
			temp.is_d = false;
			temp.is_g = false;
			temp.is_knownp = false;
			temp.is_effect_by_knownp = false;
			temp.is_big_error = false;
			data.push_back(temp);
		}
		data.pop_back();
		fid.close();
		return data;
	}
}
