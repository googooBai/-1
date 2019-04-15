#include"Function.h"

namespace ClassInter {

	void GetOrbitPosition(vector<Point>& orbit_a, vector<Point>&orbit_b, string filename)
	{
		ifstream fid;
		fid.open(filename);
		if (!fid.is_open())
			exit(EXIT_FAILURE);
		
		orbit_a.clear();
		orbit_b.clear();
		
		string temp; Time t;
		Point a, b;
		while (!fid.eof())
		{
			fid >> temp;
			fid >> t;
			a.time = b.time = t;
			getline(fid, temp);
			getline(fid, temp);
			getline(fid, temp);
			getline(fid, temp);
			fid >> temp;
			fid >> a.x>> a.y>> a.z;
			getline(fid, temp);
			fid >> temp;
			fid >> b.x >> b.y >> b.z;
			getline(fid, temp);
			orbit_a.push_back(a);
			orbit_b.push_back(b);
		}
		orbit_a.pop_back();
		orbit_b.pop_back();
		fid.close();
	}
}