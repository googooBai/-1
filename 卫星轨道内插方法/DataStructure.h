#ifndef DATASTUCTURE_H
#define DATASTUCTURE_H

namespace ClassInter {

	using Time = double;
	const Time minute = 1.0 / (24 * 60);

	struct Point {
		//Coordinate
		double x;
		double y;
		double z;
		//Time
		Time time;
	};

	struct SourceData
	{
		//Coordinate
		double x;
		double y;
		double z;
		//mm
		double xx;
		double yy;
		double zz;
		double xy;
		double xz;
		double yz;
		//TIME
		Time time;
	};
}
#endif // !DATASTUCTURE_H

