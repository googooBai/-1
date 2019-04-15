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
		//Insert Corrdinate
		double ix;
		double iy;
		double iz;
		//mm
		double xx;
		double yy;
		double zz;
		double xy;
		double xz;
		double yz;
		//TIME
		Time time_g;
		Time time_d;
		//Error
		double error_x;
		double error_y;
		double error_z;
		//d_Error
		double de_x;
		double de_y;
		double de_z;
		//Flag
		bool is_g;	//几何学轨道是否超限
		bool is_d;	//动力学轨道是否超限
	};
}
#endif // !DATASTUCTURE_H

