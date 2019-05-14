#ifndef DATASTUCTURE_H
#define DATASTUCTURE_H
namespace ClassInter {

	using Time = double;
	const Time minute = 1.0 / (24 * 60);	//分钟转换成天数
	const int m_C = 45;	//切比雪夫采样点数
	const int m_P = 15;	//N阶多项式采样点数

	struct Point {	//节点数据类型
		//历元
		Time time;
		//坐标
		double x;
		double y;
		double z;
		//协方差
		double xx;
		double yy;
		double zz;
	};

	struct SourceData
	{
		//坐标
		double x;
		double y;
		double z;
		//内插值
		double ix;
		double iy;
		double iz;
		//协方差
		double xx;
		double yy;
		double zz;
		double xy;
		double xz;
		double yz;
		//历元
		Time time_g;	//几何学轨道历元（单位：天）
		Time time_d;	//动力学轨道历元（单位：秒）
		//内插值与几何学轨道较差
		double error_x;	
		double error_y;	
		double error_z;
		//几何学轨道与动力学轨道较差
		double de_x;
		double de_y;
		double de_z;
		//Flag
		bool is_g;	//几何学轨道是否超限
		bool is_d;	//动力学轨道是否超限
		bool is_knownp;	//是否设为已知点
		bool is_effect_by_knownp;	//是否受粗差已知点影响
		bool is_big_error;	//是否大粗差
	};
	
}
#endif // !DATASTUCTURE_H

