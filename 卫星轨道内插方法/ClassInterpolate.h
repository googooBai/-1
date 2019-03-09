#pragma once
#include"DataStructure.h"
#include"include.h"



class Interpolate
{
public:
	Interpolate(){}
	virtual ~Interpolate(){}
	virtual  Point interpolate(vector<Point> data, Time t, int staff) = 0;
};


class LarrangeInterpolate :public Interpolate
{
public:
	LarrangeInterpolate(){}
	virtual ~LarrangeInterpolate(){}
	virtual Point interpolate(vector<Point> data, Time t, int staff);
};

class ChebyshevInterpolate :public Interpolate
{
public:
	ChebyshevInterpolate(){}
	virtual ~ChebyshevInterpolate(){}
	virtual Point interpolate(vector<Point> data, Time t, int staff);
};

class Spline3 :public Interpolate
{
public:
	Spline3(){}
	virtual ~Spline3(){}
	virtual Point interpolate(vector<Point> data, Time t, int staff);
};