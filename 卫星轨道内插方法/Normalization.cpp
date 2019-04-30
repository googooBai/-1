#include"Function.h"
namespace ClassInter {
	Time Normalization(vector<Point>& data,int begin,int staff, Time t)
	{
		Time start, end;
		start = data[begin].time;
		end = data[begin+staff].time;
		t = 2*(t-start)/(end-start)-1;
		return t;
	}
}