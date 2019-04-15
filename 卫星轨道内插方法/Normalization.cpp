#include"Function.h"
namespace ClassInter {
	Time Normalization(vector<Point>& data,int index,int staff, Time t)
	{
		Time start, end;
		start = data[index].time;
		end = data[index+staff].time;
		t = (t - start) / (end - start);
		return t;
	}
}