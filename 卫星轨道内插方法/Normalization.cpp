#include"Function.h"

Time Normalization(vector<Point>& data,int index,int staff, Time t)
{
	Time start, end;
	start = data[index].time;
	end = data[index+staff].time;
	t = (t - start) / (end - start);
	for (auto i = index; i<=index+staff; i++)
		data[i].time = (data[i].time - start) / (end - start);
	return t;
}