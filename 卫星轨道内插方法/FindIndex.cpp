#include"Function.h"

namespace ClassInter {

	int FindIndex(int fp, int ep,Time t, vector<Point> data)
	{
		int index = (fp + ep) >> 1;
		while (fp < ep)
		{
			if (t > data[index].time&&t <= data[index + 1].time)
				break;
			if (t <= data[index].time)
			{
				ep = index;
				index = (fp + ep) >>1;
			}
			if (t > data[index + 1].time)
			{
				fp = index + 1;
				index = (fp + ep) >>1;
			}
		}
		return index;
	}
}
