#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {
	
	Point LarrangeInterpolate::interpolate(Time t)
	{
		//查找内插点的位置
		int posfirst = 0, posend = this->knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, this->knownp);
		//处理边界
		if (index < this->Slope / 2)
			index = 0;
		else if (index >(this->knownp.size() - (this->Slope - this->Slope / 2) - 1))
			index = this->knownp.size() - this->Slope - 1;
		else
			index = index - this->Slope / 2;
		Point result{ 0.0,0.0,0.0,0.0 };
		result.time = t;
		//内插过程
		for (auto k = index; k <= index + this->Slope; k++)
		{
			double tempx = 1, tempy = 1, tempz = 1;
			for (auto i = index; i <= index + this->Slope; i++)
			{
				if (k != i)
				{
					tempx *= (t - this->knownp[i].time) / (this->knownp[k].time - this->knownp[i].time);
					tempy *= (t - this->knownp[i].time) / (this->knownp[k].time - this->knownp[i].time);
					tempz *= (t - this->knownp[i].time) / (this->knownp[k].time - this->knownp[i].time);
				}
			}
			result.x += tempx*this->knownp[k].x;
			result.y += tempy*this->knownp[k].y;
			result.z += tempz*this->knownp[k].z;
		}
		return result;
	}
}