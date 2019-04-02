#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {

	void CmpError(Interpolate& inter, int staff, vector<SourceData> data,vector<Point> knownp,ostream& os)
	{
		double error_x=0.0, error_y=0.0, error_z=0.0, error_3d=0.0;
		double RMSE_X=0.0, RMSE_Y=0.0, RMSE_Z=0.0, RMSE_3D=0.0;
		os 
			<< setw(22) << "Time"
			<< setw(15) <<"error x"
			<< setw(15) <<"error y"
			<< setw(15) <<"error z"
			<< setw(15) <<"error 3d"
			<<endl;
		for (auto i = data.begin(); i != data.end(); i++)
		{
			int index = FindIndex(0, knownp.size() - 1, i->time, knownp);
			if (index >= staff / 2 && index <= (knownp.size() - (staff - staff / 2) - 1))
			{
				Point temp=inter.interpolate(knownp, i->time, staff);
				//x
				error_x = (temp.x - i->x);
				RMSE_X += error_x*error_x;
				//y
				error_y = (temp.y - i->y);
				RMSE_Y += error_y*error_y;
				//z
				error_z = (temp.z - i->z);
				RMSE_Z += error_z*error_z;
				//3d
				error_3d = sqrt(error_x*error_x+error_y*error_y+error_z*error_z);
				RMSE_3D += error_3d*error_3d;

				os.setf(ios_base::fixed);
				os.precision(15);
				os <<setw(22)<< i->time;
				os.precision(3);
				os
					<< setw(15) << error_x
					<< setw(15) << error_y
					<< setw(15) << error_z
					<< setw(15) << error_3d
					<< endl;
			}
		}

		os
			<< setw(22) << "Total"
			<< setw(15) << "RMSE_X"
			<< setw(15) << "RMSE_Y"
			<< setw(15) << "RMSE_Z"
			<< setw(15) << "RMSE_3D"
			<< endl;

		RMSE_X = sqrt(RMSE_X / data.size());
		RMSE_Y = sqrt(RMSE_Y / data.size());
		RMSE_Z = sqrt(RMSE_Z / data.size());
		RMSE_3D = sqrt(RMSE_3D / data.size());

		os.precision(3);
		os
			<< setw(22) << " "
			<< setw(15) << RMSE_X
			<< setw(15) << RMSE_Y
			<< setw(15) << RMSE_Z
			<< setw(15) << RMSE_3D
			<< endl;
	}
}

