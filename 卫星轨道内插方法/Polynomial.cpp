#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {

	//

	Polynomial::Polynomial(const vector<Point>& data, int slope) :knownp(data), Slope(slope)
	{
		//Caculate V=TA-L
		A = new MatrixXd*[knownp.size() - m_P + 1];
		error[0] = new double[knownp.size() - m_P + 1];
		error[1] = new double[knownp.size() - m_P + 1];
		error[2] = new double[knownp.size() - m_P + 1];
		is_knownp_error = new bool[knownp.size()]{ false };
		for (int i = 0; i <= knownp.size() - m_P; i++)
		{
			A[i] = new MatrixXd(Slope + 1, 3);
			MatrixXd newT(m_P, Slope+1);
			MatrixXd newPx = MatrixXd::Zero(m_P, m_P);
			MatrixXd newPy = MatrixXd::Zero(m_P, m_P);
			MatrixXd newPz = MatrixXd::Zero(m_P, m_P);
			MatrixXd newQx = MatrixXd::Zero(m_P, m_P);
			MatrixXd newQy = MatrixXd::Zero(m_P, m_P);
			MatrixXd newQz = MatrixXd::Zero(m_P, m_P);
			MatrixXd newL(m_P, 3);
			MatrixXd V;

			vector<int> knownp_index(m_P);
			for (int j = 0; j < m_P; j++)
				knownp_index[j] = i + j;

			for (int j = 0; j < m_P; j++)
			{
				//Set T
				newT(j, 0) = 1;
				for (int k = 1; k <= Slope; k++)
					newT(j, k) = newT(j, k - 1) *
					(2 * (knownp[i + j].time - knownp[i].time) / (knownp[i + m_P - 1].time - knownp[i].time) - 1);
				//Set L
				newL(j, 0) = this->knownp[i + j].x;
				newL(j, 1) = this->knownp[i + j].y;
				newL(j, 2) = this->knownp[i + j].z;
				//Set P
				newPx(j, j) = knownp[i + j].xx;
				newPy(j, j) = knownp[i + j].yy;
				newPz(j, j) = knownp[i + j].zz;
				//Set Q
				newQx(j, j) = 1 / knownp[i + j].xx;
				newQy(j, j) = 1 / knownp[i + j].yy;
				newQz(j, j) = 1 / knownp[i + j].zz;
			}

			//Test
			/*(*A[i]).col(0) = (newT.transpose()*newPx*newT).inverse()*newT.transpose()*newPx*newL.col(0);
			(*A[i]).col(1) = (newT.transpose()*newPy*newT).inverse()*newT.transpose()*newPy*newL.col(1);
			(*A[i]).col(2) = (newT.transpose()*newPz*newT).inverse()*newT.transpose()*newPz*newL.col(2);*/

			bool check;	//粗差探测
			do
			{
				check = false;
				MatrixXd T(newT), L(newL),
					Px(newPx), Py(newPy), Pz(newPz),
					Qx(newQx), Qy(newQy), Qz(newQz);
				
				(*A[i]).col(0) = (T.transpose()*Px*T).inverse()*T.transpose()*Px*L.col(0);
				(*A[i]).col(1) = (T.transpose()*Py*T).inverse()*T.transpose()*Py*L.col(1);
				(*A[i]).col(2) = (T.transpose()*Pz*T).inverse()*T.transpose()*Pz*L.col(2);
				V = T*(*A[i]) - L;
				//Caculate Mean Error
				MatrixXd e;
				e = (V.col(0)).transpose()*Px*V.col(0);
				error[0][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
				e = (V.col(1)).transpose()*Py*V.col(1);
				error[1][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
				e = (V.col(2)).transpose()*Pz*V.col(2);
				error[2][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
				//Caculate Qvv
				MatrixXd Qvv_x, Qvv_y, Qvv_z;
				Qvv_x = Qx - T*(T.transpose()*Px*T).inverse()*T.transpose();
				Qvv_y = Qy - T*(T.transpose()*Py*T).inverse()*T.transpose();
				Qvv_z = Qz - T*(T.transpose()*Pz*T).inverse()*T.transpose();
				//Check Error
				int index_of_error = -1;
				double max_error=0;
				double ex, ey, ez;
				for (int j = 0; j < T.rows(); j++)
				{
					ex = abs(V(j, 0) / error[0][i] / sqrt(Qvv_x(j, j)));
					ey = abs(V(j, 1) / error[1][i] / sqrt(Qvv_y(j, j)));
					ez = abs(V(j, 2) / error[2][i] / sqrt(Qvv_z(j, j)));
					if (ex > 3 || ey > 3 || ez > 3)					{
						/*if (ex > max_error)
						{
							max_error = ex;
							index_of_error = j;
						}
						if (ey > max_error)
						{
							max_error = ey;
							index_of_error = j;
						}
						if (ez > max_error)
						{
							max_error = ez;
							index_of_error = j;
						}*/
						if (ex*ex + ey*ey + ez*ez > max_error)
						{
							max_error = ex*ex + ey*ey + ez*ez;
							index_of_error = j;
						}
						check = true;
					}
				}
				if (check)
				{
					is_knownp_error[knownp_index[index_of_error]] = true;
					for (int j = index_of_error; j != knownp_index.size() - 1; j++)
					{
						knownp_index[j] = knownp_index[j + 1];
					}
					knownp_index.pop_back();
					newT = MatrixXd(T.rows() - 1, Slope + 1);
					newL = MatrixXd(T.rows() - 1, 3);
					newPx = MatrixXd(T.rows() - 1, T.rows()-1);
					newPy = MatrixXd(T.rows() - 1, T.rows()-1);
					newPz = MatrixXd(T.rows() - 1, T.rows()-1);
					newQx = MatrixXd(T.rows() - 1, T.rows() - 1);
					newQy = MatrixXd(T.rows() - 1, T.rows() - 1);
					newQz = MatrixXd(T.rows() - 1, T.rows() - 1);
					for (int j = 0; j != index_of_error; j++)
					{
						newT.row(j) = T.row(j);
						newL.row(j) = L.row(j);
						newPx(j, j) = Px(j, j);
						newPy(j, j) = Py(j, j);
						newPz(j, j) = Pz(j, j);
					}
					for (int j = index_of_error + 1; j < T.rows(); j++)
					{
						newT.row(j - 1) = T.row(j);
						newL.row(j - 1) = L.row(j);
						newPx(j - 1, j - 1) = Px(j - 1, j - 1);
						newPy(j - 1, j - 1) = Py(j - 1, j - 1);
						newPz(j - 1, j - 1) = Pz(j - 1, j - 1);
					}
				}
			} while (check);
		}
	}

	Point Polynomial::interpolate(Time t)
	{
		int posfirst = 0, posend = knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, knownp);
		//处理边界
		if (index < m_P / 2)
			index = 0;
		else if (index >this->knownp.size()-1 - (m_P- m_P/2-1))
			index = this->knownp.size() - m_P;
		else
			index = index - m_P / 2;
		Point res{ 0.0,0.0,0.0,0.0 };
		res.time = t;
		t = 2 * (t - knownp[index].time) / (knownp[index + m_P-1].time - knownp[index].time) - 1;
		for (int i = Slope; i >= 0; i--)
		{
			res.x *= t;
			res.x += (*A[index])(i, 0);
			res.y *= t;
			res.y += (*A[index])(i, 1);
			res.z *= t;
			res.z += (*A[index])(i, 2);
		}
		return res;
	}

	Point Polynomial::interpolate(Time t,ostream& os)
	{
		int posfirst = 0, posend = knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, knownp);
		//处理边界
		if (index < m_P / 2)
			index = 0;
		else if (index >this->knownp.size() - 1 - (m_P - m_P / 2 - 1))
			index = this->knownp.size() - m_P;
		else
			index = index - m_P / 2;
		Point res{ 0.0,0.0,0.0,0.0 };
		res.time = t;
		t = 2 * (t - knownp[index].time) / (knownp[index + m_P - 1].time - knownp[index].time) - 1;
		for (int i = Slope; i >= 0; i--)
		{
			res.x *= t;
			res.x += (*A[index])(i, 0);
			res.y *= t;
			res.y += (*A[index])(i, 1);
			res.z *= t;
			res.z += (*A[index])(i, 2);
		}
		//test
		os << "inter point: " << endl;
		os
			<< setw(22) << "Time"
			<< setw(15) << "x"
			<< setw(15) << "y"
			<< setw(15) << "z"
			<< endl;
		os.setf(ios_base::fixed);
		os.precision(15);
		os << setw(22) << res.time;
		os.precision(3);
		os
			<< setw(15) << res.x
			<< setw(15) << res.y
			<< setw(15) << res.z
			<< endl;

		os << "known point: " << endl;
		os
			<< setw(22) << "Time"
			<< setw(15) << "x"
			<< setw(15) << "y"
			<< setw(15) << "z"
			<< setw(25) << "xx"
			<< setw(25) << "yy"
			<< setw(25) << "zz"
			<< endl;
		for (int j = index; j <= index + m_C; j++)
		{
			os.setf(ios_base::fixed);
			os.precision(15);
			os << setw(22) << knownp[j].time;
			os.precision(3);
			os
				<< setw(15) << knownp[j].x
				<< setw(15) << knownp[j].y
				<< setw(15) << knownp[j].z
				;
			os.precision(15);
			os
				<< setw(25) << knownp[j].xx
				<< setw(25) << knownp[j].yy
				<< setw(25) << knownp[j].zz
				<< endl;
		}
		return res;
	}

	void Polynomial::CmpError(
		vector<SourceData>& data,
		const vector<Point>& ob,
		ostream& os,
		ostream& os2,
		ostream& os3,
		ostream& os4
	)
	{
		double error_x = 0.0, error_y = 0.0, error_z = 0.0;
		double RMSE_X = 0.0, RMSE_Y = 0.0, RMSE_Z = 0.0;
		int not_knownp = 0;

		os
			<< setw(22) << "Time"
			<< setw(15) << "error x"
			<< setw(15) << "error y"
			<< setw(15) << "error z"
			<< endl;
		os2
			<< setw(22) << "Time"
			<< setw(15) << "error x"
			<< setw(15) << "error y"
			<< setw(15) << "error z"
			<< endl;

		os3
			<< setw(22) << "Time"
			<< setw(15) << "g_error x"
			<< setw(15) << "g_error y"
			<< setw(15) << "g_error z"
			<< setw(15) << "d_error x"
			<< setw(15) << "d_error y"
			<< setw(15) << "d_error z"
			<< endl;

		for (int i = 0; i<data.size(); i++)
		{
			int ob_index;
			{
				double t = data[i].time_d;
				int fp = 0, ep = ob.size() - 1;
				ob_index = (fp + ep) >> 1;
				while (fp < ep)
				{
					if (abs(t - ob[ob_index].time) <= 1e-6)
						break;
					else if (t - ob[ob_index].time > 1e-6)
					{
						fp = ob_index;
						ob_index = (fp + ep) >> 1;
					}
					else if (t - ob[ob_index].time < -1e-6)
					{
						ep = ob_index;
						ob_index = (fp + ep) >> 1;
					}
				}
			}

			data[i].de_x = ob[ob_index].x - data[i].x;
			data[i].de_y = ob[ob_index].y - data[i].y;
			data[i].de_z = ob[ob_index].z - data[i].z;
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= m_P / 2 && index <= this->knownp.size()-1 - (m_P - m_P / 2 - 1))
			{
				if (!data[i].is_knownp)
				{
					Point temp = this->interpolate(data[i].time_d); 
					data[i].ix = temp.x;
					data[i].iy = temp.y;
					data[i].iz = temp.z;
				}
				else
				{
					data[i].ix = data[i].x;
					data[i].iy = data[i].y;
					data[i].iz = data[i].z;
				}


				//x
				data[i].error_x = data[i].ix-data[i].x;
				if (!data[i].is_knownp)
				{
					RMSE_X += data[i].error_x*data[i].error_x;
					not_knownp++;
				}
				//y
				data[i].error_y = data[i].iy - data[i].y;
				if(!data[i].is_knownp)
					RMSE_Y += data[i].error_y*data[i].error_y;
				//z
				data[i].error_z = data[i].iz - data[i].z;
				if(!data[i].is_knownp)
					RMSE_Z += data[i].error_z*data[i].error_z;

				os.setf(ios_base::fixed);
				os.precision(15);
				os << setw(22) << data[i].time_d;
				os.precision(3);
				os
					<< setw(15) << data[i].error_x
					<< setw(15) << data[i].error_y
					<< setw(15) << data[i].error_z
					<< endl;

				os2.setf(ios_base::fixed);
				os2.precision(15);
				os2 << setw(22) << data[i].time_d;
				os2.precision(3);
				os2
					<< setw(15) << data[i].de_x
					<< setw(15) << data[i].de_y
					<< setw(15) << data[i].de_z
					<< endl;
			}
		}

		RMSE_X = sqrt(RMSE_X / not_knownp);
		RMSE_Y = sqrt(RMSE_Y / not_knownp);
		RMSE_Z = sqrt(RMSE_Z / not_knownp);

		//粗差处理
		int total{ 0 };
		int mistake{ 0 };
		int hit{ 0 };
		int knownp_index = 0;
		int big_error = 0;
		int hit_big_error = 0;
		for (int i = 0; i < data.size(); i++)
		{
			if (data[i].is_knownp)
			{
				if (this->is_knownp_error[knownp_index++])
					data[i].is_g = true;
				else
					data[i].is_g = false;
			}
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= m_P / 2 && index <= this->knownp.size()-1 - (m_P - m_P / 2 - 1))
			{
				index= index - m_P / 2;
				if (abs(data[i].de_x) > 0.06 ||
					abs(data[i].de_y) > 0.06 ||
					abs(data[i].de_z) > 0.06)
				{
					data[i].is_d = true;
					total++;
				}
				else
					data[i].is_d = false;

				if(!data[i].is_knownp)
				{
					if (abs(data[i].error_x) > 2 * RMSE_X ||
						abs(data[i].error_y) > 2 * RMSE_Y ||
						abs(data[i].error_z) > 2 * RMSE_Z)
					/*if (abs(data[i].error_x) >  this->error[0][index] * 2 ||
						abs(data[i].error_y) >  this->error[1][index] * 2 ||
						abs(data[i].error_z) >  this->error[2][index] * 2)*/
						data[i].is_g = true;
					else
						data[i].is_g = false;


					if (data[i].is_g&&data[i].is_d)
						hit++;
					else if (data[i].is_g && !data[i].is_d)
						mistake++;
				}

				if (abs(data[i].de_x) > 0.1 ||
					abs(data[i].de_y) > 0.1 ||
					abs(data[i].de_z) > 0.1)
				{
					big_error++;
					if (data[i].is_g)
						hit_big_error++;
				}

				if (!data[i].is_g)
				{
					os3.setf(ios_base::fixed);
					os3.precision(15);
					os3 << setw(22) << data[i].time_d;
					os3.precision(3);
					os3
						<< setw(15) << data[i].error_x
						<< setw(15) << data[i].error_y
						<< setw(15) << data[i].error_z
						<< setw(15) << data[i].de_x
						<< setw(15) << data[i].de_y
						<< setw(15) << data[i].de_z
						<< endl;
				}
			}
		}
		std::cout << "Amount of error: " << total << endl;
		std::cout << "Amount of big error(more than 1cm): " << big_error << "	  eliminate: " << hit_big_error << endl;
		std::cout.setf(ios_base::fixed);
		std::cout.precision(1);
		std::cout << "Rate of sucessful: %" << (double)hit / total * 100 << endl;
		std::cout << "Rate of mistake: %" << (double)mistake / data.size() * 100 << endl;
		std::cout.precision(3);
		std::cout
			<< setw(8) << RMSE_X
			<< setw(8) << RMSE_Y
			<< setw(8) << RMSE_Z
			<< endl;
		//os4 << "Amount of error: " << total << endl;
		//os4 << "Amount of big error(more than 1cm): " << big_error << "	  eliminate: " << hit_big_error << endl;
		//os4.setf(ios_base::fixed);
		//os4.precision(1);
		//os4 << "Rate of sucessful: %" << (double)hit / total * 100 << endl;
		//os4 << "Rate of mistake: %" << (double)mistake / data.size() * 100 << endl;
		//os4.precision(3);
		//os4
		//	<< setw(8) << RMSE_X
		//	<< setw(8) << RMSE_Y
		//	<< setw(8) << RMSE_Z
		//	<< endl;
		total = 0;
		hit = 0;
		mistake = 0;
		for (int i = 0; i != data.size(); i++)
		{
			if (data[i].is_knownp)
			{
				os4.setf(ios_base::fixed);
				os4.precision(15);
				os4 << setw(22) << data[i].time_d;
				
				if (data[i].is_g)
					os4 << setw(15) << 1;
				else
					os4 << setw(15) << 0;

				if (data[i].is_d)
				{
					total++;
					os4 << setw(15) << 1 << endl;
				}
				else
					os4 << setw(15) << 0 << endl;

				if (data[i].is_g&&data[i].is_d)
					hit++;
				if (data[i].is_g&&!data[i].is_d)
					mistake++;
			}
		}
		std::cout << "Amount of known point: " << total << endl;
		std::cout << "Amount of hit known point: " << hit << endl;
		std::cout << "Amount of mistake konwn point: " << mistake << endl;
	}
}