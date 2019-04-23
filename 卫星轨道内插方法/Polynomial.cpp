#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {

	Polynomial::Polynomial(const vector<Point>& data, int slope) :knownp(data), Slope(slope)
	{
		//Caculate V=TA-L
		T = new MatrixXd*[data.size() - m_P + 1];
		L = new MatrixXd*[data.size() - m_P + 1];
		A = new MatrixXd*[data.size() - m_P + 1];
		V = new MatrixXd*[data.size() - m_P + 1];
		P[0] = new MatrixXd*[data.size() - m_P + 1];
		P[1] = new MatrixXd*[data.size() - m_P + 1];
		P[2] = new MatrixXd*[data.size() - m_P + 1];
		error[0] = new double[data.size() - m_P + 1];
		error[1] = new double[data.size() - m_P + 1];
		error[2] = new double[data.size() - m_P + 1];
		is_knownp_error = new bool[knownp.size()]{ false };
		time_knownp = new int[knownp.size()]{ 0 };
		time_knownp_error = new int[knownp.size()]{ 0 };
		for (int i = 0; i <= knownp.size() - m_P; i++)
		{
			T[i] = new MatrixXd(m_P, Slope + 1);
			A[i] = new MatrixXd(Slope + 1, 3);
			L[i] = new MatrixXd(m_P, 3);
			V[i] = new MatrixXd(m_P, 3);

			MatrixXd Pinv[3];
			for (int j = 0; j < 3; j++)
			{
				(*P[j]) = new MatrixXd(m_P, m_P);
				*(*P[j]) = MatrixXd::Zero(m_P, m_P);
				Pinv[j] = MatrixXd::Zero(m_P, m_P);
			}

			for (int j = 0; j < m_P; j++)
			{
				(*T[i])(j, 0) = 1;
				for (int k = 1; k <= Slope; k++)
					(*T[i])(j, k) = (*T[i])(j, k - 1) *
					(2 * (knownp[i + j].time - knownp[i].time) / (knownp[i + m_P - 1].time - knownp[i].time) - 1);
				(*L[i])(j, 0) = this->knownp[i + j].x;
				(*L[i])(j, 1) = this->knownp[i + j].y;
				(*L[i])(j, 2) = this->knownp[i + j].z;

				//Set P
				(*(*P[0]))(j, j) = 1 / knownp[i + j].xx;
				(*(*P[1]))(j, j) = 1 / knownp[i + j].yy;
				(*(*P[2]))(j, j) = 1 / knownp[i + j].zz;
				//Set P.inverse()
				Pinv[0](j, j) = knownp[i + j].xx;
				Pinv[1](j, j) = knownp[i + j].yy;
				Pinv[2](j, j) = knownp[i + j].zz;
			}


			(*A[i]).col(0) = ((*T[i]).transpose()*(*(*P[0]))*(*T[i])).inverse()*(*T[i]).transpose()*(*(*P[0]))*(*L[i]).col(0);
			(*A[i]).col(1) = ((*T[i]).transpose()*(*(*P[1]))*(*T[i])).inverse()*(*T[i]).transpose()*(*(*P[1]))*(*L[i]).col(1);
			(*A[i]).col(2) = ((*T[i]).transpose()*(*(*P[2]))*(*T[i])).inverse()*(*T[i]).transpose()*(*(*P[2]))*(*L[i]).col(2);
			(*V[i]) = (*T[i])*(*A[i]) - (*L[i]);

			for (int j = 0; j < 3; j++)
			{
				MatrixXd e = ((*V[i]).col(j)).transpose()*(*(*P[j]))*(*V[i]).col(j);
				error[j][i] = sqrt(e(0, 0) / (m_P - Slope - 1));
			}

			MatrixXd Qv[3];
			for (int j = 0; j < 3; j++)
			{
				Qv[j] =
					Pinv[j] -
					(*T[i])*((*T[i]).transpose()*(*(*P[j]))*(*T[i])).inverse()*(*T[i]).transpose();
			}

			int not_error = 0;
			for (int j = 0; j < m_P; j++)
			{
				time_knownp[i + j]++;
				if (abs((*V[i])(j, 0) / error[0][i] / sqrt(Qv[0](j, j))) > 3 ||
					abs((*V[i])(j, 1) / error[1][i] / sqrt(Qv[1](j, j))) > 3 ||
					abs((*V[i])(j, 2) / error[2][i] / sqrt(Qv[2](j, j))) > 3)
					time_knownp_error[i + j]++;
				else
					not_error++;
			}
			if (not_error < Slope + 1)
				continue;
			MatrixXd newT(not_error, Slope + 1);
			MatrixXd newP[3];
			MatrixXd newPinv[3];
			MatrixXd newL(not_error, 3);
			MatrixXd newV(not_error, 3);
			for (int j = 0; j < 3; j++)
			{
				newP[j] = MatrixXd::Zero(not_error, not_error);
				newPinv[j] = MatrixXd::Zero(not_error, not_error);
			}
			int newTptr = 0;
			for (int j = 0; j < m_P; j++)
			{
				if (!(abs((*V[i])(j, 0) / error[0][i] / sqrt(Qv[0](j, j))) > 3 ||
					abs((*V[i])(j, 1) / error[1][i] / sqrt(Qv[1](j, j))) > 3 ||
					abs((*V[i])(j, 2) / error[2][i] / sqrt(Qv[2](j, j))) > 3))
				{
					newT.row(newTptr) = (*T[i]).row(j);

					newP[0](newTptr, newTptr) = (*(*P[0]))(j, j);
					newP[1](newTptr, newTptr) = (*(*P[1]))(j, j);
					newP[2](newTptr, newTptr) = (*(*P[2]))(j, j);

					newPinv[0](newTptr, newTptr) = Pinv[0](j, j);
					newPinv[1](newTptr, newTptr) = Pinv[1](j, j);
					newPinv[2](newTptr, newTptr) = Pinv[2](j, j);

					newL(newTptr, 0) = (*L[i])(j, 0);
					newL(newTptr, 1) = (*L[i])(j, 1);
					newL(newTptr, 2) = (*L[i])(j, 2);

					newTptr++;
				}
			}
			(*A[i]).col(0) = (newT.transpose()*newP[0] * newT).inverse()*newT.transpose()*newP[0] * newL.col(0);
			(*A[i]).col(1) = (newT.transpose()*newP[1] * newT).inverse()*newT.transpose()*newP[1] * newL.col(1);
			(*A[i]).col(2) = (newT.transpose()*newP[2] * newT).inverse()*newT.transpose()*newP[2] * newL.col(2);
		}
		for (int i = 0; i < knownp.size(); i++)
			if ((double)time_knownp_error[i] / time_knownp[i] > 0.2)
				is_knownp_error[i] = true;
			else
				is_knownp_error[i] = false;
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
		Point result{ 0.0,0.0,0.0,0.0 };
		result.time = t;
		t = 2 * (t - knownp[index].time) / (knownp[index + m_P-1].time - knownp[index].time) - 1;
		for (int i = Slope; i >= 0; i--)
		{
			result.x *= t;
			result.x += (*A[index])(i, 0);
			result.y *= t;
			result.y += (*A[index])(i, 1);
			result.z *= t;
			result.z += (*A[index])(i, 2);
		}
		return result;
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

				if (data[i].is_knownp)
				{
					if (this->is_knownp_error[knownp_index])
						data[i].is_g = true;
					else
						data[i].is_g = false;
				}
				else
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
	}
}