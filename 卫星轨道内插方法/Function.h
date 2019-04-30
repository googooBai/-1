#pragma once

#include"DataStructure.h"
#include"ClassInterpolate.h"
#include"include.h"

namespace ClassInter {

	vector<SourceData> GetData(string filename);	//��ȡ���ι������
	vector<Point> GetKnownPoint(vector<SourceData>& data,int interval);	//��ȡ��֪��
	int FindIndex(int fp,	//Ѱ������t���ڲ�λ��
		int ep,
		Time t, 
		vector<Point> data);
	void GetOrbitPosition(vector<Point>& orbit_a,	//��ȡ����ѧ�������
		vector<Point>&orbit_b,
		string filename);
	Time Normalization(vector<Point>& data,		//��һ������
		int begin, 
		int staff, 
		Time t);
}

	/*void CmpError(Interpolate& inter,	
		int staff,
		vector<SourceData>& data, 
		const vector<Point>& knownp,
		const vector<Point>& ob,
		ostream& os,
		ostream& os2, 
		ostream& os3);*/

