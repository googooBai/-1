#pragma once
#include<vector>
#include<string>
#include"DataStructure.h"

using std::string;
using std::vector;


vector<SourceData> GetData(string filename);
vector<Point> GetKnownPoint(vector<SourceData> data);
int FindIndex(int fp, int ep,Time t, vector<Point> data);
Point LagrangeInterpolating(vector<Point> data, Time t, int staff);//��֪�� �ڲ�ʱ��� ���� 
Point ChebyshevInterpolating(vector<Point> data, Time t,int staff);//��֪�� �ڲ�ʱ��� ����