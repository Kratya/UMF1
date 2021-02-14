#pragma once
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <list>

using namespace std;

// myType = double, float
 
template <typename myType>
struct area
{
	vector<myType> knots, coefs, step0;
	vector<int> counts, globalIndexKnots;
}; 

template <class myType>
class edge {
public:
	int startI, endI, step;
	myType val1, val2;
	vector<int> getVertexes();
};

template <class myType>
class Grid {
public:
	vector<myType> ox, oy;
	area<myType> areaX, areaY;
	vector<int> boundary;
	void readFirst(const char* fl);
	void readSecond(const char* fl);
	void readX(string f);
	void readY(string f);
	void initGrid();
	int size();
	int Xsize();
	int Ysize();
	myType stepX(int i, int j);
	myType stepY(int i, int j);
	vector<edge<myType>> getFirst();
	vector<edge<myType>> getSecond();

private:
	vector<edge<myType>> first;
	vector<edge<myType>> second;
	vector<edge<myType>> third;


	void read(string f, area<myType>& a);
	void initGrid_XY(vector<myType>& o, area<myType>& a);
};

template <typename myType>
vector<int> edge<myType>::getVertexes()
{
	vector<int> rs;
	int pos = startI;
	int end = endI;
	while (pos <= end)
	{
		rs.push_back(pos);
		pos += step;
	}
	return rs;
}