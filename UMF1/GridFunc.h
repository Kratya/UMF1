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
 
template <typename mytype>
struct area
{
	vector<mytype> knots, coefs, step0;
	vector<int> counts, globalIndexKnots;
}; 

template <class mytype>
class edge {
public:
	int startI, endI, step;
	mytype val1, val2, normX, normY;
	vector<int> getVertexes();
};

template <class mytype>
class Grid {
public:
	vector<mytype> ox, oy;
	area<mytype> areaX, areaY;
	vector<int> boundary;
	void readFirst(string fl);
	void readSecond(string fl);
	void readX(string f);
	void readY(string f);
	void initGrid();
	int size();
	int Xsize();
	int Ysize();
	mytype stepX(int i, int j);
	mytype stepY(int i, int j);
	vector<edge<mytype>> getFirst();
	vector<edge<mytype>> getSecond();

private:
	vector<edge<mytype>> first;
	vector<edge<mytype>> second;
	vector<edge<mytype>> third;


	void read(string f, area<mytype>& a);
	void initGrid_XY(vector<mytype>& o, area<mytype>& a);
};

template <typename mytype>
vector<int> edge<mytype>::getVertexes()
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

template <typename mytype>
void Grid<mytype>::readFirst(string fl)
{
	ifstream fin(fl);
	vector<int> roadMap;
	mytype tVal;
	int n = 0, xI1, xI2, yI1, yI2;
	fin >> n;
	int rowsz = areaX.globalIndexKnots.back() + 1;
	for (int i = 0; i < n; ++i)
	{
		fin >> xI1 >> yI1 >> xI2 >> yI2;
		fin >> tVal;
		int xGI1 = areaX.globalIndexKnots[xI1 - 1];
		int yGI1 = areaY.globalIndexKnots[yI1 - 1];
		int xGI2 = areaX.globalIndexKnots[xI2 - 1];
		int yGI2 = areaY.globalIndexKnots[yI2 - 1];
		int startI = xGI1 + yGI1 * rowsz;
		int endI = xGI2 + yGI2 * rowsz;
		int step = 1;
		if (xI1 == xI2)
			step = rowsz;
		first.emplace_back();
		first.back().startI = startI;
		first.back().endI = endI;
		first.back().step = step;
		first.back().val1 = tVal;
		for (int th : first.back().getVertexes())
			boundary.push_back(th);
	}
	sort(boundary.begin(), boundary.end());
}

template <typename mytype>
void Grid<mytype>::readSecond(string fl)
{
	ifstream fin(fl);
	mytype nX, nY;
	mytype tVal;
	vector<int> roadMap;
	int n = 0, xI1, xI2, yI1, yI2;
	fin >> n;
	int rowsz = areaX.globalIndexKnots.back() + 1;
	for (int i = 0; i < n; ++i)
	{
		fin >> xI1 >> yI1 >> xI2 >> yI2 >> nX >> nY;
		int xGI1 = areaX.globalIndexKnots[xI1 - 1];
		int yGI1 = areaY.globalIndexKnots[yI1 - 1];
		int xGI2 = areaX.globalIndexKnots[xI2 - 1];
		int yGI2 = areaY.globalIndexKnots[yI2 - 1];
		int startI = xGI1 + yGI1 * rowsz;
		int endI = xGI2 + yGI2 * rowsz;
		int step = 1;
		if (xI1 == xI2)
			step = rowsz;
		second.emplace_back();
		second.back().startI = startI;
		second.back().endI = endI;
		second.back().step = step;
		second.back().normX = nX;
		second.back().normY = nY;
		fin >> tVal;
		first.back().val1 = tVal;
		for (int th : second.back().getVertexes())
			boundary.push_back(th);
	}
	sort(boundary.begin(), boundary.end());
}

template <typename mytype>
void  Grid<mytype>::readX(string f)
{
	read(f, areaX);
}

template <typename mytype>
void  Grid<mytype>::readY(string f)
{
	read(f, areaY);
}

template <typename mytype>
void  Grid<mytype>::read(string f, area<mytype>& a)
{
	ifstream fin(f);
	int n;
	fin >> n;
	a.knots.resize(n);
	a.counts.resize(n - 1);
	a.coefs.resize(n - 1);
	a.step0.resize(n - 1);
	a.globalIndexKnots.resize(n);
	for (int i = 0; i < n; ++i)
		fin >> a.knots[i];
	a.globalIndexKnots[0] = 0;
	for (int i = 0; i < n - 1; ++i)
	{
		fin >> a.counts[i];
		a.globalIndexKnots[i + 1] = a.globalIndexKnots[i] + a.counts[i];
	}
	for (int i = 0; i < n - 1; ++i)
	{
		fin >> a.coefs[i];
		if (a.coefs[i] != 1)
		{
			a.step0[i] = (a.knots[i + 1] - a.knots[i]) * (1 - a.coefs[i]) / (1 - pow(a.coefs[i], a.counts[i]));
		}
		else
		{
			a.step0[i] = (a.knots[i + 1] - a.knots[i]) / a.counts[i];
		}
	}
}

template <typename mytype>
void Grid<mytype>::initGrid()
{
	initGrid_XY(ox, areaX);
	initGrid_XY(oy, areaY);
}

template <typename mytype>
void Grid<mytype>::initGrid_XY(vector <mytype>& o, area<mytype>& a)
{
	int n = a.counts.size();
	int m = 1;
	for (int i = 0; i < n; ++i)
	{
		m += a.counts[i];
	}
	o.resize(m, 0);
	for (int i = 0, k = 0; i < n; ++i)
	{
		mytype h = a.step0[i];
		o[k] = a.knots[i];
		for (int j = 0; j < a.counts[i]; ++j)
		{
			o[k + 1] = o[k] + h;
			h *= a.coefs[i];
			++k;
		}
	}
}

template <typename mytype>
int Grid<mytype>::size()
{
	return ox.size() * oy.size();
}

template <typename mytype>
int Grid<mytype>::Xsize()
{
	return ox.size();
}

template <typename mytype>
int Grid<mytype>::Ysize()
{
	return oy.size();
}

template <typename mytype>
mytype Grid<mytype>::stepX(int i, int j)
{
	return ox[i] - ox[j];
}

template <typename mytype>
mytype Grid<mytype>::stepY(int i, int j)
{
	return oy[i] - oy[j];
}

template <typename mytype>
vector<edge<mytype>> Grid<mytype>::getFirst()
{
	return first;
}

template <typename mytype>
vector<edge<mytype>> Grid<mytype>::getSecond()
{
	return second;
}
