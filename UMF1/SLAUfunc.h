#pragma once
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <list>
#include "GridFunc.h"

#define test2

using namespace std;

template <class mytype>
class SLAU
{
public:
	mytype eps = 1e-10;
	int maxInter = 100000000, n, m;
	mytype* diag, * l1, * l2, * r1, * r2, * f, * b_g, * x;
	void buildMatrix(Grid<mytype>& grid);
	SLAU<mytype>(Grid<mytype>& grid);
	void applyFirst(Grid<mytype>&);
	void applySecond(Grid<mytype>&);
	mytype iter(mytype*, mytype*, mytype);
	void iterYak();
	void iterZeid();
	mytype gamma(mytype, mytype);
	mytype lambda(mytype, mytype);
	mytype F(mytype, mytype);
	mytype kr1(mytype x, mytype y, int k);
	void printResult(Grid<mytype>& grid, string out);
};

template <typename mytype>
mytype SLAU<mytype>::F(mytype x, mytype y)
{
#ifdef test1
	return 2;
#endif
#ifdef test2
	return x + y;
#endif
#ifdef test3
	return x * x + y * y - 4;
#endif
#ifdef test4
	return x * x * x + y * y * y - 6 * x - 6 * y;
#endif
}


template <typename mytype>
void SLAU<mytype>::buildMatrix(Grid<mytype>& grid)
{
	int xsize = grid.Xsize();
	int ysize = grid.Ysize();
	int n = xsize * ysize;
	list<int> stack;
	unordered_set<int> old;
	auto begin = grid.boundary.begin();
	auto end = grid.boundary.end();
	int startvalue = grid.boundary[2] + xsize;
	while (find(begin, end, startvalue) != end)
		startvalue += xsize;
	stack.push_back(startvalue);
	old.insert(startvalue);
	while (stack.size() > 0)
	{
		int th = stack.front();
		stack.pop_front();
		int j = th % xsize;
		int i = th / xsize;
		mytype hy_j = grid.stepY(i + 1, i);
		mytype hy_jm1 = grid.stepY(i, i - 1);
		mytype hx_im1 = grid.stepX(j, j - 1);
		mytype hx_i = grid.stepX(j + 1, j);
		mytype lambdaL = -lambda(grid.ox[j], grid.oy[i]);
		diag[th] = lambdaL * (-2 / (hx_im1 * hx_i) - 2 / (hy_jm1 * hy_j)) + gamma(grid.ox[j], grid.oy[i]);
		l1[th - 1] = lambdaL * 2 / (hx_im1 * (hx_i + hx_im1));
		r1[th] = lambdaL * 2 / (hx_i * (hx_i + hx_im1));
		l2[th - xsize] = lambdaL * 2 / (hy_jm1 * (hy_j + hy_jm1));
		r2[th] = lambdaL * 2 / (hy_j * (hy_j + hy_jm1));


		f[th] = F(grid.ox[j], grid.oy[i]);

		if (find(begin, end, th + 1) == end
			&& find(old.begin(), old.end(), th + 1) == old.end())
		{
			stack.push_front(th + 1);
			old.insert(th + 1);
		}
		if (find(begin, end, th - 1) == end
			&& find(old.begin(), old.end(), th - 1) == old.end())
		{
			stack.push_front(th - 1);
			old.insert(th - 1);
		}
		if (find(begin, end, th + xsize) == end
			&& find(old.begin(), old.end(), th + xsize) == old.end())
		{
			stack.push_front(th + xsize);
			old.insert(th + xsize);
		}
		if (find(begin, end, th - xsize) == end
			&& find(old.begin(), old.end(), th - xsize) == old.end())
		{
			stack.push_front(th - xsize);
			old.insert(th - xsize);
		}
	}
}

template <typename mytype>
SLAU<mytype>::SLAU<mytype>(Grid<mytype>& grid)
{
	n = grid.size();
	m = grid.Xsize();
	diag = new mytype[n];
	f = new mytype[n];
	b_g = new mytype[n];
	l1 = new mytype[n - 1];
	r1 = new mytype[n - 1];
	for (int i = 0; i < n - 1; ++i)
		l1[i] = r1[i] = 0;
	l2 = new mytype[n - m];
	r2 = new mytype[n - m];
	for (int i = 0; i < n - m; ++i)
		l2[i] = r2[i] = 0;
	for (int i = 0; i < n; ++i)
	{
		diag[i] = 1;
		b_g[i] = f[i] = 0;
	}
}

template <typename mytype>
mytype SLAU<mytype>::kr1(mytype x, mytype y, int k)
{
	switch (k)
	{
	case 0: return 2;
	case 1: return x + y;
	case 2: return x * x + y * y;
	case 3: return x * x * x + y * y * y;
	default: return 0;
	}
}

template <typename mytype>
void SLAU<mytype>::applyFirst(Grid<mytype>& grid)
{
	for (edge<mytype> th : grid.getFirst())
		for (int ith : th.getVertexes())
		{
			int j = ith % m;
			int i = ith / m;

			diag[ith] = 1.0;
			l1[ith] = 0.0;
			r1[ith] = 0.0;
			l2[ith] = 0.0;
			r2[ith] = 0.0;

			f[ith] = kr1(grid.ox[j], grid.oy[i], th.val1);
		}
}

template <typename mytype>
void SLAU<mytype>::applySecond(Grid<mytype>& grid)
{
	int num, i, xsize = grid.Xsize(), j;
	for (edge<mytype> th : grid.getSecond())
	{
		auto vertx = th.getVertexes();
		mytype h = 0;
		for (i = 0; i < vertx.size(); ++i)
		{
			num = vertx[i];
			int ii = num % xsize;
			int jj = num / xsize;
			//f[num] = F(grid.ox[ii], grid.oy[jj]);
			f[num] = th.val1;

			if (th.normX > 0)
			{
				j = num % xsize;
				h = grid.stepX(j, j - 1);
				l1[num - 1] = -1.0 / h;
			}
			else if (th.normX < 0)
			{
				j = num % xsize;
				h = grid.stepX(j + 1, j);
				r1[num] = -1.0 / h;
			}
			else if (th.normY > 0)
			{
				j = num / xsize;
				h = grid.stepY(j, j - 1);
				l2[num - m] = -1.0 / h;
			}
			else if (th.normY < 0)
			{
				j = num / xsize;
				h = grid.stepY(j + 1, j);
				r2[num] = -1.0 / h;
			}
			diag[num] = 1.0 / h;
		}

	}
	for (int i = grid.Ysize() - 1; i >= 0; --i)
	{
		for (int j = 0; j < xsize; ++j)
			cout << diag[i * xsize + j] << " ";
		cout << endl;
	}
}

template <typename mytype>
mytype SLAU<mytype>::gamma(mytype x, mytype y)
{
	return 1;
}

template <typename mytype>
mytype SLAU<mytype>::lambda(mytype x, mytype y)
{
	return 1;
}

template <typename mytype>
void SLAU<mytype>::printResult(Grid<mytype>& grid, string out)
{
	ofstream fout(out);
	int sy = grid.oy.size() - 1;
	int p = n / m;
	for (int i = n - 1; i >= 0; --i)
	{
		if (i % m == m - 1)
		{
			fout << endl;
			fout << "oy[" << sy << "] = " << grid.oy[sy] << endl;
			--sy;
		}
		fout << x[p * m - i % m - 1] << ",";
		if (i % m == 0)
		{
			--p;
			fout << endl;
		}
	}
	fout << ";" << endl;
	for (int i = 0; i < grid.ox.size(); ++i)
	{
		fout << "ox[" << i << "] = " << grid.ox[i] << endl;
	}
	fout << endl;

}

template <typename mytype>
mytype  SLAU<mytype>::iter(mytype* xn, mytype* x, mytype fnorm)
{
	mytype* du2 = l2;
	mytype* du1 = l1;
	mytype* d0 = diag;
	mytype* dl1 = r1;
	mytype* dl2 = r2;
	mytype* b = f;
	int k1 = 1;
	int km = m;
	mytype descrep = 0;
	for (int i = 0; i < n; ++i)
	{
		mytype sum = d0[i] * x[i];
		if (i - k1 >= 0)
		{
			sum += du1[i - k1] * x[i - k1];
		}
		if (i - km >= 0)
		{
			sum += du2[i - km] * x[i - km];
		}
		if (i + k1 < n)
		{
			sum += dl1[i] * x[i + k1];
		}
		if (i + km < n)
		{
			sum += dl2[i] * x[i + km];
		}
		descrep += (b[i] - sum) * (b[i] - sum);
		xn[i] = x[i] + (b[i] - sum) / d0[i];
	}
	descrep = sqrt(descrep);
	return descrep / fnorm;
}

template <typename mytype>
void  SLAU<mytype>::iterYak()
{
	mytype fnorm = 0;
	int countIter;
	x = new mytype[n];
	mytype* xn = new mytype[n];
	for (int i = 0; i < n; ++i)
	{
		x[i] = xn[i] = 0;
		fnorm += f[i];
	}
	fnorm = sqrt(fnorm);
	mytype descrep = iter(xn, x, fnorm);
	for (countIter = 0; countIter <= maxInter && descrep > eps; ++countIter)
	{
		descrep = iter(xn, x, fnorm);
		mytype* temp = xn;
		xn = x;
		x = temp;
	}
	delete[] xn;
	return;
}

template <typename mytype>
void SLAU<mytype>::iterZeid()
{
	int countIter = 0;
	x = new mytype[n];
	mytype* xn = new mytype[n];
	mytype fnorm = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = xn[i] = 0;
		fnorm += f[i];
	}
	fnorm = sqrt(fnorm);

	for (int i = 0; i < n; ++i)
		x[i] = xn[i] = 0;
	mytype descrep = iter(xn, x, fnorm);
	for (countIter = 0; countIter <= maxInter && descrep > eps; ++countIter)
	{
		descrep = iter(x, x, fnorm);
	}
	delete[] xn;
	return;
}