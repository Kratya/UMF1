#include "SLAUfunc.h"
using namespace std;

int main()
{
	Grid<double> grid;
	grid.readX("gX.txt");
	grid.readY("gY.txt");
	grid.initGrid();
	grid.readFirst("FirstCondition.txt");
	grid.readSecond("SecondCondition.txt");
	SLAU<double> slae(grid);
	slae.buildMatrix(grid);
	slae.applyFirst(grid);
	slae.applySecond(grid);
	slae.iterZeid();
	slae.printResult(grid, "out.txt");
	return 0;
}