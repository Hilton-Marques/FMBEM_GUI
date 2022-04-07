#ifndef FMM_H
#define FMM_H

#include <vector>
#include "../topoSolver/Face.h"
#include "../topoSolver/Solid.h"
#include "../math/gmres.h"
#include "../math/matrix.h"
#include "Solver.h"


class FMM
{
public:

	FMM();
	~FMM();
	FMM(topoSolver::Solid* solid, int n, int ng);
	std::vector<topoSolver::Vertex*> getLeafNodesFromElements(std::vector<topoSolver::Face*> adjacentElements);
	std::vector<topoSolver::Face*> getFarElements(std::vector<topoSolver::Face*> elements, topoSolver::Face* element);
	std::vector< std::vector<topoSolver::Face*> > *m_elementsLevel;
	void checkWithConventionalBem(topoSolver::Solid* solid, std::vector<topoSolver::Face*> elements, Matrix &fmmAx);
	void checkWithConventionalBemH(topoSolver::Solid* solid, std::vector<topoSolver::Face*> elements);
	Matrix computeBvector();
	Matrix matrixVectorMulti(Matrix &x);
	int m_nLMax;
	int m_nTotalVerts;
	int m_N;
	int m_nDb;
	int m_tot;
	double m_resid;
	Solver* m_solver;
	topoSolver::Solid* m_solid;
	Matrix m_b{ m_nTotalVerts , 1 };
	Matrix m_x{ m_nTotalVerts , 1 };
	Matrix m_Hdb{ m_nTotalVerts , m_nDb }; // prescribed displacements
	Matrix m_db{ m_nDb , 0 }; // prescribed displacements
	void prepareForGMRES();
	void showB();
	void showX();
	void setFMM_ME_sizes();

	};



#endif // !FMM_H
