#ifndef SOLVER_H
#define SOLVER_H

#include "../topoSolver/Face.h"
#include "../math/RTable.h"
#include "../math/GaussQuadrature.h"
#include "../math/Matrix.h"

#include <complex>
#include <vector>

enum class SolverType {
	bVector, GMRES
};

class Solver
{
public:
	Solver(const int &N, const int &nGFMM, const int &nNodes);
	~Solver();
	void CalculateME(topoSolver::Face *element);
	void CalculateMEForGMRES(topoSolver::Face* element);
	friend std::complex<double> operator*(std::vector<std::complex<double>>& dR , Point& x);
	void CalculateNearInt(std::vector<topoSolver::Vertex*> sourcesNodes, topoSolver::Face* element);
	void CalculateFarInt(std::vector<topoSolver::Vertex*> sourcesNodes, topoSolver::Face* element);
	void CalculateMMT(topoSolver::Face* element);
	static inline Point Multi(double p[] , topoSolver::Vertex* v[3]) {
		return Point(p[0] * v[0]->m_coord + p[1] * v[1]->m_coord + p[2] * v[2]->m_coord);
	}
	double Dot(std::complex<double> R[], std::vector<std::complex<double>>* Sb);
	void reset() { m_Hd = Matrix(m_nNodes, 1); m_Gt = Matrix(m_nNodes, 1); };
	int getPos(const int& N, const int& M);
  void prepareElementMatrix(std::vector<topoSolver::Vertex*> sourceNodes, topoSolver::Face* element);
	void CalculateNearIntMatrix(std::vector<topoSolver::Vertex*> sourceNodes, topoSolver::Face*  element);
	const int m_nNodes;
	const int m_N; // truncation term
	const int m_nGFMM; // gauss point for FMM
	const int m_nGBEM = 40; // gauss point for BEM
	int m_NTot; // total of truncation terms 
	Matrix m_Gt{ m_nNodes,1 };
	Matrix m_Hd{ m_nNodes,1 };
	RTable m_RTable;
	GaussQuad m_gauss;
	SolverType m_type = SolverType::bVector;

	

};

#endif // !SOLVER_H

