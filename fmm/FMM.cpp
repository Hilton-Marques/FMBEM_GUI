#include<iostream>
#include "Solver.h"
#include "FMM.h"
#include "../math/GaussQuadrature.h"
#include <algorithm> // set_difference
#include <set> // erase duplicates

using namespace geomUtils;
using namespace topoSolver;

FMM::FMM()
{
}
FMM::FMM(topoSolver::Solid* solid, int N, int NG) {


	m_nLMax = solid->m_nL;
	//m_elementsLevel = & solid->m_elementsLevel;
	m_nTotalVerts = solid->m_nPts;
	m_N = N;
	m_solver = new Solver(N, NG, m_nTotalVerts);
	m_solid = solid;
	m_nDb = m_solid->m_vertsBd.size(); // temperature prescribed
	m_Hdb = Matrix(m_nTotalVerts, m_nDb); 

	std::vector<Face*> elements;
	elements = solid->m_elementsLevel[m_nLMax - 1];
	////Initialize  leaf elements
	m_tot = (N + 1) * (N + 1);
	for (topoSolver::Face* element : elements) {
		element->m_points[0]->m_elements.push_back(element);
		element->m_points[1]->m_elements.push_back(element);
		element->m_points[2]->m_elements.push_back(element);
		element->m_adjacentFaces = element->getAdjacentElements();

		element->setNormal();

		element->m_isLeaf = true;
	}
	//checkWithConventionalBem(solid, elements,m_b);
	//m_solver->CalculateME(elements[0]);
	//elements[0]->showMEH();
	//m_solver->CalculateFarInt(std::vector<Vertex*> {m_solid->m_verts[3] }, elements[0]);
	m_b = computeBvector();
	prepareForGMRES();
	GMRES gmres(this);
	m_x = gmres.Solver(m_b);
	m_resid = gmres.m_resid;
	//Matrix A(4, 4, { 8, 7, 6, 1, 6, 7, 2, 5, 5, 9, 3, 2, 10, 9, 9, 10 });
	//Matrix b(4, 1, { 1,0,1,0 });
	//Matrix res = gmres.BiCGM(A,b);
	//for (int i = 0; i < 4; i++)
	//{
	//	std::cout << res[{i, 0}] << "\n";
	//}




		//checkWithConventionalBemH(solid, solid->m_elementsLevel.back());
		//checkWithConventionalBem(solid, solid->m_elementsLevel.back(), m_b);
	}



std::vector<topoSolver::Vertex*> FMM::getLeafNodesFromElements(std::vector<topoSolver::Face*> adjacentElements)
{
	std::vector<topoSolver::Vertex*> sourceVerts;
	if (adjacentElements.empty())
	{
		return sourceVerts;
	}	//sourceVerts.reserve(adjacentElements.size());

	for (topoSolver::Face* element : adjacentElements)
	{
		// considere apenas elementos nao adjacentes ao filho (distantes)
		  std::vector<topoSolver::Face*> leafElements;
			element->getLeafElements(&leafElements);
			for (topoSolver::Face* leafElement : leafElements)
			{
				for (topoSolver::Vertex* pt : leafElement->m_points)
				{
					sourceVerts.push_back(pt);
				}
			}
	}
	// erase duplicates
	std::set<topoSolver::Vertex*> s(sourceVerts.begin(), sourceVerts.end());
	sourceVerts.assign(s.begin(), s.end());
	return sourceVerts;
}

std::vector<topoSolver::Face*> FMM::getFarElements(std::vector<topoSolver::Face*> elements, topoSolver::Face* element)
{
	std::vector<topoSolver::Face*> farElements;
	for (Face* adj : element->m_adjacentFaces)
	{
		adj->markTemp = true;
	}
	element->markTemp = true;
	for (topoSolver::Face* element : elements)
	{
		if (!element->markTemp)
		{
			farElements.push_back(element);
		}
	}
	for (topoSolver::Face* adj : element->m_adjacentFaces)
	{
		adj->markTemp = false;
	}
	element->markTemp = false;
	return farElements;
}

void FMM::checkWithConventionalBem(topoSolver::Solid* solid, std::vector<topoSolver::Face*> elements, Matrix &Ax)
{
	GaussQuad gauss(40);
	double* m_Ax;
	double* m_Bx;
	m_Ax = new double[solid->m_nPts];
	m_Bx = new double[solid->m_nPts];
	for (int i = 0; i < solid->m_nPts; i++)
	{
		m_Ax[i] = 0.0;
		m_Bx[i] = 0.0;
	}
	int count = 0;
#pragma omp parallel for
	for (int m = 0; m < elements.size(); m++)
	{
		Face* element = elements[m];
		//for (Face* element : elements) 

		Point* normal = &element->m_normal;
		double* jacobian = &element->m_jacobian;
		const double q = element->m_q;
		const double Const2 = Const * q * (*jacobian);
		for (int i = 0; i < 3; i++)
		{
			for (topoSolver::Vertex* source : solid->m_verts)
			{
				double G = 0.0;
				// Integrate Gq
				for (int j = 0; j < gauss.m_totN; j++)
				{
					double* Xi = &gauss.m_gaussPointsXi[j];
					double* Eta = &gauss.m_gaussPointsEta[j];
					double* W = &gauss.m_gaussWeights[j];
					double shapeF[3] = { 1. - *Xi - *Eta  , *Xi, *Eta, };

					Point xksi = (shapeF[0] * element->m_points[0]->m_coord) + (shapeF[1] * element->m_points[1]->m_coord) + (shapeF[2] * element->m_points[2]->m_coord);
					Point r = xksi - source->m_coord;
					double rNorm = r.getNorm();
					G = G + (1 / rNorm) * (*W) * shapeF[i];
				}

				m_Ax[source->m_id] = m_Ax[source->m_id] + Const2 * G;
				m_Bx[source->m_id] = m_Bx[source->m_id] + Const2 * G;
			}
			count++;
		}
	}
	std::cout.precision(15);
	for (int i = 0; i < solid->m_nPts; i++)
	{
		std::cout << abs(m_Ax[i] - m_Bx[i]) << "\n";
	}
	delete[] m_Ax;
}

void FMM::checkWithConventionalBemH(Solid* solid, std::vector<Face*> elements)
{
	GaussQuad gauss(40);
	double* m_Ax;
	m_Ax = new double[solid->m_nPts];
	for (int i = 0; i < solid->m_nPts; i++)
	{
		m_Ax[i] = 0.0;
	}
	int count = 0;
#pragma omp parallel for
	for (int m = 0; m < elements.size(); m++)
	{
		Face* element = elements[m];
		//for (Face* element : elements) 

		Point* normal = &element->m_normal;
		double* jacobian = &element->m_jacobian;
		const int u = 1;
		const double Const2 = Const * u;

		for (int i = 0; i < 3; i++)
		{
			topoSolver::Vertex* pt = element->m_points[i];

			for (topoSolver::Vertex* source : solid->m_verts)
			{
				if (source->m_id == pt->m_id)
				{
					if (!source->markTemp)
					{
						// Calculate diagonal of H Matrix
						m_Ax[source->m_id] = m_Ax[source->m_id] + Const * source->calculateSolidAngle();
						source->markTemp = true;
					}
				}
				else {				
				double H = 0.0;
				// Integrate Gq
				for (int j = 0; j < gauss.m_totN; j++)
				{
					double* Xi = &gauss.m_gaussPointsXi[j];
					double* Eta = &gauss.m_gaussPointsEta[j];
					double* W = &gauss.m_gaussWeights[j];
					double shapeF[3] = { 1. - *Xi - *Eta  , *Xi, *Eta, };

					Point xksi = (shapeF[0] * element->m_points[0]->m_coord) + (shapeF[1] * element->m_points[1]->m_coord) + (shapeF[2] * element->m_points[2]->m_coord);
					Point r = xksi - source->m_coord;
					double rNorm = r.getNorm();

					double rNorm3 = rNorm * rNorm * rNorm;
					H = H + (-(r * (*normal)) / rNorm3) * (*W) * shapeF[i];

					
				}
				m_Ax[source->m_id] = m_Ax[source->m_id] + Const2 * H;
				}

				
			}
			count++;
		}
	}
	std::cout.precision(17);
	for (int i = 0; i < solid->m_nPts; i++)
	{
		std::cout << std::fixed << m_Ax[i] << "\n";
	}
	delete[] m_Ax;
}

Matrix FMM::computeBvector()
{
	Matrix b = Matrix(m_nTotalVerts, 1);
	// Part 1
	std::vector<Face*> elements = m_solid->m_elementsLevel[m_nLMax - 2];

	//conventional BEM
#pragma omp parallel for
	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];
		std::vector<Face*> adjacentElements = element->getAdjacentElements();
		std::vector<topoSolver::Vertex*> sourceVertexs = getLeafNodesFromElements(adjacentElements);
		element->setFMMSizes(m_tot);
		for (Face* child : element->m_childrenElement)
		{
			m_solver->CalculateNearInt(sourceVertexs, child);

			// Integrate Expansion
			m_solver->CalculateME(child);
			// accumulate
			for (int k = 0; k < m_tot; k++)
			{
				element->m_MEG[k] +=  child->m_MEG[k];
				element->m_MEH[k] +=  child->m_MEH[k];
			}
			//child->m_collectVertexes.insert(child->m_collectVertexes.end(), sourceVertexs.begin(), sourceVertexs.end());
			child->reset();
		}

	}
	// Part2
	for (int i = m_nLMax - 3; i > 0; i--)
	{
		elements = m_solid->m_elementsLevel[i];

#pragma omp parallel for		
		for (int j = 0; j < elements.size(); j++)
		{
			Face* element = elements[j];
			
			
			//already evaluate FMM delivery downward pass
			std::vector<Face*> adjacentElementsMother = element->getAdjacentElements();
			std::vector<Face*> adjacentElementsChildren;
			adjacentElementsChildren.reserve(4 * adjacentElementsMother.size());
			for (Face* adj : adjacentElementsMother)
			{
				adjacentElementsChildren.insert(adjacentElementsChildren.end(), std::begin(adj->m_childrenElement), std::end(adj->m_childrenElement));
			}
			std::vector<topoSolver::Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsChildren);
			std::sort(leafPoints.begin(), leafPoints.end());


			for (Face* child : element->m_childrenElement)
			{
				std::vector<topoSolver::Vertex*> sourceVertexs;
				std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
				std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
				std::vector<topoSolver::Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
				std::sort(adjacentNodes.begin(), adjacentNodes.end());
				// got A \ B with poitns	
				std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));

				m_solver->CalculateFarInt(sourceVertexs, child);
			}
			element->setFMMSizes(m_tot);
			m_solver->CalculateMMT(element);
		}
	}
		//Part 3
	elements = m_solid->m_elementsLevel[0];
//#pragma omp parallel for 	
	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];

		//already evaluate FMM delivery downward pass
		std::vector<Face*> adjacentElementsMother = element->getAdjacentElements();
		std::vector<Face*> adjacentElementsChildren;
		adjacentElementsChildren.reserve(4 * adjacentElementsMother.size());
		for (Face* adj : adjacentElementsMother)
		{
			adjacentElementsChildren.insert(adjacentElementsChildren.end(), std::begin(adj->m_childrenElement), std::end(adj->m_childrenElement));
		}
		std::vector<topoSolver::Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsChildren);
		std::sort(leafPoints.begin(), leafPoints.end());

		auto farElements = getFarElements(elements, element);
		auto farNodes = getLeafNodesFromElements(farElements);

		for (Face* child : element->m_childrenElement)
		{
			std::vector<topoSolver::Vertex*> sourceVertexs;
			std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
			std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
			std::vector<topoSolver::Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
			std::sort(adjacentNodes.begin(), adjacentNodes.end());
			// got A \ B with poitns	
			std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));

			// accumulate nodes from far elements
			if (!farNodes.empty())
			{
				sourceVertexs.insert(sourceVertexs.end(), farNodes.begin(), farNodes.end());
			}
			// Far integral
			//m_solver->CalculateFarInt(sourceVertexs, child);			
			Point yc = child->getYc();
			std::vector<std::complex<double>> Sb;
			for (topoSolver::Vertex* source : sourceVertexs)
			{

				Point x = source->m_coord - yc;
				Sb = m_solver->m_RTable.evaluateRecursiveTableS(x);
				m_solver->m_Gt[{source->m_id, 0}] += Const * m_solver->Dot(child->m_MEG, &Sb);
				m_solver->m_Hd[{source->m_id, 0}] += Const * m_solver->Dot(child->m_MEH, &Sb);

			}
			child->reset();
		}
	}
	b = m_solver->m_Gt - m_solver->m_Hd;
	m_solver->reset();
	return b;
}

Matrix FMM::matrixVectorMulti(Matrix &x)
{

	//#pragma omp for
	for (int i = 0; i < m_nTotalVerts; i++)
	{
		topoSolver::Vertex* pt = m_solid->m_verts[i];
		if (pt->m_bd)
		{
			pt->m_u = x[{pt->m_id, 0}];
		}
	}

	// the same with elements
	std::vector<Face*> elements = m_solid->m_elementsLevel[m_nLMax - 2];
	//conventional BEM

//#pragma omp parallel for 	

	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];
		std::vector<topoSolver::Face*> adjacentElements = element->getAdjacentElements();
		std::vector<topoSolver::Vertex*> sourceVertexs = getLeafNodesFromElements(adjacentElements);
		for (topoSolver::Face* child : element->m_childrenElement)
		{
			m_solver->CalculateNearIntMatrix(sourceVertexs, child);
			// Integrate Expansion
			//m_solver->CalculateME(child);

			for (int l = 0; l < 3; l++)
			{
				if (child->m_points[l]->m_bd)
				{
					for (int k = 0; k < m_tot; k++)
					{
						element->m_MEH[k] += child->m_points[l]->m_u * child->m_MEH_p[l][k];
					}
				}
				if (child->m_bd) // have to change to take care of three types of gradiente
				{
					for (int k = 0; k < m_tot; k++)
					{
						element->m_MEG[k] += child->m_q * child->m_MEG_p[l][k];
					}
				}
			}

			//child->m_collectVertexes.insert(child->m_collectVertexes.end(), sourceVertexs.begin(), sourceVertexs.end());
			//child->reset();
		}

		// Accumulate moments

		//element->accumulate();

	}

	// Part2
	for (int i = m_nLMax - 3; i > 0; i--)
	{

		elements = m_solid->m_elementsLevel[i];
#pragma omp parallel for		
		for (int j = 0; j < elements.size(); j++)
		{
			Face* element = elements[j];

			//already evaluate FMM delivery downward pass
			std::vector<Face*> adjacentElementsMother = element->getAdjacentElements();
			std::vector<Face*> adjacentElementsChildren;
			adjacentElementsChildren.reserve(4 * adjacentElementsMother.size());
			for (Face* adj : adjacentElementsMother)
			{
				adjacentElementsChildren.insert(adjacentElementsChildren.end(), std::begin(adj->m_childrenElement), std::end(adj->m_childrenElement));
			}
			std::vector<topoSolver::Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsChildren);
			std::sort(leafPoints.begin(), leafPoints.end());


			for (Face* child : element->m_childrenElement)
			{
				std::vector<topoSolver::Vertex*> sourceVertexs;
				std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
				std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
				std::vector<topoSolver::Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
				std::sort(adjacentNodes.begin(), adjacentNodes.end());
				// got A \ B with poitns	
				std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));

				//far integral
				m_solver->CalculateFarInt(sourceVertexs, child);

			}
			m_solver->CalculateMMT(element);
		}
	}
	// Part3 ///////////////////
	elements = m_solid->m_elementsLevel[0];
//#pragma omp parallel for 	
	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];

		//already evaluate FMM delivery downward pass
		std::vector<Face*> adjacentElementsMother = element->m_adjacentFaces;
		std::vector<Face*> adjacentElementsChildren;
		adjacentElementsChildren.reserve(4 * adjacentElementsMother.size());
		for (Face* adj : adjacentElementsMother)
		{
			adjacentElementsChildren.insert(adjacentElementsChildren.end(), std::begin(adj->m_childrenElement), std::end(adj->m_childrenElement));
		}
		std::vector<topoSolver::Vertex*> leafPoints = getLeafNodesFromElements(adjacentElementsChildren);
		std::sort(leafPoints.begin(), leafPoints.end());

		auto farElements = getFarElements(elements, element);
		auto farNodes = getLeafNodesFromElements(farElements);

		for (Face* child : element->m_childrenElement)
		{
			std::vector<topoSolver::Vertex*> sourceVertexs;
			std::vector<Face*> adjacentElementsChild = child->getAdjacentElements();
			std::sort(adjacentElementsChild.begin(), adjacentElementsChild.end());
			std::vector<topoSolver::Vertex*> adjacentNodes = getLeafNodesFromElements(adjacentElementsChild);
			std::sort(adjacentNodes.begin(), adjacentNodes.end());
			// got A \ B with poitns	
			std::set_difference(leafPoints.begin(), leafPoints.end(), adjacentNodes.begin(), adjacentNodes.end(), std::inserter(sourceVertexs, sourceVertexs.begin()));

			// accumulate nodes from far elements
			if (!farNodes.empty())
			{
				sourceVertexs.insert(sourceVertexs.end(), farNodes.begin(), farNodes.end());
			}

			// Far integral
			//m_solver->CalculateFarInt(sourceVertexs, child );

			Point yc = child->getYc();
			std::vector<std::complex<double>> Sb;
			for (topoSolver::Vertex* source : sourceVertexs)
			{

				Point x = source->m_coord - yc;
				Sb = m_solver->m_RTable.evaluateRecursiveTableS(x);
				m_solver->m_Gt[{source->m_id, 0}] += Const * m_solver->Dot(child->m_MEG, &Sb);
				m_solver->m_Hd[{source->m_id, 0}] += Const * m_solver->Dot(child->m_MEH, &Sb);

			}
			child->reset();
		}
	}
	for (int i = 0; i < m_nDb; i++)
	{
		topoSolver::Vertex* pt = m_solid->m_vertsBd[i];
		m_solver->m_Hd[{pt->m_id, 0}] = x[{pt->m_id, 0}];
	}
	Matrix res = (m_solver->m_Hd - m_solver->m_Gt);
	m_solver->reset();
	return res;
}

void FMM::prepareForGMRES()
{
	m_solver->m_type = SolverType::GMRES;
	// inverte boundary conditions
	for (topoSolver::Vertex* pt : m_solid->m_verts)
	{
		if (pt->m_bd)
		{
			pt->m_bd = false;
		}
		else
		{
			pt->m_bd = true;
			pt->m_u = 1.0; // calculate ME temporary
		}
		pt->markTemp = false;
	}

	std::vector<Face*> elements = m_solid->m_elementsLevel[m_nLMax - 1];

	for (Face* element : elements)
	{
		//initialize vectors 
		element->set_FMM_ME_sizes();
		if (element->m_bd)
		{
			element->m_bd = false;
		}
		else
		{
			element->m_bd = true;
			element->m_q = 1.0; // calculate ME temporary
		}
	}

	elements = m_solid->m_elementsLevel[m_nLMax - 2];

//#pragma omp parallel for		
	for (int i = 0; i < elements.size(); i++)
	{
		Face* element = elements[i];
		std::vector<Face*> adjacentElements = element->getAdjacentElements();
		std::vector<topoSolver::Vertex*> sourceVertexs = getLeafNodesFromElements(adjacentElements);
		for (Face* child : element->m_childrenElement)
		{	
			m_solver->prepareElementMatrix(sourceVertexs, child);		
			m_solver->CalculateMEForGMRES(child);
		}
	}


}

void FMM::showB()
{
	std::cout.precision(15);
	for (int i = 0; i < m_nTotalVerts; i++)
	{
		std::cout << m_b[{i, 0}] << "\n";
	}
}

void FMM::showX()
{
	std::cout.precision(15);
	for (int i = 0; i < m_nTotalVerts; i++)
	{
		std::cout << m_x[{i, 0}] << "\n";
	}
}





FMM::~FMM()
{

}