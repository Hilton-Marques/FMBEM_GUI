#include <cmath>
#include "Solver.h"
#include "../math/GaussQuadrature.h"
#include "../math/RTable.h"
#include <omp.h>
#include <iostream>

using namespace geomUtils;

Solver::Solver(const int &N, const int &nGFMM, const int &nNodes): m_N(N), m_nGFMM(nGFMM), m_nNodes(nNodes)
{
    m_NTot = (N + 1) * (N + 1);
    m_RTable.setSizeRTable(N);
    m_RTable.setSizeRecTableR(N);
    m_RTable.setSizeRecTableS(N);
    m_gauss.init(m_nGBEM);
}

Solver::~Solver()
{
}

void Solver::CalculateME(topoSolver::Face* element)
{
    element->setFMMSizes(m_NTot);
    GaussQuad gauss(m_nGFMM);
    Point yc = element->m_fatherElement->getYc();
    Point* normal = &element->m_normal;
    const std::complex<double> CNormal(normal->m_x, -normal->m_y);
    const double ZNormal = normal->m_z;
    double* jacobian = & element->m_jacobian;


    for (int i = 0; i < 3; i++)
    {
        topoSolver::Vertex* field = element->m_points[i];
        if (field->m_bd)
        {
            const double u = field->m_u; // temperature

            for (int j = 0; j < gauss.m_totN; j++)
            {
                double* Xi = &gauss.m_gaussPointsXi[j];
                double* Eta = &gauss.m_gaussPointsEta[j];
                double* W = &gauss.m_gaussWeights[j];
                double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
                Point xksi = Multi(shapeF, element->m_points);
                Point r = xksi - yc;
                //RTable Ry(m_N, r);
                double Const2 = u * shapeF[i] * (*W);
                int zdRpos;
                std::complex<double> adR, bdR, zdR;
                std::complex<double> xdR, ydR;
                std::complex<double> I(0, 1);
                int count = 0;  // carry the index of array Ry
                std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);
               
               // m_RTable.getPosCdR(-1, -1, R);
                for (int n = 0; n < m_N + 1; n++)
                {

                  int posMean = n * n + n; // this is the position of middle elements on pascal triangle

                  zdRpos = m_RTable.getPosZdR(n, 0);
                  zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

                  bdR = m_RTable.getPosCdR(n, 0, R);
                  adR = m_RTable.getProjdR(n, 0, R);

                  xdR = -0.5 * (adR - bdR);
                  ydR = -0.5 * I*(adR + bdR);

                  element->m_MEH[posMean] +=  + Const2 * (xdR*normal->m_x + ydR*normal->m_y + zdR*normal->m_z);

                   

                  for (int m = 1; m < n + 1; m++)
                  {                                       
                    int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
                    int posMNeg = posMean - m;                  
                    
                    zdRpos = m_RTable.getPosZdR(n, m);

                    short int minus1 = (m % 2) == 0 ? 1 : -1;

                    zdR = zdRpos != -1 ? minus1 *1.0* std::conj(R[zdRpos]) : 0.0;

                    
                    bdR = m_RTable.getPosCdR(n, -m , R);
                    
                    adR = m_RTable.getProjdR(n, -m, R);


                    xdR = -0.5 * (adR - bdR);
                    ydR = -0.5 * I*(adR + bdR);

                    element->m_MEH[posMNeg] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);


                    zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

                    bdR = m_RTable.getPosCdR(n, m , R);
                    adR = m_RTable.getProjdR(n, m, R);


                    xdR = -0.5 * (adR - bdR);
                    ydR = -0.5 * I *(adR + bdR);

                    element->m_MEH[posMPos] +=  Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);



                    count++;
                  }
                }

        
            } 


        }
        if (element->m_bd)
        {       
            const double q = element->m_q; // gradient
            const double Const = (*jacobian) * q;

                for (int j = 0; j < gauss.m_totN; j++)
                {
                    double* Xi = &gauss.m_gaussPointsXi[j];
                    double* Eta = &gauss.m_gaussPointsEta[j];
                    double* W = &gauss.m_gaussWeights[j];
                    double shapeF[3] = {  1. - *Xi - *Eta, *Xi, *Eta };
                    Point xksi = Multi(shapeF, element->m_points);
                    const Point r = xksi - yc;
                    double Const2 = *W * Const * shapeF[i];
                    int count = 0 ;  // carry the index of array Ry
                    std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);
                    for (int n = 0; n < m_N + 1; n++) 
                    {
                        int posMean = n * n + n; // this is the position of middle elements on pascal triangle

                        element->m_MEG[posMean] = element->m_MEG[posMean] + Const2 *( R [count]);
                        count++;
                        for (int m = 1; m < n + 1; m++)
                        {
                           int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
                           int posMNeg = posMean - m;
                           element->m_MEG[posMPos] = element->m_MEG[posMPos] + Const2 * (R[count]);
                           
                           short int minus1 = (m % 2) == 0 ? 1 : -1; 
                           element->m_MEG[posMNeg] = element->m_MEG[posMNeg] + Const2 * (minus1*1.* std::conj((R[count])) );
                           count++;
                        }
                    }

                }
            

            
        }
    }



}

void Solver::CalculateNearInt(std::vector<topoSolver::Vertex*> sourcesNodes, topoSolver::Face* element)
{
    
    Point* normal = &element->m_normal;
    double* jacobian =  &element->m_jacobian ;
    
    
    for (int i = 0; i < 3; i++)
    {
      topoSolver::Vertex* pt = element->m_points[i];

      if (pt->m_bd)
      {
        const double u = pt->m_u;
        const double Const2 = Const * u;

        for (topoSolver::Vertex* source : sourcesNodes)
        {

          if (source->m_id == pt->m_id)
          {
            if (!source->markTemp)
            
            {
              source->markTemp = true;
              #pragma omp atomic
              m_Hd[{source->m_id, 0}] +=  Const2* pt->calculateSolidAngle();
             
            }
            // Calculate diagonal of H Matrix
          }
          else
          {
            // Integrate Hu
            double H = 0.0;
            for (int j = 0; j < m_gauss.m_totN; j++)
            {
              double* Xi = &m_gauss.m_gaussPointsXi[j];
              double* Eta = &m_gauss.m_gaussPointsEta[j];
              double* W = &m_gauss.m_gaussWeights[j];
              double shapeF[3] = { 1. - *Xi - *Eta, *Xi, *Eta, };
              Point xksi = Multi(shapeF, element->m_points);
              Point r = xksi - source->m_coord;
              double rNorm = r.getNorm();
              double rNorm3 = rNorm * rNorm * rNorm;
              H = H + (- (r * (*normal) ) / rNorm3 ) * (*W) * shapeF[i];
            }

            #pragma omp atomic
            m_Hd[{source->m_id, 0}] +=  Const2 * H;

          }
        }

      }
      if (element->m_bd)
      {
        const double q = element->m_q;
        const double Const2 = Const * q * (*jacobian);
        for (topoSolver::Vertex* source : sourcesNodes)
        {          

          // Integrate Gq
          double G = 0.0;
          for (int j = 0; j < m_gauss.m_totN; j++)
          {
            double* Xi = &m_gauss.m_gaussPointsXi[j];
            double* Eta = &m_gauss.m_gaussPointsEta[j];
            double* W = &m_gauss.m_gaussWeights[j];
            double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
            Point xksi = Multi(shapeF, element->m_points);
            Point r = xksi - source->m_coord;
            double rNorm = r.getNorm();
            G = G + (1 / rNorm) * (*W) * shapeF[i];
          }
          #pragma omp atomic
          m_Gt[{source->m_id, 0}] += Const2 * G;
        }

      }
    }


}

void Solver::CalculateFarInt(std::vector<topoSolver::Vertex*> sourcesNodes, topoSolver::Face* element)
{
    Point yc = element->getYc();
    std::vector<std::complex<double>> Sb;
    for (topoSolver::Vertex* source : sourcesNodes)
    {

        Point x = source->m_coord - yc;
        Sb = m_RTable.evaluateRecursiveTableS(x);   
        #pragma omp atomic
        m_Gt[{source->m_id, 0}] +=  Const * Dot(element->m_MEG, &Sb);
        #pragma omp atomic
        m_Hd[{source->m_id, 0}] +=  Const * Dot(element->m_MEH, &Sb);

        //if (value != 0)
        //{
        //  if (source->m_id == 13)
        //  {
        //    Face* el = element->m_childrenElement[0];
        //    element->showMEH();
        //  }

        //}

      
    }
}

void Solver::CalculateMMT(topoSolver::Face* element)
{
    
    Point ycL = element->getYc();
    for (topoSolver::Face* child : element->m_childrenElement)
    {
        Point yc = child->getYc();
        Point r = yc - ycL;
        int id;
        std::vector<std::complex<double>> R = m_RTable.evaluateRecursiveTableR(r);
        int count = 0;
        for (int n = 0; n < m_N + 1; n++)
        {
            for (int m = -n; m < n+1; m++)
            {
                for (int nc = 0; nc < m_N + 1; nc++)
                {
                    for (int mc = -nc; mc < nc + 1; mc++)
                    {
                        id = getPos(n - nc, m - mc);
                        if ( id > -1)
                        {
                            element->m_MEG[count] = element->m_MEG[count] + R[getPos(nc, mc)] * child->m_MEG[id];
                            element->m_MEH[count] = element->m_MEH[count] + R[getPos(nc, mc)] * child->m_MEH[id];
                            
                        }
                    }
                }
                count++;
            }
        }
        child->reset();
    }
}



double Solver::Dot(std::complex<double> R[], std::vector<std::complex<double>>* Sb)
{
  std::complex <double> dot = 0.0;
  for (int i = 0; i < m_NTot; i++)
  {
    dot = dot + (R[i] * (*Sb)[i]);
  }
  return std::real(dot);
}

int Solver::getPos(const int& N, const int& M)
{

        const int mTot = 2 * N + 1;
        const int mVirtual = M + (mTot - 1) * 0.5 ;
        if (mVirtual >= mTot || mVirtual < 0)
        {
            return -1;
        }
        const int nTot = (N ) * (N );
        return (nTot + mVirtual);

 
}




std::complex<double> operator*(const std::vector<std::complex<double>>& dR,const Point& x)
{
    return std::complex<double>(dR[0]*x.m_x + dR[1]*x.m_y + dR[2]*x.m_z);
}

void Solver::prepareElementMatrix(std::vector<topoSolver::Vertex*> sourceNodes , topoSolver::Face* element)
{
  Point* normal = &element->m_normal;
  double* jacobian = &element->m_jacobian;
  Matrix GG(sourceNodes.size(), 3);
  Matrix HH(sourceNodes.size(), 3);

  for (int j = 0; j < 3; j++)
  {
    topoSolver::Vertex* pt = element->m_points[j];

    if (pt->m_bd)
    {
      const double Const2 = Const;

     for (size_t i = 0; i < sourceNodes.size(); i++)
       {
       topoSolver::Vertex* source = sourceNodes[i];
        if (source->m_id == pt->m_id)
        {
          if (!source->markTemp)
          {
            HH[{i,j}] = Const2 * pt->calculateSolidAngle();
            source->markTemp = true;
          }
          // Calculate diagonal of H Matrix
        }
        else
        {
          // Integrate Hu
          double H = 0.0;
          for (int k = 0; k < m_gauss.m_totN; k++)
          {
            double* Xi = &m_gauss.m_gaussPointsXi[k];
            double* Eta = &m_gauss.m_gaussPointsEta[k];
            double* W = &m_gauss.m_gaussWeights[k];
            double shapeF[3] = { 1. - *Xi - *Eta, *Xi, *Eta, };
            Point xksi = Multi(shapeF, element->m_points);
            Point r = xksi - source->m_coord;
            double rNorm = r.getNorm();
            double rNorm3 = rNorm * rNorm * rNorm;
            H = H + (-(r * (*normal)) / rNorm3) * (*W) * shapeF[j];
          }
          HH[{i, j}] = Const2 * H;
        }
      }
    }
    if (element->m_bd)
    {
      const double Const2 = Const * (*jacobian);
      for (size_t i = 0; i < sourceNodes.size(); i++)
      {
        topoSolver::Vertex* source = sourceNodes[i];
        // Integrate Gq
        double G = 0.0;
        for (int k = 0; k < m_gauss.m_totN; k++)
        {
          double* Xi = &m_gauss.m_gaussPointsXi[k];
          double* Eta = &m_gauss.m_gaussPointsEta[k];
          double* W = &m_gauss.m_gaussWeights[k];
          double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
          Point xksi = Multi(shapeF, element->m_points);
          Point r = xksi - source->m_coord;
          double rNorm = r.getNorm();
          G = G + (1 / rNorm) * (*W) * shapeF[j];
        }
        GG[{i, j}] =  Const2 * G;
      }

    }
  }
  element->m_H = HH;
  element->m_G = GG;

}
void Solver::CalculateNearIntMatrix(std::vector<topoSolver::Vertex*> sourceNodes , topoSolver::Face* element)
{
  Matrix d(3, 1);
  Matrix t(3, 1);
  for (int i = 0; i < 3; i++)
  {
    topoSolver::Vertex* pt = element->m_points[i];
    if (pt->m_bd)
    {
      d[{i, 0}] = pt->m_u;
    }
    if (element->m_bd)
    {
      t[{i, 0}] = element->m_q;
    }
  }
  for (int i = 0; i < sourceNodes.size(); i++)
  {
    topoSolver::Vertex* pt = sourceNodes[i];
    Matrix rowH = element->m_H.getSubMatrix(i + 1, i + 1, 1, 3);
    Matrix rowG = element->m_G.getSubMatrix(i+1, i+1, 1, 3);
    #pragma omp atomic
    m_Hd[{pt->m_id, 0}] +=  rowH.dot(d);
    #pragma omp atomic
    m_Gt[{pt->m_id, 0}] += rowG.dot(t);
  }

}
void Solver::CalculateMEForGMRES(topoSolver::Face* element)
{
  GaussQuad gauss(m_nGFMM);
  Point yc = element->m_fatherElement->getYc();
  Point* normal = &element->m_normal;
  const std::complex<double> CNormal(normal->m_x, -normal->m_y);
  const double ZNormal = normal->m_z;
  double* jacobian = &element->m_jacobian;


  for (int i = 0; i < 3; i++)
  {
    topoSolver::Vertex* field = element->m_points[i];
    std::complex<double>* MEH_i = element->m_MEH_p[i];
    std::complex<double>* MEG_i = element->m_MEG_p[i];

    if (field->m_bd)
    {
      const double u = field->m_u; // temperature

      for (int j = 0; j < gauss.m_totN; j++)
      {
        double* Xi = &gauss.m_gaussPointsXi[j];
        double* Eta = &gauss.m_gaussPointsEta[j];
        double* W = &gauss.m_gaussWeights[j];
        double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
        Point xksi = Multi(shapeF, element->m_points);
        Point r = xksi - yc;
        //RTable Ry(m_N, r);
        double Const2 = u * shapeF[i] * (*W);
        int zdRpos;
        std::complex<double> adR, bdR, zdR;
        std::complex<double> xdR, ydR;
        std::complex<double> I(0, 1);
        int count = 0;  // carry the index of array Ry
        std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);

        // m_RTable.getPosCdR(-1, -1, R);
        for (int n = 0; n < m_N + 1; n++)
        {

          int posMean = n * n + n; // this is the position of middle elements on pascal triangle

          zdRpos = m_RTable.getPosZdR(n, 0);
          zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

          bdR = m_RTable.getPosCdR(n, 0, R);
          adR = m_RTable.getProjdR(n, 0, R);

          xdR = -0.5 * (adR - bdR);
          ydR = -0.5 * I * (adR + bdR);

          MEH_i[posMean] += +Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);



          for (int m = 1; m < n + 1; m++)
          {
            int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
            int posMNeg = posMean - m;

            zdRpos = m_RTable.getPosZdR(n, m);

            short int minus1 = (m % 2) == 0 ? 1 : -1;

            zdR = zdRpos != -1 ? minus1 * 1.0 * std::conj(R[zdRpos]) : 0.0;


            bdR = m_RTable.getPosCdR(n, -m, R);

            adR = m_RTable.getProjdR(n, -m, R);


            xdR = -0.5 * (adR - bdR);
            ydR = -0.5 * I * (adR + bdR);

            MEH_i[posMNeg] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);


            zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

            bdR = m_RTable.getPosCdR(n, m, R);
            adR = m_RTable.getProjdR(n, m, R);


            xdR = -0.5 * (adR - bdR);
            ydR = -0.5 * I * (adR + bdR);

            MEH_i[posMPos] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);



            count++;
          }
        }


      }


    }
    if (element->m_bd)
    {
      const double q = element->m_q; // gradient
      const double Const = (*jacobian) * q;

      for (int j = 0; j < gauss.m_totN; j++)
      {
        double* Xi = &gauss.m_gaussPointsXi[j];
        double* Eta = &gauss.m_gaussPointsEta[j];
        double* W = &gauss.m_gaussWeights[j];
        double shapeF[3] = { 1. - *Xi - *Eta, *Xi, *Eta };
        Point xksi = Multi(shapeF, element->m_points);
        const Point r = xksi - yc;
        double Const2 = *W * Const * shapeF[i];
        int count = 0;  // carry the index of array Ry
        std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);
        for (int n = 0; n < m_N + 1; n++)
        {
          int posMean = n * n + n; // this is the position of middle elements on pascal triangle

          MEG_i[posMean] = element->m_MEG[posMean] + Const2 * (R[count]);
          count++;
          for (int m = 1; m < n + 1; m++)
          {
            int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
            int posMNeg = posMean - m;
            MEG_i[posMPos] = element->m_MEG[posMPos] + Const2 * (R[count]);

            short int minus1 = (m % 2) == 0 ? 1 : -1;
            MEG_i[posMNeg] = element->m_MEG[posMNeg] + Const2 * (minus1 * 1. * std::conj((R[count])));
            count++;
          }
        }

      }



    }
  }



}
