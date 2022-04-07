#include "handleselectedobjects.h"

HandleSO::HandleSO()
{
}

HandleSO::HandleSO(Model* model, QTableWidget* logPoints, QTableWidget* logEdges, QTableWidget* logFaces)
{
	m_model = model;
	m_logPoints = logPoints;
	m_logEdges = logEdges;
	m_logFaces = logFaces;
}

void HandleSO::selectObjects(QVector2D point, QMatrix4x4* transform)
{
	switch (m_selectionType)
	{
	case SelectionType::SELECTION:
		selectVertices(point, transform);
		if (!m_isPointSelected)
			selectEdges(point, transform);
		if (!m_isPointSelected && !m_isEdgeSelected)
			selectFaces(point, transform);
		break;
	case SelectionType::SELECTIONBC:
		selectVertices(point, transform);
		if (!m_isPointSelected && !m_isEdgeSelected)
			selectFaces(point, transform);
		break;
	}



}

void HandleSO::selectVertices(QVector2D point, QMatrix4x4* transform)
{
	m_isPointSelected = false;
	
	qreal minDepth = 1.0;
	int idSelected = -1;

	for (int i = 0; i < m_model->m_vertexes.size(); i++)
	{
		Vertex* v = m_model->m_vertexes[i];
		if (v->isScreenPointClose(point, transform))
		{
			if (v->getDepth() < minDepth)
			{
				minDepth = v->getDepth();
				idSelected = i;
			}
		}
	}

	if (!(idSelected == -1))
	{
		Vertex* selectedVertex = m_model->m_vertexes[idSelected];
		m_isPointSelected = true;
		selectedVertex->setIsSelected(true);
		// se ele já foi selecionado nao precisa colocar na pilha
		if (!(selectedVertex->getColor() == selectedVertex->getSelectedColor()))
		{
			m_selectedPoints.push_back(selectedVertex);
		}
	}
}

void HandleSO::selectEdges(QVector2D point, QMatrix4x4* transform)
{ 
	m_isEdgeSelected = false;
	int idSelected = -1;
	qreal minDepth = 1.0;

	// Edges
	for (int i = 0; i < m_model->m_edges.size(); i++)
		{
			Edge* v = m_model->m_edges[i];
			if (v->isScreenPointClose(point, transform))
			{
				if (v->getDepth() < minDepth)
				{
					minDepth = v->getDepth();
					idSelected = i;
				}
			}
		}
		if (!(idSelected == -1))
		{
			Edge* selectedEdge = m_model->m_edges[idSelected];
			m_isEdgeSelected = true;
			selectedEdge->setIsSelected(true);
			// se nao for selecionado , nao insira
			if (!(selectedEdge->getColorEdge() == selectedEdge->getSelectedColorEdge()))
			{
				m_selectedEdges.push_back(selectedEdge);

			}

		}

}

void HandleSO::selectFaces(QVector2D point, QMatrix4x4* transform)
{
	m_isFaceSelected = false;
	int idSelected = -1;
	qreal minDepth = 1.0;

	//Faces

	for (int i = 0; i < m_model->m_faces.size(); i++)
	{
		Face* v = m_model->m_faces[i];
		if (v->isScreenPointClose(point, transform))
		{
			if (v->getDepth() < minDepth)
			{
				minDepth = v->getDepth();
				idSelected = i;
			}
		}
	}

	if (!(idSelected == -1))
	{
		Face* selectedFace= m_model->m_faces[idSelected];
		m_isFaceSelected = true;
		selectedFace->setIsSelected(true);
		if ( !(selectedFace->getCurrentColor() == selectedFace->getSelectedColor()) )
		{
			m_selectedFaces.push_back(selectedFace);
		}
		

		
	}
}

std::vector<Face*> HandleSO::getSelectedFaces()
{
	return m_selectedFaces;
}

void HandleSO::reset()
{
	m_isPointSelected = false;
	m_isEdgeSelected = false;
	m_isFaceSelected = false;

	for (Vertex* v : m_model->m_vertexes)
	{
		v->setColor(v->getNormalColor()); //red
		v->setIsSelected(false);
	}

	for (Edge* v : m_model->m_edges)
	{	
		v->setColorEdge(v->getNormalColorEdge()); //green
		v->setColorCone(v->getNormalColorCone()); //yellow
		v->setIsSelected(false);
	}

	for (Face* v : m_model->m_faces)
	{
		v->setColor(v->getNormalColor());
		v->setIsSelected(false);
	}

	// reset stacks

	m_selectedPoints.clear();
	m_selectedEdges.clear();
	m_selectedFaces.clear();

	// reset tables
	m_logPoints->setRowCount(0);
	m_logPointsInit = false;
	m_logEdges->setRowCount(0);
	m_logEdgesInit = false;
	m_logFaces->setRowCount(0);
	m_logFacesInit = false;
}

void HandleSO::changeColorSO()
{
	for (Vertex* v : m_selectedPoints)
	{
		v->setColor(v->getSelectedColor()); //magenta
	}
	for (Edge* v : m_selectedEdges)
	{
		v->setColorEdge(v->getSelectedColorEdge()); //magenta
		v->setColorCone(v->getSelectedColorCone()); //magenta
	}
	for (Face* v : m_selectedFaces)
	{
		v->setColor(v->getSelectedColor());
	}

	updateTable();
}

bool HandleSO::isAnySelected()
{
	switch (m_selectionType)
	{
	case SelectionType::SELECTION:
	{
		if ((m_isPointSelected == false) && (m_isEdgeSelected == false) && (m_isFaceSelected == false))
		{
			return false;
		}
		return true;
		break;
	}
	case SelectionType::SELECTIONBC:
	{
		if ((m_isPointSelected == false) && (m_isFaceSelected == false))
		{
			return false;
		}
		return true;
		break;
	}
	}

}

void HandleSO::updateTable()
{
	QString str;
	m_logPoints->setRowCount(0);
	m_logEdges->setRowCount(0);
	m_logFaces->setRowCount(0);
	switch (m_selectionType)
	{
	case SelectionType::SELECTION:
		str = "P";
		for (Vertex* v : m_selectedPoints)
		{
			QString idStr = str + QString::number(v->getId());
			m_logPoints->setRowCount(m_logPoints->rowCount() + 1); //new line
			QTableWidgetItem* idItem = new QTableWidgetItem(idStr);
			QTableWidgetItem* xItem = new QTableWidgetItem(QString::number(v->X()));
			QTableWidgetItem* yItem = new QTableWidgetItem(QString::number(v->Y()));
			QTableWidgetItem* zItem = new QTableWidgetItem(QString::number(v->Z()));
			m_logPoints->setItem(m_logPoints->rowCount() - 1, 0, idItem);
			m_logPoints->setItem(m_logPoints->rowCount() - 1, 1, xItem);
			m_logPoints->setItem(m_logPoints->rowCount() - 1, 2, yItem);
			m_logPoints->setItem(m_logPoints->rowCount() - 1, 3, zItem);

		}
		str = "A";
		for (Edge* e : m_selectedEdges)
		{
			int p1 = e->getIdP1();
			int p2 = e->getIdP2();

			QString idStr = str + QString::number(e->getId());
			m_logEdges->setRowCount(m_logEdges->rowCount() + 1); //new line
			QTableWidgetItem* idItem = new QTableWidgetItem(idStr);
			QTableWidgetItem* p1Id = new QTableWidgetItem(QString::number(p1));
			QTableWidgetItem* p2Id = new QTableWidgetItem(QString::number(p2));

			m_logEdges->setItem(m_logEdges->rowCount() - 1, 0, idItem);
			m_logEdges->setItem(m_logEdges->rowCount() - 1, 1, p1Id);
			m_logEdges->setItem(m_logEdges->rowCount() - 1, 2, p2Id);


		}
		str = "F";
		for (Face* f : m_selectedFaces)
		{
			int e1Id = f->getE1()->getId();
			int e2Id = f->getE2()->getId();
			int e3Id = f->getE3()->getId();

			QString idStr = str + QString::number(f->getId());
			m_logFaces->setRowCount(m_logFaces->rowCount() + 1); //new line
			QTableWidgetItem* idItem = new QTableWidgetItem(idStr);
			QTableWidgetItem* e1IdItem = new QTableWidgetItem(QString::number(e1Id));
			QTableWidgetItem* e2IdItem = new QTableWidgetItem(QString::number(e2Id));
			QTableWidgetItem* e3IdItem = new QTableWidgetItem(QString::number(e3Id));

			m_logFaces->setItem(m_logFaces->rowCount() - 1, 0, idItem);
			m_logFaces->setItem(m_logFaces->rowCount() - 1, 1, e1IdItem);
			m_logFaces->setItem(m_logFaces->rowCount() - 1, 2, e2IdItem);
			m_logFaces->setItem(m_logFaces->rowCount() - 1, 3, e3IdItem);
		}
		break;
	case SelectionType::SELECTIONBC:
		str = "P";
		for (Vertex* v : m_selectedPoints)
		{
			QString idStr = str + QString::number(v->getId());
			QString temperature;
			bool bc = v->getBC();
			if (bc)
			{
				temperature = QString::number(v->getTemperature());
			}
			m_logPoints->setRowCount(m_logPoints->rowCount() + 1); //new line
			QTableWidgetItem* idItem = new QTableWidgetItem(idStr);
			QTableWidgetItem* BC = new QTableWidgetItem(QString::number(bc));
			QTableWidgetItem* value = new QTableWidgetItem(temperature);
			m_logPoints->setItem(m_logPoints->rowCount() - 1, 0, idItem);
			m_logPoints->setItem(m_logPoints->rowCount() - 1, 1, BC);
			m_logPoints->setItem(m_logPoints->rowCount() - 1, 2, value);

		}
		str = "F";
		for (Face* f : m_selectedFaces)
		{
			QString gradient;
			bool bc = f->getBC();
			if (bc)
			{
				gradient = QString::number(f->getGradient());
			}

			int e1Id = f->getE1()->getId();
			int e2Id = f->getE2()->getId();
			int e3Id = f->getE3()->getId();

			QString idStr = str + QString::number(f->getId());
			m_logFaces->setRowCount(m_logFaces->rowCount() + 1); //new line
			QTableWidgetItem* idItem = new QTableWidgetItem(idStr);
			QTableWidgetItem* BC = new QTableWidgetItem(QString::number(bc));
			QTableWidgetItem* value = new QTableWidgetItem(gradient);

			m_logFaces->setItem(m_logFaces->rowCount() - 1, 0, idItem);
			m_logFaces->setItem(m_logFaces->rowCount() - 1, 1, BC);
			m_logFaces->setItem(m_logFaces->rowCount() - 1, 2, value);

		}
		break;

	}



}

void HandleSO::deleteSelectedObjects()
{
	m_model->deleteSO();
	m_model->updateIds();
	m_selectedFaces.clear();
	m_selectedPoints.clear();
	m_selectedEdges.clear();	
	updateTable();
}
