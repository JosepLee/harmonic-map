/*!
*      \file GaussBonnetCurvatureMesh.h
*      \brief Mesh for total GaussBonnetian curvature
*	   \author David Gu
*      \date Documented 02/13/2014
*
*/

#ifndef  _GAUSS_BONNET_MESH_H_
#define  _GAUSS_BONNET_MESH_H_

#include <map>
#include <vector>
/*
#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
*/

#include "BaseMesh.h"
#include "Vertex.h"
#include "HalfEdge.h"
#include "Edge.h"
#include "Face.h"
#include "iterators.h"
#include "boundary.h"
#include "parser.h"
#define ON_FIRE 1
#define IN_FIRE 2
#define OUT_FIRE 3
namespace MeshLib
{
	/*!
	*	\brief CGaussBonnetVertex class
	*
	*	Vertex class for GaussBonnet Theorem
	*/
	class CGaussBonnetVertex : public CVertex
	{

	public:
		/*!
		*	CHarmonicVertex constructor
		*/
		CGaussBonnetVertex() { m_k = 0; };
		/*!
		*	CHarmonicVertex destructor
		*/
		~CGaussBonnetVertex() {};
		/*
		*	GaussBonnet curvature
		*/
		double & k() { return m_k; };
		/*
		*	Vertex normal
		*/
		CPoint & normal() { return m_normal; };

	protected:
		/*
		*	GaussBonnet curvature
		*/
		double m_k;
		/*!
		*	Vertex normal
		*/
		CPoint m_normal;
	};

	/*!
	*	\brief CGaussBonnetEdge class
	*
	*	Edge class for GaussBonnet curvature
	*/
	class CGaussBonnetEdge : public  CEdge
	{
	public:
		/*!	CGaussBonnetEdge constructor
		*/
		CGaussBonnetEdge() { m_length = 0; };
		/*!	CGaussBonnetEdge destructor
		*/
		~CGaussBonnetEdge(){};
		/*! edge length trait
		*/
		double & length() { return m_length; };

	protected:
		/*! edge length trait */
		double   m_length;
	};



	/*!
	*	\brief CGaussBonnetEdge class
	*
	*	Face class for GaussBonnet curvature
	*/
	class CGaussBonnetFace : public  CFace
	{
	public:
		/*!	CGaussBonnetFace constructor
		*/
		CGaussBonnetFace() {};
		/*!	CGaussBonnetFace destructor
		*/
		~CGaussBonnetFace(){};
		/*!
		*	face area
		*/
		double & area() { return m_area; };
		/*!
		*	face normal
		*/
		CPoint & normal(){ return m_normal; };

	protected:
		double m_area;
		CPoint m_normal;
	};


	/*!
	*	\brief CGaussBonnetHalfEdge class
	*
	*	HalfEdge class for GaussBonnet curvature
	*/
	class CGaussBonnetHalfEdge : public  CHalfEdge
	{
	public:
		/*!	CHarmonicHalfEdge constructor
		*/
		CGaussBonnetHalfEdge() {};
		/*!	CHarmonicHalfEdge destructor
		*/
		~CGaussBonnetHalfEdge(){};
		/*!	Corner angle trait
		*/
		double & angle() { return m_angle; };

	protected:
		/*! Corner angle trait */
		double m_angle;
	};

	/*-------------------------------------------------------------------------------------------------------------------------------------

	GaussBonnet Mesh Class

	--------------------------------------------------------------------------------------------------------------------------------------*/
	/*!
	*	\brief CGaussBonnetMesh class
	*
	*	Mesh class for GaussBonnet curvature
	*/
	template<typename V, typename E, typename F, typename H>
	class CGaussBonnetMesh : public CBaseMesh<V, E, F, H>
	{
	public:

		typedef V CVertex;
		typedef E CEdge;
		typedef F CFace;
		typedef H CHalfEdge;

		typedef CGaussBonnetMesh<V, E, F, H> M;

		typedef CBoundary<M> CBoundary;
		typedef CLoop<M> CLoop;

		typedef MeshVertexIterator<M> MeshVertexIterator;
		typedef MeshEdgeIterator<M> MeshEdgeIterator;
		typedef VertexVertexIterator<M> VertexVertexIterator;
		typedef VertexEdgeIterator<M> VertexEdgeIterator;
		typedef MeshFaceIterator<M> MeshFaceIterator;
		typedef FaceVertexIterator<M> FaceVertexIterator;
		typedef VertexFaceIterator<M> VertexFaceIterator;
		typedef FaceHalfedgeIterator<M> FaceHalfedgeIterator;
		typedef VertexInHalfedgeIterator<M> VertexInHalfedgeIterator;
		typedef FaceEdgeIterator<M> FaceEdgeIterator;

		int Total_C;
	};

	/*! Mesh class for CHarmonicMapper class, Abbreviated as 'CGCMesh'
	*/
	typedef CGaussBonnetMesh<CGaussBonnetVertex, CGaussBonnetEdge, CGaussBonnetFace, CGaussBonnetHalfEdge> CGBMesh;

};
#endif  _GAUSS_BONNET_MESH_H_