/*!
*      \file GaussCurvature.h
*      \brief Algorithm for Gauss Curvature
*	   \author David Gu
*      \date Document 02/13/2014
*
*/


#ifndef _GAUSS_CURVATURE_H_
#define _GAUSS_CURVATURE_H_

#include <vector>
#include "Mesh/iterators.h"
#include "GaussCurvatureMesh.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif
namespace MeshLib
{

	template<typename M>
	class CGaussCurvature
	{
	public:
		/*!	CGaussCurvature constructor
		 *	\param pMesh the input mesh
		 */
		CGaussCurvature( M* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CGaussCurvature();
		/*!  Compute the harmonic map using direct method
		 */
		void _calculate_curvature();
		/*!	Compute face normal
		 */
		void _calculate_face_normal();
		/*!	Compute vertex normal
		 */
		void _calculate_vertex_normal();
		/*! Compute Euler number
		 */
		void _calculate_Euler_characteristics();

	protected:
		/*!	The input surface mesh
		 */
		M* m_pMesh;
		/*	boundary
		 */
		typename M::CBoundary m_boundary;
	};

double _cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
template<typename M>
CGaussCurvature<M>::CGaussCurvature( M* pMesh ): m_pMesh( pMesh ), m_boundary( pMesh )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->k() = 0;
	}
};

/*!
 *	CGaussCurvature destructor
 */
template<typename M>
CGaussCurvature<M>::~CGaussCurvature()
{
};


//Set the boundary vertices to the unit circle
/*!
 *	Fix the boundary using arc length parameter
 */
template<typename M>
void CGaussCurvature<M>::_calculate_curvature()
{
	//insert your code here
	//calculate edge length
	//calculate corner angle
	//caculate curvature

	double sum = 0;
	std::cout << "Total Curvature is " <<  sum/PI << " PI" << std::endl;;
}

/*!
 *	Compute face normal
 */
template<typename M>
void CGaussCurvature<M>::_calculate_face_normal()
{
	//insert your code here
}

/*!
 *	Compute face normal
 */
template<typename M>
void CGaussCurvature<M>::_calculate_vertex_normal()
{
	//insert your code here
}

/*!
 *  Calculate Euler characteristics number
 */
template<typename M>
void CGaussCurvature<M>::_calculate_Euler_characteristics()
{
	//insert your code in this function
	int V, E, F;
	std::cout << "Vertices: " << V << " Faces: " << F << " Edges: " << E << std::endl;
	int euler = V + F - E;
	std::cout << "Euler Characteristic Number " << euler << std::endl;	
	int b;
	int g;
	std::cout << "Number of boundaries " << b << std::endl;
	std::cout << "Genus " << g << std::endl;
}
};
#endif

