/*!
*      \file GaussCurvature.h
*      \brief Algorithm for Gauss Curvature
*	   \author David Gu
*      \date Document 02/13/2014
*
*/


#ifndef _GAUSS_CURVATURE_H_
#define _GAUSS_CURVATURE_H_
#define ON_FIRE 1
#define IN_FIRE 2
#define OUT_FIRE 3
#include <vector>
//#include "Mesh/iterators.h"
#include "CGaussBonnetMesh.h"
#include "iterators.h"
#include<math.h>
#include<Windows.h>
#include <iostream>
#include "CGaussBonnetMesh.h"
#include <string>
#include<stack>
#include<time.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif
using namespace std;
namespace MeshLib
{
	using namespace std;
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
		void _calculate_K();
		void init_mapping();
		void mapping();
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
	//	CPoint m_vertex;
		typename M::CBoundary m_boundary;
	};

double _cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


/*!	constructor 
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
void CGaussCurvature<M>::_calculate_curvature()//计算曲率
{
	//insert your code here
	//calculate edge length


	//caculate curvature曲率
	double sum = 0;

	for (M::MeshVertexIterator mviter(m_pMesh); !mviter.end(); ++mviter)
	{
		M::CVertex* vertex = mviter.value();
		double temp_sum = 0;
		for (M::VertexInHalfedgeIterator vihiter(m_pMesh, vertex); !vihiter.end(); ++vihiter)
		{
			
			temp_sum += vihiter.value()->angle();
		}
		if (vertex->boundary())
		{
			vertex->k() = PI - temp_sum;
		}
		else
		{
			vertex->k() = 2 * PI - temp_sum;
		}
		sum += vertex->k();
	}

	std::cout << "Total Curvature is " << sum / PI << " PI" << std::endl;;



}
/*
计算每条边的权值
*/
template<typename M>
void CGaussCurvature<M>::_calculate_K()
{
	double sum = 0;
	double E = 0;
	for (M::MeshEdgeIterator me_iter(m_pMesh); !me_iter.end(); ++me_iter)
	{
		
		CEdge* M_Edge = me_iter.value();
		if (M_Edge->boundary()) {
		CPoint* Vi= &M_Edge->halfedge(0)->source()->point();
		CPoint* Vj = &M_Edge->halfedge(0)->vertex()->point();
		CPoint* Vk = &M_Edge->halfedge(0)->he_prev()->source()->point();
		//CPoint* Vl = &M_Edge->halfedge(1)->he_prev()->source()->point();
		double K = 0;
		K = (((*Vi) - (*Vk))*((*Vj) - (*Vk))) / (((*Vi) - (*Vk)) ^ ((*Vj) - (*Vk))).norm();
		M_Edge->m_k = K;
		M_Edge->length = ((*Vi) - (*Vj)).norm();
		cout <<"Boundary:"<<M_Edge->m_k<<"length"<< M_Edge->length <<endl;
		E += M_Edge->m_k*M_Edge->length*M_Edge->length;
		}

		if (!M_Edge->boundary()) {
			CPoint* Vi = &M_Edge->halfedge(0)->source()->point();
			CPoint* Vj = &M_Edge->halfedge(0)->vertex()->point();
			CPoint* Vk = &M_Edge->halfedge(0)->he_prev()->source()->point();
			CPoint* Vl = &M_Edge->halfedge(1)->he_prev()->source()->point();
			double K = 0;
			K = (((*Vi) - (*Vk))*((*Vj) - (*Vk))) / (((*Vi) - (*Vk)) ^ ((*Vj) - (*Vk))).norm()+ (((*Vi) - (*Vl))*((*Vj) - (*Vl))) / (((*Vi) - (*Vl)) ^ ((*Vj) - (*Vl))).norm();
			M_Edge->m_k = K;
			M_Edge->length = ((*Vi) - (*Vj)).norm();
			cout << "Not Boundary:" << M_Edge->m_k << "length" << M_Edge->length << endl;
			E += M_Edge->m_k*M_Edge->length*M_Edge->length;
		}

	}
	cout << "Energy:" << E << endl;

}
/*
初始映射
*/
template<typename M>
void CGaussCurvature<M>::init_mapping()//初始映射
{
	//calculate edge length
	double sum_length1 = 0;
	double new_length = 0;
	CVertex* begin_vertex;
	double theit = 0;
	double R = 0;
	for (M::MeshVertexIterator me_iter(m_pMesh); !me_iter.end(); ++me_iter)
	{
		
		
			if (me_iter.value()->boundary()) {
				begin_vertex = me_iter.value();
				break;
			}
	}

	for (M::MeshEdgeIterator me_iter(m_pMesh); !me_iter.end(); ++me_iter)
	{
		
		CEdge* M_Edge = me_iter.value();
		if (M_Edge->boundary()) {
			CPoint* Vi = &M_Edge->halfedge(0)->source()->point();
			CPoint* Vj = &M_Edge->halfedge(0)->vertex()->point();
			sum_length1 += M_Edge->length;
			//CPoint* Vl = &M_Edge->halfedge(1)->he_prev()->source()->point();
			cout << "Length:" << M_Edge->length << endl;
		}
		if (!M_Edge->boundary()) {
			CPoint* Vi = &M_Edge->halfedge(0)->source()->point();
			CPoint* Vj = &M_Edge->halfedge(0)->vertex()->point();
			Vi->changePoint(0, 0, 0);
			Vj->changePoint(0, 0, 0);

			cout << "Zero" << endl;
		}


	}

	cout << "Sum_Length:" << sum_length1 << endl;
	//init mapping
	R = sum_length1 / 2 / PI;
	begin_vertex->point() = CPoint(R, 0, 0);
	new_length += begin_vertex->most_clw_out_halfedge()->edge()->length;
	for (CVertex *cur_vert= begin_vertex->most_clw_out_halfedge()->target();cur_vert!=begin_vertex; cur_vert=cur_vert->most_clw_out_halfedge()->target())
	{
		
		if (cur_vert->boundary()) {
			CPoint* NewP = &cur_vert->point();
			new_length += cur_vert->most_clw_out_halfedge()->edge()->length;

			theit = new_length / sum_length1 * PI * 2;
			R = sum_length1 / 2 / PI;
			NewP->changePoint(cos(theit)*R, sin(theit)*R, 0);
			cout << "Sucsess:"<<"X:"<<(*NewP)(0)<<"Y:"<< (*NewP)(1)<< endl;
		}
		if (!cur_vert->boundary()) {
			cout << "something wrong" << endl;
		}

	}
	cout << new_length << "and" << sum_length1 << endl;

}
/*
迭代调和映照
*/
template<typename M>
void CGaussCurvature<M>::mapping()
{

	double k_sum=0;
	
	double Ener_value = 87838.2;
	int count = 0;
	while (1) {
			double EE = 0;
			for (M::MeshVertexIterator mviter(m_pMesh); !mviter.end(); ++mviter)
			{

				
				//计算新坐标
				if (!mviter.value()->boundary()) {
					double k_sum = 0;
					
					for (M::VertexInHalfedgeIterator vi_hiter(m_pMesh, mviter.value()); !vi_hiter.end(); ++vi_hiter)
					{
						M::CHalfEdge* M_Edge = vi_hiter.value();
						k_sum += M_Edge->edge()->m_k;
					}
					CPoint New_Pos=CPoint(0,0,0);
					for (M::VertexInHalfedgeIterator vi_hiter(m_pMesh, mviter.value()); !vi_hiter.end(); ++vi_hiter)
					{
						M::CHalfEdge* M_Edge = vi_hiter.value();
						double t_k = M_Edge->edge()->m_k;
						CPoint t_p = M_Edge->source()->point();
						New_Pos += t_p * t_k;
					}
					mviter.value()->point() = New_Pos / k_sum;
					//New_Pos /= k_sum;
					//cout << "Sucsess:" << "X:" << New_Pos(0) << "Y:" << New_Pos(1) << "Z:" << New_Pos(2) << endl;
					//vertex->point().changePoint(New_Pos(0), New_Pos(1), New_Pos(2));
					//cout << "all is well" << k_sum<<endl;
					//CPoint* Check = &vertex->point();
					//cout << "Sucsess:" << "X:" << (*Check)[0] << "Y:" << (*Check)[1] << "Z:" << (*Check)[2] <<endl;				
				}




			}
			//计算迭代之后的能量
			for (M::MeshEdgeIterator me_iter(m_pMesh); !me_iter.end(); ++me_iter)
			{
				CPoint c1, c2;
				c1 = me_iter.value()->halfedge(0)->source()->point();
				c2 = me_iter.value()->halfedge(0)->target()->point();
				me_iter.value()->length() = (c1 - c2).norm();
				EE += me_iter.value()->m_k * me_iter.value()->length() * me_iter.value()->length();
			}

			//Ener_value -= EE;
			cout << "Energy" << EE << endl;
			count++;
			cout << "迭代次数" << count << "能量差值" << fabs(Ener_value - EE) <<endl;
			if (fabs(Ener_value - EE) < 1e-7) {
				break;
			}

			Ener_value = EE;

			if (count == 5000) {
				//break;
			}
			
	}
	

}
//brain.m 第 52362 次迭代  E = 9.99157e-08	
//oldman.m 第 70985 次迭代 E = 9.99999e-08  

void iterate_globe(char *fileName, char *outName)
{

	time_t start, end;
	time(&start);
	CGBMesh mesh;
	mesh.read_m(fileName);
	CGaussCurvature<CGBMesh> GC(&mesh);
	double a = 0.01; //阿尔法
	//double e = 0.001;

	std::cout << mesh.numVertices() << " " << mesh.numEdges() << " " << mesh.numFaces() << std::endl;
	std::cout << mesh.vertices().size() << std::endl;

	//求K值
	int b_edge_sum = 0;
	for (MeshEdgeIterator<CGBMesh> edgeIter(&mesh); !edgeIter.end(); edgeIter++)
	{
		CHalfEdge *halfEdge_0 = edgeIter.value()->halfedge(0);
		CHalfEdge *halfEdge_1 = edgeIter.value()->halfedge(1);
		CPoint vi, vj, vk, vl;
		vi = halfEdge_0->source()->point();
		vj = halfEdge_0->target()->point();
		vk = halfEdge_0->he_next()->target()->point();
		vl = halfEdge_1->he_next()->target()->point();
		CPoint px = vi - vk;
		CPoint py = vj - vk;
		CPoint pm = vi - vl;
		CPoint pn = vj - vl;
		double k1 = (px * py) / (px^py).norm();
		double k2 = (pm * pn) / (pm^pn).norm();
		edgeIter.value()->m_k = k1 + k2;
		edgeIter.value()->m_k = k1 + k2;
		//cout << edgeIter.value()->m_k << std::endl;
	}

	//构造初始映射
	//遍历所有点，求向量x乘和
	for (MeshVertexIterator<CGBMesh> vertexIter(&mesh); !vertexIter.end(); vertexIter++)
	{
		CVertex *t_vert = vertexIter.value();
		VertexOutHalfedgeIterator<CGBMesh> halfIter(&mesh, vertexIter.value());
		CHalfEdge *h_curr, *h_pre = halfIter.value();
		for (halfIter++; !halfIter.end(); halfIter++)
		{
			h_curr = halfIter.value();
			CPoint vec_curr = h_curr->target()->point() - h_curr->source()->point();
			CPoint vec_pre = h_pre->target()->point() - h_pre->source()->point();
			t_vert->normal() = vec_pre ^ vec_curr;
			vertexIter.value()->normal() = vec_pre ^ vec_curr;
			//cout << t_vert->normal().norm() << endl;
			
			h_pre = h_curr;
		}
		//cout << t_vert->normal().norm() << endl;

	}
	//将所有点映射到单位球
	for (MeshVertexIterator<CGBMesh> vertexIter(&mesh); !vertexIter.end(); vertexIter++) {
		//cout << vertexIter.value()->normal().norm() << endl;
		vertexIter.value()->point() = vertexIter.value()->normal() / vertexIter.value()->normal().norm();
//标准化
	}

	mesh.write_m("START.m");
	//迭代
	double k1 = 0;
	for (MeshEdgeIterator<CGBMesh> edgeIter(&mesh); !edgeIter.end(); edgeIter++)
	{
		CPoint c1, c2;
		c1 = edgeIter.value()->halfedge(0)->source()->point();
		c2 = edgeIter.value()->halfedge(0)->target()->point();
		edgeIter.value()->length() = (c1 - c2).norm();
		//验证
	   // cout << edgeIter.value()->length()<< edgeIter.value()->m_k << endl;
		k1 += edgeIter.value()->m_k * edgeIter.value()->length() * edgeIter.value()->length();
		//cout << k1<<endl;
	}
	for (int i = 0; true; i++)
	{
		//遍历所有点，向切向量方向位移一小段
		for (MeshVertexIterator<CGBMesh> vertexIter(&mesh); !vertexIter.end(); vertexIter++)
		{
			CVertex *t_vert = vertexIter.value();
			CHalfEdge *h_curr;
			CPoint v_inter; //半边向量和
			for (VertexOutHalfedgeIterator<CGBMesh> halfIter(&mesh, vertexIter.value()); !halfIter.end(); halfIter++)
			{
				h_curr = halfIter.value();
				CPoint t_v = h_curr->target()->point() - h_curr->source()->point();
				v_inter += t_v;
			}
			CPoint v_tangent = v_inter - vertexIter.value()->normal() * (v_inter * vertexIter.value()->normal()); //切向量
			vertexIter.value()->point() += v_tangent * a;
			vertexIter.value()->point() /= vertexIter.value()->point().norm();  //规格化
		}
		//规格化法向量
		CPoint p_center; //质心
		for (MeshVertexIterator<CGBMesh> vertexIter(&mesh); !vertexIter.end(); vertexIter++)
			p_center += vertexIter.value()->point();
		p_center /= mesh.numVertices();
		for (MeshVertexIterator<CGBMesh> vertexIter(&mesh); !vertexIter.end(); vertexIter++)
		{
			CPoint v_norm = vertexIter.value()->point() - p_center;
			vertexIter.value()->normal() = v_norm / v_norm.norm();
			vertexIter.value()->point() = v_norm / v_norm.norm();
		}
		std::cout << "第 " << i << " 次迭代" << std::endl;
		//if (i == 1000)
		//{
		//	mesh.write_m("brain_1000.m");
		//}
		//if (i == 10000)
		//{
		//	mesh.write_m("brain_10000.m");
		//}
		//if(i == 15000)
		//{
		//	mesh.write_m("brain_15000.m");
		//}
		//if (i == 20000)
		//{
		//	mesh.write_m("brain_20000.m");
		//}
		//调和能量
		double k2 = 0;
		for (MeshEdgeIterator<CGBMesh> edgeIter(&mesh); !edgeIter.end(); edgeIter++)
		{
			CPoint c1, c2;
			c1 = edgeIter.value()->halfedge(0)->source()->point();
			c2 = edgeIter.value()->halfedge(0)->target()->point();
			edgeIter.value()->length() = (c1 - c2).norm();
			k2 += edgeIter.value()->m_k * edgeIter.value()->length() * edgeIter.value()->length();
		}
		//std::cout << "k1=" << k1 << std::endl;
		//std::cout << "k2=" << k2 << std::endl;
		cout << k1 << k2 << endl;
		std::cout << "E=" << fabs(k1 - k2) << std::endl;
		if (fabs(k1 - k2) < 1e-7)
			break;
		k1 = k2;
	}


	//输出
	mesh.write_m(outName);
	std::cout << "输出.m" << std::endl;
	time(&end);
	cout << difftime(end, start) << endl;
}


//同伦群基底
void homotopy(char *fileName, char *outName)
{
	time_t start, end;
	time(&start);
	CGBMesh mesh;
	mesh.read_m(fileName);

	std::cout << mesh.numVertices() << " " << mesh.numEdges() << " " << mesh.numFaces() << std::endl;

	std::vector<CFace*> f_face;
	//起始面
	MeshFaceIterator<CGBMesh> it(&mesh);
	it.value()->state() = ON_FIRE;
	//it.value()->string() = "sharp = true";
	f_face.push_back(it.value());
	//BFS
	for (int m = 0; true; m++)
	{
		//遍历所有火内的面
		int t_size = f_face.size();
		std::cout << t_size << endl;
		if (t_size >= mesh.numFaces())
		{
			break;
		}
		cout << mesh.numFaces() << endl;
		for (int i = 0; i < t_size; i++) //如果在火边缘
		{
			CFace *t_face = f_face[i];
			if (t_face->state() == ON_FIRE)
			{
				cout << "GOGO" << endl;
				t_face->state() = IN_FIRE; //改变face状态至火内
				//每个Face三个halfEdge， 一圈， 三个halfEdge与一个Face相对应
				for (FaceHalfedgeIterator<CGBMesh> it((MeshLib::CGBMesh::CFace*)t_face); !it.end(); it++)
				{
					CHalfEdge *t_h1 = it.value();
					CHalfEdge *t_h2 = t_h1->edge()->other(t_h1); //半边相对的另一条半边
					CFace *t_f = t_h2->face();
					if (t_f->state() == OUT_FIRE)
					{
						t_f->state() = ON_FIRE;
						f_face.push_back(t_f);
						//判断是否烧到了另一边
						for (FaceHalfedgeIterator<CGBMesh> itt((MeshLib::CGBMesh::CFace*)t_f); !itt.end(); itt++)
						{
							CHalfEdge *o_half = itt.value()->edge()->other(itt.value());
							if (o_half->face()->state() != OUT_FIRE && o_half->face() != t_face)
							{
								o_half->edge()->string() = "sharp = true";
							}
						}
					}
				}
			}
		}
	}
	std::cout << "Finish iter" << endl;
	mesh.write_m("fenzhi.m");
	//去分支
	while (true)
	{
		bool flag = false;
		for (MeshVertexIterator<CGBMesh> it(&mesh); !it.end(); it++)
		{
			CVertex *t_vert = it.value();
			int e_sum = 0;
			CEdge *r_edge;
			for (VertexOutHalfedgeIterator<CGBMesh> halfIt(&mesh, it.value()); !halfIt.end(); halfIt++)
			{
				CEdge *t_e = halfIt.value()->edge();
				if (t_e->string() == "sharp = true")
				{
					r_edge = t_e;
					e_sum++;
				}
			}
			if (e_sum == 1) //顶点周围只有一条sharp边，去掉
			{
				flag = true;
				r_edge->string() = "";
			}
		}
		if (flag == false)
			break;
	}
	mesh.write_m(outName);
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//创建二叉树储存Cut-Graph
	CEdge *root_edge = NULL; //根边
	for (MeshEdgeIterator<CGBMesh> edgeIt(&mesh); !edgeIt.end(); edgeIt++)
	{
		if (edgeIt.value()->string() == "sharp = true")
		{
			root_edge = edgeIt.value();
			//root_edge->string() = "";
			break;
		}
	}
	root_edge->isVisited() = true;
	CVertex *root_vert = root_edge->halfedge(0)->source(); //根顶点
	root_vert->isVisited() = true;
	root_vert->isNode() = true;
	vector<CEdge*> res_edges; //叶子节点之间的余下的边
	//DFS
	stack<CEdge*> e_stack;
	e_stack.push(root_edge);
	while (!e_stack.empty())
	{
		CEdge *cur_edge = e_stack.top(); //当前边
		cur_edge->isVisited() = true;
		e_stack.pop();
		CVertex *cur_vert = NULL; //当前边的out点
		if (cur_edge->halfedge(0)->target()->isVisited())
		{
			cur_vert = cur_edge->halfedge(0)->source();
		}
		else
		{
			cur_vert = cur_edge->halfedge(0)->target();
		}
		cur_vert->isVisited() = true;
		//cur_vert->string() = "sharp = true";
		vector<CEdge*> next_edges; //储存下一边
		for (VertexOutHalfedgeIterator<CGBMesh> halfIt(&mesh, (CGaussBonnetVertex*)cur_vert); !halfIt.end(); halfIt++)
		{
			CHalfEdge *next_half = halfIt.value();
			if (next_half->edge()->isVisited() == false && next_half->edge()->string() == "sharp = true") //未被访问的G山的吓一跳边
			{
				next_edges.push_back(next_half->edge());
			}
		}
		if (next_edges.size() == 0) //没有下一条边，叶子节点
		{
			res_edges.push_back(cur_edge);
		}
		if (next_edges.size() >= 2) //分支节点
		{
			cur_vert->isNode() = true;
		}
		for (int i = 0; i < next_edges.size(); i++)
		{
			next_edges[i]->preEdge() = cur_edge; //赋值”前“边
			e_stack.push(next_edges[i]);
		}
		//////////////////////
		//mesh.write_m(outName);
	}
	//删除res_edges中重复的边
	vector<CEdge*> result_edges;
	for (int i = 0; i < res_edges.size(); i++)
	{
		bool flag = false;
		for (int j = 0; j < result_edges.size(); j++)
		{
			if (result_edges[j] == res_edges[i])
				flag = true;
		}
		if (flag == false)
			result_edges.push_back(res_edges[i]);
	}
	//从叶子节点往回找 求4个基底
	vector<vector<CEdge*>> bases; //储存所有基底
	for (int i = 0; i < result_edges.size(); i++)
	{
		CEdge* r_edge = result_edges[i]; //叶子边
		CVertex *r_v1 = r_edge->halfedge(0)->target();
		CVertex *r_v2 = r_edge->halfedge(0)->source();
		CVertex *beg_v = NULL; //叶子顶点
		if (r_v1->isNode() == true)
			beg_v = r_v2;
		else
			beg_v = r_v1;
		vector<CEdge*> b_edges; //储存一个基底
		b_edges.push_back(r_edge);
		CEdge *cur_edge = r_edge;
		CVertex *cur_vt = beg_v;
		CVertex *pre_vt = NULL;
		while (true)
		{
			cur_edge = cur_edge->preEdge(); //往回找
			b_edges.push_back(cur_edge);
			if (cur_edge->halfedge(0)->target() == cur_vt)
				pre_vt = cur_edge->halfedge(0)->source();
			else
				pre_vt = cur_edge->halfedge(0)->target();
			cur_vt = pre_vt;
			if (pre_vt == r_v1 || pre_vt == r_v2)
			{
				break;
			}
		}
		bases.push_back(b_edges);
	}
	//输出四个基底
	for (int i = 0; i < bases.size(); i++)
	{
		//清理mesh
		for (MeshEdgeIterator<CGBMesh> it(&mesh); !it.end(); it++)
			it.value()->string() = "";
		//高亮一个同伦群基底
		vector<CEdge*> t_base = bases[i];
		for (int j = 0; j < t_base.size(); j++)
		{
			t_base[j]->string() = "sharp = true";
		}
		//写
		string s = "BASE_";
		char c = i + 49;
		s = s + c;
		s = s + '.';
		s = s + 'm';
		mesh.write_m(s.c_str());
	}

	time(&end);
	cout << difftime(end, start) << endl;

	////高亮 分支节点 res边
	//for (MeshVertexIterator<CGBMesh> i(&mesh); !i.end(); i++)
	//{
	//	if (i.value()->isNode() == true)
	//		i.value()->string() = "sharp = true";
	//}
	//for (int i = 0; i < res_edges.size(); i++)
	//{
	//	res_edges[i]->string() = "";
	//}

	//mesh.write_m(outName);
}
void hull(char *fileName, char *outName)
{
	time_t start, end;
	time(&start);
	CGBMesh mesh;
	mesh.read_vertices_m(fileName);

	cout << mesh.numVertices() << endl;

	//提升到3D
	for (list<CGaussBonnetVertex*>::iterator it = mesh.vertices().begin(); it != mesh.vertices().end(); it++)
	{
		CVertex *t_vertex = *it;
		t_vertex->string() = "sharp = true";
		CPoint t_pot = t_vertex->point();
		t_pot = CPoint(t_pot[0], t_pot[1], t_pot[0] * t_pot[0] + t_pot[1] * t_pot[1]); //(x, y, x^2+y^2)
		t_vertex->point() = t_pot;
	}
	//mesh.write_m("meshes/VVVVVV.m");

	//初始四个点
	CGaussBonnetVertex *ver1 = mesh.idVertex(31);
	//ver1->string() = "sharp = true";
	CGaussBonnetVertex *ver2 = mesh.idVertex(0);
	//ver2->string() = "sharp = true";
	CGaussBonnetVertex *ver3 = mesh.idVertex(1);
	//ver3->string() = "sharp = true";
	CGaussBonnetVertex *ver4 = mesh.idVertex(9);
	//ver4->string() = "sharp = true";
	//创建四面体
	int faceId = 0;
	CGaussBonnetVertex *face1[3] = { ver1, ver2, ver3 };  //注意点顺序
	CGaussBonnetVertex *face2[3] = { ver1, ver3, ver4 };
	CGaussBonnetVertex *face3[3] = { ver1, ver4, ver2 };
	CGaussBonnetVertex *face4[3] = { ver4, ver3, ver2 };
	mesh.createHullFace(face1, faceId++);
	mesh.createHullFace(face2, faceId++);
	mesh.createHullFace(face3, faceId++);
	mesh.createHullFace(face4, faceId++);
	mesh.writeHull_m(outName);
	//Smesh.write_m(outName);
	//建凸包
	//初始化四个点
	ver1->state() = IN_FIRE;
	ver2->state() = IN_FIRE;
	ver3->state() = IN_FIRE;
	ver4->state() = IN_FIRE;
	//开始遍历其他点
	for (list<CGaussBonnetVertex*>::iterator vertexIt = mesh.vertices().begin(); vertexIt != mesh.vertices().end(); vertexIt++)
	{
		CVertex *t_ver = *vertexIt;
		//t_ver->string() = "";
		///////////////////////////////
		//mesh.writeHull_m(outName);
		///////////////////////////////
		if (t_ver->state() != IN_FIRE) //未被访问
		{
			t_ver->state() = IN_FIRE; //标记为访问过
			vector<CFace*> det_faces;//用于储存det为负的面
			for (list<CGaussBonnetFace*>::iterator faceIt = mesh.faces().begin(); faceIt != mesh.faces().end(); faceIt++) //遍历所有面 判断是否朝向当前点 
			{
				CFace *t_face = *faceIt;
				CVertex *t_vers[3];//当前面上的三个点
				for (int i = 0; i < 3; i++)
				{
					t_vers[i] = t_face->faceVertices(i);
				}
				int x = 0, y = 1, z = 2;
				int j = 0, k = 1, l = 2;
				//行列式！！！
				double det = (t_vers[j]->point()[x] - t_ver->point()[x])*
					(t_vers[k]->point()[y] - t_ver->point()[y])*
					(t_vers[l]->point()[z] - t_ver->point()[z]) +
					(t_vers[l]->point()[x] - t_ver->point()[x])*
					(t_vers[j]->point()[y] - t_ver->point()[y])*
					(t_vers[k]->point()[z] - t_ver->point()[z]) +
					(t_vers[j]->point()[z] - t_ver->point()[z])*
					(t_vers[k]->point()[x] - t_ver->point()[x])*
					(t_vers[l]->point()[y] - t_ver->point()[y]) -
					(t_vers[j]->point()[z] - t_ver->point()[z])*
					(t_vers[k]->point()[y] - t_ver->point()[y])*
					(t_vers[l]->point()[x] - t_ver->point()[x]) -
					(t_vers[j]->point()[x] - t_ver->point()[x])*
					(t_vers[k]->point()[z] - t_ver->point()[z])*
					(t_vers[l]->point()[y] - t_ver->point()[y]) -
					(t_vers[l]->point()[z] - t_ver->point()[z])*
					(t_vers[j]->point()[y] - t_ver->point()[y])*
					(t_vers[k]->point()[x] - t_ver->point()[x]);
				if (det < 0) //面积为负 建面
				{
					//mesh.deleteFace((CGaussBonnetFace*)t_face); //不能直接删，否则空指针异常！！！
					det_faces.push_back(t_face); //先存起来
					//t_face->string() = "sharp = true";
				}
			}
			//便利所有det为负的面
			vector<vector<CVertex*>> edges_line; //存储建面序列
			for (int d = 0; d < det_faces.size(); d++)
			{
				CFace *det_t_face = det_faces[d];
				for (int i = 0; i < 3; i++) //两个点两个点地加入
				{
					vector<CVertex*> two_v;
					two_v.push_back(det_t_face->faceVertices(i));
					two_v.push_back(det_t_face->faceVertices((i + 1) % 3));
					edges_line.push_back(two_v);
				}
				mesh.deleteHullFace((CGaussBonnetFace*)det_t_face); //删除原有面
			}
			///////////////////////////////
			//mesh.writeHull_m(outName);
			///////////////////////////////
			//遍历所有边界边 （建面序列中非重复边）
			for (int i = 0; i < edges_line.size(); i++)
			{
				vector<CVertex*> create_face_line = edges_line[i];
				bool flag = false; //不是重复的
				for (int j = 0; j < edges_line.size(); j++)
				{
					vector<CVertex*> com_line = edges_line[j];
					if (com_line[0]->id() == create_face_line[1]->id() && com_line[1]->id() == create_face_line[0]->id())
					{
						flag = true;  //重复边
					}
				}
				if (flag == false)
				{
					CGaussBonnetVertex *new_face_vers[3] = { (CGaussBonnetVertex*)create_face_line[0], (CGaussBonnetVertex*)create_face_line[1], (CGaussBonnetVertex*)t_ver };
					mesh.createHullFace(new_face_vers, faceId++);
				}
			}
			///////////////////////////////
			//mesh.writeHull_m(outName);
			///////////////////////////////
		}
	}
	//mesh.writeHull_m("meshes/VV.m");
	//删掉向上的面
	vector<CGaussBonnetFace*> nor_up_faces;
	for (list<CGaussBonnetFace*>::iterator faceIt = mesh.faces().begin(); faceIt != mesh.faces().end(); faceIt++)
	{
		CFace *t_f = *faceIt;
		CPoint vi = t_f->faceVertices(0)->point() - t_f->faceVertices(1)->point();
		CPoint vj = t_f->faceVertices(1)->point() - t_f->faceVertices(2)->point();
		CPoint f_n = vi ^ vj;
		if (f_n * CPoint(0, 0, 1) > 0)
		{
			nor_up_faces.push_back((CGaussBonnetFace*)t_f);
		}
	}
	for (int i = 0; i < nor_up_faces.size(); i++)
	{
		mesh.deleteHullFace(nor_up_faces[i]);
	}
	mesh.writeHull_m("RAND_R.m");
	//投影回平面
	for (list<CGaussBonnetVertex*>::iterator vertexIt = mesh.vertices().begin(); vertexIt != mesh.vertices().end(); vertexIt++)
	{
		CVertex *t_v = *vertexIt;
		CPoint t_p = t_v->point();
		t_v->point() = CPoint(t_p[0], t_p[1], 0);
	}

	//mesh.writeHull_m(outName);
	mesh.writeHull_m(outName);
	time(&end);
	cout << difftime(end, start) << endl;
	system("pause");
}
void iterate_circle(char *fileName, char *outName)
{
	CGBMesh mesh;
	mesh.read_m(fileName);
	time_t start, end;
	time(&start);
	//CGaussCurvature<CGBMesh> GC(&mesh);

	//double e = 1e-7; //调和能量差值

	//GC._calculate_curvature();

	//GC._calculate_Euler_characteristics();

	//验证调和映照
	std::cout << mesh.numVertices() << " " << mesh.numEdges() << " " << mesh.numFaces() << std::endl;
	std::cout << mesh.vertices().size() << std::endl;

	//求K值
	int b_edge_sum = 0;
	for (MeshEdgeIterator<CGBMesh> edgeIter(&mesh); !edgeIter.end(); edgeIter++)
	{
		CHalfEdge *halfEdge_0 = edgeIter.value()->halfedge(0);
		CHalfEdge *halfEdge_1 = edgeIter.value()->halfedge(1);
		CPoint vi, vj, vk, vl;
		vi = halfEdge_0->source()->point();
		vj = halfEdge_0->target()->point();
		vk = halfEdge_0->he_next()->target()->point();
		if (halfEdge_1) //非边界的边 半边不为空
		{
			vl = halfEdge_1->he_next()->target()->point();
		}
		if (edgeIter.value()->boundary())
		{
			CPoint px = vi - vk;
			CPoint py = vj - vk;
			double k = (px * py) / (px^py).norm();
			//std::cout << k << std::endl;
			edgeIter.value()->m_k = k;
			edgeIter.value()->m_k = k;
			b_edge_sum++;
			//std::cout << edgeIter.value()->m_k << std::endl;
		}
		else
		{
			CPoint px = vi - vk;
			CPoint py = vj - vk;
			CPoint pm = vi - vl;
			CPoint pn = vj - vl;
			double k1 = (px * py) / (px^py).norm();
			double k2 = (pm * pn) / (pm^pn).norm();
			edgeIter.value()->m_k = k1 + k2;
			edgeIter.value()->m_k = k1 + k2;
			//std::cout << edgeIter.value()->m_k << std::endl;
		}
	}
	//std::cout << "b_edge_sum= " << b_edge_sum << std::endl;

	//遍历所有边界点 562
	CVertex *begin_vert;
	int b_vert_sum = 0;
	for (MeshVertexIterator<CGBMesh> meshIter(&mesh); !meshIter.end(); meshIter++)
	{
		if (meshIter.value()->boundary())
		{
			//std::cout << meshIter.value()->point()[0] <<","<< meshIter.value()->point()[1]<<","<< meshIter.value()->point()[2]<< std::endl;
			begin_vert = meshIter.value();
			//b_vert_sum++;
			break;
		}
	}
	//std::cout << "b_vert_sum= " << b_vert_sum << std::endl;
	//std::cout << begin_vert->point()[0] << "," << begin_vert->point()[1] << "," << begin_vert->point()[2] << std::endl;

	//遍历边界点，计算各边长
	std::list<double> lens;
	double len_sum = 0;
	CVertex *pre_vert = begin_vert;
	CVertex *next_vert;
	for (next_vert = begin_vert->most_clw_out_halfedge()->target(); next_vert != begin_vert; //逆时针100次 互相直来直去？？？
		next_vert = next_vert->most_clw_out_halfedge()->target())
	{
		double len = (next_vert->point() - pre_vert->point()).norm();
		lens.push_back(len);
		len_sum += len;
		pre_vert = next_vert;
		std::cout << len << std::endl;
	}
	//最后一次
	double len = (next_vert->point() - pre_vert->point()).norm();
	lens.push_back(len);
	len_sum += len;
	pre_vert = next_vert;
	std::cout << len << std::endl;

	std::cout << lens.size();
	std::cout << "len_sum= " << len_sum << std::endl; //609.703

	//将边界点移至拓扑圆盘
	double c_cir = len_sum;
	double r_cir = c_cir / (2 * PI);
	double p_len = 0;
	begin_vert->point() = CPoint(r_cir, 0, 0);
	for (CVertex *cur_vert = begin_vert->most_clw_out_halfedge()->target(); cur_vert != begin_vert;
		cur_vert = cur_vert->most_clw_out_halfedge()->target())
	{
		p_len += lens.front();
		double p_x = cos(p_len / len_sum * 2 * PI) * r_cir;
		double p_y = sin(p_len / len_sum * 2 * PI) * r_cir;
		lens.pop_front();
		cur_vert->point() = CPoint(p_x, p_y, 0);
	}
	//将非边界点归零
	for (MeshVertexIterator<CGBMesh> meshIter(&mesh); !meshIter.end(); meshIter++)
	{
		if (!meshIter.value()->boundary())
			meshIter.value()->point() = CPoint(0, 0, 0);
	}
	//////////////////////
	//mesh.write_m("meshes/1111.m");
	//开始迭代
	double k1 = 0;
	for (MeshEdgeIterator<CGBMesh> edgeIter(&mesh); !edgeIter.end(); edgeIter++)
	{
		CPoint c1, c2;
		c1 = edgeIter.value()->halfedge(0)->source()->point();
		c2 = edgeIter.value()->halfedge(0)->target()->point();
		edgeIter.value()->length() = (c1 - c2).norm();
		k1 += edgeIter.value()->m_k * edgeIter.value()->length() * edgeIter.value()->length();
	}
	for (int i = 0; true; i++)
	{
		//计算点一次迭代之后的坐标
		for (MeshVertexIterator<CGBMesh> verIter(&mesh); !verIter.end(); verIter++)
		{
			if (!verIter.value()->boundary())
			{
				double k_sum = 0;
				for (VertexOutHalfedgeIterator<CGBMesh> edgeIter(&mesh, verIter.value()); !edgeIter.end(); edgeIter++)
					k_sum += edgeIter.value()->edge()->m_k;
				CPoint k_pot = CPoint(0, 0, 0);
				for (VertexOutHalfedgeIterator<CGBMesh> edgeIter(&mesh, verIter.value()); !edgeIter.end(); edgeIter++)
				{
					double t_k = edgeIter.value()->edge()->m_k;
					CPoint t_p = edgeIter.value()->target()->point();
					k_pot += t_p * t_k;
				}
				verIter.value()->point() = k_pot / k_sum;
			}
		}
		std::cout << "第" << i << "次迭代" << std::endl;
		//调和能量
		double k2 = 0;
		for (MeshEdgeIterator<CGBMesh> edgeIter(&mesh); !edgeIter.end(); edgeIter++)
		{
			CPoint c1, c2;
			c1 = edgeIter.value()->halfedge(0)->source()->point();
			c2 = edgeIter.value()->halfedge(0)->target()->point();
			edgeIter.value()->length() = (c1 - c2).norm();
			k2 += edgeIter.value()->m_k * edgeIter.value()->length() * edgeIter.value()->length();
		}
		//std::cout << "k1=" << k1 << std::endl;
		//std::cout << "k2=" << k2 << std::endl;
		std::cout << "E=" << fabs(k1 - k2) << std::endl;
		if (fabs(k1 - k2) < 1e-7)
			break;
		k1 = k2;
		//if (i == 1000)
		//{
		//	mesh.write_m("meshes/2222.m");
		//}
	}

	//UV
	for (MeshVertexIterator<CGBMesh> it(&mesh); !it.end(); it++)
	{
		CVertex *v = it.value();
		CPoint p = v->point();
		v->uv()[0] = p[0];
		v->uv()[1] = p[1];
		v->_to_string();
	}

	mesh.write_m(outName);
	std::cout << "输出.m" << std::endl;
	time(&end);
	cout << difftime(end, start) << endl;
}

/*!
 *	Compute face normal
 */
template<typename M>
void CGaussCurvature<M>::_calculate_face_normal()//计算面的法向量
{
	//insert your code here
}

/*!
 *	Compute face normal
 */
template<typename M>
void CGaussCurvature<M>::_calculate_vertex_normal()//计算面的法向量
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
	V = m_pMesh->numVertices();
	E =m_pMesh->numEdges();
	F = m_pMesh->numFaces();
	std::cout << "Vertices: " << V << " Faces: " << F << " Edges: " << E << std::endl;
	int euler = V + F - E;//欧拉示性数
	std::cout << "Euler Characteristic Number " << euler << std::endl;	
	int b;

	/*CBoundary boundary(m_pMesh);
	b = boundary.loops().size();
	int g=(2-b-euler)/2;
	std::cout << "Number of boundaries " << b << std::endl;
	std::cout << "Genus " << g << std::endl;*/
	//insert your code in this function
	

	CBoundary<CGBMesh> boundary(m_pMesh);
	b = boundary.loops().size();

	int g = (2 - b - euler) / 2;
	std::cout << "Number of boundaries " << b << std::endl;
	std::cout << "Genus " << g << std::endl;

}


};
#endif

