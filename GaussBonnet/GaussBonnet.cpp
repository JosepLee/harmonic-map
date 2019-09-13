// GaussBonnet.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>

#include "GaussCurvature.h"
#include "CGaussBonnetMesh.h"
using namespace MeshLib;

void iterate_circle() {
		CGBMesh mesh;
		mesh.read_m("sophie.m");
		CGaussCurvature<CGBMesh> GC(&mesh);
		GC._calculate_K();
		GC.init_mapping();
		mesh.write_m("sophie_init.m");
		GC.mapping();
		mesh.write_m("sophie_out.m");
}
void GaussBonnet() {
	CGBMesh mesh;
	mesh.read_m("eight.m");
	CGaussCurvature<CGBMesh> GC(&mesh);
	GC._calculate_curvature();
	GC._calculate_Euler_characteristics();

}

	int main(int argc, char ** argv)
	{
		CGBMesh mesh;
		if (argc < 2)
		{
			std::cout << "Usage: GaussBonnet.exe inputMesh.m" << std::endl;
			//system("PAUSE");
			//return 0;
		}
		//高斯波涅定理
		GaussBonnet();
		//调和映照拓扑圆盘
		iterate_circle("sophie.m", "sophie_final_2.m");
		//球体调和映照
		//iterate_globe("brain.m", "brain_final.m");   //拓扑球面

		//同伦群基底

		//homotopy("eight.m", "eight_cutGraph.m");

		//homotopy("meshes/loveme.m", "meshes/loveme_cutGraph.m");

		//homotopy("meshes/BoyAndGirl.m", "meshes/BoyAndGirl_cutGraph.m");

		//凸包映射

		//hull("rand_1.m", "RAND_1.m");



		//std::cout << "Check :" << GC.total_C - 2 * GC.lambda*PI << endl;
		system("PAUSE");
		return 0;
	}

