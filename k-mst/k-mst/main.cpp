#include <iostream>
#include <stdio.h>
#include <math.h> 
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream> 
#include <vector> 
#include "helper.h" 
#include "simplify.h"
#include <pcl/point_types.h>
#include <pcl/point_cloud.h> 
#include "EMst.h"
#include <filesystem>
//函数声明
//int order_element_in_vector(const vector<int> &v, int element);
//没有用到main.cpp中的内容，可以将实现粒子滤波的过程写在EMst.h中，这样在EMst.h中的函数就可以直接调用了
int main(int argc, char* argv[])
{
	EMst g("G://Code//data//1_AIM//AIM.pcd", false);
	//EMst g("G://Code//data//2_castle//castle.pcd", false);
	//EMst g("G://Code//data//3_boundary.pcd", false);
	//EMst g("G://Code//data//7_old//old1.pcd", false);
	//EMst g("G://Code//data//5_jiduchurch//jiduchurch.pcd", false);
	int K = 30;
	pcl::PointCloud<RawPointType>::Ptr m_knnpoint(new pcl::PointCloud<RawPointType>);
	int T = 500;//粒子扩增T次
	int Pa = 300;//粒子的个数
	vector<trace> paths;
	m_knnpoint = (g.raw).makeShared();
	g.get_new_widght(m_knnpoint, Pa, T, K, paths);
	g.outputPLFile(m_knnpoint, paths, "//./result/1_AIM/1.4-300-500-RD.pl");// //./result/7_OLD/2000-20-OLD
	// "//./result/3_boundary/3_boundary-3000-50-k30.pl
}

////void mains(EMst &g, int Orderitreator) {//接受一个EMst对象和一个整数作为参数
////
////	//PCL点云类型
////	int K = 20;
////	pcl::PointCloud<MyPointType>::Ptr m_knnpoint(new pcl::PointCloud<MyPointType>);//创建一个新的点云对象，并将其存储在m_knnpoint指针，使用new关键字在堆上分配了内存，并通过调用pcl::PointCloud<MyPointType>的构造函数进行初始化。
////	if (pcl::io::loadPCDFile("castle.pcd", *m_knnpoint) == -1) {
////		PCL_ERROR("Couldn't read file input.pcd\n");
////		return(-1);
////	}
////	//pcl::PointCloud<MyPointType>::Ptr c_knnpoint(new pcl::PointCloud<RawPointType>);
////	//pcl::PointCloud<MyPointType>::Ptr b_knnpoint(new pcl::PointCloud<RawPointType>);
////	//深度拷贝
////	m_knnpoint = (g.mypoint).makeShared();//在不影响 g.knnPoint 的情况下，创建了 m_knnpoint 的一个共享指针，以便进行后续的操作，而不复制整个点云数据
////	/*c_knnpoint = (g.cornerknnPoint).makeShared();
////	b_knnpoint = (g.boundaryknnPoint).makeShared();*/
////	//构图
////	g.CreateALGraph(m_knnpoint, K);//调用了g对象的CreateALGraph函数，传递了m_knnpoint,c_knnpoint和20作为参数。
////
////	//调用EMst.h中的函数，进行粒子滤波处理
////	g.get_new_widght(m_knnpoint, K);
////}
//void txtToPCD(string txtname)
//{
//	int num_txt;
//	cout << "Now is translating txt file to PCD and the point type is XYZ8D" << endl;
//	//定义一种类型表示TXT中的点云格式xyz
//	typedef struct TXT_Point_XYZ
//	{
//		double x;
//		double y;
//		double z;
//		float intensity;
//		float gvalue;
//		float p0;
//		float p1;
//		float ptype;
//
//	}TOPOINT_XYZ8D;
//
//	//读取txt文件
//	FILE *fp_txt;
//	TXT_Point_XYZ txt_points;
//	vector<TXT_Point_XYZ> my_vTxtPoints;
//	fp_txt = fopen(txtname.c_str(), "r");
//
//	if (fp_txt)
//	{
//		while (fscanf(fp_txt, "%lf %lf %lf %f %f %f %f %f", &txt_points.x, &txt_points.y, &txt_points.z, &txt_points.intensity,
//			&txt_points.gvalue, &txt_points.p0, &txt_points.p1, &txt_points.ptype) != EOF)
//			/*&txt_points.x, &txt_points.y, &txt_points.z,&txt_points.Tlambda2, &txt_points.Tlambda1,
//				&txt_points.Tlambda0, &txt_points.gz, &txt_points.gy, &txt_points.gx, &txt_points.gvalue, &txt_points.pca_lamda2,
//				&txt_points.pca_lamda1, &txt_points.pca_lamda0, &txt_points.intensity*/
//		{//将点存入容器尾部
//			my_vTxtPoints.push_back(txt_points);
//		}
//	}
//	else cout << "读取txt文件失败" << endl;
//
//	num_txt = my_vTxtPoints.size();
//
//	//写入点云数据
//	pcl::PointCloud<MyPointType> ::Ptr Bcloud(new pcl::PointCloud<MyPointType>);
//	Bcloud->width = num_txt;
//	Bcloud->height = 1;
//	Bcloud->is_dense = false;
//	Bcloud->points.resize(Bcloud->width * Bcloud->height);
//	for (int i = 0; i < Bcloud->points.size(); ++i)
//	{
//		Bcloud->points[i].x = my_vTxtPoints[i].x;
//		Bcloud->points[i].y = my_vTxtPoints[i].y;
//		Bcloud->points[i].z = my_vTxtPoints[i].z;
//		Bcloud->points[i].intensity = my_vTxtPoints[i].intensity;
//		Bcloud->points[i].gvalue = my_vTxtPoints[i].gvalue;
//		Bcloud->points[i].p0 = my_vTxtPoints[i].p0;
//		Bcloud->points[i].p1 = my_vTxtPoints[i].p1;
//		Bcloud->points[i].ptype = my_vTxtPoints[i].ptype;
//	}
//	string tname = txtname.substr(txtname.length() - 10, txtname.length() - 4);
//	string suffix = ".pcd";
//	string pcdname = tname + suffix;
//	pcl::io::savePCDFileASCII("1.pcd", *Bcloud);
//	cout << "从 txt_pcd.txt读取" << Bcloud->points.size() << "点写入txtTO.pcd" << endl;
//}