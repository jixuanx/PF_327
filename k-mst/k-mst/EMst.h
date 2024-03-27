#pragma once
#include<iostream>
#include<string>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common_headers.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <numeric>//提供了一些数值算法和操作的函数,以下用到的是accumulate()：对指定范围内的元素进行累积操作，返回累积结果。
#include <Eigen/Dense>
#include <array>
#include <cstdlib>
#include <ctime>
#include <random>
#include <filesystem>
//#include<omp.h>
#define path "G:/Code/new/k-mst/k-mst/lookprocess/"
typedef int VertexType;
typedef float EdgeType;
#define MAXVEX 5000
#define Inf 255800.000//Inf 是一个用于表示无穷大的常量。通过将 Inf 定义为 255800.000，你可以在代码中使用 Inf 来表示一个较大的值。
using namespace std;
//the structure of data in MST
struct trace {
	vector<int> index;
	float weight = 0.0;
};
typedef struct EdgeNode {
	int adjvex;
	EdgeType weight;
	struct EdgeNode *next;
}EdgeNode;
typedef struct VertexNode {
	VertexType data;
	int degree;
	EdgeNode *firstedge;
	bool visit;
	bool is_cor;
}AdjList;
typedef struct {
	AdjList *adjList;
	int numVertsxes, numEdges;
}GraphAdjList;
struct EndNode {
	int data;
	int degree;
	vector<int> lineorder;
};
typedef struct degreeMore2Node
{
	int data;
	int degree;
	vector<EndNode> endnode;
	int maxlength;

}crticalNode;
struct save {
	int depthnum;
	vector<int> dirdata;
};
struct closedge {
	float lowcost; //closedge[i].lowcost表示当前与第i个顶点相连的的权值最小的边
	int vex; //closedeg[i].vex表示与第i个顶点相连的顶点的位置 
};
struct MyPointType
{
	PCL_ADD_POINT4D;  //该点类型有4个元素      

	/*尝试新增一个自定义*/
	float intensity;
	float gvalue;
	float p0;
	float p1;
	float ptype;
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW   //确保new操作符对齐操作 

}EIGEN_ALIGN16;   //强制SSE 对齐
POINT_CLOUD_REGISTER_POINT_STRUCT(MyPointType,    //注册点类型宏
(float, x, x)
(float, y, y)
(float, z, z)
(float, intensity, intensity)
(float, gvalue, gvalue)
(float, p0, p0)
(float, p1, p1)
(float, ptype, ptype)
)
struct RawPointType    //定义点类型结构
{
	PCL_ADD_POINT4D;  //该点类型有4个元素      

	/*尝试新增一个自定义*/
	float intensity;
	float gvalue;
	float gx;
	float gy;
	float gz;
	float Tlambda0;
	float Tlambda1;
	float Tlambda2;
	int iscorner;
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW   //确保new操作符对齐操作 

}EIGEN_ALIGN16;   //强制SSE 对齐
POINT_CLOUD_REGISTER_POINT_STRUCT(RawPointType, //注册点类型宏
(float, x, x)
(float, y, y)
(float, z, z)
(float, intensity, intensity)
(float, gvalue, gvalue)
(float, gx, gx)
(float, gy, gy)
(float, gz, gz)
(float, Tlambda0, Tlambda0)
(float, Tlambda1, Tlambda1)
(float, Tlambda2, Tlambda2)
(int, iscorner, iscorner)
)
typedef pcl::PointCloud<pcl::PointXYZ> PointXYZ;
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr p_PointXYZ;
typedef pcl::PointCloud<MyPointType> Mypoint;//用更短的名称 Mypoint 来代替 pcl::PointCloud<MyPointType>
typedef pcl::PointCloud<MyPointType>::Ptr p_Mypoint;
typedef pcl::PointCloud<RawPointType> rawPoint;
typedef pcl::PointCloud<RawPointType>::Ptr p_rawPoint;
void pruning4_plus(GraphAdjList* tree1, int index, int* depth, vector<int> &a);//根据树的特定条件进行修剪，从而移除或保留某些节点或子树。
static EdgeNode * make_node(const int pos, const float distance);//创建一个新的边节点（EdgeNode）。
static void initial_graph(GraphAdjList * graph, GraphAdjList * kruskal_tree);//初始化一个图结构的副本（kruskal_tree）,并没有复制原始图中的边信息，只是将顶点的数据信息和边表初始化为相同的结构。
template<typename T, typename... U>
void logger(T t, U... ts);
static EdgeNode * make_node(const int pos, const float distance)
{
	EdgeNode * new_node = (EdgeNode *)malloc(sizeof(EdgeNode));
	if (new_node == NULL)
		exit(1);

	new_node->next = NULL;
	new_node->weight = distance;
	new_node->adjvex = pos;

	return new_node;
}
static void initial_graph(GraphAdjList * graph, GraphAdjList * kruskal_tree)
{
	int i;
	kruskal_tree->numVertsxes = graph->numVertsxes;
	kruskal_tree->numEdges = graph->numEdges;

	for (i = 0; i < graph->numVertsxes; i++)
	{
		kruskal_tree->adjList[i].data = graph->adjList[i].data;
		kruskal_tree->adjList[i].firstedge = NULL;
	}
}
void pruning4_plus(GraphAdjList* tree1, int index, int* depth, vector<int> &a)
{
	EdgeNode *p = tree1->adjList[index].firstedge;
	if (tree1->adjList[index].degree == 1)
	{
		tree1->adjList[index].visit = true;
		a.push_back(index);
		(*depth)++;
	}

	if (tree1->adjList[index].degree == 2)
	{
		if (!tree1->adjList[index].visit)
		{
			tree1->adjList[index].visit = true;
			a.push_back(index);
			if (tree1->adjList[p->adjvex].visit)
			{
				if (p->next) {
					p = p->next;
				}
				if (!tree1->adjList[p->adjvex].visit && tree1->adjList[p->adjvex].degree <= 2)
				{
					pruning4_plus(tree1, p->adjvex, depth, a);
				}
				else
				{
					if (tree1->adjList[p->adjvex].degree > 2) {
						a.push_back(p->adjvex);
					}
					else
					{
						p = tree1->adjList[index].firstedge;
						a.push_back(p->adjvex);
					}

				}
			}
			else
			{
				(*depth)++;
				pruning4_plus(tree1, p->adjvex, depth, a);
			}
		}
	}

	if (tree1->adjList[index].degree >= 3)
	{
		tree1->adjList[index].visit = true;
		a.push_back(index);
		(*depth)++;
	}

}
struct EdgeInfo
{
	PCL_ADD_POINT4D;  //该点类型有4个元素      
	/*尝试新增自定义*/
	int start;
	int end;
	long double wight_value;
	EdgeInfo()
	{
		memset(this, 0, sizeof(EdgeInfo));
	}
};
class EMst
{//指定了以下成员的访问权限为公共
public:
	GraphAdjList Graph;
	GraphAdjList CornerGraph;
	GraphAdjList BoundaryGraph;
	//Mypoint点云数据类型为5维<I,g,p0,p1,ptype>
	Mypoint mypoint;
	Mypoint knnPoint;
	//rawPoint点云类型为8维数据<I,g,gx,gy,gz,lamb1,lamb2,lamb3>
	rawPoint raw;
	rawPoint corraw;
	rawPoint bouraw;
	rawPoint rawknnPoint;
	rawPoint cornerknnPoint;
	rawPoint boundaryknnPoint;
public:
	EMst(string allPointCloudPath, bool flag);//定义了 EMst 类的构造函数，该函数接受四个参数。
	void CreateALGraph(p_rawPoint P, int K)//该函数接受两个参数，用于创建无向图结构。
	{
		if (this->Graph.adjList == nullptr) {
			// 可能需要在这里进行 Graph 对象的初始化，具体根据需要修改
		}
		//get_new_widght(P, Pa, T, K, paths);//作用是什么++++++++++++++++++++++++++++++++++++++++++++++++
		return;//直接返回，结束了 CreateALGraph 函数的执行，因此函数将不会执行后续的代码。
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::copyPointCloud(*P, *cloud);
		pcl::KdTreeFLANN<pcl::PointXYZ> searchtree;
		searchtree.setInputCloud(cloud);
		pcl::PointXYZ searchpoint, pointj;
		std::vector<int> Neighbors;
		std::vector<float> SquaredDistance;
		int i, j, k;
		EdgeNode* e;
		cout << "knn点数：" << this->Graph.numVertsxes << endl;
		//\\后续要将计算好的权值赋值到这里来
		for (i = 0; i < this->Graph.numVertsxes; i++) {
			//ofstream fp("po.txt", ios::app);
			//fp << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << " " << P->points[i].p0 << " " << P->points[i].p1 << " " << endl;
			searchpoint.x = P->points[i].x;
			searchpoint.y = P->points[i].y;
			searchpoint.z = P->points[i].z;
			searchtree.nearestKSearch(searchpoint, K, Neighbors, SquaredDistance);
			float p0 = 0.5f *(P->points[i].intensity + P->points[i].gvalue);
			float pi = p0;  //当前点的概率
			float pj = 0.0;
			float wij = 0.0;
			// 权值函数（粒子滤波）
			for (k = 1; k < Neighbors.size(); k++) {
				j = Neighbors[k]; //全局索引号，而k是局部的近邻索引
				float p0j = 0.5f *(P->points[j].intensity + P->points[j].gvalue);
				pj = p0j; //第k个邻域点的概率
				wij = 0;
				j = Neighbors[k]; // 全局索引号，而 k 是局部的近邻索引
				pj =p0j; // 第 k 个邻域点的概率
				//////pro = P->points[i].p1 /(P->points[i].p1* SquaredDistance[k]); // 计算 pro 值
				//////weight= P->points[i].p1 *  P->points[j].p1*exp(-1 * SquaredDistance[k] * SquaredDistance[k]);
				e = (EdgeNode*)malloc(sizeof(EdgeNode));
				e->adjvex = j;
				e->weight = wij;//
				e->next = this->Graph.adjList[i].firstedge;
				this->Graph.adjList[i].firstedge = e;

				e = (EdgeNode*)malloc(sizeof(EdgeNode));
				e->adjvex = i;
				e->weight = wij;
				e->next = this->Graph.adjList[j].firstedge;
				this->Graph.adjList[j].firstedge = e;
			}
		}
	}
	void get_new_widght(p_rawPoint P, int Pnum, int time, int k, vector<trace>& paths);
	void outputPLFile(p_rawPoint P, const vector<trace>& paths, string name);
	//std::vector<std::vector<long double>> Score(pcl::KdTreeFLANN<pcl::PointXYZ> kdtree,int last_start_index, int candidate_start_index, pcl::PointXYZ candidate_start, p_Mypoint P, p_rawPoint C, p_rawPoint B, int iter, int k, std::vector<int> index, std::vector<float> dist);
	/*std::vector<int> N1Score(p_rawPoint& P, int Pa, long double threshold);
	int Sampling(std::vector<std::vector<long double>>& scores_array);*/
	//void graph_to_Graph(std::vector<std::vector<float>>&graph);
	//void writeToCSV(const p_rawPoint& P, const std::vector<EdgeInfo>& edges, string name);
	~EMst()
	{}
};
float ProLine = 0.0;
float threshold_p = 0.0;
EMst::EMst(string allPointCloudPath, bool flag)
{
	if (flag) {
		pcl::io::loadPCDFile(allPointCloudPath, this->mypoint);
		cout << "点云数据输入完成." << endl;
	}
	else{
		pcl::io::loadPCDFile(allPointCloudPath, this->raw);
		logger("点云数据输入完成."); 
		int Point = this->raw.size();
		cout << "全部点的数量" << Point << "：" << endl;
		for (int i = 0; i < Point; i++) {
			float p0 = 0.5f *(this->raw.points[i].intensity + this->raw.points[i].gvalue);
			float p1 = 1 - p0;
			if (p0 * p1 < 0)logger("错误发生在第", i, "行");
			else {
				RawPointType p;
				p.x = this->raw.points[i].x;
				p.y = this->raw.points[i].y;
				p.z = this->raw.points[i].z;
				p.intensity = this->raw.points[i].intensity;
				p.gvalue = this->raw.points[i].gvalue;
				p.Tlambda0 = this->raw.points[i].Tlambda0;
				p.Tlambda1 = this->raw.points[i].Tlambda1;
				p.Tlambda2 = this->raw.points[i].Tlambda2;
				//this->raw.points.push_back(p);
			}
		}
	}
	//全部点云数据构图
	int Point = this->raw.size();
	//cout << "点云个数Point" << Point <<": " << endl;
	//define the points in KNN graph
	for (int i = 0; i < Point; i++) {
		float p0 = 0.5f *(this->raw.points[i].intensity + this->raw.points[i].gvalue);
		ProLine = p0;
		if (ProLine > threshold_p) this->rawknnPoint.push_back(this->raw.points[i]);
	}
	int rawknnPointSIZE = this->rawknnPoint.size();
	this->Graph.adjList = new AdjList[rawknnPointSIZE];
	this->Graph.numVertsxes = rawknnPointSIZE;
	//cout << "点云个数rawknnPoint" << ": " << rawknnPointSIZE  << endl;//点云个数Point81366 int NumknnPoint = this->rawknnPoint.size();
	GraphAdjList *G = &this->Graph;
	G->adjList = new AdjList[rawknnPointSIZE];
	G->numVertsxes = rawknnPointSIZE;
	for (int i = 0; i < rawknnPointSIZE; i++) {
		G->adjList[i].data = i;
		G->adjList[i].degree = 0;
		G->adjList[i].firstedge = NULL;
	}
}
//-------------------------------------------------------||||||||||粒子滤波核心部分||||||||||--------------------------------------------------------------------------------------
template <typename T> void releaseVector(vector<T>& vec) {
	vec.clear();
	vector<T>(vec).swap(vec);
}
template <typename T> void normalizeVector(std::vector<T>& vec) {
	// 计算向量中所有元素的和
	float sum = 0.0f;
	for (float element : vec) {
		sum += element;
	}
	// 归一化向量
	for (float& element : vec) {
		element /= sum;
	}
}
vector<int> RandomGenerateParticle(const int range, const int Pa)
{
	//1. 随机生成Pnum个粒子,存到容器中
	// 创建包含所有可能整数的向量
	std::vector<int> allIntegers;
	if (range <= 0) {
        // 处理 range 为非正数的情况
		cout << "range 为非正数" << endl;
        return allIntegers; // 返回一个空的向量
    }
	for (int i = 0; i < range; ++i) {
		allIntegers.push_back(i);
	}
	// 检查 allIntegers 是否为空
	if (allIntegers.empty()) {
		// 处理 allIntegers 为空的情况
		cout << "allIntegers 为空" << endl;
		return allIntegers; // 返回一个空的向量
	}
	// 使用随机设备生成器
	std::random_device rd;
	std::mt19937 g(rd());
	// 打乱序列
	std::shuffle(allIntegers.begin(), allIntegers.end(), g);
	// 检查 allIntegers 的大小是否足够
	if (Pa > allIntegers.size()) {
		// 处理 Pa 大于 allIntegers 大小的情况
		return allIntegers; // 返回全部元素
	}
	// 从打乱后的序列中选择前 numElements 个元素
	std::vector<int> randomIntegers(allIntegers.begin(), allIntegers.begin() + Pa);
	// 输出结果
	std::cout << "随机生成的" << Pa << "个不相同的整数：" << std::endl;
	return randomIntegers;
}
//****************************|||||||||||||粒子滤波的核心函数||||||||||||||*********************************
void EMst::get_new_widght(p_rawPoint P, int Pa, int T, int k, vector<trace>& paths)
{
	if (this->Graph.adjList == nullptr) {
		// 可能需要在这里进行 Graph 对象的初始化，具体根据需要修改
	}
	//1--定义粒子存储结构
	int Psize = P->points.size();
	cout << "P->points.size:  " << Psize << endl;
	vector<int> randomStart = RandomGenerateParticle(P->points.size(), Pa);//从所有点中随机选Pa个
	paths.resize(randomStart.size());//vector 是一种动态数组
	for (int i = 0; i < randomStart.size(); i++) {
		paths[i].index.push_back(randomStart[i]);
		paths[i].weight = 1.0 / randomStart.size();
		//cout << "定义粒子存储结构paths[i].weight: " << paths[i].weight << endl;
		//paths[i].weight = 1.0f / randomStart.size();
	}
	//输出随机选择的Pa个粒子
	outputPLFile(P, paths, "//./result/1_AIM/1.4-RD.pl");
	releaseVector(randomStart);
	//2--计算领域内分布律，先对KNN邻域进行构建
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::copyPointCloud(*P, *cloud);
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud(cloud);
	vector<int> kneighbor;
	vector<float> dist; 
	//3--粒子扩增
	for (int iter = 0; iter < T; iter++)
	{
		cout << iter << endl;
		for (int j = 0; j < paths.size(); j++)
		{
			//cout << "开始扩增" << endl;
			int curindex = *paths[j].index.rbegin();
			pcl::PointXYZ curPt;
			curPt.x = P->points[curindex].x;
			curPt.y = P->points[curindex].y;
			curPt.z = P->points[curindex].z;
			//cout << "建立Kdtree" << endl;
			kdtree.nearestKSearch(curPt, k, kneighbor, dist); //计算邻域内分布律-------const int neiborknei = 30;
			int samplingIndex = 0;
			for (int t = 0; t < 5; ) {//5是一个自定义的范围阈值，在最近的五个点中找到不包括在之前的路径中的扩增点
				//cout << "采样：" << samplingIndex << endl;
				auto it = find(paths[j].index.begin(), paths[j].index.end(), kneighbor[samplingIndex + 1]);//find函数是标准库 <algorithm> 中的函数，用于在指定范围内查找特定值
				if (it != paths[j].index.end()) {
					samplingIndex++;
					t++;
				}
				else break;
			}
			paths[j].index.push_back(kneighbor[samplingIndex + 1]);
			paths[j].weight += 1.0 / dist[samplingIndex + 1];
			releaseVector(kneighbor);
			releaseVector(dist);
		}
		cout << "扩增完成，计算权" << endl;
		vector<float> patheweight;
		for (int j = 0; j < paths.size(); j++) {
			patheweight.push_back(paths[j].weight);
		}
		normalizeVector(patheweight);
		float weightsqr = 0.0;
		for (int j = 0; j < paths.size(); j++) {
			paths[j].weight = patheweight[j];
			weightsqr += pow(paths[j].weight, 2.0);
		}
		//cout << "计算有效样本：" << endl;
		int Neffective = ceil(1.0f / weightsqr);
		cout << "计算有效样本：" << Neffective << endl;
		cout << "权计算完成" << endl;
	}
}

//输出pl文件的函数（√）
void EMst::outputPLFile(p_rawPoint P, const vector<trace>& paths, string name)
{
	string filename = path + name; // 使用预定义的路径信息拼接完整的文件路径
	ofstream PLfile1(filename, ios::app); // 打开文件流
	cout << "共有路径Path.size =" << paths.size() << endl;
	for (int m = 0; m < paths.size(); m++)
	{
		PLfile1 << "GOCAD PLine \n";
		PLfile1 << "HEADER{\n";
		PLfile1 << "name:" << "1-random-ponints" << m << "\n";
		PLfile1 << "value:" << paths[m].weight << "\n";
		PLfile1 << "}\n";
		PLfile1 << "ILINE\n";
		for (int i = 0; i < paths[m].index.size(); i++)
		{
			int num = paths[m].index[i];
			PLfile1 << "VRTX " << i << " " << P->points[num].x << " " << P->points[num].y << " " << P->points[num].z << "\n";
		}
		PLfile1 << "END";
	}
	PLfile1.close(); // 关闭文件流
}
////void EMst::outputPLFile(p_rawPoint P, const vector<trace>& paths, string name)
////{
////	string filename = path + name; // 使用预定义的路径信息拼接完整的文件路径
////	ofstream PLfile1(filename, ios::app); // 打开文件流
////	cout << "共有路径Path.size =" << paths.size() << endl;
////	for (int m = 0; m < paths.size(); m++)
////	{
////		PLfile1 << "GOCAD PLine \n";
////		PLfile1 << "HEADER{\n";
////		PLfile1 << "name:" << "1-random-ponints" << m << "\n";
////		//PLfile1 << "value:" << paths[m].weight << "\n";
////		PLfile1 << "}\n";
////		PLfile1 << "IPOINT\n";
////		for (int i = 0; i < paths[m].index.size(); i++)
////		{
////			int num = paths[m].index[i];
////			PLfile1 << "VRTX " << i << " " << P->points[num].x << " " << P->points[num].y << " " << P->points[num].z << "\n";
////		}
////		PLfile1 << "END";
////	}
////	PLfile1.close(); // 关闭文件流
////}

////////////void EMst::outputPLFile(p_rawPoint& P,  std::vector<int> startIndices, string name)
////////////{
////////////	string filename = path + name + "3c04-30-30-vw0.005.pl";
////////////	ofstream PLfile1(filename, ios::app);
////////////	PLfile1 << "GOCAD PLine 0.02sim\n";
////////////	PLfile1 << "HEADER{\n";
////////////	PLfile1 << "name:" << "FP" << name << "\n";
////////////	PLfile1 << "}\n";
////////////	PLfile1 << "ILINE\n";
////////////	//std::cout << "selected_edges size: " << selected_edges.size() << std::endl;
////////////	//std::cout << "selected_edges[1] size: " << selected_edges[1].size() << std::endl;
////////////	for (int m = 0; m < startIndices.size(); m++){
////////////		int i = startIndices[m];
////////////		PLfile1 << "VRTX " << m << " " << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << "\n";
////////////	}
////////////	for (int m = 0; m < startIndices.size() - 1; m++) {
////////////		PLfile1 << "SEG " << m << " " << m + 1 << "\n";
////////////	}
////////////	PLfile1 << "ENDGOCAD PLine 0.02sim\n";
////////////	PLfile1 << "\n"; // 在内层循环结束后添加换行符
////////////	//PLfile1 << "END\n";
////////////	PLfile1.close();
////////////}





//////////void EMst::outputPLFile(p_rawPoint P, const vector<trace>& paths, string name)
//////////{
//////////	string filename = path + name;
//////////	ofstream PLfile1(filename, ios::app); // 打开文件流
//////////	if (!PLfile1.is_open()) {
//////////		cerr << "无法打开文件：" << filename << endl;
//////////		return;
//////////	}
//////////	cout << "共有路径Path.size =" << paths.size() << endl;
//////////	for (int m = 0; m < paths.size(); m++)
//////////	{
//////////		PLfile1 << "GOCAD PLine \n";
//////////		PLfile1 << "HEADER{\n";
//////////		PLfile1 << "name:" << "FP" << m << "\n";
//////////		PLfile1 << "value:" << paths[m].weight << "\n";
//////////		PLfile1 << "}\n";
//////////		PLfile1 << "ILINE\n";
//////////		for (int i = 0; i < paths[m].index.size(); i++)
//////////		{
//////////			int num = paths[m].index[i];
//////////			PLfile1 << "VRTX " << i << " " << P->points[num].x << " " << P->points[num].y << " " << P->points[num].z << "\n";
//////////		}
//////////		PLfile1 << "END\n";
//////////	}
//////////	PLfile1.close(); // 关闭文件流
//////////}

////void EMst::outputPLFile(p_rawPoint P, const vector<trace>& paths, string name)
////{
////	string filename = name;
////	ofstream PLfile1(filename, ios::app);
////	cout << "共有路径Path.size =" << paths.size() << endl;
////	for (int m = 0; m < paths.size(); m++)
////	{
////		PLfile1 << "GOCAD PLine \n";
////		PLfile1 << "HEADER{\n";
////		PLfile1 << "name:" << "FP" << m << "\n";
////		PLfile1 << "value:" << paths[m].weight << "\n";
////		PLfile1 << "}\n";
////		PLfile1 << "ILINE\n";
////		for (int i = 0; i < paths[m].index.size(); i++)
////		{
////			int num = paths[m].index[i];
////			PLfile1 << "VRTX " << i << " " << P->points[num].x << " " << P->points[num].y << " " << P->points[num].z << "\n";
////		}
////		PLfile1 << "END\n";
////	}
////	PLfile1.close();
////}



//之前备份的代码，勿删

////////////EMst::EMst(string allPointCloudPath, bool flag)
////////////{
////////////	if (flag) {
////////////		pcl::io::loadPCDFile(allPointCloudPath, this->mypoint);
////////////		cout << "点云数据输入完成." << endl;
////////////	}
////////////	else {
////////////		// 在构造函数中加载两个点云文件
////////////		pcl::io::loadPCDFile(allPointCloudPath, this->raw);
////////////		logger("点云数据输入完成.");
////////////		int Point = this->raw.size();
////////////		cout << "全部点的数量" << Point << "：" << endl;
////////////		for (int i = 0; i < Point; i++) {
////////////			float p0 = 0.5f *(this->raw.points[i].intensity + this->raw.points[i].gvalue);
////////////			float p1 = 1 - p0;
////////////			if (p0 * p1 < 0)logger("错误发生在第", i, "行");
////////////			else {
////////////				RawPointType p;
////////////				p.x = this->raw.points[i].x;
////////////				p.y = this->raw.points[i].y;
////////////				p.z = this->raw.points[i].z;
////////////				p.intensity = this->raw.points[i].intensity;
////////////				p.gvalue = this->raw.points[i].gvalue;
////////////				p.Tlambda0 = this->raw.points[i].Tlambda0;
////////////				p.Tlambda1 = this->raw.points[i].Tlambda1;
////////////				p.Tlambda2 = this->raw.points[i].Tlambda2;
////////////				//this->raw.points.push_back(p);
////////////			}
////////////		}
////////////	}
////////////	//全部点云数据构图
////////////	int Point = this->raw.size();
////////////	//cout << "点云个数Point" << Point <<": " << endl;
////////////	//define the points in KNN graph
////////////	for (int i = 0; i < Point; i++) {
////////////		float p0 = 0.5f *(this->raw.points[i].intensity + this->raw.points[i].gvalue);
////////////		ProLine = p0;
////////////		if (ProLine > threshold_p) this->rawknnPoint.push_back(this->raw.points[i]);
////////////	}
////////////	int rawknnPointSIZE = this->rawknnPoint.size();
////////////
////////////	this->Graph.adjList = new AdjList[rawknnPointSIZE];
////////////	this->Graph.numVertsxes = rawknnPointSIZE;
////////////
////////////	//cout << "点云个数rawknnPoint" << ": " << rawknnPointSIZE  << endl;//点云个数Point81366 int NumknnPoint = this->rawknnPoint.size();
////////////	GraphAdjList *G = &this->Graph;
////////////	G->adjList = new AdjList[rawknnPointSIZE];
////////////	G->numVertsxes = rawknnPointSIZE;
////////////	for (int i = 0; i < rawknnPointSIZE; i++) {
////////////		G->adjList[i].data = i;
////////////		G->adjList[i].degree = 0;
////////////		G->adjList[i].firstedge = NULL;
////////////	}
////////////}

////////////void MiniTree_Prim(GraphAdjList G, int v, GraphAdjList * result)
////////////{
////////////	int i, j;
////////////	bool *visit = new bool[G.numVertsxes];
////////////	closedge *closedges = new closedge[G.numVertsxes];
////////////	initial_graph(&G, result);
////////////	//初始化closedge
////////////	for (i = 0; i < G.numVertsxes; i++)
////////////	{
////////////		visit[i] = false;
////////////		closedges[i].vex = v;
////////////		closedges[i].lowcost = Inf;
////////////	}
////////////	visit[v] = true;
////////////	EdgeNode *p, *tmp;
////////////	p = G.adjList[v].firstedge;
////////////	while (p)
////////////	{
////////////		closedges[p->adjvex].lowcost = p->weight;
////////////		p = p->next;
////////////	}
////////////	//找出closedge中最小的边，输出并进行更新
////////////	for (j = 0; j < G.numVertsxes; j++)
////////////	{
////////////		double min = Inf;
////////////		int t = Inf;
////////////		for (i = 0; i < G.numVertsxes; i++)
////////////		{
////////////			if (closedges[i].lowcost < min && visit[i] == false)
////////////			{
////////////				min = closedges[i].lowcost;
////////////				t = i;
////////////			}
////////////		}
////////////		if (t != Inf) {
////////////			//printf("%d,%d,%f\n", closedges[t].vex, t,closedges[t].lowcost);
////////////			int star = closedges[t].vex;
////////////			int to = t;
////////////			visit[t] = true;
////////////			if (result->adjList[t].firstedge == NULL) {
////////////				result->adjList[t].firstedge = make_node(closedges[t].vex, closedges[t].lowcost);
////////////				if (result->adjList[closedges[t].vex].firstedge == NULL) {
////////////					result->adjList[closedges[t].vex].firstedge = make_node(t, closedges[t].lowcost);
////////////				}
////////////				else {
////////////					tmp = result->adjList[closedges[t].vex].firstedge;
////////////					while (tmp->next != NULL)
////////////						tmp = tmp->next;
////////////					tmp->next = make_node(t, closedges[t].lowcost);
////////////				}
////////////			}
////////////			else {
////////////				tmp = result->adjList[t].firstedge;
////////////				while (tmp->next != NULL)
////////////					tmp = tmp->next;
////////////				tmp->next = make_node(closedges[t].vex, closedges[t].lowcost);
////////////				if (result->adjList[closedges[t].vex].firstedge == NULL) {
////////////					result->adjList[closedges[t].vex].firstedge = make_node(t, closedges[t].lowcost);
////////////				}
////////////				else {
////////////					tmp = result->adjList[closedges[t].vex].firstedge;
////////////					while (tmp->next != NULL)
////////////						tmp = tmp->next;
////////////					tmp->next = make_node(t, closedges[t].lowcost);
////////////				}
////////////			}
////////////			//更新 
////////////			p = G.adjList[t].firstedge;
////////////			while (p)
////////////			{
////////////				if (closedges[p->adjvex].lowcost > p->weight)
////////////				{
////////////					closedges[p->adjvex].lowcost = p->weight;
////////////					closedges[p->adjvex].vex = t;
////////////				}
////////////				p = p->next;
////////////			}
////////////
////////////		}
////////////		else {
////////////
////////////			logger("The graph whose initial point is ", v, " has been over!");
////////////			//cout << "The graph whose initial point is " << v << " has been over!" << endl;
////////////			break;
////////////		}
////////////	}
////////////	delete[] closedges;
////////////	delete[] visit;
////////////}
////////////void findline(GraphAdjList * tree1, int index, vector<EndNode> &line) {
////////////	save *b = new save[tree1->adjList[index].degree];
////////////	if (tree1->adjList[index].degree >= 3)
////////////	{
////////////		int i = 0;
////////////		EdgeNode* p;
////////////		p = tree1->adjList[index].firstedge;
////////////		while (p)
////////////		{
////////////			int depth = 0;
////////////			b[i].dirdata.push_back(index);
////////////			tree1->adjList[index].visit = true;
////////////			pruning4_plus(tree1, p->adjvex, &depth, b[i].dirdata);
////////////			b[i].depthnum = depth;
////////////			++i;
////////////			p = p->next;
////////////		}
////////////
////////////	}
////////////	for (int de = 0; de < tree1->adjList[index].degree; de++)
////////////	{
////////////		EndNode temp;
////////////		if (b[de].dirdata.size())
////////////		{
////////////			temp.lineorder = b[de].dirdata;
////////////			//尾结点：度为1 或者 度为大于等于3的结点，因此需要记录下来；
////////////			int end = b[de].dirdata.back();
////////////			temp.degree = tree1->adjList[end].degree;
////////////			temp.data = end;
////////////			line.push_back(temp);
////////////		}
////////////
////////////	}
////////////}

	//////std::vector<std::vector<float>>graph;
	////////无向图分配内存
	//////graph.resize(P->points.size());
	//////for (auto &p : graph)
	//////{
	//////	p.resize(P->points.size());
	//////}
	//////pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	//////pcl::copyPointCloud(*P, *cloud);
	//////pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	//////kdtree.setInputCloud(cloud);
	//////kdtree.setSortedResults(true);
	//////int point_size = P->points.size();//获取点云中的数量
	//////std::vector<ProInfo>vProInfo;
	//////vProInfo.reserve(point_size * 20);//20是寻找20个邻域
	//////for (int i = 0; i < point_size; i++)
	//////{
	//////	pcl::PointXYZ curPt;
	//////	curPt.x = P->points[i].x;
	//////	curPt.y = P->points[i].y;
	//////	curPt.z = P->points[i].z;
	//////	vector<int>index(k + 1);
	//////	vector<float>dist(k + 1);
	//////	//20个邻域全局索引在index
	//////	kdtree.nearestKSearch(curPt, k + 1, index, dist);
	//////	for (int j = 1; j < k + 1; j++)
	//////	{
	//////		ProInfo info;
	//////		info.start = i;
	//////		info.end = index[j];
	//////		info.dist = dist[j];
	//////		info.pro_value = P->points[i].p1 / (P->points[index[j]].p0 * sqrt(dist[j]));
	//////		graph[i][index[j]] = P->points[i].p1 * P->points[index[j]].p1 * exp(-dist[j]);
	//////		vProInfo.emplace_back(info);
	//////	}
	//////}
	////////排序
	//////std::sort(vProInfo.begin(), vProInfo.end(), [](ProInfo info1, ProInfo info2) {return info1.pro_value > info2.pro_value; });//降序排序
	////////求初始权重
	//////std::vector<EdgeInfo>vEdge;
	//////for (int i = 0; i < Pa; i++)
	//////{
	//////	EdgeInfo cur_edge;
	//////	cur_edge.start = vProInfo[i].start;
	//////	cur_edge.end = vProInfo[i].end;
	//////	cur_edge.wight_value = P->points[vProInfo[i].start].p1 * P->points[vProInfo[i].end].p1 * exp(-vProInfo[i].dist);
	//////	vEdge.push_back(cur_edge);
	//////}
	//粒子搜索后加入的新边

	//////std::vector<std::vector<EdgeInfo>>all_edges;//T次递增后所有的粒子集合
	////////all_edges.push_back(vEdge);//加入初始粒子
	//////for (int i = 0; i < Pa; i++)
	//////{
	//////	std::vector<EdgeInfo> new_edge;
	//////	new_edge.push_back(vEdge[i]); // 将 vEdge 中的每条边作为一个元素添加到 new_edge 中
	//////	all_edges.push_back(new_edge); // 将 new_edge 添加到 all_edges 中
	//////}
	//////for (int iter = 0; iter < T; iter++)
	//////{
	//////	std::vector<EdgeInfo> new_edges;
	//////	for (int i = 0; i < Pa; i++)
	//////	{
	//////		pcl::PointXYZ cur_end;
	//////		cur_end.x = P->points[all_edges[i].back().end].x;
	//////		cur_end.y = P->points[all_edges[i].back().end].y;
	//////		cur_end.z = P->points[all_edges[i].back().end].z;
	//////		std::vector<int> index(k + 1);
	//////		std::vector<float> dist(k + 1);
	//////		// 20个邻域全局索引在index
	//////		kdtree.nearestKSearch(cur_end, k + 1, index, dist);
	//////		EdgeInfo new_edge;
	//////		for (int j = 1; j < k + 1; j++)
	//////		{
	//////			float cur_widht = P->points[all_edges[i].back().end].p1 * P->points[index[j]].p1 * exp(-dist[j]);
	//////			if (cur_widht > new_edge.wight_value)
	//////			{
	//////				new_edge.start = all_edges[i].back().end;
	//////				new_edge.end = index[j];
	//////				new_edge.wight_value = cur_widht;
	//////			}
	//////		}
	//////		all_edges[i].push_back(new_edge); // 将新的边对象添加到对应的 all_edges[i] 中
	//////	}
	//////}
	//////std::cout << "all_edges size: " << all_edges.size() << std::endl;
	//////std::cout << "all_edges[1] size: " << all_edges[1].size() << std::endl;
	////////存储Pa个权重和
	//////std::vector<float>all_wights;
	//////for (int i = 0; i < Pa; i++)
	//////{
	//////	float cur_sum_widght = 0.0;
	//////	for (auto& edge : all_edges[i])
	//////	{
	//////		cur_sum_widght += edge.wight_value;
	//////	}
	//////	all_wights.push_back(cur_sum_widght);
	//////}
	//////std::cout << "all_edges size: " << all_edges.size() << std::endl;
	////////计算所有权重和
	//////float sum_widght = std::accumulate(all_wights.begin(), all_wights.end(), 0.0);
	//////// 权重归一化并赋值
	//////for (int i = 0; i < Pa; i++) {
	//////	all_edges[i][0].wight_value = all_wights[i] / sum_widght;
	//////}
	//////// 1. 比较 all_edges[i][0].wight_value 的值，返回最大的120个 all_edges
	//////std::partial_sort(all_edges.begin(), all_edges.begin() + 500, all_edges.end(), [](const std::vector<EdgeInfo>& a, const std::vector<EdgeInfo>& b) {
	//////	return a[0].wight_value > b[0].wight_value;
	//////});
	//////// 2.获取前500个最大的 all_edges
	//////std::vector<std::vector<EdgeInfo>> selected_edges(all_edges.begin(), all_edges.begin() + 500);
	//////// 3.调用输出函数，输出pl文件
	//////std::cout << "selected_edges size: " << selected_edges.size() << std::endl;
	//////for (int p = 0; p < selected_edges.size(); p++) {
	//////	outputPLFile(P, selected_edges, p, "selected_edges_10_try_300_500_321_");
	//////}

//////////void EMst::outputPLFile(p_Mypoint& P, const vector<trace>& paths, string name)
//////////{
//////////	string filename = name;
//////////	ofstream PLfile1(filename, ios::app);
//////////	PLfile1 << "GOCAD PLine 0.02sim\n";
//////////	PLfile1 << "HEADER{\n";
//////////	PLfile1 << "name:" << "FP" << name << "\n";
//////////	PLfile1 << "}\n";
//////////	PLfile1 << "ILINE\n";
//////////	std::cout << "selected_edges[particleIndex] size: " << selected_edges[particleIndex].size() << std::endl;
//////////	for (int m = 0; m < 301; m++) {//selected_edges[particleIndex].size()
//////////		int i = selected_edges[particleIndex][m].start;
//////////		int j = selected_edges[particleIndex][m].end;
//////////		PLfile1 << "VRTX " << 2 * m << " " << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << "\n";
//////////		PLfile1 << "VRTX " << 2 * m + 1 << " " << P->points[j].x << " " << P->points[j].y << " " << P->points[j].z << "\n";
//////////	}
//////////	for (int m = 0; m < 300; m++) {
//////////		PLfile1 << "SEG " << 2 * m << " " << 2 * m + 1 << "\n";
//////////	}
//////////	PLfile1 << "\n"; // 在内层循环结束后添加换行符
//////////	//PLfile1 << "ENDGOCAD PLine 0.02sim\n";
//////////	PLfile1 << "END\n";
//////////	PLfile1.close();
//////////}

////////////void EMst::outputPLFile(p_rawPoint& P,  std::vector<int> startIndices, string name)
////////////{
////////////	string filename = path + name + "3c04-30-30-vw0.005.pl";
////////////	ofstream PLfile1(filename, ios::app);
////////////	PLfile1 << "GOCAD PLine 0.02sim\n";
////////////	PLfile1 << "HEADER{\n";
////////////	PLfile1 << "name:" << "FP" << name << "\n";
////////////	PLfile1 << "}\n";
////////////	PLfile1 << "ILINE\n";
////////////	//std::cout << "selected_edges size: " << selected_edges.size() << std::endl;
////////////	//std::cout << "selected_edges[1] size: " << selected_edges[1].size() << std::endl;
////////////	for (int m = 0; m < startIndices.size(); m++){
////////////		int i = startIndices[m];
////////////		PLfile1 << "VRTX " << m << " " << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << "\n";
////////////	}
////////////	for (int m = 0; m < startIndices.size() - 1; m++) {
////////////		PLfile1 << "SEG " << m << " " << m + 1 << "\n";
////////////	}
////////////	PLfile1 << "ENDGOCAD PLine 0.02sim\n";
////////////	PLfile1 << "\n"; // 在内层循环结束后添加换行符
////////////	//PLfile1 << "END\n";
////////////	PLfile1.close();
////////////}

//// 求点云数据的五维向量平均值（√）
//std::array<long double, 5> compute5DVectorMean(const p_rawPoint& cloud) {
//	std::array<long double, 5> mean = { 0.0 };
//	if (cloud->empty()) {
//		std::cout << "点云数据异常，为空值" << std::endl;
//		return mean;  // 或者返回一个默认的全0向量
//	}
//	for (const auto& point : cloud->points) {
//		mean[0] += point.intensity;
//		mean[1] += point.gvalue;
//		mean[2] += point.Tlambda0;
//		mean[3] += point.Tlambda1;
//		mean[4] += point.Tlambda2;
//	}
//	const long double size = static_cast<long double>(cloud->size());
//	for (auto& val : mean) {
//		val /= size;
//	}
//	return mean;
//}
//// 计算协方差矩阵（√）
//Eigen::Matrix<long double, 5, 5> covariance(const p_rawPoint& Points, const std::array<long double, 5>& mean) {
//	Eigen::Matrix<long double, 5, 5> covariance = Eigen::Matrix<long double, 5, 5>::Zero();
//	if (!Points || Points->points.empty()) {
//		std::cout << "计算协方差矩阵时数据为空" << std::endl;
//		return covariance;
//	}
//	else {
//		const long double size = static_cast<long double>(Points->points.size());
//		for (const auto& point : Points->points) {
//			Eigen::Matrix<long double, 5, 1> X;
//			X << static_cast<long double>(point.intensity), static_cast<long double>(point.gvalue),
//				static_cast<long double>(point.Tlambda0), static_cast<long double>(point.Tlambda1),
//				static_cast<long double>(point.Tlambda2);
//			Eigen::Matrix<long double, 5, 1> diff = X - Eigen::Map<const Eigen::Matrix<long double, 5, 1>>(mean.data());
//			covariance += diff * diff.transpose();
//		}
//		covariance /= size; // 计算协方差矩阵
//		return covariance;
//	}
//}

// 初始化粒子的起点：应随机选择，计算每个点的角点的高斯概率密度分布的得分（√）
//std::vector<int> EMst::N1Score(p_rawPoint& P, int Pa, long double threshold) {
//	if (!P || P->points.empty()) {
//		std::cout << "点云数据为空" << std::endl;
//		return {}; // Return an empty vector
//	}
//	int point_size = P->points.size();
//	std::vector<std::vector<long double>> scores_of_first(point_size, std::vector<long double>(2, 0.0f));
//	std::vector<int> selected_points_indices;
//	array<long double, 5> meanPVector = compute5DVectorMean(P);
//	Eigen::Matrix<long double, 5, 5> covariancePMatrix = covariance(P, meanPVector);
//	for (int j = 0; j < point_size; j++) {
//		std::array<long double, 5> firstVector = {
//			P->points[j].intensity,
//			P->points[j].gvalue,
//			P->points[j].Tlambda0,
//			P->points[j].Tlambda1,
//			P->points[j].Tlambda2
//		};
//		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> X(firstVector.data());
//		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> PMean(meanPVector.data());
//		Eigen::Matrix<long double, 5, 1> diff = X - PMean;
//		long double exponent = -0.5f * diff.transpose() * covariancePMatrix.inverse() * diff;
//		long double determinant = covariancePMatrix.determinant(); 		
//		long double constant = 1.0 / (pow(2 * 3.14159265359f, 5 / 2.0) * sqrt(covariancePMatrix.determinant()));		
//		long double P_density = constant * std::exp(exponent);
//		scores_of_first[j][0] = j; // 索引号
//		scores_of_first[j][1] = P_density;
//		//cout << "exponent:   " << exponent << endl;
//		//cout << "covariancePMatrix:\n" << covariancePMatrix << endl;
//		//cout << "covariancePMatrix.determinant():\n" << determinant << endl;
//		//cout << "P_density:   " << P_density << endl;
//		//cout << "constant:   " << constant << endl;
//		if (P_density > threshold) {
//			selected_points_indices.push_back(j);
//		}
//	}
//	// 随机打乱大于阈值的点索引
//	std::random_device rd;
//	std::mt19937 g(rd());
//	std::shuffle(selected_points_indices.begin(), selected_points_indices.end(), g);
//	// 截取前 Pa 个点作为初始点
//	std::vector<int> selected_initial_points(selected_points_indices.begin(), selected_points_indices.begin() + Pa);
//	return selected_initial_points;
//}

//|||||||||||||||||||||||+++++++++++++++++++++计算运动模型扩增的总得分++++++++++++++++++++++++++++++||||||||||||||||||||||||
//std::vector<std::vector<long double>> EMst::Score(pcl::KdTreeFLANN<pcl::PointXYZ> kdtree,int last_start_index, int candidate_start_index, pcl::PointXYZ candidate_start, p_rawPoint P, p_rawPoint C, p_rawPoint B, int iter, int k, std::vector<int> index, std::vector<float> dist)
//{//------------------------0（√）------------------------输入原始数据、构建K近邻搜索---------------------------------
//	std::vector<std::vector<long double>> scores_array(20, std::vector<long double>(2));// 20行2列的二维数组，存储0索引号和1对应的 score 值
//	long double w1 = 0.7;
//	long double a1 = exp(-1 * w1 * iter);
//	long double a2 = 1 - a1; 
//	constexpr long double pi = 3.14159265359f;
//	int dimensions = 5; // 维度数
//	std::array<long double, 5> meanCVector = compute5DVectorMean(C);
//	std::array<long double, 5> meanBVector = compute5DVectorMean(B);
//	Eigen::Matrix<long double, 5, 5> covarianceCMatrix = covariance(C, meanCVector);
//	long double determinantC = covarianceCMatrix.determinant();//直接计算答案是不对的->是对的，在matlab中计算的反而不对，控制台上只输出了前几位，不是完整的元素
//	Eigen::Matrix<long double, 5, 5> covarianceBMatrix = covariance(B, meanBVector);
//	long double determinantB = covarianceBMatrix.determinant();
//	//cout << "covarianceCMatrix:\n" << covarianceCMatrix << endl;
//	//cout << "determinantC:  " << determinantC << endl;
//	//cout << "covarianceBMatrix:\n" << covarianceBMatrix << endl;
//	//cout << "determinantB:  " << determinantB << endl;
//	if (determinantC <= 0) {
//		cout << "协方差矩阵的行列式小于零，存在两个或者若干个随机变量之间有完全的函数关系" << endl;
//		covarianceCMatrix += Eigen::Matrix<long double, 5, 5>::Identity() * 0.001; // 正则化协方差矩阵
//		determinantC = covarianceCMatrix.determinant(); // 重新计算行列式
//	}
//	int k_index = index.size();
////	cout << "k:  " << k << endl;
////	cout << "k_index:  " << k_index << endl;
//	for (int j = 0; j < k - 1 ; j++) {//多次（超过20次）计算出来的X向量没有一个重合的，不合理，合理的，每次初始化的点都不一样
//		int candidate_end_index = index[j];
//		long double distance = dist[j];
//		//cout << "索引为： " << j << "距离为：	" << distance << endl;
//		std::array<long double, 5> candidate_end = {
//			P->points[candidate_end_index].intensity,
//			P->points[candidate_end_index].gvalue,
//			P->points[candidate_end_index].Tlambda0,
//			P->points[candidate_end_index].Tlambda1,
//			P->points[candidate_end_index].Tlambda2
//		};
//		//------------------------1（√）----------------------------计算Pzx-------------------------------------------------
//		//------------------------1-1（√）-------N1计算++++++++++++++++++++++++++++++++++高斯分布中，μ和协方差矩阵由角点数据确定，只有X是邻域内的点
//		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> X(candidate_end.data());
//		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> CMean(meanCVector.data());
//		Eigen::Matrix<long double, 5, 1> diff = X - CMean;//total_covarianceC /= k; // 取所有邻域协方差的均值,不需要均值，均值的意义是什么！！！！！！
//		long double Cexponent = -0.5f * diff.transpose() * covarianceCMatrix.inverse() * diff;
//		long double QUexponent = diff.transpose() * covarianceCMatrix.inverse() * diff;
//		long double constant = 1.0 / (pow(2 * pi, dimensions / 2.0) * sqrt(determinantC));
//		long double C_density = constant * std::exp(Cexponent);
//		//------------------------（√）N2计算---------------------
//		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> XB(candidate_end.data());
//		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> BMean(meanBVector.data());
//		Eigen::Matrix<long double, 5, 1> diffB = XB - BMean;
//		long double Bexponent = -0.5f * diff.transpose() * covarianceBMatrix.inverse() * diff;
//		long double QUBexponent = diff.transpose() * covarianceBMatrix.inverse() * diff;
//		long double constantB = 1.0 / (pow(2 * pi, dimensions / 2.0) * sqrt(determinantB));
//		long double B_density = constant * std::exp(Bexponent);
//		if (j==1){
//			//cout << "CMean：\n" << CMean << endl;
//			//cout << "BMean：\n" << BMean << endl;
//			//cout << "determinantC:  " << determinantC << endl;
//			//cout << "determinantB:  " << determinantB << endl;
//			//cout << "constant:  " << constant << endl;
//			//cout << "constantB:  " << constantB << endl;
//		}
////		cout << "X：\n" << X << endl;
//		//cout << "diff：\n" << diff << endl;
//		//cout << "exponent:  " << Cexponent << endl;
//		//cout << "C_density:  " << C_density << endl;
//		//cout << "             B           " << endl;
//		//cout << "diffB：\n" << diffB << endl;
//		//cout << "Bexponent:  " << Bexponent << endl;
//		//cout << "B_density:  " << C_density << endl;
//
//		// 计算总的 Pzx
//		long double Pzx = a1 * C_density + a2 * B_density;
////		cout << "<<<<<<<<<<<<<<<<索引号:	" << j <<"		" << ": Total Pzx:   " << Pzx << endl;
//
//		//-------------------------2（√）---------------------------计算Pxx-------------------------------------------------// 声明一个变量来存储 PXX 的总和
//		long double sum_Pxx = 0.0;
//		long double start_x = P->points[candidate_start_index].x;
//		long double start_y = P->points[candidate_start_index].y;
//		long double start_z = P->points[candidate_start_index].z;
//		long double last_x = P->points[last_start_index].x;
//		long double last_y = P->points[last_start_index].y;
//		long double last_z = P->points[last_start_index].z;
//		long double cur_x = P->points[index[j]].x;
//		long double cur_y = P->points[index[j]].y;
//		long double cur_z = P->points[index[j]].z;
//		long double A1 = start_x - last_x;
//		long double A2 = start_y - last_y;
//		long double A3 = start_z - last_z;
//		long double D1 = cur_x - start_x;
//		long double D2 = cur_y - start_y;
//		long double D3 = cur_z - start_z;
//		////cout << "start_x：   " << start_x << endl;
//		////cout << "start_y：   " << start_y << endl;
//		////cout << "start_z：   " << start_z << endl;
//		////cout << "last_x：    " << last_x << endl;
//		////cout << "last_y：    " << last_y << endl;
//		////cout << "last_z：    " << last_z << endl;
//		////cout << "cur_x：     " << cur_x << endl;
//		////cout << "cur_y：     " << cur_y << endl;
//		////cout << "cur_z：     " << cur_z << endl;
//		////cout << "A1：	     " << A1 << endl;
//		////cout << "A2：	     " << A2 << endl;
//		////cout << "A3：	     " << A3 << endl;
//		////cout << "D1：	     " << D1 << endl;
//		////cout << "D2：	     " << D2 << endl;
//		////cout << "D3：	     " << D3 << endl;
//		//long double mo = std::sqrt(std::pow((A2 * D3 - A3 * D2), 2) + std::pow((A3 * D1 - A1 * D3), 2) + std::pow((A1 * D2 - A2 * D1), 2));// 计算 mo 值
//		long double dotProduct = A2 * D3 + A3 * D1 + A1 * D2;
////		cout << "	dotProduct:		" << dotProduct << endl;
//		long double w2 = 0.5;
//		for (int n = 1; n < k + 1; n++) {
//			long double dotProduct_sum = A2 * D3 + A3 * D1 + A1 * D2;
//			//long double mo = std::sqrt(std::pow((A2 * D3 - A3 * D2), 2) + std::pow((A3 * D1 - A1 * D3), 2) + std::pow((A1 * D2 - A2 * D1), 2));
//			long double pxx = exp(w2 * dotProduct_sum);//去掉负号后比较符合整体的结果
//			sum_Pxx += pxx;
//			//cout << "dotProduct_sum：   " << dotProduct_sum << endl;
//			//cout << "pxx：   " << pxx << endl;
//			//cout << "sum_Pxx：  " << sum_Pxx << endl;
//		}
//		long double Pxx = exp(w2 * dotProduct) / sum_Pxx;
//		long double Score = Pzx * Pxx;
////		cout << "	最后计算结果Pxx：  " << Pxx << endl;
////		cout << "		索引号：	" << j <<"		 "<< "最后计算结果Score：  " << Score << endl;
//		/*scores_array[j - 1][0] = index[j];*/
//		scores_array[j][0] = j;
//		scores_array[j][1] = Score;
//	}
////	cout << "第 " << iter << " 次全部计算完毕" << endl;
//	return scores_array;
//}

//++++++++++++++++++++++|||||||||||||||||||||||根据离散分布律的随机采样函数|||||||||||||||||||||||||||||||||||++++++++++++++++++++++++++++++++++++++
//int EMst::Sampling(std::vector<std::vector<long double>>& scores_array) {
//	// 计算累积概率
//	std::vector<long double> cumulativeProbabilities(20);
//	long double totalScore = 0.0;
//	int emptyCount = 0; // 记录空值个数
//	for (size_t i = 0; i < 20; ++i) {
//		if (scores_array[i][1] == 0.0) { // 假设空值为0
//			emptyCount++;
//		}
//		else {
//			totalScore += scores_array[i][1];
//		}
//	}
//	if (emptyCount == 20) { // 如果全部为空值，生成一个随机索引，初始化后第一次扩增的过程是全为空值的，希望可以找最近的
//		int randomIndex = 0;
////		std::cout << "索引号为:	 "<< randomIndex <<"得分值为:	"<< scores_array[randomIndex][0] << std::endl;
//		return scores_array[randomIndex][0];
//	}
//	for (size_t i = 0; i < 20; ++i) {
//		if (scores_array[i][1] == 0.0) {
//			cumulativeProbabilities[i] = 1.0 / 20; // 平均分布
//		}
//		else {
//			cumulativeProbabilities[i] = scores_array[i][1] / totalScore;
//			if (i > 0) {
//				cumulativeProbabilities[i] += cumulativeProbabilities[i - 1];
//			}
//		}
//	}
//	// 生成随机数
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_real_distribution<> dis(0.0, 1.0);
//	long double randomValue = dis(gen); // 生成0到1之间的随机数
//	int sampled_index = -1;
//	for (size_t i = 0; i < 20; ++i) {
//		if (randomValue <= cumulativeProbabilities[i]) {
//			sampled_index = i;
//			break;
//		}
//	}
//	if (sampled_index != -1) {
////		std::cout << "Sampling随机生成的局部索引号为: " << scores_array[sampled_index][0] << std::endl;
//		return scores_array[sampled_index][0];
//	}
//	else {
//		std::cerr << "Error: Unable to determine sampled index." << std::endl;
//		return -1; // 或者其他适合的错误代码
//	}
//}

//输出csv格式的函数（√）
//void EMst::writeToCSV(const p_rawPoint& P, const std::vector<EdgeInfo>& edges, string name) {
//	string filename = path + name + "6_lowbuild_point.csv";
//	ofstream csvFile(filename, ios::app);
//	if (csvFile.is_open()) {
//		csvFile << "x,y,z,weight\n";  // 写入 CSV 文件的标题行
//		for (const EdgeInfo& edge : edges) {
//			int startIndex = edge.start;
//			int endIndex = edge.end;
//			long double weight = edge.wight_value;
//			// 写入起点坐标和权值
//			csvFile << P->points[startIndex].x << "," << P->points[startIndex].y << "," << P->points[startIndex].z << "," << weight << "\n";
//			// 写入终点坐标和权值
//			csvFile << P->points[endIndex].x << "," << P->points[endIndex].y << "," << P->points[endIndex].z << "," << weight << "\n";
//		}
//		csvFile.close();
//	}
//	else {
//		std::cerr << "Unable to open or create CSV file: " << filename << std::endl;
//	}
//}

////计算行列式的自定义·函数
//static long double CalcDeterminant(vector<vector<long double>> Matrix)
//{
//	long double det = 0; // the determinant value will be stored here
//	if (Matrix.size() == 1)
//	{
//		return Matrix[0][0]; // no calculation needed
//	}
//	else if (Matrix.size() == 2)
//	{
//		det = (Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0]);
//		return det;
//	}
//	else
//	{
//		for (long double p = 0; p < Matrix[0].size(); p++)
//		{
//			vector<vector<long double>> TempMatrix; // to hold the shaped matrix;
//			for (long double i = 1; i < Matrix.size(); i++)
//			{
//				vector<long double> TempRow;
//				for (long double j = 0; j < Matrix[i].size(); j++)
//				{
//					if (j != p)
//					{
//						TempRow.push_back(Matrix[i][j]);//add current cell to TempRow 
//					}
//				}
//				if (TempRow.size() > 0)
//					TempMatrix.push_back(TempRow);
//			}
//			det = det + Matrix[0][p] * pow(-1, p) * CalcDeterminant(TempMatrix);
//		}
//		//cout << "det:  " << det << endl;
//		return det;
//	}
//}
//// 检查两个矩阵是否相同
//bool areMatricesEqual(const std::vector<std::vector<long double>>& matrix1, const std::vector<std::vector<long double>>& matrix2) {
//	if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
//		return false;  // 矩阵尺寸不同，直接返回 false
//	}
//
//	for (size_t i = 0; i < matrix1.size(); ++i) {
//		for (size_t j = 0; j < matrix1[0].size(); ++j) {
//			if (matrix1[i][j] != matrix2[i][j]) {
//				cout << "matrix1[i][j] != matrix2[i][j]" << i << j << endl;
//			}
//		}
//		return false;  // 元素不同，返回 false
//	}
//
//	return true;  // 矩阵相同
//}

//计算协方差矩阵，是由角点数据确定的，和邻域内的点无关
//////Eigen::Matrix<long double, 5, 5> total_covarianceC = Eigen::Matrix<long double, 5, 5>::Zero(); // 总协方差矩阵
//////for (int n = 1; n < k + 1; n++) {
//////	int neighbour_index = index[n];
//////	std::array<long double, 5> CneighbourVector = {
//////		P->points[neighbour_index].intensity,
//////		P->points[neighbour_index].gvalue,
//////		P->points[neighbour_index].Tlambda0,
//////		P->points[neighbour_index].Tlambda1,
//////		P->points[neighbour_index].Tlambda2
//////	};
//////	Eigen::Map<const Eigen::Matrix<long double, 5, 1>> Cneighbour(CneighbourVector.data());
//////	Eigen::Matrix<long double, 5, 1> Cneighbour_diff = Cneighbour - CMean;
//////	total_covarianceC += Cneighbour_diff * Cneighbour_diff.transpose();
//////}

////	int thread_ = 2 * omp_get_num_procs();
////	std::cout << thread_ << "s" << endl;
////#pragma omp parallel for num_threads(2* omp_get_num_procs() - 1)

////////////////std::vector<int> startIndices;//是一个存储了起点索引号的容器
////////////////for (int p = 0; p < Pa; p++) {
////////////////	const std::vector<EdgeInfo>& edges = all_edges[p];
////////////////	startIndices.clear(); // 清空 startIndices 向量
////////////////	for (const EdgeInfo& edge : edges) {
////////////////		int startIndex = edge.start; // 获取边的起点索引号
////////////////		startIndices.push_back(startIndex); // 将起点索引号添加到当前粒子的起点索引号向量中
////////////////	}
////////////////	outputPLFile(P, startIndices, "selected_edges_935_1000_corner_19_3_");
////////////////}

//////std::vector<std::vector<long double>> EMst::Score(pcl::KdTreeFLANN<pcl::PointXYZ> kdtree, int last, int candidate_start, pcl::PointXYZ cur_end, p_rawPoint P, p_rawPoint C, p_rawPoint B, int iter, int k)
//////{
//////	//------------------------0------------------------输入原始数据、构建K近邻搜索---------------------------------
//////	std::vector<int> index;
//////	std::vector<float> dist;
//////	kdtree.nearestKSearch(cur_end, k + 1, index, dist);
//////	std::vector<std::vector<long double>> scores_array(20, std::vector<long double>(2)); // 20行2列的二维数组，存储索引号和对应的 score 值
//////	for (int j = 1; j < k + 1; j++) {
//////		int candidate_end = index[j];
//////		long double dis = dist[j];
//////		std::vector<long double> scores_of_neighbours;//存储 20 个邻域的得分值
//////		std::array<long double, 5> meanCVector = corner5mean(C);
//////		std::array<long double, 5> meanBVector = boundary5mean(B);
//////		std::array<long double, 5> curVector = {
//////			P->points[candidate_end].intensity,
//////			P->points[candidate_end].gvalue,
//////			P->points[candidate_end].Tlambda0,
//////			P->points[candidate_end].Tlambda1,
//////			P->points[candidate_end].Tlambda2
//////		};
//////		//------------------------1----------------------------计算Pzx-------------------------------------------------
//////		long double w1 = 0.7;
//////		//long double a1 = exp(-1 * w1 * (iter + 0.01));//为了防止第一次
//////		long double a1 = exp(-1 * w1 * iter);
//////		long double a2 = 1 - a1;
//////		//-------------------------------N1计算+++++++++++++++++++++++++++++++++++++
//////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> X(curVector.data());
//////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> CMean(meanCVector.data());
//////		Eigen::Matrix<long double, 5, 1> diff = X - CMean;
//////		Eigen::Matrix<long double, 5, 5> total_covarianceC = Eigen::Matrix<long double, 5, 5>::Zero(); // 总协方差矩阵
//////		long double total_determinant = 0.0; // 总协方差行列式之和
//////		long double total_exponent = 0.0; // 总指数之和
//////		for (int n = 1; n < k + 1; n++) {
//////			int neighbour_index = index[n];
//////			std::array<long double, 5> neighbourVector = {
//////				P->points[neighbour_index].intensity,
//////				P->points[neighbour_index].gvalue,
//////				P->points[neighbour_index].Tlambda0,
//////				P->points[neighbour_index].Tlambda1,
//////				P->points[neighbour_index].Tlambda2
//////			};
//////			Eigen::Map<const Eigen::Matrix<long double, 5, 1>> neighbour(neighbourVector.data());
//////			Eigen::Matrix<long double, 5, 1> neighbour_diff = neighbour - CMean;
//////			total_covarianceC += neighbour_diff * neighbour_diff.transpose();
//////		}
//////		total_covarianceC /= k; // 取所有邻域协方差的均值
//////		long double total_determinantC = total_covarianceC.determinant(); // 计算总协方差行列式
//////		if (total_determinantC <= 0) {
//////			cout << "GaussianDensityC 协方差矩阵可能不可逆或奇异，处理异常情况" << endl;
//////			total_covarianceC += Eigen::Matrix<long double, 5, 5>::Identity() * 0.001; // 正则化协方差矩阵
//////			total_determinantC = total_covarianceC.determinant(); // 重新计算行列式
//////		}
//////		constexpr long double pi = 3.14159265359f;
//////		long double exponent = -0.5f * diff.transpose() * total_covarianceC.inverse() * diff;
//////		int dimensions = 5; // 维度数
//////		long double constant = 1.0f / (std::sqrt(std::pow(2 * pi, dimensions) * total_determinantC));
//////		long double C_density = constant * std::exp(exponent);
//////		//////////////////////////////////////cout << "C_density:  " << C_density << endl;
//////		//------------------------N2计算---------------------
//////		std::array<long double, 5> curVectorB = {
//////			P->points[candidate_end].intensity,
//////			P->points[candidate_end].gvalue,
//////			P->points[candidate_end].Tlambda0,
//////			P->points[candidate_end].Tlambda1,
//////			P->points[candidate_end].Tlambda2
//////		};
//////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> XB(curVectorB.data());
//////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> BMean(meanBVector.data());
//////		Eigen::Matrix<long double, 5, 1> diffB = XB - BMean;
//////		Eigen::Matrix<long double, 5, 5> total_covarianceB = Eigen::Matrix<long double, 5, 5>::Zero(); // 总协方差矩阵
//////		for (int n = 1; n < k + 1; n++) {
//////			int neighbour_index = index[n];
//////			std::array<long double, 5> neighbourVector = {
//////				P->points[neighbour_index].intensity,
//////				P->points[neighbour_index].gvalue,
//////				P->points[neighbour_index].Tlambda0,
//////				P->points[neighbour_index].Tlambda1,
//////				P->points[neighbour_index].Tlambda2
//////			};
//////			Eigen::Map<const Eigen::Matrix<long double, 5, 1>> neighbourB(neighbourVector.data());
//////			Eigen::Matrix<long double, 5, 1> neighbour_diffB = neighbourB - BMean;
//////			total_covarianceB += neighbour_diffB * neighbour_diffB.transpose();
//////		}
//////		total_covarianceB /= k; // 取所有邻域协方差的均值
//////		long double determinantB = total_covarianceB.determinant();
//////		////////cout << "Total Covariance Matrix_B: " << endl;
//////		////////for (int i = 0; i < total_covarianceB.rows(); ++i) {
//////		////////	for (int j = 0; j < total_covarianceB.cols(); ++j) {
//////		////////		cout << total_covarianceB(i, j) << " ";
//////		////////	}
//////		////////	cout << endl;
//////		////////}
//////		Eigen::JacobiSVD<Eigen::Matrix<long double, 5, 5>> svdB(total_covarianceB, Eigen::ComputeThinU | Eigen::ComputeThinV);
//////		long double thresholdC = 1e-6; // 设置一个阈值，判断奇异性
//////		if (svdB.singularValues()(4) < thresholdC) {
//////			cout << "GaussianDensityBB 协方差矩阵可能奇异，处理异常情况" << endl;
//////			// 此时可能需要对矩阵进行调整或使用其他方法处理奇异性
//////		}
//////		determinantB = svdB.singularValues().prod(); // 重新计算行列式
//////		//////////////////////////////////////cout << "DeterminantBBBB of the covariance matrix: " << determinantB << endl;
//////		long double exponentB = -0.5f * diffB.transpose() * total_covarianceB.inverse() * diffB;
//////		long double constantB = 1.0f / (std::sqrt(std::pow(2 * pi, dimensions) * determinantB));
//////		long double B_density = constantB * std::exp(exponentB);
//////		//////////////////////////////////////cout << "B_density:   " << B_density << endl;
//////		// 计算总的 Pzx
//////		long double Pzx = a1 * C_density + a2 * B_density;
//////		//////////////////////////////////////cout << "Total Pzx:   " << Pzx << endl;
//////		//-------------------------2---------------------------计算Pxx-------------------------------------------------
//////		long double start_x = P->points[candidate_start].x;
//////		long double start_y = P->points[candidate_start].y;
//////		long double start_z = P->points[candidate_start].z;
//////		long double last_x = P->points[last].x;
//////		long double last_y = P->points[last].y;
//////		long double last_z = P->points[last].z;
//////		long double cur_x = P->points[index[j]].x;
//////		long double cur_y = P->points[index[j]].y;
//////		long double cur_z = P->points[index[j]].z;
//////		long double A1 = last_x - start_x;
//////		long double A2 = last_y - start_y;
//////		long double A3 = last_z - start_z;
//////		long double D1 = start_x - cur_x;
//////		long double D2 = start_y - cur_y;
//////		long double D3 = start_z - cur_z;
//////		//两个三维向量
//////		std::array<long double, 5> A = { A1, A2, A3 };
//////		std::array<long double, 5> D = { D1, D2, D3 };
//////		double mo = std::sqrt(std::pow((A2 * D3 - A3 * D2), 2) + std::pow((A3 * D1 - A1 * D3), 2) + std::pow((A1 * D2 - A2 * D1), 2));
//////		long double w2 = 0.5;
//////		long double sum_pxx = 0.0; // 初始化 px 总和为 0
//////		// 计算每个邻域内的 exp(w * mo) 并求和
//////		for (int n = 1; n < k + 1; n++) {
//////			double mo = std::sqrt(std::pow((A2 * D3 - A3 * D2), 2) + std::pow((A3 * D1 - A1 * D3), 2) + std::pow((A1 * D2 - A2 * D1), 2));
//////			long double pxx = exp(-w2 * mo);
//////			sum_pxx += pxx; // 将 pxx 加到总和中
//////		}
//////		long double Pxx = exp(w2 * mo) / sum_pxx; // 计算 Pxx
//////		//////////////////////////////////////cout << "Pxx:   " << Pxx << endl;
//////		long double Score = Pzx * Pxx;
//////		//////////////////////////////////////cout << "Score:   " << Score << endl;
//////		//-------------------------3---------------------将得分存入数组或向量----------------------------------
//////		scores_array[j - 1][0] = index[j]; // 存储索引号
//////		scores_array[j - 1][1] = Score;   // 存储对应的 score 值
//////	}
//////	return scores_array; // 返回二维数组
//////}
//////
////////根据离散分布律的随机采样函数
//////int EMst::Sampling(std::vector<std::vector<long double>>& scores_array) {
//////	// 计算累积概率
//////	std::vector<long double> cumulativeProbabilities(20);
//////	long double totalScore = 0.0;
//////	int emptyCount = 0; // 记录空值个数
//////	for (size_t i = 0; i < 20; ++i) {
//////		if (scores_array[i][1] == 0.0) { // 假设空值为0
//////			emptyCount++;
//////		}
//////		else {
//////			totalScore += scores_array[i][1];
//////		}
//////	}
//////	if (emptyCount == 20) { // 如果全部为空值，生成一个随机索引
//////		std::random_device rd;
//////		std::mt19937 gen(rd());
//////		std::uniform_int_distribution<> dis(0, 19);
//////		int randomIndex = dis(gen);
//////		////////////////////////std::cout << "随机生成的局部索引号为: " << scores_array[randomIndex][0] << std::endl;
//////		return scores_array[randomIndex][0];
//////	}
//////	for (size_t i = 0; i < 20; ++i) {
//////		if (scores_array[i][1] == 0.0) {
//////			cumulativeProbabilities[i] = 1.0 / 20; // 平均分布
//////		}
//////		else {
//////			cumulativeProbabilities[i] = scores_array[i][1] / totalScore;
//////			if (i > 0) {
//////				cumulativeProbabilities[i] += cumulativeProbabilities[i - 1];
//////			}
//////		}
//////	}
//////	// 生成随机数
//////	std::random_device rd;
//////	std::mt19937 gen(rd());
//////	std::uniform_real_distribution<> dis(0.0, 1.0);
//////	long double randomValue = dis(gen); // 生成0到1之间的随机数
//////	int sampled_index = -1;
//////	for (size_t i = 0; i < 20; ++i) {
//////		if (randomValue <= cumulativeProbabilities[i]) {
//////			sampled_index = i;
//////			break;
//////		}
//////	}
//////	if (sampled_index != -1) {
//////		std::cout << "随机生成的局部索引号为: " << scores_array[sampled_index][0] << std::endl;
//////		return scores_array[sampled_index][0];
//////	}
//////	else {
//////		std::cerr << "Error: Unable to determine sampled index." << std::endl;
//////		return -1; // 或者其他适合的错误代码
//////	}
//////}

////之前以角点开始的代码
	////for (int i = 0; i < Pa; i++){
	////	for (int j = 0; j < point_size; j++){
	////		pcl::PointXYZ first;
	////		first.x = P->points[j].x;
	////		first.y = P->points[j].y;
	////		first.z = P->points[j].z;
	////		std::vector<std::vector<float>> scores_first_array = N1Score(P, C, B, k);
	////		int sampled_first_index = Sampling(scores_first_array);//返回一个candidate_end的全局索引号
	////		cout << "sampled_first_index:      " << sampled_first_index << endl;
	////	}
	////	////if (P->points[i].iscorner == 1){
	////	////	EdgeInfo first_edge;
	////	////	first_edge.start = i;//全局索引
	////	////	first_edge.end = i;
	////	////	first_edge.wight_value = NULL;
	////	////	cEdge.push_back(first_edge);
	////	////}
	////}
	////for (int i = 0; i < corner_size; i++){
	////	std::vector<EdgeInfo> temp_edge;
	////	temp_edge.push_back(cEdge[i]);
	////	all_edges.push_back(temp_edge);
	////}

//////std::vector<std::vector<float>> EMst::N1Score(p_rawPoint P, p_rawPoint C, int k, int Pa)
//////{//输入PC，输出全部点云的高斯分布的概率密度值，二维数组
//////	int point_size = P->points.size();
//////	int corner_size = C->points.size();
//////	std::vector<std::vector<float>>scores_of_first(point_size, vector<float>(2, 0.0f));//创建了一个二维的 vector，有 point_size 行，每行包含两个元素，初始值都是 0.0f
//////	array<float, 5> meanCVector = corner5mean(C);
//////	Eigen::Matrix<float, 5, 5> covarianceCMatrix = corner5covariance(C, meanCVector);
//////	for (int j = 0; j < point_size; j++){
//////		std::array<float, 5> firstVector = {
//////			P->points[j].intensity,
//////			P->points[j].gvalue,
//////			P->points[j].Tlambda0,
//////			P->points[j].Tlambda1,
//////			P->points[j].Tlambda2
//////		};
//////		Eigen::Map<const Eigen::Matrix<float, 5, 1>> X(firstVector.data());
//////		Eigen::Map<const Eigen::Matrix<float, 5, 1>> CMean(meanCVector.data());
//////		Eigen::Matrix<float, 5, 1> diff = X - CMean;
//////		float exponent = -0.5f * diff.transpose() * covarianceCMatrix.inverse() * diff;
//////		float constant = 1.0 / (pow(2 * 3.14159265359f, 5 / 2.0) * sqrt(covarianceCMatrix.determinant()));
//////		float C_density = constant * std::exp(exponent);
//////		scores_of_first[j][0] = j; // 索引号
//////		scores_of_first[j][1] = C_density;
//////	}
//////	return std::vector<std::vector<float>>();
//////}
//扩增过程中的总得分

//void EMst::graph_to_Graph(std::vector<std::vector<float>>&graph)
//{
//	EdgeNode* e;
//	for (int i = 0; i < graph.size(); i++)
//	{
//		for (int j = 0; j < graph[i].size(); j++)
//		{
//			e = (EdgeNode*)malloc(sizeof(EdgeNode));
//			e->adjvex = j;
//			e->weight = graph[i][j];//
//			e->next = this->Graph.adjList[i].firstedge;
//			this->Graph.adjList[i].firstedge = e;
//
//			e = (EdgeNode*)malloc(sizeof(EdgeNode));
//			e->adjvex = i;
//			e->weight = graph[i][j];
//			e->next = this->Graph.adjList[j].firstedge;
//			this->Graph.adjList[j].firstedge = e;
//		}
//	}
//}

//有重复点检查的粒子滤波函数1219
////void EMst::get_new_widght(p_rawPoint P, p_rawPoint C, p_rawPoint B, int k)
////{
////	iscorner(P, C);
////	std::vector<std::vector<float>>graph;//存储二维动态数组
////	graph.resize(P->points.size());	//无向图分配内存
////	for (auto &p : graph) {
////		p.resize(P->points.size());
////	}
////	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);//这行代码创建了一个名为 cloud 的指针，指向一个 pcl::PointCloud 类型的点云对象，该对象的点类型为 pcl::PointXYZ，用于存储三维空间中的点数据。
////	pcl::copyPointCloud(*P, *cloud);//将 P 指向的点云数据复制到 cloud 指向的点云对象中
////	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;//构建 Kd 树数据结构。
////	kdtree.setInputCloud(cloud);
////	kdtree.setSortedResults(true);//后续用到kdtree.nearestKSearch(curPt, k + 1, index, dist);会返回距离当前点 curPt 最近的 k + 1 个点的索引和距离。
////	//之后可以通过P->points[candidate_start].p1 来调用P,C中的五维向量
////	int point_size = P->points.size();//获取点云中的数量 3c-02.pcd  ：40683
////	int corner_size = C->points.size();
////	cout << "0-点云中的数量" << "：" << point_size << endl;
////	cout << "0-角点点云中的数量/粒子个数" << "：" << corner_size << endl;
////	//--1--用角点定义粒子的第一条边
////	std::vector<std::vector<EdgeInfo>>all_edges;//T次递增后所有的粒子集合
////	std::vector<EdgeInfo>cEdge;//cEdge 是存储以角点为起始点的935条初始粒子（边）
////	//从角点开始存储粒子结构
////	for (int i = 0; i < point_size; i++)
////	{
////		if (P->points[i].iscorner == 1)
////		{
////			EdgeInfo first_edge;
////			first_edge.start = i;//全局索引
////			first_edge.end = i;
////			first_edge.wight_value = NULL;
////			cEdge.push_back(first_edge);
////		}
////	}
////	for (int i = 0; i < corner_size; i++)
////	{
////		std::vector<EdgeInfo> temp_edge;
////		temp_edge.push_back(cEdge[i]);
////		all_edges.push_back(temp_edge);
////	}
////	//--2--粒子正式搜索后加入的新边
////	for (int iter = 0; iter < T; iter++) {
////		std::vector<EdgeInfo> new_edges;
////		std::vector<float> scores_of_neighbours;//存储 20 个邻域的得分值，离散分布律
////		for (int p = 0; p < corner_size; p++) {//p就是角点个数、粒子个数
////			if (!all_edges[p].empty()) {
////				pcl::PointXYZ cur_end;
////				cur_end.x = P->points[all_edges[p].back().end].x;
////				cur_end.y = P->points[all_edges[p].back().end].y;
////				cur_end.z = P->points[all_edges[p].back().end].z;
////				std::vector<int> index;
////				std::vector<float> dist;
////				kdtree.nearestKSearch(cur_end, k + 1, index, dist);
////				EdgeInfo new_edge;
////				new_edge.wight_value = 0;
////				float min_widht = std::numeric_limits<float>::max(); //初始化最小cur_widht为一个较大的值
////				int du = 0;
////				for (int j = 1; j < k + 1; j++) {
////					int candidate_start = all_edges[p].back().end;//上一条边终点的索引号
////					int candidate_end = index[j];
////					float dis = dist[j];
////					//-------------------------------得分函数-----------------------------------
////					//1220目前需要：得分函数输出正确的scores_of_neighbours，之后传输给sampling函数，输出一个全局索引号作为下一个点的索引号
////					Score(candidate_start, candidate_end, new_edge, P, C, B, iter, p, index, dis);
////					Sampling(scores_of_neighbours);//返回一个candidate_end的全局索引号
////
////					//在邻域20个都遍历完之后，检查是否与之前选点相同
////					//bool isDuplicate = false;// 检查新选的边的终点是否与已选择的边的起点或终点相同
////					//for (const auto& edge : all_edges[p]){
////					//	if (edge.end == candidate_end || edge.start == candidate_end){
////					//		du++;
////					//		isDuplicate = true;
////					//		break;
////					//	}
////					float cur_widht = dist[j];
////					min_widht = cur_widht;
////					new_edge.start = candidate_start;
////					new_edge.end = candidate_end;
////					new_edge.wight_value = cur_widht;
////				}
////				////if (du == k){
////				////	new_edge.start = all_edges[p].back().start;
////				////	new_edge.end = all_edges[p].back().end;
////				////}
////				////else{
////				////	for (int j = 1; j < k + 1; j++){
////				////		int candidate_start = all_edges[p].back().end;
////				////		int candidate_end = index[j];
////				////		float cur_widht = dist[j]; //从角点开始扩增的第一条边是得分函数值，之后是上一次粒子的权重*得分函数
////				////		bool isDuplicate = false;
////				////		for (const auto& edge : all_edges[p]){
////				////			if (edge.end == candidate_end || edge.start == candidate_end){
////				////				isDuplicate = true;
////				////				break;
////				////			}
////				////		}
////				////		if (cur_widht < min_widht && !isDuplicate){
////				////			min_widht = cur_widht;
////				////			new_edge.start = candidate_start;
////				////			new_edge.end = candidate_end;
////				////			new_edge.wight_value = cur_widht;
////				////		}
////				////	}
////				////}
////
////				all_edges[p].push_back(new_edge);
////			}
////		}
////	}

//////////void EMst::Score(pcl::KdTreeFLANN<pcl::PointXYZ> kdtree, int last, int candidate_start, pcl::PointXYZ cur_end, p_rawPoint P, p_rawPoint C, p_rawPoint B, int iter, int i, int k)
//////////{
//////////	//------------------------0------------------------输入原始数据、构建K近邻搜索---------------------------------
//////////	std::vector<int> index;
//////////	std::vector<float> dist;
//////////	kdtree.nearestKSearch(cur_end, k + 1, index, dist);
//////////	EdgeInfo new_edge;
//////////	new_edge.wight_value = 0;
//////////	float min_widht = std::numeric_limits<float>::max(); //初始化最小cur_widht为一个较大的值
//////////	for (int j = 1; j < k + 1; j++) {
//////////		int candidate_end = index[j];
//////////		float dis = dist[j];
//////////		std::array<float, 5> meanCVector = corner5mean(C);
//////////		std::array<float, 5> meanBVector = boundary5mean(B);
//////////		std::vector<float> scores_of_neighbours;//存储 20 个邻域的得分值
//////////		float cur_int = P->points[candidate_end].intensity;//candidate_end是调用了全部索引index[j]
//////////		float cur_gva = P->points[candidate_end].gvalue;
//////////		float cur_t0 = P->points[candidate_end].Tlambda0;
//////////		float cur_t1 = P->points[candidate_end].Tlambda1;
//////////		float cur_t2 = P->points[candidate_end].Tlambda2;
//////////		float cur_x = P->points[candidate_end].x;
//////////		float cur_y = P->points[candidate_end].y;
//////////		float cur_z = P->points[candidate_end].z;
//////////		std::array<float, 5> curVector = { cur_int, cur_gva, cur_t0, cur_t1, cur_t2 };
//////////		//------------------------1----------------------------计算Pzx-------------------------------------------------
//////////		float w1 = 0.7;
//////////		float a1 = exp(-1 * w1 * (iter + 0.01));
//////////		float a2 = 1 - a1;
//////////		//float C_density = GaussianDensity(curVector, meanCVector);
//////////		//float B_density = GaussianDensity(curVector, meanBVector);
//////////		//N1计算
//////////		Eigen::Map<const Eigen::Matrix<float, 5, 1>> CX(curVector.data());
//////////		Eigen::Map<const Eigen::Matrix<float, 5, 1>> CMean(meanCVector.data());
//////////		Eigen::Matrix<float, 5, 1> diff = CX - CMean;
//////////
//////////		Eigen::Matrix<float, 5, 5> covariance = (diff * diff.transpose()) / (curVector.size() - 1);
//////////		float determinant = covariance.determinant();
//////////		if (determinant <= 0) {
//////////			cout << "GaussianDensity  协方差矩阵可能不可逆或奇异，处理异常情况" << endl;
//////////		}
//////////		float exponent = -0.5f * diff.transpose() * covariance.inverse() * diff;
//////////		constexpr float pi = 3.14159265359f;
//////////		int dimensions = 5; // 维度数
//////////		float constant = 1.0f / (std::sqrt(std::pow(2 * pi, dimensions) * determinant));
//////////		float C_density = constant * std::exp(exponent);
//////////
//////////
//////////
//////////		//float Pzx = a1 * C_density + a2 * B_density;
//////////		//cout << "iter:   " << iter << endl;
//////////		//cout << "a1:   " << a1 << endl;
//////////		//cout << "C_density" << C_density << ":" << endl;//nan值
//////////		//cout << "B_density" << B_density << ":" << endl;//0
//////////		//cout << "Pzx" << Pzx << ":" << endl;
//////////		//-------------------------2---------------------------计算Pxx-------------------------------------------------
//////////		float start_x = P->points[candidate_start].x;
//////////		float start_y = P->points[candidate_start].y;
//////////		float start_z = P->points[candidate_start].z;
//////////		float last_x = P->points[last].x;
//////////		float last_y = P->points[last].y;
//////////		float last_z = P->points[last].z;
//////////		float A1 = last_x - start_x;
//////////		float A2 = last_y - start_y;
//////////		float A3 = last_z - start_z;
//////////		float D1 = start_x - cur_x;
//////////		float D2 = start_y - cur_y;
//////////		float D3 = start_z - cur_z;
//////////		//两个三维向量
//////////		std::array<float, 5> A = { A1, A2, A3 };
//////////		std::array<float, 5> D = { D1, D2, D3 };
//////////		double magnitude = std::sqrt(std::pow((A2 * D3 - A3 * D2), 2) + std::pow((A3 * D1 - A1 * D3), 2) + std::pow((A1 * D2 - A2 * D1), 2));
//////////		float w2 = 0.5;
//////////		float pxx = exp(-1 * w2 * magnitude);
//////////		float sum_pxx += pxx;
//////////		float Pxx;
//////////		if (i == k) {
//////////			Pxx = pxx / sum_pxx;
//////////			float Score = Pzx * Pxx;
//////////			//-------------------------3---------------------将得分存入数组或向量----------------------------------
//////////			float score = Score;
//////////			scores_of_neighbours.push_back(score);
//////////		}
//////////	}
//////////}

////////int EMst::Sampling(std::vector<float>& scores_of_neighbours) {
////////	// 计算累积概率
////////	std::vector<double> cumulativeProbabilities(scores_of_neighbours.size());
////////	double totalScore = 0.0;
////////	for (size_t i = 0; i < scores_of_neighbours.size(); ++i) {
////////		totalScore += scores_of_neighbours[i];
////////	}
////////	//cumulativeProbabilities[0] = scores_of_neighbours[0] / totalScore;
////////	for (size_t i = 1; i < scores_of_neighbours.size(); ++i) {
////////		cumulativeProbabilities[i] = cumulativeProbabilities[i - 1] + scores_of_neighbours[i] / totalScore;
////////	}
////////	// 生成随机数
////////	std::random_device rd;
////////	std::mt19937 gen(rd());
////////	std::uniform_real_distribution<> dis(0.0, 1.0);
////////	double randomValue = dis(gen); // 生成0到1之间的随机数
////////	int j = 0;
////////	for (size_t i = 0; i < cumulativeProbabilities.size(); ++i) {
////////		if (randomValue <= cumulativeProbabilities[i]) {
////////			j = i; // 找到对应的索引号
////////			break;
////////		}
////////	}
////////	std::cout << "随机生成的局部索引号为: " << j << std::endl;
////////	return j;
////////}

//std::vector<std::vector<float>> scores_array
//////////////int EMst::Sampling(std::vector<std::vector<float>>& scores_array){
//////////////	// 计算累积概率
//////////////	std::vector<double> cumulativeProbabilities(20);
//////////////	double totalScore = 0.0;
//////////////	for (size_t i = 0; i < 20; ++i) {
//////////////		totalScore += scores_array[i][1];
//////////////	}
//////////////	for (size_t i = 0; i < 20; ++i) {
//////////////		cumulativeProbabilities[i] = scores_array[i][1] / totalScore;
//////////////		if (i > 0) {
//////////////			cumulativeProbabilities[i] += cumulativeProbabilities[i - 1];
//////////////		}
//////////////	}
//////////////	// 生成随机数
//////////////	std::random_device rd;
//////////////	std::mt19937 gen(rd());
//////////////	std::uniform_real_distribution<> dis(0.0, 1.0);
//////////////	double randomValue = dis(gen); // 生成0到1之间的随机数
//////////////	int sampled_index = -1;
//////////////	for (size_t i = 0; i < 20; ++i) {
//////////////		if (randomValue <= cumulativeProbabilities[i]) {
//////////////			sampled_index = i;
//////////////			break;
//////////////		}
//////////////	}
//////////////	std::cout << "随机生成的局部索引号为: " << scores_array[sampled_index][0] << std::endl;
//////////////	return scores_array[sampled_index][0];
//////////////}

////从角点开始，判断是不是角点的函数
////void EMst::iscorner(p_rawPoint P, p_rawPoint C) {
////	int point_size = P->points.size();
////	int corner_size = C->points.size();
////	float threshold = 0.001; // 设置阈值
////	for (int i = 0; i < point_size; i++) {
////		pcl::PointXYZ cur_point;
////		cur_point.x = P->points[i].x;
////		cur_point.y = P->points[i].y;
////		cur_point.z = P->points[i].z;
////		for (int j = 0; j < corner_size; j++) {
////			pcl::PointXYZ corner_point;
////			corner_point.x = C->points[j].x;
////			corner_point.y = C->points[j].y;
////			corner_point.z = C->points[j].z;
////			// 检查坐标之间的差异是否小于阈值
////			float diff_x = std::abs(cur_point.x - corner_point.x);
////			float diff_y = std::abs(cur_point.y - corner_point.y);
////			float diff_z = std::abs(cur_point.z - corner_point.z);
////			if (diff_x < threshold && diff_y < threshold && diff_z < threshold) {
////				P->points[i].iscorner = 1;
////				break; // 已经找到匹配的角点，可以退出循环
////			}
////		}
////	}
////}

////std::vector<std::vector<long double>> EMst::Score(pcl::KdTreeFLANN<pcl::PointXYZ> kdtree, int last, int candidate_start, pcl::PointXYZ cur_end, p_rawPoint P, p_rawPoint C, p_rawPoint B, int iter, int k)
////{
////	//------------------------0------------------------输入原始数据、构建K近邻搜索---------------------------------
////	std::vector<int> index;
////	std::vector<float> dist;
////	kdtree.nearestKSearch(cur_end, k + 1, index, dist);
////	std::vector<std::vector<long double>> scores_array(20, std::vector<long double>(2)); // 20行2列的二维数组，存储索引号和对应的 score 值
////	for (int j = 1; j < k + 1; j++) {
////		int candidate_end = index[j];
////		long double dis = dist[j];
////		std::vector<long double> scores_of_neighbours;//存储 20 个邻域的得分值
////		std::array<long double, 5> meanCVector = corner5mean(C);
////		std::array<long double, 5> meanBVector = boundary5mean(B);
////		std::array<long double, 5> curVector = {
////			P->points[candidate_end].intensity,
////			P->points[candidate_end].gvalue,
////			P->points[candidate_end].Tlambda0,
////			P->points[candidate_end].Tlambda1,
////			P->points[candidate_end].Tlambda2
////		};
////		//------------------------1----------------------------计算Pzx-------------------------------------------------
////		long double w1 = 0.7;
////		long double a1 = exp(-1 * w1 * (iter + 0.01));//为了防止第一次
////		long double a2 = 1 - a1;
////		//-------------------------------N1计算+++++++++++++++++++++++++++++++++++++
////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> X(curVector.data());
////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> CMean(meanCVector.data());
////		Eigen::Matrix<long double, 5, 1> diff = X - CMean;
////		Eigen::Matrix<long double, 5, 5> total_covarianceC = Eigen::Matrix<long double, 5, 5>::Zero(); // 总协方差矩阵
////		long double total_determinant = 0.0; // 总协方差行列式之和
////		long double total_exponent = 0.0; // 总指数之和
////		for (int n = 1; n < k + 1; n++) {
////			int neighbour_index = index[n];
////			std::array<long double, 5> neighbourVector = {
////				P->points[neighbour_index].intensity,
////				P->points[neighbour_index].gvalue,
////				P->points[neighbour_index].Tlambda0,
////				P->points[neighbour_index].Tlambda1,
////				P->points[neighbour_index].Tlambda2
////			};
////
////			Eigen::Map<const Eigen::Matrix<long double, 5, 1>> neighbour(neighbourVector.data());
////			Eigen::Matrix<long double, 5, 1> neighbour_diff = neighbour - CMean;
////			total_covarianceC += neighbour_diff * neighbour_diff.transpose();
////		}
////		total_covarianceC /= k; // 取所有邻域协方差的均值
////		long double total_determinantC = total_covarianceC.determinant(); // 计算总协方差行列式
////		if (total_determinantC <= 0) {
////			cout << "GaussianDensityC 协方差矩阵可能不可逆或奇异，处理异常情况" << endl;
////			total_covarianceC += Eigen::Matrix<long double, 5, 5>::Identity() * 0.001; // 正则化协方差矩阵
////			total_determinantC = total_covarianceC.determinant(); // 重新计算行列式
////		}
////		constexpr long double pi = 3.14159265359f;
////		long double exponent = -0.5f * diff.transpose() * total_covarianceC.inverse() * diff;
////		int dimensions = 5; // 维度数
////		long double constant = 1.0f / (std::sqrt(std::pow(2 * pi, dimensions) * total_determinantC));
////		long double C_density = constant * std::exp(exponent);
////		//////////////////////////////////////cout << "C_density:  " << C_density << endl;
////		//------------------------N2计算---------------------
////		std::array<long double, 5> curVectorB = {
////			P->points[candidate_end].intensity,
////			P->points[candidate_end].gvalue,
////			P->points[candidate_end].Tlambda0,
////			P->points[candidate_end].Tlambda1,
////			P->points[candidate_end].Tlambda2
////		};
////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> XB(curVectorB.data());
////		Eigen::Map<const Eigen::Matrix<long double, 5, 1>> BMean(meanBVector.data());
////		Eigen::Matrix<long double, 5, 1> diffB = XB - BMean;
////		Eigen::Matrix<long double, 5, 5> total_covarianceB = Eigen::Matrix<long double, 5, 5>::Zero(); // 总协方差矩阵
////		for (int n = 1; n < k + 1; n++) {
////			int neighbour_index = index[n];
////			std::array<long double, 5> neighbourVector = {
////				P->points[neighbour_index].intensity,
////				P->points[neighbour_index].gvalue,
////				P->points[neighbour_index].Tlambda0,
////				P->points[neighbour_index].Tlambda1,
////				P->points[neighbour_index].Tlambda2
////			};
////			Eigen::Map<const Eigen::Matrix<long double, 5, 1>> neighbourB(neighbourVector.data());
////			Eigen::Matrix<long double, 5, 1> neighbour_diffB = neighbourB - BMean;
////			total_covarianceB += neighbour_diffB * neighbour_diffB.transpose();
////		}
////		total_covarianceB /= k; // 取所有邻域协方差的均值
////		long double determinantB = total_covarianceB.determinant();
////		////////cout << "Total Covariance Matrix_B: " << endl;
////		////////for (int i = 0; i < total_covarianceB.rows(); ++i) {
////		////////	for (int j = 0; j < total_covarianceB.cols(); ++j) {
////		////////		cout << total_covarianceB(i, j) << " ";
////		////////	}
////		////////	cout << endl;
////		////////}
////		Eigen::JacobiSVD<Eigen::Matrix<long double, 5, 5>> svdB(total_covarianceB, Eigen::ComputeThinU | Eigen::ComputeThinV);
////		long double thresholdC = 1e-6; // 设置一个阈值，判断奇异性
////		if (svdB.singularValues()(4) < thresholdC) {
////			cout << "GaussianDensityBB 协方差矩阵可能奇异，处理异常情况" << endl;
////			// 此时可能需要对矩阵进行调整或使用其他方法处理奇异性
////		}
////		determinantB = svdB.singularValues().prod(); // 重新计算行列式
////		//////////////////////////////////////cout << "DeterminantBBBB of the covariance matrix: " << determinantB << endl;
////		long double exponentB = -0.5f * diffB.transpose() * total_covarianceB.inverse() * diffB;
////		long double constantB = 1.0f / (std::sqrt(std::pow(2 * pi, dimensions) * determinantB));
////		long double B_density = constantB * std::exp(exponentB);
////		//////////////////////////////////////cout << "B_density:   " << B_density << endl;
////		// 计算总的 Pzx
////		long double Pzx = a1 * C_density + a2 * B_density;
////		//////////////////////////////////////cout << "Total Pzx:   " << Pzx << endl;
////		//-------------------------2---------------------------计算Pxx-------------------------------------------------
////		long double start_x = P->points[candidate_start].x;
////		long double start_y = P->points[candidate_start].y;
////		long double start_z = P->points[candidate_start].z;
////		long double last_x = P->points[last].x;
////		long double last_y = P->points[last].y;
////		long double last_z = P->points[last].z;
////		long double cur_x = P->points[index[j]].x;
////		long double cur_y = P->points[index[j]].y;
////		long double cur_z = P->points[index[j]].z;
////		long double A1 = last_x - start_x;
////		long double A2 = last_y - start_y;
////		long double A3 = last_z - start_z;
////		long double D1 = start_x - cur_x;
////		long double D2 = start_y - cur_y;
////		long double D3 = start_z - cur_z;
////		//两个三维向量
////		std::array<long double, 5> A = { A1, A2, A3 };
////		std::array<long double, 5> D = { D1, D2, D3 };
////		double mo = std::sqrt(std::pow((A2 * D3 - A3 * D2), 2) + std::pow((A3 * D1 - A1 * D3), 2) + std::pow((A1 * D2 - A2 * D1), 2));
////		long double w2 = 0.5;
////		long double sum_pxx = 0.0; // 初始化 px 总和为 0
////		// 计算每个邻域内的 exp(w * mo) 并求和
////		for (int n = 1; n < k + 1; n++) {
////			double mo = std::sqrt(std::pow((A2 * D3 - A3 * D2), 2) + std::pow((A3 * D1 - A1 * D3), 2) + std::pow((A1 * D2 - A2 * D1), 2));
////			long double pxx = exp(-w2 * mo);
////			sum_pxx += pxx; // 将 pxx 加到总和中
////		}
////		long double Pxx = exp(w2 * mo) / sum_pxx; // 计算 Pxx
////		//////////////////////////////////////cout << "Pxx:   " << Pxx << endl;
////		long double Score = Pzx * Pxx;
////		//////////////////////////////////////cout << "Score:   " << Score << endl;
////		//-------------------------3---------------------将得分存入数组或向量----------------------------------
////		scores_array[j - 1][0] = index[j]; // 存储索引号
////		scores_array[j - 1][1] = Score;   // 存储对应的 score 值
////	}
////	return scores_array; // 返回二维数组
////}
//////根据离散分布律的随机采样函数
////int EMst::Sampling(std::vector<std::vector<long double>>& scores_array) {
////	// 计算累积概率
////	std::vector<long double> cumulativeProbabilities(20);
////	long double totalScore = 0.0;
////	int emptyCount = 0; // 记录空值个数
////	for (size_t i = 0; i < 20; ++i) {
////		if (scores_array[i][1] == 0.0) { // 假设空值为0
////			emptyCount++;
////		}
////		else {
////			totalScore += scores_array[i][1];
////		}
////	}
////	if (emptyCount == 20) { // 如果全部为空值，生成一个随机索引
////		std::random_device rd;
////		std::mt19937 gen(rd());
////		std::uniform_int_distribution<> dis(0, 19);
////		int randomIndex = dis(gen);
////		////////////////////////std::cout << "随机生成的局部索引号为: " << scores_array[randomIndex][0] << std::endl;
////		return scores_array[randomIndex][0];
////	}
////	for (size_t i = 0; i < 20; ++i) {
////		if (scores_array[i][1] == 0.0) {
////			cumulativeProbabilities[i] = 1.0 / 20; // 平均分布
////		}
////		else {
////			cumulativeProbabilities[i] = scores_array[i][1] / totalScore;
////			if (i > 0) {
////				cumulativeProbabilities[i] += cumulativeProbabilities[i - 1];
////			}
////		}
////	}
////	// 生成随机数
////	std::random_device rd;
////	std::mt19937 gen(rd());
////	std::uniform_real_distribution<> dis(0.0, 1.0);
////	long double randomValue = dis(gen); // 生成0到1之间的随机数
////	int sampled_index = -1;
////	for (size_t i = 0; i < 20; ++i) {
////		if (randomValue <= cumulativeProbabilities[i]) {
////			sampled_index = i;
////			break;
////		}
////	}
////	if (sampled_index != -1) {
////		std::cout << "随机生成的局部索引号为: " << scores_array[sampled_index][0] << std::endl;
////		return scores_array[sampled_index][0];
////	}
////	else {
////		std::cerr << "Error: Unable to determine sampled index." << std::endl;
////		return -1; // 或者其他适合的错误代码
////	}
////}