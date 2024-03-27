#pragma once

#include "helper.h"
#include "EMst.h"

using namespace std;
//����û��ʹ��ָ����������������Ϊ���ǵ���������ַ�����������к�ᶪʧ�������Ⱥ�ӳ�����Ϣ,����ʹ�õݹ鷽�����Ա�����Щ����
//�����ر���Ҫע�⣬������push_back�൱���������һ��������ȥ
typedef struct PathTree {
	vector<int> self;
	vector<PathTree> childChains;
}PathTree;

vector<int> deprecatedList;
//������ֵ
int PathMinLimit = 30;
//����С��ֵ
int TreeMinLimit = 30;
//����ֵ
float HeightLimit = 0.005f;

//vector<int> DPSimplify(p_Mypoint &P, vector<int> &nodeNumList, vector<int> singlePath, float &tolerance);
//void outputSimplePolylineFile(p_Mypoint &P, PathTree &pt, float tolerance);
//
//template <typename PointT>
//float getHeight(PointT ps, PointT pe, PointT p0);
//
//ostream & operator<<(ostream &cout, PathTree &p);
//int searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos);
//PathTree searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos, int flag, p_Mypoint &P);
//p0��ʾ��ߵĵ�

template <typename PointT>
float getHeight(PointT ps, PointT pe, PointT p0)
{
	float line_S0 = sqrt(pow(ps.x - p0.x, 2) + pow(ps.y - p0.y, 2) + pow(ps.z - p0.z, 2));
	float line_E0 = sqrt(pow(pe.x - p0.x, 2) + pow(pe.y - p0.y, 2) + pow(pe.z - p0.z, 2));
	float line_SE = sqrt(pow(ps.x - pe.x, 2) + pow(ps.y - pe.y, 2) + pow(ps.z - pe.z, 2));
	if (line_S0 + line_E0 > line_SE && line_S0 + line_SE > line_E0 && line_SE + line_E0 > line_S0) {
		{
			//�����ܳ�
			float s = (line_S0 + line_E0 + line_SE) / 2;
			//���׹�ʽ
			float area = sqrt(s*(s - line_S0)*(s - line_E0)*(s - line_SE));
			//����p0��Ӧ��
			return area / line_SE * 2;
		}
	}
	else return min(line_S0, line_E0);
}

//���룺ȫ���б������б�������ֵ
//�������룺�����������꣩������ű�ӳ��ʹ�ã�
//���أ���·��
vector<int> DPSimplify(p_Mypoint &P, vector<int> &nodeNumList, vector<int> singlePath, float &tolerance)
{
	if (singlePath.size() > 2)
	{
		int ps = nodeNumList[singlePath[0]];
		int pe = nodeNumList[singlePath[singlePath.size() - 1]];

		int indexMax = 0;
		float hMax = 0;
		for (int i = 1; i < singlePath.size() - 1; i++)
		{
			int pi = nodeNumList[singlePath[i]];
			float dH = getHeight(P->points[ps], P->points[pe], P->points[pi]);
			if (dH > tolerance && dH > hMax)
			{
				hMax = dH;
				indexMax = i;
			}
		}
		//���Ϊ0����;��ȫ����ȥ
		if (hMax == 0)
		{
			singlePath.erase(singlePath.begin() + 1, singlePath.end() - 1);
			return singlePath;
		}
		else {
			vector<int> left;

			for (int i = 0; i < indexMax + 1; i++) {
				left.push_back(singlePath[i]);
			}

			vector<int> right;
			for (int i = indexMax; i < singlePath.size(); i++) {
				right.push_back(singlePath[i]);
			}
			vector<int> all;
			connect(all, DPSimplify(P, nodeNumList, left, tolerance));
			connect(all, DPSimplify(P, nodeNumList, right, tolerance));

			return all;
		}
	}
	else return singlePath;

}

template <typename PointT>
float getArea(PointT p1, PointT p2, PointT p3)
{
	float line_13 = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2) + pow(p1.z - p3.z, 2));
	float line_23 = sqrt(pow(p2.x - p3.x, 2) + pow(p2.y - p3.y, 2) + pow(p2.z - p3.z, 2));
	float line_12 = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));

	//�����ܳ�
	float s = (line_13 + line_23 + line_12) / 2;
	//���׹�ʽ
	float area = sqrt(s*(s - line_13)*(s - line_23)*(s - line_12));
	//����p3��Ӧ��
	return area;

}

//���룺ȫ���б������б������ֵ
//�������룺�����������꣩������ű�ӳ��ʹ�ã�
//���أ���·��
//���룺ȫ���б������б������ֵ
//�������룺�����������꣩������ű�ӳ��ʹ�ã�
//���أ���·��
vector<int> VMSimplify(p_Mypoint &P, vector<int> &nodeNumList, vector<int> singlePath, float &tolerance)
{
	//���ﲻͬ������Ϊʹ���������Σ����뱣֤·���ڵ�������3
	if (singlePath.size() > 3)
	{
		//singlePath��ԭʼ���ӳ����洢����·���ϵĵ��ڽ���ű��е�λ�ã����ڷ����ô��Σ����Կ��������޸�����
		//���������Σ���singlePath��˳��洢����

		//��ʼ�����������˳�����Ҫȷ������Ϊ���ӳ���ĳ���-2
		vector<float> areaList;
		for (int i = 1; i < singlePath.size() - 1; i++) {
			int ps = nodeNumList[singlePath[i - 1]];
			int pm = nodeNumList[singlePath[i]];
			int pe = nodeNumList[singlePath[i + 1]];
			areaList.push_back(getArea(P->points[ps], P->points[pm], P->points[pe]));
		}

		//��˳�����������С����������λ�úͶ�Ӧ����
		//��������˵��������ֻ�оֲ������˱仯Ӧ������̬���£�����ʱ����������ֱ��ʹ����������Сֵ�ĺ���.
		int count = 0;
		int baseSize = singlePath.size();
		while (true)
		{
			//�����С�����ε��±�
			int min_data = min_element(areaList.begin(), areaList.end()) - areaList.begin();
			//�������Ҫ����ֹѭ�������߳����쳣ѭ��ʱҲ��Ҫ��ֹ
			if (areaList[min_data] > tolerance || singlePath.size() < 4)break;
			if (count > baseSize)break;
			//�����С�����ζ�Ӧ������ӳ����е��±�
			int pmin = min_data + 1;
			logger(singlePath[singlePath.size() - 1]);
			//logger("pos:",min_data,",size:",areaList.size(),",size2:", singlePath.size());
			//������Ҫ�ų���λ�ڵڶ����͵����ڶ������������Σ������λ���м������
			if ((pmin - 1) == 0)
			{
				int p1 = nodeNumList[singlePath[0]];
				//int pm = nodeNumList[singlePath[pmin]];
				int p3 = nodeNumList[singlePath[2]];
				int p4 = nodeNumList[singlePath[3]];
				float a1 = getArea(P->points[p1], P->points[p3], P->points[p4]);
				//��������б�
				areaList.erase(areaList.begin(), areaList.begin() + 2);
				areaList.insert(areaList.begin(), a1);
				//����·��ӳ���
				singlePath.erase(std::begin(singlePath) + 1);
			}
			else if ((pmin + 1) == (singlePath.size() - 1))
			{
				int p1 = nodeNumList[singlePath[pmin - 2]];
				int p2 = nodeNumList[singlePath[pmin - 1]];
				//int pm = nodeNumList[singlePath[pmin]];
				int p4 = nodeNumList[singlePath[pmin + 1]];
				float a1 = getArea(P->points[p1], P->points[p2], P->points[p4]);
				//��������б�
				areaList.erase(areaList.end() - 2, areaList.end());
				areaList.push_back(a1);
				//����·��ӳ���
				singlePath.erase(std::end(singlePath) - 2);
			}
			else {
				int p1 = nodeNumList[singlePath[pmin - 2]];
				int p2 = nodeNumList[singlePath[pmin - 1]];
				//int pm = nodeNumList[singlePath[pmin]];
				int p4 = nodeNumList[singlePath[pmin + 1]];
				int p5 = nodeNumList[singlePath[pmin + 2]];
				float a1 = getArea(P->points[p1], P->points[p2], P->points[p4]);
				float a2 = getArea(P->points[p2], P->points[p4], P->points[p5]);
				//��������б�
				areaList.erase(areaList.begin() + (min_data - 1), areaList.begin() + (min_data + 2));
				areaList.insert(areaList.begin() + (min_data - 1), a1);
				areaList.insert(areaList.begin() + (min_data), a2);
				//����·��ӳ���
				singlePath.erase(singlePath.begin() + pmin);
			}
			if (singlePath.size() != areaList.size() + 2)logger(singlePath.size(), "//", areaList.size());
			count++;
		}


	}
	return singlePath;
}





template <typename PointT>
float getDistance(PointT p1, PointT p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

//��������ֵ
float PointInterval = 0.01000f;
vector<MyPointType> createConsecutivePoints(const Mypoint &P, vector<int> &nodeNumList, vector<int> &pathList) {
	vector<MyPointType> pSet;
	int size = pathList.size();
	//ѹ����ʼ��
	pSet.push_back(P[nodeNumList[pathList[0]]]);
	for (int i = 0; i < size-1; i++)
	{
		int pos1 = nodeNumList[pathList[i]];
		int pos2 = nodeNumList[pathList[i+1]];

		//����һ��������
		float d = getDistance(P[pos1], P[pos2]);
		//logger(d);
		//�������м�������
		int Num = round((d / PointInterval) - 1);
		float x0 = P[pos1].x >P[pos2].x ? P[pos1].x : P[pos2].x;
		float y0 = P[pos1].y > P[pos2].y ? P[pos1].y : P[pos2].y;
		float z0 = P[pos1].z > P[pos2].z ? P[pos1].z : P[pos2].z;
		if (Num > 0)
		{
			//���·�����
			float dx = fabs((P[pos1].x - P[pos2].x) / ((float)Num));
			float dy = fabs((P[pos1].y - P[pos2].y) / ((float)Num));
			float dz = fabs((P[pos1].z - P[pos2].z) / ((float)Num));
			//logger(dx);
			for (int d=1; d <= Num; d++)
			{
				float x = x0 + d*dx;
				float y = y0 + d*dy;
				float z = z0 + d*dz;
				MyPointType dp;
				dp.x = x;
				dp.y = y;
				dp.z = z;
				pSet.push_back(dp);
			}
		}

		//���ڶ����յ�ѹ��
		pSet.push_back(P[pos2]);
	}

	return pSet;

}
void adjustProbability(EMst &g, vector<int> &nodeNumList, vector<int> &pathList)
{
	
	//logger(nodeNumList);
	//������С����
	int size = nodeNumList.size();

	vector<MyPointType> pSet = createConsecutivePoints(g.knnPoint, nodeNumList, pathList);
	int childSize = pSet.size();
	for (int i = 0; i < size; i++)
	{
		int pos = nodeNumList[i];
		float baseD = Inf;
		
		for (int j = 0; j < childSize; j++)
		{
			float d = getDistance(g.knnPoint[pos], pSet[j]);
			if (d < baseD)baseD = d;
		}
		float p0 = 0.5 * (g.knnPoint[pos].p0 + exp(-(pow(baseD, 2))));
		//cout << "����p0��" << p0 << endl;
		//if (p0 > threshold_p) 
		//{
		//	nodeNumList.front()
		//	g.knnPoint.erase(i); i--; }
		g.knnPoint[pos].p0 = p0;
		g.knnPoint[pos].p1 = 1 - p0;
		g.knnPoint[pos].ptype = p0 > 0.5 ? 0 : 1;
	}
}

//���������Ϊ�˲�ѯԪ��������������Ĵ���
int order(vector<int> v, int element) {
	for (int order = 0; order < v.size(); order++) {
		if (element == v[order])
			return order;
	}
	//Ϊ���������������
	return -1;
}
//���ӻ�1029
void outputSimplePolylineFile(EMst &g, p_Mypoint &P, PathTree &pt, float tolerance,int Orderitreator)
{
	vector<int> all = pt.self;
	if (pt.childChains.size() > 0)
	{
		//��Ҳ��һ��ӳ�䣬�������path����ĵ���ʵ�������б��
		vector<int> rootNodeList;

		for (int i = 0; i < pt.childChains.size(); i++)
		{

			vector<int> c = pt.childChains[i].self;
			int pos = order(all, c[c.size() - 1]);
			if (pos != -1) {
				rootNodeList.push_back(pos);
			}
			else {
				logger("�쳣" + s(all[0]));
			}
			outputSimplePolylineFile(g, P, pt.childChains[i], tolerance, Orderitreator);
		}
		//����˳�򣬷�������
		quickSort(rootNodeList, 0, rootNodeList.size() - 1);
		vector<int> singlePath;
		//outputPolylineFile(P, pt.self, createSequentialList(0, pt.self.size()), s(pt.self[0]));
		for (int i = 0; i < all.size(); i++)
		{
			singlePath.push_back(all[i]);
			if ((rootNodeList.size() > 0 && i == rootNodeList[0]) ||
				i == all.size() - 1)
			{
				PathTree pt0 = *new PathTree();
				pt0.self = singlePath;
				//pt0.root = vectorToList(singlePath);
				outputSimplePolylineFile(g, P, pt0, tolerance, Orderitreator);

				if (i < all.size() - 1) {
					rootNodeList.erase(rootNodeList.begin());
					singlePath.clear();
					singlePath.push_back(all[i]);
				}

			}
		}
	}
	else
	{
		vector<int> pathMap = createSequentialList(0, all.size());

		//debugList(pt.self);
		if (all.size() > 0)
		{
			//�����ٴ�һ��δ�򻯰汾��ȥ�򻯼���
			//outputPolylineFile(P, pt.self, pathMap, s(pt.self[0]));
			//vector<int> resultPath = DPSimplify(P, all, pathMap, tolerance);
			vector<int> resultPath = VMSimplify(P, all, pathMap, tolerance);
			outputPolylineFile(P, all, resultPath, s(Orderitreator));
			adjustProbability(g, all, resultPath);
		}
	}
}


ostream & operator<<(ostream &cout, PathTree &p)
{
	cout << "[";
	int i = 0;
	for (; i < p.self.size() - 1; i++)cout << p.self[i] << ",";
	if (i == p.self.size() - 1)cout << p.self[i] << "]";
	if (p.childChains.size() > 0)
	{
		cout << "[";
		for (int j = 0; j < p.childChains.size(); j++) {
			PathTree pt = p.childChains[j];
			cout << pt << ",";
		}
		cout << "]";
	}
	return cout;
}

/*
����ѡ�������ȵ�ԭ������Ϊ������Σ����Ƕ�Ҫ�������еķ�֧�Ա���©
��������㷨��

*/
//startPos����Ӧ��ָ������node�б������λ��
//nodeNumList�洢��Ӧλ�õĵ���ȫͼ�е��ڽӱ�λ��
//����ֵΪҲ����nodeNumList��λ��
void dfs(GraphAdjList* tree, vector<int> *nodeNumList, vector<int> *visitedList, int root, int dep, vector<int> *depthList) {
	int childDep = dep + 1;
	//cout << root << endl;
	EdgeNode *first = tree->adjList[nodeNumList->at(root)].firstedge;
	while (first)
	{
		if (order_element_in_vector(*(vector<int>*) visitedList, first->adjvex) == -1)
		{
			int childPositon = order_element_in_vector(*(vector<int>*) nodeNumList, first->adjvex);
			if (childPositon != -1)
			{
				visitedList->push_back(first->adjvex);
				dfs(tree, nodeNumList, visitedList, childPositon, childDep, depthList);
				//Ϊʲô�������д���±�0����depthList[0][childPositon]
				(*depthList)[childPositon] = childDep;
				//cout << first->adjvex << "�ӽڵ�" <<","<< childDep << depthList->at(childPositon) << endl;
			}
		}
		first = first->next;
	}
}
//������·����ʵ�ȼ�¼��ȣ�Ȼ������Ѱ����ȵ͵ĵط����ɣ���Ⱥ�������ȶ��ܱ�ʾ
//������Ȼ��ݲ�ȷ��Ҫ��Ҫ�Ż���ʼ��ѡ���������Ϊ���ھӵ���ʼ�㣩
//�������Ϊ����������
vector<int> backtracking(GraphAdjList* tree, vector<int> &nodeNumList, vector<int> &depthList, int maxDepthPos, int root) {

	//cout << nodeNumList[maxDepthPos] << endl;
	//���
	int dep = depthList[maxDepthPos];
	//�����ʾ����Ľ��ĸ����
	//�����и�������û�а�δ�������������ڽ�
	EdgeNode *endFather = tree->adjList[nodeNumList[maxDepthPos]].firstedge;
	while (endFather) {
		int nPos = order_element_in_vector(nodeNumList, endFather->adjvex);
		if (depthList[nPos] < depthList[maxDepthPos])
		{
			dep = depthList[nPos];

			break;
		}
		endFather = endFather->next;
	}

	//·���б��ɺ���ǰ
	vector<int> pathList;
	//����һ����ѹ��·���б�
	pathList.push_back(maxDepthPos);
	//�������������Ϊ����0
	while (dep > 0)
	{
		//���Ƚ��ø��ڵ�ѹ��·���б�
		pathList.push_back(order_element_in_vector(nodeNumList, endFather->adjvex));
		//�����Ǳ������ڵ���ھӣ����縸�ڵ���һ����ȱȸ��ڵ�С�ĵ㣬����Ϊ���ڵ�ĸ��ڵ�
		vector<EdgeNode *> neighborList;
		//ѡ�񸸽ڵ���ھ�
		EdgeNode *friendNode = tree->adjList[endFather->adjvex].firstedge;
		while (friendNode)
		{
			//��¼�����ھ�
			neighborList.push_back(friendNode);
			//ǰ����һ���ھ�
			friendNode = friendNode->next;
		}
		//Ѱ�����С���ھ�
		for (int i = 0; i < neighborList.size(); i++)
		{
			int neighborPos = order_element_in_vector(nodeNumList, neighborList[i]->adjvex);
			if (depthList[neighborPos] < dep)
			{
				//order_element_in_vector(nodeNumList, friendNode->adjvex)
				dep = depthList[neighborPos];
				endFather = neighborList[i];
				break;
			}
		}
		//�����ھӱ�ռ�
		neighborList.clear();
	}
	//��������������⣬����������Ҫ�������
	pathList.push_back(root);
	return pathList;
}




//����flag��ʾ�Ƿ����·��
PathTree searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos, int flag) {
	logger(s(flag));
	int totalAmount = nodeNumList.size();

	vector<int> visitedList;

	//depthList��nodeNumList��һ��ӳ���ϵ��ͬλԪ���Ƕ�ĳ����Ĳ�ͬ���Եı��������Ҳ���Կ����ýṹ�崦��
	//��˽�ֹ�������е���������м򵥲��������Ҫ��������Ҫ������������
	vector<int> depthList;
	for (int vec = 0; vec < nodeNumList.size(); vec++)depthList.push_back(-1);

	visitedList.push_back(nodeNumList[startPos]);
	//tagList[nodeNumList[startPos]] = 1;
	depthList[startPos] = 0;
	dfs(tree, &nodeNumList, &visitedList, startPos, 0, &depthList);


	//����ֻ��Ҫ��ȡ��Զ�㣬��ô������ȷ��
	int maxDepthPos = 0;
	for (int i = 0; i < depthList.size(); i++)
	{
		if (depthList[i] > depthList[maxDepthPos])maxDepthPos = i;
	}
	visitedList.swap(vector<int>());

	//if (flag == 0)return depthList.size() > 0 ? maxDepthPos : -1;
	if (flag == 0) {
		PathTree tmp0;
		vector<int> tmp; tmp.push_back(maxDepthPos);
		tmp0.self = tmp;
		return tmp0;
	}

	//ע��·������Ҫ��Ͻ��ӳ���ʹ��
	vector<int> pathList = backtracking(tree, nodeNumList, depthList, maxDepthPos, startPos);
	depthList.swap(vector<int>());
	//��������Ļ������Ľ�������flag�������ص�����
	//if (flag == 0)return *new vector<int>->push_back(maxDepthPos);

	//������ô��������Ҫ����·������������������·�����������ȫ����ע��·������Ҫ��Ͻ��ӳ���ʹ��
	PathTree childPath;
	vector<int> realPath = pathList;
	if (pathList.size() < PathMinLimit)
	{
		connect(deprecatedList, nodeNumList);
		//connectByMapping(deprecatedList, pathList, nodeNumList);

		return *new PathTree();
	}
	else
	{
		inverseMapping(realPath, nodeNumList);

		childPath.self = realPath;

		//debugList(realPath);
		connectByMapping(deprecatedList, pathList, nodeNumList);
	}


	//���·�����ļ���
	string filename = to_string(nodeNumList[startPos]);
	//outputPolylineFile(P, nodeNumList, pathList, filename);
	//logger(filename);

	/*����pathList�ĵ��ʾ�Ѿ������Ϊ���ɵĵ㣬��������б�����ֱ�ӱ�����ʹ�ã���Ϊ����Ŀ������Ϊǰ������������
	������Ҫ��nodeNumListʣ�µĵ��н����Ų飬��Щ�㶼������ȣ������������ϢϢ���
	��Щ������µĽ���ű�
	�����ϵĽ����Ƿ�֧����㣬Ҳ�Ƿ�֧�����͵ĵط�
	*/

	//ֻ�е�flagΪ2ʱ��ִ�з�֧����
	/*
	���Ҫ�����֧�ķ�֧�����������ô�������ɸѡ��������ڿ�ʼ���ۣ�����������������޶�Լ��
	��һ���������ĸ��ڵ��������������ɽ�㣬���������ɵ���㣻
	�ڶ���������������Ľڵ㲻��Ҫ�����ǣ�
		��һ�����Ѿ�����������������̫�٣���ô�����������е�����������ϣ������������������ͬ��ȫ������ü��ϣ����������Ͽ��Կ��Ǵ�ȫ���ӳ���ű��Ͻ����Ƴ������ǻ���Ҫͬʱ������ӳ����в�ѯ��λ��ʱ����Ϊ-1�����
		�ڶ������Ѿ���֤�������ɵĵ㣨ע�����������Ҫ�ر���������

	������ѯ�Ĺ�����һ���ݹ���̣���������
	1�����ɲ�ѯ��㣬����������ӳ���ű�����̫С����ɣ�����̽����������ڵ㣬���ݲ�ѯ����Ƿ��꣬��������·��������ȫ�����ɣ�
		n�����������е�ǰ���ɵĽڵ�˳���ͬ��һ�����������㹻��������������n�����ɵĲ�ѯ
	ǿ����ֹ������û�е���δ���༯��


	*/
	if (flag > 0)
		for (int i = 0; i < pathList.size(); i++)
		{
			int adjvexMainBranch = nodeNumList[pathList[i]];
			//�ȴ���2�ſ�����ǿ׳��֧�����
			if (tree->adjList[adjvexMainBranch].degree > 2)
			{

				//��ѯ����ھ�
				EdgeNode *friendNode = tree->adjList[adjvexMainBranch].firstedge;
				while (friendNode)
				{

					//�����ӽڵ��ھӶ��ԣ�����ǰ����ȿ��ܲ��䣨�����ͬһ�����ӽ����ڶ������������Ҫ������һ��������ӽڵ�
					//һ����ʱ�ԵĽ��������������Ҫ�ṩ�����Ľ���ű������Ƴ������ɣ���Ҫ���¹������������ֱ���ù�����ȴ���
					int friendAddressGlobal = friendNode->adjvex;
					//ֻ�е��ھӲ�������ʱ������Ϊ�ӽڵ�
					if (order_element_in_vector(pathList, order_element_in_vector(nodeNumList, friendAddressGlobal)) == -1
						&& order_element_in_vector(deprecatedList, friendAddressGlobal) == -1)
					{

						queue<int> nodeBranchQueque;
						//�����ڼ������������Ľ���ű����������Ҫʹ��ȫ���±꣬�����ǽ���ű��ӳ���±�
						nodeBranchQueque.push(friendAddressGlobal);
						vector<int> childNodeNumList;
						childNodeNumList.push_back(friendAddressGlobal);
						while (nodeBranchQueque.size()) {

							//̽���ӽڵ�
							queue<int> childLayerQueque;
							while (nodeBranchQueque.size())
							{
								EdgeNode *first = tree->adjList[nodeBranchQueque.front()].firstedge;
								while (first)
								{
									if (first->adjvex != adjvexMainBranch
										&& order_element_in_vector(childNodeNumList, first->adjvex) == -1
										)
									{
										childLayerQueque.push(first->adjvex);
										childNodeNumList.push_back(first->adjvex);
									}
									first = first->next;
								}
								nodeBranchQueque.pop();
							}
							nodeBranchQueque = childLayerQueque;
							childLayerQueque.empty();
						}
						//�����������Ϊȫ����ֵ
						if (childNodeNumList.size() > TreeMinLimit)
						{	//������Ӧ����0����Ϊ�˱��������һ��λ�ã�������˳���ң���ܿ죩
							childNodeNumList.insert(childNodeNumList.begin(), adjvexMainBranch);

							PathTree reslut = searchByDistance(tree, childNodeNumList, order_element_in_vector(childNodeNumList, adjvexMainBranch), 2);
							if (reslut.self.size() > 0)
								childPath.childChains.push_back(reslut);

							childNodeNumList.swap(vector<int>());
							//cout << childNodeNumList[B] << endl;
						}
						else
						{
							connect(deprecatedList, childNodeNumList);
						}
					}
					//ǰ����һ���ھ�
					friendNode = friendNode->next;
				}
			}
		}
	if (flag > 0)
	{

		return childPath;
	}
	return  *new PathTree();

}

int searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos) {
	vector<int> result = searchByDistance(tree, nodeNumList, startPos, 0).self;
	return result.size() > 0 ? result[0] : -1;
}