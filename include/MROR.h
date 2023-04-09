#pragma once
#include <iostream>
#include <cmath>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/pcd_io.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <string.h>
using namespace std;

class MROR{
public:
	MROR() = default;
	~MROR() = default;
	
	void SetRadius(double radius_multiplier);
    double GetRadius();
    void SetAngle(double azimuth_angle);
    double GetAngle();
    void SetMinNeighbors(double min_neighbors);
    double GetMinNeighbors();
    void SetMinRadius(double min_search_radius);
    double GetMinSearchRadius();
	void SetKSearch(double KSearch);
	double GetKSearch();
	void SetMaDistance(double min_madistance);

	template<typename T>
	void Filter(typename pcl::PointCloud<T>::Ptr& input_cloud,
		    	typename pcl::PointCloud<T>& filtered_cloud);
	template<typename T>
	float getMahalanobisDistance(const T& centerpoint,
		                         vector<int>& pointIndex,
								 typename pcl::PointCloud<T>& input_cloud);

private:
    	double radius_multiplier_{3};
    	double azimuth_angle_{0.04};
    	double min_neighbors_{3};
    	double min_search_radius_{0.04};
    	double KSearch_{4};
		double min_madistance_{0.4};
};

template<typename T>
void MROR::Filter(typename pcl::PointCloud<T>::Ptr& input_cloud, typename pcl::PointCloud<T>& filtered_cloud){
   	//using 实现类型的重新定义
    using KdTreePtr = typename pcl::KdTreeFLANN<T>::Ptr;
    filtered_cloud.clear();

    // kd树初始化
    KdTreePtr kd_tree_(new pcl::KdTreeFLANN<T>());
    kd_tree_->setInputCloud(input_cloud);
    for(typename pcl::PointCloud<T>::iterator it = input_cloud->begin(); it != input_cloud->end(); ++it){
            float x_i = it->x;
            float y_i = it->y;
            float range_i = sqrt(pow(x_i, 2) + pow(y_i, 2));
            float search_radius_dynamic = radius_multiplier_ * azimuth_angle_ * M_PI /180 * pow(range_i / 5, 3);
	    	vector<int> pointIndxKNNSearch;			//nearstSearch的下标集合
	    	vector<float> pointSquareDistance;		//nearstSearch的距离集合
	    	vector<int> pointIndxRadiusSearch;		//RadiusNNSearch的下标集合
	    	vector<float> pointRadiusDistance;		//RadiusNNSearch的距离集合

	    	//基于马氏距离的判定
        	if(search_radius_dynamic < min_search_radius_){  
				kd_tree_->nearestKSearch(*it, KSearch_, pointIndxKNNSearch, pointSquareDistance);
				float m_distance = getMahalanobisDistance<T>(*it, pointIndxKNNSearch, *input_cloud);
				if(m_distance < min_madistance_){
					filtered_cloud.push_back(*it);
				}
            }else{	//基于RadiusNN的判定
	    		int neighbors = kd_tree_->radiusSearch(*it, search_radius_dynamic, pointIndxRadiusSearch, pointRadiusDistance);
				if(neighbors >= min_neighbors_){
					filtered_cloud.push_back(*it);
				}
	    }
    }
}

template<typename T>
float MROR::getMahalanobisDistance(const T& centerpoint, vector<int>& pointIndex, typename pcl::PointCloud<T>& input_cloud){
	//计算均值
	float sum_x = 0.0, sum_y = 0.0;
	float x_mean = 0.0, y_mean = 0.0;
	vector<float> x_arr, y_arr;

	for(vector<int>::iterator it = pointIndex.begin(); it != pointIndex.end(); ++it){
		sum_x += input_cloud.points[*it].x;
		sum_y += input_cloud.points[*it].y;
		x_arr.push_back(input_cloud.points[*it].x);
		y_arr.push_back(input_cloud.points[*it].y);
	}
	x_mean = sum_x / pointIndex.size();
	y_mean = sum_y / pointIndex.size();
	
	//计算协方差
	float covxx = 0.0, covxy = 0.0, covyy = 0.0;
	for(int i = 0; i < x_arr.size(); i++){
		covxx += (x_arr[i] - x_mean)*(x_arr[i] - x_mean);
		covxy += (x_arr[i] - x_mean)*(y_arr[i] - y_mean);
		covyy += (y_arr[i] - y_mean)*(y_arr[i] - y_mean);
	}

	//计算协方差矩阵
	Eigen::MatrixXd cov(2,2);
	float length = pointIndex.size() - 1;
	float det = covxx*covyy - covxy*covxy;
	cov << covxx/length, covxy/length, covxy/length, covyy/length;
	// cout << "协方差矩阵：" << endl << cov << endl << endl;

	//计算马氏距离
	Eigen::MatrixXd mm(2,1);
	Eigen::MatrixXd tmp(0,0);
	float ans = 0.0;
	if(fabs(det) > 1e-5){
		// cout << "逆矩阵：" << endl << cov.inverse() << endl;
		mm << (centerpoint.x - x_mean), (centerpoint.y - y_mean);
		tmp = mm.transpose() * cov.inverse() * mm;
		ans = sqrt(sqrt(tmp(0, 0)));
		// cout << "马氏距离为：" << ans << endl;
	}
	return ans;
}