#pragma once
#include <fstream>
#include <sstream>
#include <limits>
#include <cstdlib>
#include<stdio.h>
#include<string>
#include <exception>
#include <algorithm>
#include <boost/functional/hash.hpp> 
#include <tuple>
#include <unordered_map>
#include <queue>
#include <ranges>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/common/io.h>
#include <pcl/common/common.h>
#include <pcl/common/distances.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/filters/voxel_grid_occlusion_estimation.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/segmentation/progressive_morphological_filter.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/features/normal_3d.h>


//PCL TypeDefinations
typedef pcl::PointXYZ PointXYZ;
typedef pcl::PointCloud<PointXYZ>::Ptr PointCloudPtr;
typedef const PointCloudPtr CPointCloudPtr;
typedef std::string str;
typedef const str cStr;
typedef std::vector<double> Distances;

class PointCloudPreprocessing {

private:
    // attributes of PointCloudPreprocessing class
    str m_pointCloudFilePath;
    PointCloudPtr m_pointCloud;
    PointCloudPtr m_filteredPointCloud;
    struct ClusterWithDensity { PointCloudPtr cloud;double density;};

public:
    typedef std::vector<ClusterWithDensity> ClustersWithDensity;

protected:
    // Function to perform DBSCAN-like clustering and return clusters with densities
    void performDBSCAN(CPointCloudPtr cloud, ClustersWithDensity& clusters, double epsilon, int minPts);

    // Function to remove extracted clusters from the original point cloud
    void removeExtractedClusters(CPointCloudPtr cloud, PointCloudPtr remainingCloud, ClustersWithDensity& clusters, double epsilon);

    // Function to compute derivative using finite differences
    void computeSecondDerivative(const Distances& minDistances, Distances& secondDerivative);

    // compute the second order derivative in given ranfe 
    int computeSecondDerivativeInRange(const Distances& minDistances,int first, int last );

    //find the index having close value to given value
    int  findNearest(Distances mindist, double dist);

    //Function to Compute Hausdorff Distance
    std::tuple<float, float, float> ComputeHausdorffDistance(CPointCloudPtr source, CPointCloudPtr target, bool identicalPointClouds);


public:
    //default constructor
  //  PointCloudPreprocessing()=default;

    //parametric constructor
    PointCloudPreprocessing(str& filepath);

    //destructor
    virtual ~PointCloudPreprocessing();

    // copy constructor
    PointCloudPreprocessing(const PointCloudPreprocessing&) = delete;

    // move constructor
    PointCloudPreprocessing(PointCloudPreprocessing&&) = delete;

    // assignment operator for copy construcotr
    PointCloudPreprocessing& operator=(const PointCloudPreprocessing&) = delete;

    // assignment operator for move constructor
    PointCloudPreprocessing& operator=(PointCloudPreprocessing&&) = delete;
    
    //to check the uniformity of point cloud
    bool isPointCloudUniformWithDensity();

    // getter for input point cloud
    PointCloudPtr getPointCloud();

    //Function that computes the minimum distance of other point to each point in the point cloud 
    void minDistances(Distances&);

    // process of filteration of Point Cloud
    void filterCloud();

    // getter for Filtered point cloud
    PointCloudPtr getFilteredPointCloud();

    // basic i/o functions
    void readPTSFile(cStr& filename, CPointCloudPtr input_cloud);
    void readXYZFile(cStr& filename, CPointCloudPtr input_cloud);
    void readPointCloud(str& filename, CPointCloudPtr input_cloud);
};