#pragma once
#include <fstream>
#include <limits>
#include <cstdlib>
#include<stdio.h>
#include<string>
#include <exception>
#include <algorithm>
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
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/extract_clusters.h>

#define DEBUG_MODE 1

//PCL TypeDefinations
typedef pcl::PointXYZ PointXYZ;
typedef pcl::PointCloud<PointXYZ>::Ptr PointCloudPtr;
typedef const PointCloudPtr CPointCloudPtr;
typedef std::string str;
typedef const str cStr;
typedef std::vector<float> Distances;

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

    // compute the second order derivative in given range 
    int computeSecondDerivativeInRange(const Distances& minDistances,int first, int last );

    // find the index which is close to given value
    int  findNearest(Distances mindist, double dist);

    // to calculate the centroid of each cluster
   void centroidOfCluster(std::vector<ClusterWithDensity> clusters, PointCloudPtr extractedNewCloud);


public:
    // parametric constructor
    PointCloudPreprocessing(str& filepath);

    // destructor
    virtual ~PointCloudPreprocessing();

    // copy constructor
    PointCloudPreprocessing(const PointCloudPreprocessing&) = delete;

    // move constructor
    PointCloudPreprocessing(PointCloudPreprocessing&&) = delete;

    // assignment operator for copy construcotr
    PointCloudPreprocessing& operator=(const PointCloudPreprocessing&) = delete;

    // assignment operator for move constructor
    PointCloudPreprocessing& operator=(PointCloudPreprocessing&&) = delete;
    
    // to check the density of point cloud is uniform or not
    bool isPointCloudDensityUniform();
    bool isPointCloudDensityUniform(Distances minDis);

    // getter for input point cloud
    PointCloudPtr getPointCloud();

    // computes the minimum distances(sphere of influence) of each point in the point cloud 
    void minDistances(Distances&);

    

    // process of filteration of Point Cloud
    bool filterCloud();

    // getter for Filtered point cloud
    PointCloudPtr getFilteredPointCloud();

    // basic i/o functions
    bool readPTSFile(cStr& filename, PointCloudPtr input_cloud);
    bool readXYZFile(cStr& filename, PointCloudPtr input_cloud);
    bool readPointCloud(str& filename, PointCloudPtr input_cloud);
};