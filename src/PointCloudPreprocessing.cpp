#include "PointCloudPreprocessing.h"

namespace 
{
    str getFileExtension(cStr& filename) {
        size_t dotIndex = filename.find_last_of('.');
        if (dotIndex != str::npos) {
            return filename.substr(dotIndex + 1);
        }
        return ""; //  return empty string if no extension found
    }
}


bool PointCloudPreprocessing::readPTSFile(cStr& filename, PointCloudPtr input_cloud)
{
    std::vector<PointXYZ> points;

    try
    {
        std::ifstream file(filename);
        if (!file.is_open()) 
        {
            if (DEBUG_MODE == 1)
                std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        str line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue; // Skip comments and empty lines
            PointXYZ pt;
            std::istringstream iss(line);
            if (!(iss >> pt.x >> pt.y >> pt.z)) 
            {
                if (DEBUG_MODE == 1)
                    std::cerr << "Error: Invalid point data in file " << filename << std::endl;
            }
            input_cloud->push_back(pt);
            points.push_back(pt);
        }
        return true;

    }
    catch (...) 
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }
    return true;
}


bool PointCloudPreprocessing::readXYZFile(cStr& filename, PointCloudPtr input_cloud)
{
    str line;
    std::ifstream file(filename);
    Distances floatValues;

    try
    {
        if (file.is_open())
        {
            while (getline(file, line))
            {
                std::istringstream ss(line);
                double floatValue;
                int c = 0;
                double x{};
                double y{};
                pcl::PointXYZ pts;
                while (ss >> floatValue)
                {
                    // Successfully read a float value
                    floatValues.push_back(floatValue);
                    pcl::PointXYZ p;
                    if (c == 0)
                    {
                        // x
                        x = floatValue;
                        pts.x = floatValue;
                    }
                    else if (c == 1)
                    {
                        // y
                        y = floatValue;
                        pts.y = floatValue;
                    }
                    else if (c == 2)
                    {
                        // z
                        pts.z = floatValue;
                    }
                    else
                    {
                        // other attributes of points like rgb
                    }

                    ++c;
                }
                input_cloud->push_back(pts);
            }
            file.close();
        }
        return true;
    }
    catch (...) 
    { 
        if (DEBUG_MODE == 1)
        std::cerr << "Error: Could not open file " << filename << std::endl; 
        return false;
    }
    return true;
}


bool PointCloudPreprocessing::readPointCloud(str& filename, PointCloudPtr input_cloud)
{
    bool errorFlag = false;
    str extension = getFileExtension(filename);

    if (extension == "pts") 
    {
        errorFlag=readPTSFile(filename, input_cloud);
    }
    else if (extension == "pcd") 
    {
        pcl::io::loadPCDFile(filename, *input_cloud);
    }
    else if (extension == "xyz")
    {
        errorFlag = readXYZFile(filename, input_cloud);
        if (DEBUG_MODE == 1 && errorFlag)
            std::cout << "Error: Could not open file" << std::endl;
    }
    else
    {
        std::cout << "Unknown file format" << std::endl;
    }

    return errorFlag;
}


void PointCloudPreprocessing::minDistances(Distances& minDis)
{
        std::vector<float> distances;
        pcl::search::KdTree<PointXYZ, pcl::KdTreeFLANN<PointXYZ>> kdTree;

        kdTree.setInputCloud(this->m_pointCloud);

        for (const auto& itr : *this->m_pointCloud)
        {
            // O(logn)
            pcl::Indices index;
            std::vector<float> squaredDistance;

            auto  minDist = std::numeric_limits<float>::min();
            minDist = 100000;
            kdTree.nearestKSearch(itr, 3, index, squaredDistance);
            for (int8_t i = 0; i < 3; ++i)
            {
                if (squaredDistance.at(i) > 0)
                {
                    if (squaredDistance.at(i) < minDist) 
                    {
                        minDist = squaredDistance.at(i);
                    }
                }
            }

            minDis.push_back(minDist);
        }

    std::sort(minDis.begin(), minDis.end());
}


bool  PointCloudPreprocessing::isPointCloudDensityUniform()
{
    if (this->m_pointCloud->empty())
        return false;
    
    bool uniformDensityFlag = false;
    Distances minDis;
    
    minDistances(minDis);
   
    double standardDevitation = 0, mean=0;
    pcl::getMeanStdDev(minDis, mean, standardDevitation);
    
    // uniformity check
    if ((standardDevitation <=DBL_EPSILON) ||(standardDevitation <= minDis.at(0)))
        uniformDensityFlag = true;
    
    return uniformDensityFlag;
}


bool  PointCloudPreprocessing::isPointCloudDensityUniform(Distances minDis)
{
    if (this->m_pointCloud->empty())
        return false;

    bool uniformDensityFlag = false;

    double standardDevitation = 0, mean = 0;
    pcl::getMeanStdDev(minDis, mean, standardDevitation);

    // uniformity check
    if ((standardDevitation <= DBL_EPSILON) || (standardDevitation <= minDis.at(0)))
        uniformDensityFlag = true;

    return uniformDensityFlag;
}


 PointCloudPreprocessing::PointCloudPreprocessing(str& filepath):
     m_pointCloud(new pcl::PointCloud<PointXYZ>),
     m_filteredPointCloud(new pcl::PointCloud<PointXYZ>),
     m_pointCloudFilePath(filepath)
 {
     readPointCloud(this->m_pointCloudFilePath, this->m_pointCloud);
 }


PointCloudPreprocessing::~PointCloudPreprocessing() 
{
}


PointCloudPtr PointCloudPreprocessing::getPointCloud()
{
    return this->m_pointCloud;
}


void PointCloudPreprocessing::performDBSCAN(CPointCloudPtr cloud, ClustersWithDensity& clusters, double epsilon, int minPts) 
{
    // Create cluster extraction object
    pcl::EuclideanClusterExtraction<PointXYZ> extractedClusters;
    extractedClusters.setClusterTolerance(epsilon);
    extractedClusters.setMinClusterSize(minPts);
    extractedClusters.setInputCloud(cloud);

    // Extract clusters
    std::vector<pcl::PointIndices> clusterIndices;
    extractedClusters.extract(clusterIndices);

    // Calculate density for each cluster and store clusters
    for (const auto& indices : clusterIndices)
    {
        ClusterWithDensity cluster;
        cluster.cloud.reset(new pcl::PointCloud<PointXYZ>);
        cluster.density = indices.indices.size() / extractedClusters.getClusterTolerance();

        for (const auto& idx : indices.indices) 
        {
            cluster.cloud->emplace_back((*cloud)[idx]);
        }

        clusters.emplace_back(std::move(cluster));	
    }
}


void PointCloudPreprocessing::removeExtractedClusters(CPointCloudPtr cloud, PointCloudPtr remainingCloud, ClustersWithDensity& clusters, double epsilon) 
{
    for (const auto& point : cloud->points)
    {
        bool inExtractedCluster = false;
        for (const auto& cluster : clusters)
        {
            for (const auto& extractedPoint : cluster.cloud->points) 
            {
                if (pcl::euclideanDistance(point, extractedPoint) < epsilon)
                {
                    inExtractedCluster = true;
                    break;
                }
            }

            if (inExtractedCluster)
            {
                break;
            }
        }

        if (!inExtractedCluster)
        {
            remainingCloud->emplace_back(point);
        }
    }
}


void PointCloudPreprocessing::computeSecondDerivative(const Distances& minDistances, Distances& secondDerivative)
{
    for (size_t i = 1; i < minDistances.size() - 1; ++i)
    {
        secondDerivative[i] = minDistances[i + 1] - 2 * minDistances[i] + minDistances[i - 1];
    }
}


int PointCloudPreprocessing::computeSecondDerivativeInRange(const Distances& minDistances, int first, int last)
{
    if (last- first<0 && minDistances.size()<last) {
        if(DEBUG_MODE)
        std::cout << "Give valid Range!!!!" << std::endl;
        return ;
    }
    else {
        std::vector<double> secondDerivative;
        secondDerivative.reserve(minDistances.size());
        double tempDiffer = 0;

        for (size_t i = 0; i < minDistances.size(); ++i)
        {
            if (i >= first && i <= last)
            {
                auto temp = minDistances[i + 1] - 2 * minDistances[i] + minDistances[i - 1];
                if (tempDiffer < temp)
                {
                    tempDiffer = temp;
                }

                secondDerivative.push_back(temp);
            }
        }

        auto maxelementit = std::max_element(std::begin(secondDerivative), std::end(secondDerivative));
        return std::distance(secondDerivative.begin(), maxelementit);
    }
}


int  PointCloudPreprocessing::findNearest(Distances mindist, double dist) 
{
    std::vector<double>diffMindistToMean;
    for (int i = 0; i < mindist.size(); ++i)
    {
        diffMindistToMean.emplace_back(std::abs(mindist.at(i) - dist));
    }

    auto minelementit = std::min_element(std::begin(diffMindistToMean), std::end(diffMindistToMean));
    return static_cast<int>(std::distance(diffMindistToMean.begin(), minelementit));
}


void PointCloudPreprocessing::centroidOfCluster(ClustersWithDensity clusters, PointCloudPtr extractedNewCloud)
{
    for (const auto& cluster : clusters)
    {
        float cx = 0, cy = 0, cz = 0;

        for (auto point : *cluster.cloud)
        {
            cx = cx + point.x;
            cy = cy + point.y;
            cz = cz + point.z;
        }
        cx = cx / cluster.cloud->size();
        cy = cy / cluster.cloud->size();
        cz = cz / cluster.cloud->size();

        extractedNewCloud->emplace_back(PointXYZ(cx, cy, 0));
    }
}


bool PointCloudPreprocessing::filterCloud()
{
    Distances mindist;
    minDistances(mindist);

    if (isPointCloudDensityUniform(mindist))
    {
        if(DEBUG_MODE)
        std::cout << "Point Cloud having Uniform Density!!!" << std::endl;
        pcl::copyPointCloud(*this->m_pointCloud, *this->m_filteredPointCloud);
        return false;
    }
    else 
    { 
        //median
        auto m = mindist.begin() + mindist.size() / 2;
        std::nth_element(mindist.begin(), m, mindist.end());
        double median = mindist[mindist.size() / 2];
        auto idxMedian = findNearest(mindist, mindist[mindist.size() / 2]);

        //mean
        auto mean = std::accumulate(std::begin(mindist), std::end(mindist), 0.0) / std::size(mindist);
        auto idxMean = findNearest(mindist, mean);

        //index of threshold value for dbscan mehtod
        auto idxThreshold = computeSecondDerivativeInRange(mindist, std::min(idxMean, idxMedian), std::max(idxMean, idxMedian));

        // Perform DBSCAN-like clustering with specified parameters
        double epsilon = mindist.at(idxThreshold); // neighborhood radius
        int minPts = 2;      // minimum number of points to form a cluster
        std::vector<ClusterWithDensity> clusters;
        performDBSCAN(this->m_pointCloud, clusters, epsilon, minPts);


        // Combine extracted clusters into one point cloud
        PointCloudPtr extractedNewCloud(new pcl::PointCloud<PointXYZ>);
        extractedNewCloud->reserve(m_pointCloud->size());

        //centroid of each cluster
        centroidOfCluster(clusters, extractedNewCloud);

        // Remove the extracted clusters from the original point cloud
        PointCloudPtr remainingCloud(new pcl::PointCloud<PointXYZ>);
        removeExtractedClusters(this->m_pointCloud, remainingCloud, clusters, epsilon);

        // Save the remaining point cloud
        //pcl::io::savePCDFile("..\\..\\remaining_cloud.pcd", *remainingCloud);
        *remainingCloud += *extractedNewCloud;

        // another method to concatenate centroid of clusters and remaining point cloud
        //pcl::concatenateFields(*remainingCloud, *extractedNewCloud, *remainingCloud);

        pcl::copyPointCloud(*remainingCloud, *this->m_filteredPointCloud);
       
        return true;
    }
}


PointCloudPtr PointCloudPreprocessing::getFilteredPointCloud() 
{
    if (this->m_filteredPointCloud->size() < 2)
    {
        if(DEBUG_MODE)
        std::cout << "ERROR!!! Size of  point cloud: " << this->m_filteredPointCloud->size() << std::endl;
        return nullptr;
    }
    else 
    {
        return this->m_filteredPointCloud;
    }
};