#include "PointCloudPreprocessing.h"


str getFileExtension(cStr& filename) {
    size_t dotIndex = filename.find_last_of('.');
    if (dotIndex != str::npos) {
        return filename.substr(dotIndex + 1);
    }
    return ""; // No extension found
}

void PointCloudPreprocessing::readPTSFile(cStr& filename, CPointCloudPtr input_cloud) {
    std::vector<PointXYZ> points;

    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
        }

        str line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue; // Skip comments and empty lines
            PointXYZ pt;
            std::istringstream iss(line);
            if (!(iss >> pt.x >> pt.y >> pt.z)) {
                std::cerr << "Error: Invalid point data in file " << filename << std::endl;
                // return points;
            }
            input_cloud->push_back(pt);
            points.push_back(pt);
        }

    }
    catch (...) { std::cerr << "Error: Could not open file " << filename << std::endl; }
}

void PointCloudPreprocessing::readXYZFile(cStr& filename, CPointCloudPtr input_cloud) {
    str line;
    std::ifstream file(filename);
    Distances floatValues;
    try {
        if (file.is_open()) {
            while (getline(file, line)) {
                std::istringstream ss(line);
                double floatValue;
                int c = 0;
                double x{};
                double y{};
                pcl::PointXYZ pts;
                while (ss >> floatValue) {
                    // Successfully read a float value
                    floatValues.push_back(floatValue);

                    pcl::PointXYZ p;
                    if (c == 0) {
                        // std::cout <<"x:" << floatValue;
                        x = floatValue;
                        pts.x = floatValue;


                    }
                    else if (c == 1) {
                        //std::cout << " y:" << floatValue ;
                        y = floatValue;
                        pts.y = floatValue;
                    }
                    else if (c == 2) {
                        //std::cout << " z:" << floatValue ;
                        pts.z = floatValue;
                    }
                    else {
                        // std::cout << " z:" << floatValue << "\n";
                    }

                    ++c;
                }
                input_cloud->push_back(pts);

            }
            file.close();
        }
    }
    catch (...) { std::cerr << "Error: Could not open file " << filename << std::endl; }
}

void PointCloudPreprocessing::readPointCloud(str& filename, CPointCloudPtr input_cloud) {
    str extension = getFileExtension(filename);
    if (extension == "pts") {
        readPTSFile(filename, input_cloud);
    }
    else if (extension == "pcd") {
        pcl::io::loadPCDFile(filename, *input_cloud);
    }
    else if (extension == "xyz") {
        readXYZFile(filename, input_cloud);
    }
    else {
        std::cout << "Unknown file format" << std::endl;
    }
}

void PointCloudPreprocessing::minDistances(Distances& minDis)
{

    std::vector<pcl::PointXYZ> points;
    double minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
    double minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
    for (int i = 0; i < this->m_pointCloud->size(); i++)
    {
        auto p = this->m_pointCloud->at(i);

        if (p.x > maxX)
        {
            maxX = p.x;
        }
        if (p.y > maxY)
        {
            maxY = p.y;
        }
        if (p.x < minX)
        {
            minX = p.x;
        }
        if (p.y < minY)
        {
            minY = p.y;
        }

        points.push_back(p);
    }


    std::vector<std::pair<double, int>> freqOfMinDist;
    int cc = 0;
    auto ddy = maxY - minY;
    ddy = ddy * ddy;
    auto ddx = maxX - minX;
    ddx = ddx * ddx;
    auto dimtr = std::sqrt(ddx + ddy);
    //std::vector<double> minDis;
    for (int i = 0; i < points.size(); ++i) {
        double m = dimtr;
        for (int j = 0; j < points.size(); ++j) {
            auto d = pcl::euclideanDistance(points.at(i), points.at(j));
            if (!d == 0 && d < m) {
                m = d;
            }
        }
        minDis.push_back(m);

        ++cc;
    }
    std::sort(minDis.begin(), minDis.end());


    // auto source = this->m_pointCloud;

    auto DistanceComputation = [&](const CPointCloudPtr& source) -> std::vector<float>
    {
        std::vector<float> distances;
        pcl::search::KdTree<PointXYZ, pcl::KdTreeFLANN<PointXYZ>> kdTree;

        kdTree.setInputCloud(source);
        for (const auto& itr : *source)
        {
            // O(logn)
            pcl::Indices index;
            std::vector<float> squaredDistance;

            // computing minimum distance between any 2 points in the same point cloud
            // 2 points are used (same point and another point) because nearest point
            // will be itself in which case distance = 0
            //index.resize(2);
            //squaredDistance.resize(2);
            auto  minDist = std::numeric_limits<float>::min();
            minDist = 100000;
            kdTree.nearestKSearch(itr, 3, index, squaredDistance);
            for (int8_t i = 0; i < 3; ++i)
            {
                if (squaredDistance.at(i) > 0)
                {


                    if (squaredDistance.at(i) < minDist) {
                        minDist = squaredDistance.at(i);
                    }

                }
            }
            distances.push_back(minDist);


        }
        return distances;
    };




    auto m = DistanceComputation(this->m_pointCloud);



    std::cout << m.size() << std::endl;
}

bool  PointCloudPreprocessing::isPointCloudUniformWithDensity() {
    if (this->m_pointCloud->size() == 0)
    return - 1;
    bool f = false;
    Distances minDis;
    minDistances(minDis);
    double mean = 0;
    for (auto m : minDis) {
        mean = mean + m;
    }
    mean = mean / minDis.size();

    //std deviation form mean for each point and whole
    Distances deviations(minDis.size());
    double standardDevitation = 0;
    for (auto m : minDis) {
        deviations.emplace_back(std::abs(mean - m));
        standardDevitation = standardDevitation + (std::abs(mean - m) * std::abs(mean - m));
      //  std::cout <<"dist:"<< m<< " Devitation:" << std::abs(mean - m) << std::endl;
    }

    standardDevitation = standardDevitation / minDis.size();
    standardDevitation = std::sqrt(standardDevitation);
    std::cout << "mean:" << mean << std::endl;
    std::cout << "Standard Devitation:" << standardDevitation << std::endl;

    if ((standardDevitation == 0) ||(standardDevitation <= minDis.at(0)))f = true;

    auto m = minDis.begin() + minDis.size() / 2;
    std::nth_element(minDis.begin(), m, minDis.end());
    std::cout << "\nThe median is " << minDis[minDis.size() / 2] << '\n';

    auto mean1 = std::accumulate(std::begin(minDis), std::end(minDis), 0.0) / std::size(minDis);
    std::cout << "\nThe mean is " << mean1 << std::endl;

    return f;
}

 PointCloudPreprocessing::PointCloudPreprocessing(str& filepath): m_pointCloud(new pcl::PointCloud<PointXYZ>),m_filteredPointCloud(new pcl::PointCloud<PointXYZ>),m_pointCloudFilePath(filepath) {
     readPointCloud(this->m_pointCloudFilePath, this->m_pointCloud);
}

PointCloudPreprocessing::~PointCloudPreprocessing() {

}

PointCloudPtr PointCloudPreprocessing::getPointCloud() {

    return this->m_pointCloud;
}

void PointCloudPreprocessing::performDBSCAN(CPointCloudPtr cloud, ClustersWithDensity& clusters, double epsilon, int minPts) {
    // Function to perform DBSCAN-like clustering and return clusters with densities

    // Create cluster extraction object
    pcl::EuclideanClusterExtraction<PointXYZ> extractedClusters;
    extractedClusters.setClusterTolerance(epsilon);
    extractedClusters.setMinClusterSize(minPts);
    extractedClusters.setInputCloud(cloud);

    // Extract clusters
    std::vector<pcl::PointIndices> clusterIndices;
    extractedClusters.extract(clusterIndices);

    // Calculate density for each cluster and store clusters
    for (const auto& indices : clusterIndices) {
        ClusterWithDensity cluster;
        cluster.cloud.reset(new pcl::PointCloud<PointXYZ>);
        cluster.density = indices.indices.size() / extractedClusters.getClusterTolerance();
        for (const auto& idx : indices.indices) {
            cluster.cloud->push_back((*cloud)[idx]);
        }
        clusters.push_back(cluster);
    }
}

void PointCloudPreprocessing::removeExtractedClusters(CPointCloudPtr cloud, PointCloudPtr remainingCloud, ClustersWithDensity& clusters, double epsilon) {
    // Function to remove extracted clusters from the original point cloud

    for (const auto& point : cloud->points) {
        bool inExtractedCluster = false;
        for (const auto& cluster : clusters) {
            for (const auto& extractedPoint : cluster.cloud->points) {
                if (pcl::euclideanDistance(point, extractedPoint) < epsilon) {
                    inExtractedCluster = true;
                    break;
                }
            }
            if (inExtractedCluster) {
                break;
            }
        }
        if (!inExtractedCluster) {
            //remainingCloud->push_back(point);
            remainingCloud->emplace_back(point);
        }
    }
}

void PointCloudPreprocessing::computeSecondDerivative(const Distances& minDistances, Distances& secondDerivative) {
   // std::vector<double> secondDerivative(minDistances.size());

    // Compute second derivative using central difference method

    for (size_t i = 1; i < minDistances.size() - 1; ++i) {
        secondDerivative[i] = minDistances[i + 1] - 2 * minDistances[i] + minDistances[i - 1];
    }

}

int PointCloudPreprocessing::computeSecondDerivativeInRange(const Distances& minDistances, int first, int last) {
    std::vector<double> secondDerivative;
    double tempDiffer = 0;
    for (size_t i = 0; i < minDistances.size(); ++i) {
        if (i >= first && i <= last) {
            auto temp = minDistances[i + 1] - 2 * minDistances[i] + minDistances[i - 1];
            if (tempDiffer< temp) {
                tempDiffer = temp;
            }
            
            secondDerivative.emplace_back(temp);
        }
        else {
            secondDerivative.emplace_back(0);
        }
    }
    auto maxelementit = std::max_element(std::begin(secondDerivative), std::end(secondDerivative));
    return (int)std::distance(secondDerivative.begin(), maxelementit);
}

std::tuple<float, float, float> PointCloudPreprocessing::ComputeHausdorffDistance(CPointCloudPtr source, CPointCloudPtr target, bool identicalPointClouds)
{
    // HausdorffDistance is the greatest of the minimum distance between 2 discrete sets
    // the same code can be used for computing minimum distance in identical point clouds
    // and also the smallest of the minimum distance between 2 non-identical point clouds
    // when the point clouds are identical, minimum and maximum are equal values
    float minDist = std::numeric_limits<float>::max(), maxDist = std::numeric_limits<float>::min();
    float avgeDist = 0.0f;
    using pair = std::pair<float, float>;
    using triplet = std::array<float, 3>;
    pair minmax{ minDist, maxDist };
    triplet distances{ minDist, maxDist, avgeDist };

    auto FindMinMax = [](const float thisDistance, pair& distances) -> void
    {
        if (thisDistance < distances.first)
            distances.first = thisDistance;
        if (thisDistance > distances.second)
            distances.second = thisDistance;
    };

    auto DistanceComputation = [&](const pcl::PointCloud<PointXYZ>::Ptr& source,
        const pcl::PointCloud<PointXYZ>::Ptr& target) -> triplet
    {
        // reinitialize the values
        minmax = { minDist, maxDist };
        avgeDist = 0.0f;
        pcl::search::KdTree<pcl::PointXYZ, pcl::KdTreeFLANN<pcl::PointXYZ>> kdTree;
        kdTree.setInputCloud(source);
        for (const auto& itr : *target)
        {
            // O(logn)
            pcl::Indices index;
            std::vector<float> squaredDistance;
            if (identicalPointClouds)
            {
                // computing minimum distance between any 2 points in the same point cloud
                // 2 points are used (same point and another point) because nearest point
                // will be itself in which case distance = 0
                //index.resize(2);
                //squaredDistance.resize(2);
                kdTree.nearestKSearch(itr, 3, index, squaredDistance);
                for (int8_t i = 0; i < 3; ++i)
                {
                    if (squaredDistance.at(i) > 0)
                    {
                        avgeDist += std::sqrtf(squaredDistance.at(i));
                        FindMinMax(std::sqrtf(squaredDistance.at(i)), minmax);
                        break;
                    }
                }
            }
            else
            {
                // computing greatest of minimum distance (2 nearest points) in different point clouds
                // computing smallest of minimum distance (2 nearest points) in different point clouds
                // the former is the real Hausdorff distance
                index.resize(1);
                squaredDistance.resize(1);
                kdTree.nearestKSearch(itr, 1, index, squaredDistance);
                avgeDist += std::sqrtf(squaredDistance.at(0));
                FindMinMax(std::sqrtf(squaredDistance.at(0)), minmax);
            }
        }

        avgeDist /= target->size();
        distances[0] = minmax.first, distances[1] = minmax.second, distances[2] = avgeDist;
        return distances;
    };

    // compare point cloud 1 to point cloud 2 and then point cloud 2 to point cloud 1
    // take the average of minimum and maximum distance as their mean
    triplet distancePair[2];
    if (!identicalPointClouds)
    {
        distancePair[0] = DistanceComputation(source, target);
        distancePair[1] = DistanceComputation(target, source);
    }
    else
    {
        distancePair[0] = distancePair[1] = DistanceComputation(source, target);
    }

    minDist = 0.5f * (distancePair[0].at(0) + distancePair[1].at(0));
    maxDist = 0.5f * (distancePair[0].at(1) + distancePair[1].at(1));
    avgeDist = 0.5f * (distancePair[0].at(2) + distancePair[1].at(2));
    return { minDist, maxDist, avgeDist };
}

int  PointCloudPreprocessing::findNearest(Distances mindist, double dist) {
    std::vector<double>diffMindistToMean;
    for (int i = 0; i < mindist.size(); ++i) {
        diffMindistToMean.emplace_back(std::abs(mindist.at(i) - dist));
    }
    auto minelementit = std::min_element(std::begin(diffMindistToMean), std::end(diffMindistToMean));
    
    return (int)std::distance(diffMindistToMean.begin(), minelementit);
}

void PointCloudPreprocessing::filterCloud() {

    if (isPointCloudUniformWithDensity()) {
        std::cout << "Point Cloud having Uniform Density!!!" << std::endl;
        this->m_filteredPointCloud = this->m_pointCloud;
    }
    else{
        // Hausdorff distance is computed after downsampling because nearestKSearch or radiusSearch queries are expensive
       // float maxSampleDist = 0.0f; float minSampleDistance = 0.0f; float inlierDistance = 0.0f;
      //  std::tie(minSampleDistance, maxSampleDist, inlierDistance) = ComputeHausdorffDistance(this->m_pointCloud, this->m_pointCloud, 1);

        Distances mindist;
        minDistances(mindist);
        
        //median
        auto m = mindist.begin() + mindist.size() / 2;
        std::nth_element(mindist.begin(), m, mindist.end());
        double median = mindist[mindist.size() / 2];
        auto idxMeadian = findNearest(mindist, mindist[mindist.size() / 2]);

        //mean
        auto mean = std::accumulate(std::begin(mindist), std::end(mindist), 0.0) / std::size(mindist);
        auto idxMean = findNearest(mindist, mean);
        
        //index of threshold value for dbscan mehtod
        auto idxThreshold = computeSecondDerivativeInRange(mindist, std::min(idxMean, idxMeadian), std::max(idxMean, idxMeadian));


        // Perform DBSCAN-like clustering with specified parameters
        double epsilon = mindist.at(idxThreshold); // neighborhood radius
        int minPts = 2;      // minimum number of points to form a cluster
        std::vector<pcl::PointIndices> cluster_indices;
        std::vector<ClusterWithDensity> clusters;
        performDBSCAN(this->m_pointCloud, clusters, epsilon, minPts);


        // Combine extracted clusters into one point cloud
        PointCloudPtr extractedNewCloud(new pcl::PointCloud<PointXYZ>);
        for (const auto& cluster : clusters) {
            double cx = 0, cy = 0, cz = 0;
            for (auto p : *cluster.cloud) {
                cx = cx + p.x;
                cy = cy + p.y;
                cz = cz + p.z;
            }
            cx = cx / cluster.cloud->size();
            cy = cy / cluster.cloud->size();
            cz = cz / cluster.cloud->size();
            extractedNewCloud->push_back(PointXYZ(cx, cy, 0));
        }


        // Remove the extracted clusters from the original point cloud
        PointCloudPtr remainingCloud(new pcl::PointCloud<PointXYZ>);
        removeExtractedClusters(this->m_pointCloud, remainingCloud, clusters, epsilon);

        // Save the remaining point cloud
        //pcl::io::savePCDFile("E:\\Avinash C drive\\Desktop\\clusters\\remaining_cloud.pcd", *remainingCloud);
        *remainingCloud += *extractedNewCloud;

        this->m_filteredPointCloud = remainingCloud;
    }

}

PointCloudPtr PointCloudPreprocessing::getFilteredPointCloud() {
    if (this->m_filteredPointCloud->size() < 2){

        std::cout << "ERROR!!! Size of  point cloud: " << this->m_filteredPointCloud->size() << std::endl;
        return nullptr;
    }
    else {
        return this->m_filteredPointCloud;
    }
       

};