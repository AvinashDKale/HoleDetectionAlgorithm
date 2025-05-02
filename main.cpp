#include"PointCloudPreprocessing.h"
#include "FindHole.h"
#include <chrono>
#include <ctime>
int main() {

    // Start timing
    auto start = std::chrono::steady_clock::now();
    std::cout << "Hole Detection Algorithm Project..." << std::endl;
    std::string filename = "E:\\Avinash C drive\\Desktop\\textPCDFormatDATA1.pcd";// sample input with non-uniform density
   // filename = "E:\\Avinash C drive\Desktop\\Algorithm-Hole Detection\\HoleDetectionAlgorithmProject\\emptyPointcloud.xyz";// to test the empty point cloud 

    PointCloudPreprocessing p(filename);
    p.filterCloud();
    auto d = p.getPointCloud();
   auto filteredpoints= p.getFilteredPointCloud();
    std::cout << "size:" << d->size() << std::endl;                                       









    std::vector<Point> points;
    std::vector<std::vector<Point>> holes;
    Edges holeEdges;


    pcl::PointCloud<PointXYZ>::Ptr input_cloud(new pcl::PointCloud<PointXYZ>);

   // pcl::io::loadPCDFile(filename, *input_cloud);
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











    for (int i = 0; i < input_cloud->size(); ++i) {
        points.push_back(Point(input_cloud->at(i).x, input_cloud->at(i).y));
    }




    //for (int i = 0; i < filteredpoints->size();++i) {
    //    points.push_back(Point(filteredpoints->at(i).x, filteredpoints->at(i).y));
    //}
   

    FindHole h(points);
    auto numberofholes = h.getholes(holes, holeEdges);

    std::cout << "InputPoints:" << input_cloud->size() << std::endl;
    std::cout << "numberofholes:" << numberofholes << std::endl;
    std::cout << "holeEdges:" << holeEdges.size() << std::endl;




   
    // Stop timing
    auto end = std::chrono::steady_clock::now();

    // time in milliseconds
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << elapsed.count() << " ms" << std::endl;
    return 0;
}
