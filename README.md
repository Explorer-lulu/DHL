# DHL

## Project Description

This is the source code of our ICDE 2023 manuscript "Double Hierarchical Labeling Shortest Distance Querying in Time-dependent Road Networks". The related algorithms have been explained in the paper.  We will do our best to maintain the entire project. If our work enlightens you, please cite our paper.

## Data Description

Data folder stores the all data of DHL. We also give a sample data for the code.

sample file contain a sample graph with 10,000 vertices and 50,000 edges. Each row represents (Vi, Vj, Weight at t=0) i.e., "1, 2, 5" means the weight between V1 and V2 at 0 is 5.  

query file represents the query set which contains 1,000 query pairs (Vi, Vj)

updata file indicates the time-dependent weight. For example, if "1, 2, 10" is recorded in updata, it means that the weight between v1 and v2 will change to 10 at time 5.

dataorder_sample file contains all vertex's hierarchy level number of the graph.

## Running 

All codes is made by C++, compiled with G++ under the O3 optimization. You can run the code after make.
