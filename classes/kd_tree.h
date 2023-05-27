#ifndef KD_TREE_H
#define KD_TREE_H

#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <iostream>
#include "vector.h"

class KDNode {
public:
    KDNode() : left(NULL), right(NULL) {}
    KDNode(const Vector& p, int idx, int axis, double mean, KDNode* left, KDNode* right) : point(p), idx(idx), axis(axis), mean(mean), left(left), right(right) {}
    Vector point;
    int idx, axis;
    double mean;
    KDNode* left;
    KDNode* right;
};

class KDTree {
public:
    KDNode* root;

    KDTree() : root(NULL) {}
    
    KDTree(std::vector<std::pair<Vector, int> >& points) {
        root = buildTreeRecursive(points);
    }

    void findKNearestNeighbors(std::vector<int>& result, const Vector& target, int k) const {
        std::priority_queue<std::pair<double, int>> nearestNeighbors;
        //std::cout << "Query point: " << target[0] << ' ' << target[1] << ' ' << target[2] << '\n';
        findKNearestNeighborsRecursive(root, target, k, nearestNeighbors);
        int i = 0;
        result.resize(nearestNeighbors.size());
        while (!nearestNeighbors.empty()) {
            //std::cout << nearestNeighbors.top().first << '\n';
            result[i++] = nearestNeighbors.top().second;
            nearestNeighbors.pop();
        }
        std::reverse(result.begin(), result.end()); // To get the nearest neighbors in ascending order
    }

private:
    KDNode* buildTreeRecursive(std::vector<std::pair<Vector, int> >& points, int depth = 0, int start_index = 0, int end_index = -1)
    {
        if (end_index == -1) end_index = points.size();
        if (end_index <= start_index) return NULL;
        if (end_index - start_index == 1) return new KDNode(points[start_index].first, points[start_index].second, 0, 0, NULL, NULL);

        double
            min_x = std::numeric_limits<double>::max(), 
            max_x = std::numeric_limits<double>::min(), 
            min_y = std::numeric_limits<double>::max(), 
            max_y = std::numeric_limits<double>::min(), 
            min_z = std::numeric_limits<double>::max(), 
            max_z = std::numeric_limits<double>::min();
        for (int i = start_index; i < end_index; i++)
        {
            min_x = std::min(min_x, points[i].first[0]);
            max_x = std::max(max_x, points[i].first[0]);
            min_y = std::min(min_y, points[i].first[1]);
            max_y = std::max(max_y, points[i].first[1]);
            min_z = std::min(min_z, points[i].first[2]);
            max_z = std::max(max_z, points[i].first[2]);
        }
        double dx = max_x - min_x, dy = max_y - min_y, dz = max_z - min_z;
        int axis = (dx > dy && dx > dz) ? 0 : (dy > dz) ? 1 : 2;
        double mean = (axis == 0) ? (min_x + max_x)/2 : (axis == 1) ? (min_y + max_y)/2 : (min_z + max_z)/2;
        int pivot = start_index;
        for (int i = start_index; i < end_index; i++)
        {
            if (points[i].first[axis] < mean)
            {
                std::swap(points[i], points[pivot]);
                pivot++;
            }
        }
        KDNode* left_node = buildTreeRecursive(points, depth + 1, start_index, pivot);
        KDNode* right_node = buildTreeRecursive(points, depth + 1, pivot+1, end_index);
        return new KDNode(points[pivot].first, points[pivot].second, axis, mean, left_node, right_node);
    }

    void findKNearestNeighborsRecursive(
        const KDNode* node, const Vector& target, int k,
        std::priority_queue<std::pair<double, int>>& nearestNeighbors) const
    {
        if (node == NULL) return;
        //std::cout << "Reaching " << node->point[0] << ", " << node->point[1] << ", " << node->point[2] << ' ' << node->axis << '\n';
        double axis_distance = target[node->axis] - node->mean;
        if (axis_distance < 0) {
            findKNearestNeighborsRecursive(node->left, target, k, nearestNeighbors);
        
            double distance = (node->point - target).norm2();
            nearestNeighbors.push(std::make_pair(distance, node->idx));
            //std::cout << "Attempting to insert " << node->point[0] << ", " << node->point[1] << ", " << node->point[2] << ' ' << distance << '\n';
            if (nearestNeighbors.size() > k) nearestNeighbors.pop(); 
        
            if (nearestNeighbors.size() < k || axis_distance * axis_distance < nearestNeighbors.top().first)
                findKNearestNeighborsRecursive(node->right, target, k, nearestNeighbors);
        }
        else {
            findKNearestNeighborsRecursive(node->right, target, k, nearestNeighbors);
        
            double distance = (node->point - target).norm2();
            nearestNeighbors.push(std::make_pair(distance, node->idx));
            //std::cout << "Attempting to insert " << node->point[0] << ", " << node->point[1] << ", " << node->point[2] << ' ' << distance << '\n';
            if (nearestNeighbors.size() > k) nearestNeighbors.pop();
            
            if (nearestNeighbors.size() < k || axis_distance * axis_distance < nearestNeighbors.top().first)
                findKNearestNeighborsRecursive(node->left, target, k, nearestNeighbors);
        }
    }
};

#endif