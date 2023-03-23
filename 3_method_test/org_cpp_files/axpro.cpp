#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <limits>
#include <queue>
#include <algorithm>

//#include "ipp.h" // Intel API extension header file

using namespace std;

// Generate a random graph with a specified number of vertices and edges
void generateGraph(int n, int m, vector<vector<pair<int, double>>>& adj) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1, n);

    // Initialize adjacency list
    adj.resize(n);

    // Generate edges with random weights
    for (int i = 0; i < m; i++) {
        int u = dis(gen) - 1;
        int v = dis(gen) - 1;
        double w = dis(gen);
        adj[u].emplace_back(v, w);
    }
}

// Dijkstra's algorithm to find shortest paths from a source vertex to all other vertices
void dijkstra(int src, const vector<vector<pair<int, double>>>& adj, vector<double>& dist) {
    // Get the number of vertices in the graph
    int n = adj.size();

    // Initialize all distances to infinity
    dist.assign(n, numeric_limits<double>::infinity());

    // Set the distance of the source vertex to 0
    dist[src] = 0.0;

    // Create a priority queue to store vertices ordered by their distance from the source vertex
    // We use a min-heap here, so we store pairs of (distance, vertex) with distance as the first element
    // so that the priority queue orders them by distance
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    // Add the source vertex to the priority queue with distance 0
    pq.emplace(0.0, src);

    // Continue until all vertices have been processed
    while (!pq.empty()) {
        // Get the vertex with the smallest distance from the priority queue
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        // If we've already processed this vertex with a shorter distance, skip it
        if (d > dist[u]) continue;

        // Iterate over all neighbors of the current vertex
        for (const auto& [v, w] : adj[u]) {
            // Compute the distance from the source vertex to the neighbor through the current vertex
            double new_dist = dist[u] + w;

            // If this distance is shorter than the previous distance we computed for the neighbor, update it
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                // Add the neighbor to the priority queue with its updated distance
                pq.emplace(new_dist, v);
            }
        }
    }
}


int main() {
    int n = 1000; // number of vertices
    int m = 100000; // number of edges
    vector<vector<pair<int, double>>> adj;// adjacency list representing the weighted graph
    vector<double> dist;// output vector to store the distances from the source vertex

    generateGraph(n, m, adj);

    //int src = 0; // the source vertex

    for (int j = 0; j < n; j++){
        int src = j;
        dijkstra(src, adj, dist);

        // Print the result
        for (int i = 0; i < dist.size(); i++) {
            cout << "Shortest path from " << src << " to " << i << " is " << dist[i] << endl;
        }

    }

    

    // // Start profiling with the Intel API extension
    // ippInit();
    // Ipp64u t0 = ippGetCpuClocks();

    // // Run Dijkstra's algorithm from vertex 0
    // dijkstra(0, adj, dist);

    // // End profiling
    // Ipp64u t1 = ippGetCpuClocks();
    // cout << "Time elapsed: " << (t1 - t0) << " cycles\n";

    return 0;
}
