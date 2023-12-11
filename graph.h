// Name: Zaheer Safi
// Date: 11/26/2023
// CS_251 Project_5 : Open street maps

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <list>
#include <unordered_map>


using namespace std;

template<typename VertexT, typename WeightT>
class graph 
{
private:
  
  // adjancy list that map the vertex to all the neighbor veritices and theier weights
  unordered_map<VertexT, list<pair<VertexT, WeightT>>> AdjList;
  
  // vector to store all the vertices
  vector<VertexT> Vertices; 

  // private function to check if a vertex exist in the graph using the adjasncy list
  int _LookupVertex(VertexT v) const 
  {
    // check if the vertex existes in the map if yes return 1
    auto it = AdjList.find(v);
    if (it != AdjList.end()) 
    {
      return 1;
    }
    // else return -1
    return -1;
  }
 
public:
  
  // default constructor 
  graph(){}

  // function to return the number of vertices in our graph
  int NumVertices() const 
  {
    return static_cast<int>(this->AdjList.size());
  }

  // function to return the number of edges in the graph
  int NumEdges() const 
  {
    int count = 0;

    // go thorugh the adjancy list map and count the number of edges and return it
    for (auto i : AdjList)
    {
      count += i.second.size();
    }
    return count;
  }

  // function to add a vertex to the graph give vertex v
  bool addVertex(VertexT v) 
  {
    // check if the vertex already exists if yes then return false
    if (_LookupVertex(v) == 1) 
    {
      return false;
    }

    // otherwise add the vertex to the vector of vertices and also to the map with empty list of adjasant vertices
    list<pair<VertexT, WeightT>> new_list;
    this->Vertices.push_back(v);
    this->AdjList.emplace(v, new_list);

    return true;
  }

  // function ot add a edge from give from vertex to vertex "to" with given weight "weight"
  bool addEdge(VertexT from, VertexT to, WeightT weight) 
  {
    // Check if the vertices exist in the graph
    if (AdjList.find(from) == AdjList.end() || AdjList.find(to) == AdjList.end()) 
    {
      return false;
    }

    // Check if the edge already exists
    for (auto& neighbor : AdjList[from]) 
    {
      if (neighbor.first == to) 
      {
        // Edge already exists, overwrite the weight
        neighbor.second = weight;
        return true;
      }
    }

    // Edge does not exist, add it
    AdjList[from].emplace_back(to, weight);
    return true;
  }

  // getter function to return the weight of edge form givem from "from" to "to" vertex and update the weight parameter
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const 
  {
    // check if the edge does not exist if yes return false
    if (AdjList.find(from) == AdjList.end() || AdjList.find(to) == AdjList.end()) 
    {
      return false;
    }
    
    // search for the from vertex
    auto it = AdjList.find(from);

    // now go through the list of given vertex
    for (auto& neighbor : it->second) 
    {
      // if the given vertex is found then update the weight parameter and return true
      if (neighbor.first == to) 
      {
        weight = neighbor.second;
        return true;
      }
    }

    return false;
  }

  // funciton to return the neighbors of a given vertex by looking over the adjList
  set<VertexT> neighbors(VertexT v) const 
  {
    set<VertexT>  S; // veriable to return 

    auto it = AdjList.find(v); // find the vertex in the map

    if (it == AdjList.end()) 
    { 
      return S;
    }

    // now go through the list for the given vertex and add the elements ot the set and return it
    for (const auto& neighbor : it->second) 
    {
      S.insert(neighbor.first);
    }

    return S;
  }

  // function to return the vertices by return the vector
  vector<VertexT> getVertices() const 
  {
    return this->Vertices; 
  }

  // dump function to dispaly all the vertices of the graph and all neighbors of the vertices
  void dump(ostream& output) const 
  {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    // output the number of vertices and the number of edges
    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    // go throught the vertices vector to output all the vertices
    for (int i = 0; i < this->NumVertices(); ++i) 
    {
      output << " " << i << ". " << this->Vertices[i] << endl;
    }

    output << endl;
    output << "**Edges:" << endl;
    // go through the adjList to show all the edges
    for (auto& i : AdjList)
    {
      output << i.first << ": ";
      for (auto& neighbor : i.second) 
      {
        output << "(" << neighbor.first << ", " << neighbor.second << "), " ;
      }
      output << endl;
    }
    
    output << "**************************************************" << endl;
  }
};
