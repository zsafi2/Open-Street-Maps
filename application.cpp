// Name: Zaheer Safi
// Date: 11/26/2023
// CS_251 Project_5 : Open street maps

#include <iostream>
#include <iomanip> 
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <set>
#include <stack>

#include "tinyxml2.h"
#include "dist.h"
#include "graph.h"
#include "osm.h"

using namespace std;
using namespace tinyxml2;

const double INF = numeric_limits<double>::max(); // define infinity

// struct for the priority queue to set up the priority for the dijiksitra's algorithm
struct prioritize 
{
  // to set up the "()"" operator
  bool operator()(const pair<long long, double>& lhs, const pair<long long, double>& rhs) const 
  {
    return lhs.second > rhs.second; // the lowest value has the highest priority
  }
};

// function to search a building int the buildings vector where every element is a buildingInfo
BuildingInfo searchBuilding(string query, vector<BuildingInfo>& Buildings)
{
  // veriable to return the found building
  BuildingInfo found;
  size_t search;
  
  // go through the Buildings vector
  for (int i = 0; i < Buildings.size(); i++)
  {
    // if the Abbreviations matches the query return the building
    if (Buildings[i].Abbrev == query)
    {
      found = Buildings[i];
      return found;
    }
    
    // if the query is not a abbreviation then search the query in the full name
    search = Buildings[i].Fullname.find(query);
    
    // if the query exists in the Given building then return the builing
    if (search != std::string::npos)
    {
      found = Buildings[i];
      return found;
    } 
  }
}

// given the Building which is a buildingInfo class search through Footways and find the nearest node and return the coordinates
Coordinates nearest_Node(vector<FootwayInfo>& Footways, map<long long, Coordinates>& Nodes, BuildingInfo Building)
{
  Coordinates nearest; // veriable the return the nearest coordinate
  Coordinates point1; // temprary veriable to find the minimum coordinate
  
  double min_dist = INF;
  double reg_dist;
  
  // go through all the footway nodes
  for (int i = 0; i < Footways.size(); i++)
  {
    // go through the member veriable Nodes of the footway
    for (int j = 0; j < Footways[i].Nodes.size(); j++)
    {
      // veriable to hold the Node value
      point1 = Nodes[Footways[i].Nodes[j]];
      // find the distance between the node and the given building coordinates
      reg_dist = distBetween2Points(point1.Lat, point1.Lon, Building.Coords.Lat, Building.Coords.Lon);
      // find the minimum value
      if (reg_dist < min_dist)
      {
        min_dist = reg_dist;
        nearest = point1;
      }

    }
  }
  return nearest;
}

// This function implements an algorithm to find the nearest path to all the nodes of the given graph from the given ID vertex.
// The algorithm uses Dijkstra's shortest path algorithm with some modifications.
void algorithm(graph<long long, double>& G, long long ID, map<long long, double>& distances, map<long long, long long>& predecessors)
{
  // Priority queue to store vertices and their tentative distances from the source vertex.
  priority_queue<pair<long long, double>, vector<pair<long long, double>>, prioritize> unvisitedQueue;
  
  // Get the list of all vertices in the graph.
  vector<long long> nodes = G.getVertices();
  
  // Initialize distances and predecessors for all vertices in the graph.
  for (int i = 0; i < G.NumVertices(); i++)
  {
    unvisitedQueue.push({nodes[i], INF});
    distances.emplace(nodes[i], INF);
    predecessors.emplace(nodes[i], 0);
  }
  
  // Set the distance of the source vertex to 0 and push it to the priority queue.
  unvisitedQueue.push({ID, 0.0});
  distances[ID] = 0.0;
  
  // Variables to store current vertex, visited vertices, adjacent vertices, edge weight, and alternate path distance.
  pair<long long, double> currentV;
  set<long long> visited;
  set<long long> adjVertices;
  double edgeweight;
  double alternatePathDistance;
  
  // Main loop to explore the graph and update distances.
  while (!unvisitedQueue.empty())
  {
    // Get the vertex with the smallest tentative distance.
    currentV = unvisitedQueue.top();
    
    // Check if the vertex has already been visited.
    auto it = visited.find(currentV.first);
    
    // Remove the vertex from the priority queue.
    unvisitedQueue.pop();
    
    // If the tentative distance is INF, break the loop.
    if (currentV.second == INF)
    {
      break;
    }
    // If the vertex has already been visited, skip it.
    else if(it != visited.end())
    {
      continue;
    }
    // Otherwise, process the vertex.
    else
    {
      // Get the set of adjacent vertices.
      adjVertices = G.neighbors(currentV.first);
      
      // Iterate over adjacent vertices and update distances.
      for (auto adjV : adjVertices)
      {
        // Get the edge weight between currentV and adjV.
        G.getWeight(currentV.first, adjV, edgeweight);
        
        // Calculate the alternate path distance.
        alternatePathDistance = edgeweight + distances[currentV.first];
        
        // If the alternate path is shorter, update distances and predecessors.
        if (alternatePathDistance < distances[adjV])
        {
          distances[adjV] = alternatePathDistance;
          predecessors[adjV] = currentV.first;
          
          // Push the adjacent vertex and its updated distance to the priority queue.
          unvisitedQueue.push({adjV, alternatePathDistance});
        }
      }
    }
  }

  return;
}

// This function prints the path from the source vertex to the given destination vertex using the predecessors map.
void printpath(map<long long, long long>& predecessors, long long destination)
{
  // Variable to track the current vertex in the path.
  long long currentV = destination;
  
  // Vector to store the resulting path.
  vector<long long> resulting_path;
  
  // Stack to build the path in reverse order.
  stack<long long> path;
  
  // Build the path by backtracking through predecessors until the source vertex (0) is reached.
  while (currentV != 0)
  {
    path.push(currentV);
    currentV = predecessors[currentV];
  }
  
  // Reverse the order of vertices in the path and store it in the resulting_path vector.
  while (!path.empty())
  {
    currentV = path.top();
    path.pop();
    resulting_path.push_back(currentV);
  }
  
  // Print the resulting path.
  cout << "Path: ";
  for (int i = 0; i < resulting_path.size(); i++)
  {
    cout << resulting_path[i] << "->";
  }
  cout << endl;
}

// This function represents the main application logic that calculates and displays information about the paths between two persons in a given map.
void application(map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
                  vector<BuildingInfo>& Buildings, graph<long long, double>& G) 
{
    // Variables to store input building names, person information, midpoint, and central building.
    string person1Building, person2Building;
    
    // Prompt the user to enter the name or abbreviation of person 1's building.
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    // Main loop to continue processing until the user enters "#".
    while (person1Building != "#") 
    {
      // Prompt the user to enter the name or abbreviation of person 2's building.
      cout << "Enter person 2's building (partial name or abbreviation)> ";
      getline(cin, person2Building);

      // Search for building information based on user input.
      BuildingInfo person1 = searchBuilding(person1Building, Buildings);
      BuildingInfo person2 = searchBuilding(person2Building, Buildings);
      Coordinates midpoint;
      BuildingInfo central_Building;

      // Check if building information for person 1 is found.
      if (person1.Fullname == "")
      {
        cout << "Person 1's building not found" << endl;
      }

      // Check if building information for person 2 is found.
      if (person2.Fullname == "")
      {
        cout << "Person 2's building not found" << endl;
      }

      // If both person 1 and person 2 building information is found.
      if (person1.Fullname != "" && person2.Fullname != "")
      {
        double minimum_dist = INF;
        double regular_dist;

        // Calculate the midpoint between person 1 and person 2.
        midpoint = centerBetween2Points(person1.Coords.Lat, person1.Coords.Lon, person2.Coords.Lat, person2.Coords.Lon);

        // Find the central building nearest to the midpoint.
        for (int i = 0; i < Buildings.size(); i++)
        {
            regular_dist = distBetween2Points(Buildings[i].Coords.Lat, Buildings[i].Coords.Lon, midpoint.Lat, midpoint.Lon);
            if (regular_dist < minimum_dist)
            {
                minimum_dist = regular_dist;
                central_Building = Buildings[i];
            }
        }

        // Find the nearest nodes on the footways for person 1, person 2, and the central building.
        Coordinates nearest_Node_build1 = nearest_Node(Footways, Nodes, person1);
        Coordinates nearest_Node_build2 = nearest_Node(Footways, Nodes, person2);
        Coordinates nearest_Node_build3 = nearest_Node(Footways, Nodes, central_Building);

        // Display information about person 1's point.
        cout << "Person 1's point:" << endl;
        cout << " " << person1.Fullname << endl;
        cout << " (" << person1.Coords.Lat << ", " << person1.Coords.Lon << ")" << endl;

        // Display information about person 2's point.
        cout << "Person 2's point:" << endl;
        cout << " " << person2.Fullname << endl;
        cout << " (" << person2.Coords.Lat << ", " << person2.Coords.Lon << ")" << endl;

        // Display information about the destination building.
        cout << "Destination Building:" << endl;
        cout << " " << central_Building.Fullname << endl;
        cout << " (" << central_Building.Coords.Lat << ", " << central_Building.Coords.Lon << ")" << endl;

        // Display information about the nearest node for person 1.
        cout << "Nearest P1 node:" << endl;
        cout << " " << nearest_Node_build1.ID << endl;
        cout << " (" << nearest_Node_build1.Lat << ", " << nearest_Node_build1.Lon << ")" << endl;

        // Display information about the nearest node for person 2.
        cout << "Nearest P2 node:" << endl;
        cout << " " << nearest_Node_build2.ID << endl;
        cout << " (" << nearest_Node_build2.Lat << ", " << nearest_Node_build2.Lon << ")" << endl;

        // Display information about the nearest node for the destination building.
        cout << "Nearest destination node:" << endl;
        cout << " " << nearest_Node_build3.ID << endl;
        cout << " (" << nearest_Node_build3.Lat << ", " << nearest_Node_build3.Lon << ")" << endl;

        // Create maps to store distances and predecessors for the Dijkstra's algorithm.
        map<long long, double> distances;
        map<long long, long long> predecessors;
        map<long long, double> distances1;
        map<long long, long long> predecessors1;

        // Run Dijkstra's algorithm for person 1's nearest node.
        algorithm(G, nearest_Node_build1.ID, distances, predecessors);

        // Run Dijkstra's algorithm for person 2's nearest node.
        algorithm(G, nearest_Node_build2.ID, distances1, predecessors1);

        // Check if the destination is unreachable for either person.
        if (distances[nearest_Node_build3.ID] >= INF || distances1[nearest_Node_build3.ID] >= INF)
        {
          cout << "Sorry, destination unreachable." << endl;
        }
        else 
        {
          // Display the distance and path for person 1 to the destination.
          cout << "Person 1's distance to dest: " << distances[nearest_Node_build3.ID] << " miles" << endl;
          printpath(predecessors, nearest_Node_build3.ID);

          // Display the distance and path for person 2 to the destination.
          cout << "Person 2's distance to dest: " << distances1[nearest_Node_build3.ID] << " miles" << endl;
          printpath(predecessors1, nearest_Node_build3.ID);
        }
      }

      // Prompt the user to enter person 1's building again.
      cout << endl;
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
    }    
}

// This is the main function for the UIC open street map navigation program.
int main() 
{
  // Create a graph to represent the map.
  graph<long long, double> G;

  // Maps a Node ID to its coordinates (lat, lon).
  map<long long, Coordinates> Nodes;
  
  // Information about each footway, in no particular order.
  vector<FootwayInfo> Footways;
  
  // Information about each building, in no particular order.
  vector<BuildingInfo> Buildings;
  
  // XML document to store the map data.
  XMLDocument xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  // Set the default filename for the map file.
  string def_filename = "map.osm";
  string filename;

  // Prompt the user to enter the map filename.
  cout << "Enter map filename> ";
  getline(cin, filename);

  // Use the default filename if the user does not provide one.
  if (filename == "") {
    filename = def_filename;
  }

  // Load the XML-based map file.
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  // Read the nodes, which are the various known positions on the map.
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  // Read the footways, which are the walking paths.
  int footwayCount = ReadFootways(xmldoc, Footways);

  // Read the university buildings.
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  // Ensure consistency in counts.
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  // Display statistics about the loaded map data.
  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  // Add vertices to the graph based on the map nodes.
  for (auto &i : Nodes)
  {
    G.addVertex(i.first);
  }

  double distance = 0;
  Coordinates point1;
  Coordinates point2;

  // Iterate through footways and add edges to the graph with distances based on coordinates.
  for (int i = 0; i < Footways.size(); i++)
  {
    for (int j = 0; j < Footways[i].Nodes.size() - 1; j++)
    {
      point1 = Nodes[Footways[i].Nodes[j]];
      point2 = Nodes[Footways[i].Nodes[j+1]];
      distance = distBetween2Points(point1.Lat, point1.Lon, point2.Lat, point2.Lon);
      G.addEdge(Footways[i].Nodes[j], Footways[i].Nodes[j+1], distance);
      G.addEdge(Footways[i].Nodes[j+1], Footways[i].Nodes[j], distance);
    }
  }

  // Display statistics about the built graph.
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  // Run the application to navigate the map.
  application(Nodes, Footways, Buildings, G);

  cout << "** Done **" << endl;
  return 0;
}

