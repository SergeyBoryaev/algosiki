// Copyright 2020 Osmanov Islam
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>


std::vector<int> global_best_tour;
double global_best_dist = DBL_MAX;

double GetRandomNumberFloat(double min, double max, int precision)
{
    // Установить стартовую точку
    srand(time(NULL));

    double value;

    // получить случайное число как целое число с порядком precision
    value = rand() % (int)pow(10, precision);

    // получить вещественное число
    value = min + (value / pow(10, precision)) * (max - min);

    return value;
}


class Edge{
public:
    int a, b;
    double pheromone, weight;
    Edge(int a, int b, double pheromone, double weight) :
            a(a), b(b), pheromone(pheromone), weight(weight) {}
};
class Ant{
public:
    double alpha, beta, distance = 0.0;
    int num_coords;
    std::vector<int> tour;
    std::vector<std::vector<Edge>> edges;
    Ant(double alpha, double beta, int coords,
        std::vector<std::vector<Edge>> edges) :
            alpha(alpha), beta(beta), distance(distance), edges(edges) {}
    int select_coord() {
        double wheel = 0.0;
        double heuristic = 0.0;
        std::vector<int> unused;
        for (int coord =0; coord < num_coords; coord++){
            if (std::find(tour.begin(),tour.end(), coord) == tour.end())
                unused.push_back(coord);
        }
        for (int un_element : unused){
            heuristic += edges[tour[tour.size() - 1]][un_element].weight;
        }
        for (int un_element : unused) {
            double f_i = pow(edges[tour[tour.size() - 1]][un_element].pheromone, alpha);
            double f_N = pow((heuristic / edges[tour[tour.size() - 1]][un_element].weight), beta);
            wheel += f_N * f_i;
        }
        double rand_val = GetRandomNumberFloat(0.0, wheel, 5);
        double wheel_pos = 0.0;
        for (int un_element : unused) {
            double f_i = pow(edges[tour[tour.size() - 1]][un_element].pheromone, alpha);
            double f_N = pow((heuristic / edges[tour[tour.size() - 1]][un_element].weight), beta);
            wheel_pos += f_N * f_i;
            if (wheel_pos >= rand_val) {
                return un_element;
            }
        }
    }

    void tour_find() {
        tour.push_back(rand() % num_coords);
        while (tour.size() < num_coords)
            tour.push_back(select_coord());
    }

    void calc_dist(){
        distance = 0.0;
        for (int i = 0; i < num_coords; i++){
            distance += edges[tour[i]][tour[(i + 1) % num_coords]].weight;
        }
    }
};


class ACO {

public:
    int colony_size, steps;
    double scaling, alpha, beta, ro, pheromone_dep_weight;
    std::vector<std::pair<double, double>> coords;
    std::vector<std::vector<Edge>> edges;
    std::vector<Ant> ants;
    ACO(    std::vector<std::pair<double, double>> coords,int colony_size, int steps, double scaling=0.001, double alpha=1.0, double beta=3.0, double ro=0.1,
        double pheromone_dep_weight=1.0) : colony_size(colony_size), steps(steps), scaling(scaling),
        alpha(alpha), beta(beta), ro(ro), pheromone_dep_weight(pheromone_dep_weight), coords(coords) {
        for (int i = 0; i < coords.size(); i++) {
            for(int j = i + 1; j < coords.size(); j++) {
                double weight = sqrt(pow(coords[i].first - coords[j].first, 2.0) +
                                     pow(coords[i].second - coords[j].second, 2.0) );
                edges[i][j] = edges[j][i] = Edge(i,j,weight,1.0);
            }
        }
        for (int i=0; i < colony_size; i++) {
            ants.push_back(Ant(alpha, beta, coords.size(), edges));
        }


    }

    void add_pher(std::vector<int> tour, double distance) {
        for (int i = 0; i < coords.size(); i++) {
            edges[tour[i]][tour[(i + 1) % coords.size()]].pheromone += (pheromone_dep_weight / distance);
        }
    }

    void acs() {
        for (int step = 0; step < steps; step++) {
            for(auto ant : ants) {
                ant.tour_find();
                ant.calc_dist();
                add_pher(ant.tour, ant.distance);
                if(ant.distance < global_best_dist) {
                    global_best_tour = ant.tour;
                    global_best_dist = ant.distance;
                }
            }
        }

        for (int i = 0; i < coords.size(); i++) {
            for(int j = i + 1; j < coords.size(); j++) {
                edges[i][j].pheromone *= (1.0 - ro);
            }
        }
        std::cout<<global_best_dist<<std::endl;
    }

};

int main() {

}