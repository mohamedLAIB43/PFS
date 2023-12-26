#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>

using namespace std; 

class Instance
{
public:
    int n; 
    int m; 
    string name; 
    vector<vector<int>> P; 
    long Cmax; 
    vector<int> permutation; 
};

class JobScheduler
{
public:
    
    static int lire_instance(const string &chemin_instance, Instance &instance)
    {
        ifstream fichier(chemin_instance);
        if (!fichier.is_open())
            return 1;

        fichier >> instance.n >> instance.m;

        instance.P.resize(instance.n, vector<int>(instance.m, 0));

        for (int i = 0; i < instance.n; ++i)
            for (int j = 0; j < instance.m; ++j)
                fichier >> instance.P[i][j];

        return 0;
    }


    static void ecrire_solution(const Instance &instance)
    {
        string nom_fichier_solution = "solution_" + instance.name;

        ofstream fichier(nom_fichier_solution);
        if (!fichier.is_open())
        {
            cerr << "Erreur lors de l'écriture du fichier de solution." << endl;
            return;
        }

        fichier << "Permutation : ";
        for (auto &p : instance.permutation)
            fichier << p << " ";
        fichier << "\nMakespan : " << instance.Cmax << endl;
    }

    
    static void differential_evolution_discrete(Instance &instance)
    {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0.0, 1.0);

        const int taille_population = 100;
        const int max_generations = 100;
        const double F = 0.8; 
        const double CR = 0.9; 

        vector<vector<int>> population(taille_population);
        for (auto &perm : population)
            perm = permutation_aleatoire(instance.n);

        for (int generation = 0; generation < max_generations; ++generation)
        {
            for (int i = 0; i < taille_population; ++i)
            {
                int r1, r2, r3;
                do
                {
                    r1 = i;
                    r2 = uniform_int_distribution<>(0, taille_population - 1)(gen);
                    r3 = uniform_int_distribution<>(0, taille_population - 1)(gen);
                } while (r1 == r2 || r1 == r3 || r2 == r3);

                vector<int> solution_essai(instance.n);
                for (int j = 0; j < instance.n; ++j)
                {
                    if (dis(gen) < CR || j == uniform_int_distribution<>(0, instance.n - 1)(gen))
                        solution_essai[j] = population[r1][j] + F * (population[r2][j] - population[r3][j]);
                    else
                        solution_essai[j] = population[i][j];

                    solution_essai[j] = max(0, min(instance.n - 1, solution_essai[j])); 
                }

                long makespan_essai = evaluer_makespan(instance, solution_essai);
                long makespan_actuel = evaluer_makespan(instance, population[i]);

                if (makespan_essai < makespan_actuel)
                {
                    population[i] = solution_essai;
                }
            }
        }
        auto meilleure_solution = min_element(population.begin(), population.end(), [&](const auto &a, const auto &b) {
            return evaluer_makespan(instance, a) < evaluer_makespan(instance, b);
        });

        instance.permutation = *meilleure_solution;
        instance.Cmax = evaluer_makespan(instance, *meilleure_solution);
    }

private:

    static vector<int> permutation_aleatoire(int taille)
    {
        vector<int> perm(taille);
        iota(perm.begin(), perm.end(), 0);
        random_shuffle(perm.begin(), perm.end());
        return perm;
    }

    static long evaluer_makespan(const Instance &instance, const vector<int> &permutation)
    {
        vector<vector<int>> C(instance.n, vector<int>(instance.m, 0));

        C[0][0] = instance.P[permutation[0]][0];
        for (int i = 1; i < instance.n; i++)
            C[i][0] = C[i - 1][0] + instance.P[permutation[i]][0];

        for (int j = 1; j < instance.m; j++)
            C[0][j] = C[0][j - 1] + instance.P[permutation[0]][j];

        for (int j = 1; j < instance.m; j++)
        {
            for (int i = 1; i < instance.n; i++)
            {
                C[i][j] = instance.P[permutation[i]][j] + max(C[i - 1][j], C[i][j - 1]);
            }
        }

        long makespan = max_element(C.begin(), C.end(), [](const auto &a, const auto &b) {
            return a > b;
        })[0][instance.m - 1];

        return makespan;
    }
};

int main()
{
    
    vector<string> chemins_instances = {
        "C:\\Users\\Mohamed laib\\Desktop\\flowshop\\taillard_N10_M2_1.txt",
        
    };

    
    for (const auto &chemin_instance : chemins_instances)
    {
        
        Instance instance;
        instance.name = chemin_instance;

       
        if (JobScheduler::lire_instance(chemin_instance, instance) != 0)
        {
            cerr << "Erreur lors de la lecture du fichier d'instance : " << chemin_instance << endl;
            continue; 
        }

        auto debut_execution = chrono::high_resolution_clock::now();

        
        JobScheduler::differential_evolution_discrete(instance);

        
        auto fin_execution = chrono::high_resolution_clock::now();

        
        cout << "Instance : " << chemin_instance << endl;
        cout << "Meilleure Permutation : ";
        for (auto &p : instance.permutation)
            cout << p << " ";
        cout << endl;

        cout << "Makespan : " << instance.Cmax << endl;

        
        JobScheduler::ecrire_solution(instance);

        
        auto duree = chrono::duration_cast<chrono::milliseconds>(fin_execution - debut_execution).count();
        cout << "Temps d'exécution : " << duree << " millisecondes" << endl;
        cout << "----------------------------------------" << endl;
    }

    return 0;
}
