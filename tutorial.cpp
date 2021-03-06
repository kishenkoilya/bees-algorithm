#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "testfuncs/benchmarks.hpp"
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <math.h>

using namespace std;

using BM = Benchmark<double>;

const bool INTERSECTION_0 = true;
const int INTERSECTION_SCOUTS = 2; //0 = shift; 1 = disband worst; 2 = disband closest;
const bool INTERSECTION_SQUADS = true;
const bool DISPLAY_SQUADS = true;
const double RANGESHRINK = 0.9;
const double INITIALSHRINK = 0.05;
const double PRECISION = 0.0;

class bee
{
public:
    int dim;
    double v = numeric_limits<double>::max();
    vector<double> x;
    double shrink = 1.0;
    bee() {}
    bee(const BM &bm)
    {
        dim = bm.getDim();
        x.resize(dim);
        x.reserve(dim + 1);
        for (int i = 0; i < dim; i++)
        {
            double a = bm.getBounds()[i].first;
            double b = bm.getBounds()[i].second;
            x[i] = a + ((double)rand() / (RAND_MAX + 1.0)) * (b - a);
        }
        v = bm.calcFunc(x);
    }

    void flyto(const BM &bm, vector<double> leaderpos, double rangeshrink = 1.0)
    {
        dim = bm.getDim();
        for (int i = 0; i < dim; i++)
        {
            double a = bm.getBounds()[i].first;
            double b = bm.getBounds()[i].second;
            double minval = leaderpos[i] - ((b - a) * rangeshrink / 2);
            if (minval < a)
                minval = a;
            double maxval = leaderpos[i] + ((b - a) * rangeshrink / 2);
            if (maxval > b)
                maxval = b;
            x[i] = minval + ((double)rand() / (RAND_MAX + 1.0)) * (maxval - minval);
        }
        v = bm.calcFunc(x);
    }
    ~bee()
    {
        x.clear();
        x.shrink_to_fit();
    }
};

class membee
{
public:
    int dim;
    double v = numeric_limits<double>::max();
    vector<double> x;
    double rangeshrink;
    int funccounter;
    membee() {}
    membee(const BM &bm, bee b, int func)
    {
        dim = bm.getDim();
        x.resize(dim);
        x.reserve(dim + 1);
        x = b.x;
        v = b.v;
        funccounter = func;
        if (b.shrink == 1.0)
        {
            rangeshrink = INITIALSHRINK; //1.0 / sites;
        }
        else
        {
            rangeshrink = b.shrink;
        }
    }
    ~membee()
    {
        x.clear();
        x.shrink_to_fit();
    }
};

class hive
{
public:
    vector<bee> swarm;
    vector<membee> membees;
    int scoutbees;
    int selectedbees;
    int elitebees;
    int selectedsites;
    int elitesites;
    int beestotal;
    int maxfunccounter;
    int maxstag;
    hive(const BM &bm, int scouts, int selbees, int bestbees, int selsites, int bestsites, int maxfunccount, int maxstagnation)
    {
        scoutbees = scouts;
        selectedbees = selbees;
        elitebees = bestbees;
        selectedsites = selsites;
        elitesites = bestsites;
        maxfunccounter = maxfunccount;
        maxstag = maxstagnation;
        beestotal = selectedbees * selectedsites + elitebees * elitesites + scoutbees;
        swarm.resize(beestotal);
        swarm.reserve(beestotal + 1);
        membees.resize(selectedsites + elitesites);
        membees.reserve(selectedsites + elitesites + 1);
        srand(time(NULL));
        #pragma omp parallel for
        for (int i = 0; i < beestotal; i++) //hive creation
        {
            swarm[i] = bee(bm);
        }
        #pragma omp parallel for
        for (int i = 0; i < elitesites + selectedsites; i++)
        {
            membees[i] = membee(bm, swarm[i], 0);
        }
    }

    int findmin(int b, int e)
    {
        int min = b;
        for (int i = b + 1; i < e; i++) 
        {
            if (swarm[i].v < swarm[min].v) min = i;
        }
        return min;
    }

    void beeswap(int a, int b)
    {
        bee sw = swarm[a];
        swarm[a] = swarm[b];
        swarm[b] = sw;
    }

    bool intersects(const BM &bm, int m, int bee, bool mem = false)
    {
        int dim = bm.getDim();
        bool inter = true;
        for (int i = 0; i < dim; i++)
        {
            double a = bm.getBounds()[i].first;
            double b = bm.getBounds()[i].second;
            double minval = membees[m].x[i] - ((b - a) * membees[m].rangeshrink / 2);
            if (minval < a)
                minval = a;
            double maxval = membees[m].x[i] + ((b - a) * membees[m].rangeshrink / 2);
            if (maxval > b)
                maxval = b;
            if (mem)
            {
                if (membees[bee].x[i] > maxval)
                {
                    inter = false;
                    break;
                }
                if (membees[bee].x[i] < minval)
                {
                    inter = false;
                    break;
                }
            }
            else
            {
                if (swarm[bee].x[i] > maxval)
                {
                    inter = false;
                    break;
                }
                if (swarm[bee].x[i] < minval)
                {
                    inter = false;
                    break;
                }
            }
        }
        return inter;
    }

    void membeesshift(int b, int e, bool reverse = false)
    {
        if (reverse)
        {
            for (int i = b; i < e; i++)
            {
                membees[i] = membees[i + 1];
            }
        }
        else
        {
            for (int i = e; i > b; i--)
            {
                membees[i] = membees[i - 1];
            }
        }
    }

    void intsort(const BM &bm, int sortrange, bool init = false)
    {
        int intersected = 0;
        if (init)
        {
            int min = findmin(0, sortrange);
            beeswap(0, min);
            membees[0] = membee(bm, swarm[0], 0);
            int done = 1;
            while (done < elitesites + selectedsites) 
            {
                min = findmin(done, sortrange - intersected);
                bool inter = false;
                for (int i = 0; i < done; i++) 
                {
                    if (intersects(bm, i, min)) 
                    {
                        inter = true;
                        if (sortrange == intersected) 
                        {
                            while (done < elitesites + selectedsites)
                            {
                                min = findmin(done, sortrange);
                                beeswap(done, min);
                                membees[0] = membee(bm, swarm[done], 0);
                                done++;
                            }
                        }
                        else
                        {
                            beeswap(min, sortrange - intersected - 1);
                            intersected++;
                            break;
                        }
                    }
                }
                if (!inter) 
                {
                    beeswap(done, min);
                    membees[0] = membee(bm, swarm[done], 0);
                    done++;
                }
            }
        }
        else
        {
            for (int j = 0; j < sortrange;)
            {
                if (j >= sortrange - intersected) break;
                int min = findmin(j, sortrange - intersected);
                if (swarm[min].v < membees[elitesites + selectedsites - 1].v)
                {
                    int betterthan = elitesites + selectedsites - 1;
                    for (int i = 0; i < elitesites + selectedsites; i++)
                    {
                        if (swarm[min].v < membees[i].v)
                        {
                            betterthan = i;
                            break;
                        }
                    }
                    if (betterthan < elitesites + selectedsites - 1)
                    {
                        bool inter = false;
                        for (int i = 0; i < betterthan; i++)
                        {
                            if (intersects(bm, i, min))
                            {
                                inter = true;
                                beeswap(min, sortrange - intersected - 1);
                                intersected++;
                                break;
                            }
                        }
                        beeswap(j, min);
                        if (!inter)
                        {
                            bool interlist[elitesites + selectedsites - betterthan] = {false};
                            int intersections = 0;
                            for (int i = betterthan; i < elitesites + selectedsites; i++)
                            {
                                if (intersects(bm, i, j))
                                {
                                    interlist[i - betterthan] = true;
                                    intersections++;
                                }
                            }
                            if (intersections == 0 || INTERSECTION_SCOUTS == 0)
                            {
                                membeesshift(betterthan, elitesites + selectedsites - 1);
                                membees[betterthan].x = swarm[j].x;
                            }
                            else if (INTERSECTION_SCOUTS == 1)
                            {
                                int disband;
                                for (int i =  - betterthan - 1; i >= 0; i--)
                                {
                                    if (interlist[i])
                                    {
                                        disband = i;
                                        break;
                                    }
                                }
                                membeesshift(betterthan, disband + betterthan);
                                membees[betterthan].x = swarm[j].x;
                            }
                            else
                            {
                                int disband;
                                double min = numeric_limits<double>::max();
                                for (int i = 0; i < elitesites + selectedsites - betterthan; i++)
                                {
                                    if (interlist[i])
                                    {
                                        int dim = bm.getDim();
                                        double distance = 0;
                                        for (int k = 0; k < dim; k++)
                                        {
                                            distance += (membees[i + betterthan].x[k] - swarm[j].x[k]) * (membees[i + betterthan].x[k] - swarm[j].x[k]);
                                        }
                                        if (distance < min) 
                                        {
                                            min = distance;
                                            disband = i;
                                        }
                                    }
                                }
                                membeesshift(betterthan, disband + betterthan);
                                membees[betterthan].x = swarm[j].x;
                            }
                        }
                        j++;
                    }
                    else 
                    {
                        beeswap(j, min);
                        membees[betterthan].x = swarm[j].x;
                        j++;
                    }
                }
                else break;
            }
        }
    }

    void sortswarm(const BM &bm, int sortrange)
    {
        for (int i = 0; i < elitesites + selectedsites; i++)
        {
            membee bestbee = membees[i];
            bee bestswarmbee = swarm[i];
            int bestnum = i;
            bool ind = false;
            for (int j = i + 1; j < sortrange; j++)
            {
                if (swarm[j].v < bestbee.v && swarm[j].v < bestswarmbee.v)
                {
                    bestswarmbee = swarm[j];
                    bestnum = j;
                    ind = true;
                }
            }
            if (ind)
            {
                swarm[bestnum] = swarm[i];
                swarm[i] = bestswarmbee;
                for (int k = 0; k < elitesites + selectedsites; k++)
                {
                    if (bestswarmbee.v < membees[k].v)
                    {
                        membee reserve = membees[k];
                        membees[k] = membee(bm, bestswarmbee, 0);
                        for (int j = elitesites + selectedsites - 1; j > k + 1; j--)
                        {
                            membees[j] = membees[j - 1];
                        }
                        if (k < elitesites + selectedsites - 1)
                        {
                            membees[k + 1] = reserve;
                        }
                        reserve.~membee();
                        break;
                    }
                }
            }
            bestswarmbee.~bee();
            bestbee.~membee();
        }
    }

    void sortmembees(const BM &bm)
    {
        for (int i = 0; i < elitesites + selectedsites; i++)
        {
            membee best = membees[i];
            int bestnum = i;
            bool ind = false;
            for (int j = i + 1; j < elitesites + selectedsites; j++)
            {
                if (membees[j].v < best.v)
                {
                    best = membees[j];
                    bestnum = j;
                    ind = true;
                }
            }
            if (ind)
            {
                membees[bestnum] = membees[i];
                membees[i] = best;
                best.~membee();
            }
        }
        if (INTERSECTION_SQUADS) 
        {
            int b = 0;
            for (int i = elitesites + selectedsites - 1; i > 0; i--)
            {
                bool inter = false;
                for (int j = 0; j < i; j++)
                {
                    inter = intersects(bm, j, i, true);
                    if (inter) break;
                }
                if (inter)
                {
                    if (i != elitesites + selectedsites - 1) membeesshift(i, elitesites + selectedsites - 1, true);
                    int min = findmin(b, scoutbees);
                    membees[elitesites + selectedsites - 1] = membee(bm, swarm[min], 0);
                    beeswap(b, min);
                    b++;
                }     
            }
        }
    }

    void beesatwork(const BM &bm)
    {
        for (int i = 0; i < elitesites + selectedsites; i++) //squad search cycle
        {
            membee bestsbee = membees[i];
            bee bestsitebee = swarm[i];
            if (i < elitesites) //best spots squads
            {
                bool ind = false;
                #pragma omp parallel for
                for (int j = scoutbees + i * elitebees; j < scoutbees + (i + 1) * elitebees; j++)
                {
                    swarm[j].flyto(bm, membees[i].x, membees[i].rangeshrink);
                    if (swarm[j].v < bestsbee.v && swarm[j].v < bestsitebee.v)
                    {
                        bestsitebee = swarm[j];
                        ind = true;
                    }
                }
                if (ind)
                {
                    membees[i] = membee(bm, bestsitebee, membees[i].funccounter + 1);
                }
                else
                {
                    membees[i].funccounter++;
                    if (membees[i].funccounter >= maxfunccounter)
                    {
                        membees[i].funccounter = 0;
                        membees[i].rangeshrink *= RANGESHRINK;
                    }
                }
            }
            else //selected sites squads
            {
                bool ind = false;
                #pragma omp parallel for
                for (int j = scoutbees + elitesites * elitebees + (i - elitesites) * selectedbees; j < scoutbees + elitesites * elitebees + (i - elitesites + 1) * selectedbees; j++)
                {
                    swarm[j].flyto(bm, membees[i].x, membees[i].rangeshrink);
                    if (swarm[j].v < bestsbee.v && swarm[j].v < bestsitebee.v)
                    {
                        bestsitebee = swarm[j];
                        ind = true;
                    }
                }
                if (ind)
                {
                    membees[i] = membee(bm, bestsitebee, membees[i].funccounter + 1);
                }
                else
                {
                    membees[i].funccounter++;
                    if (membees[i].funccounter >= maxfunccounter)
                    {
                        membees[i].funccounter = 0;
                        membees[i].rangeshrink *= RANGESHRINK;
                    }
                }
            }
        }

        int dim = bm.getDim();
        vector<double> center(dim); //variables for scouts //can be sent to hive variables

        for (int i = 0; i < dim; i++)
        {
            double a = bm.getBounds()[i].first;
            double b = bm.getBounds()[i].second;
            center[i] = (a + b) / 2.0;
        }
        #pragma omp parallel for //num_threads(8)
        for (int i = 0; i < scoutbees; i++)
        {
            swarm[i].flyto(bm, center);
        }
    }

    double runcycle(const BM &bm, int maxiter, ofstream &fout, ofstream &fout2, ofstream &fout3)
    {
        double bestfit = numeric_limits<double>::max();
        int stagnationcount = 0;
        int z;
        for (z = 0; z < maxiter; z++)
        {
            // auto begin = chrono::high_resolution_clock::now();
            beesatwork(bm);
            // auto end = chrono::high_resolution_clock::now();
            // fout2 << chrono::duration_cast<chrono::nanoseconds>(end - begin).count() << ";";
            sortmembees(bm);
            if (INTERSECTION_0) intsort(bm, scoutbees);
            else 
            sortswarm(bm, scoutbees);

            if (membees[0].v < bestfit)
            {
                bestfit = membees[0].v;
                stagnationcount = 0;
            }
            else
            {
                stagnationcount++;
            }
            if (stagnationcount >= maxstag)
            {
                cout << "algorithm reached stagnation on iteration N" << z << endl;
                fout << z << ";";
                break;
            }
            if (z == maxiter - 1)
            {
                fout << maxiter << ";";
            }
            if (bestfit - bm.getGlobMinY() < PRECISION && PRECISION != 0.0)
            {
                cout << "algorithm reached appropriate level of precision on iteration N" << z << endl;
                fout << z << ";";
                break;
            }
        }
        if (DISPLAY_SQUADS)
        {
            fout3 << bm.getDesc() << endl;
            int dim = bm.getDim();
            double dist;
            fout3 << " " << bm.getGlobMinY() << ";" << "0;";
            for (int i = 0; i < dim; i++)
            {
                fout3 << " " <<  bm.getGlobMinX()[i] << ";";
            }
            fout3 << endl;
            for (int j = 0; j < elitesites + selectedsites; j++)
            {
                dist = 0;
                fout3 << " " << membees[j].v - bm.getGlobMinY() << ";";
                for (int i = 0; i < dim; i++)
                {
                    dist += (membees[j].x[i] - bm.getGlobMinX()[i]) * (membees[j].x[i] - bm.getGlobMinX()[i]);
                }
                fout3 << " " << sqrt(dist) << ";";
                for (int i = 0; i < dim; i++)
                {
                    fout3 << " " << membees[j].x[i] << ";";
                }
                fout3 << endl;
            }
            fout3 << endl;
        }
        return bestfit;
    }
    ~hive()
    {
        for (int i = 0; i < beestotal; i++)
        {
            swarm[i].~bee();
        }
        for (int i = 0; i < elitesites + selectedsites; i++)
        {
            membees[i].~membee();
        }
        membees.clear();
        membees.shrink_to_fit();
        swarm.clear();
        swarm.shrink_to_fit();
    }
};

double findMin(const BM &bm, ofstream &fout, int param[8], ofstream &fout2, ofstream &fout3)
{
    int scoutbees = param[0]; //my parameters
    int selectedbees = param[1];
    int elitebees = param[2];
    int selectedsites = param[3];
    int elitesites = param[4];
    int maxfunccounter = param[5];
    int maxiter = param[6];
    int maxstag = param[7];
    const int beestotal = selectedbees * selectedsites + elitebees * elitesites + scoutbees;
    // auto begin = chrono::high_resolution_clock::now();
    hive *h = new hive(bm, scoutbees, selectedbees, elitebees, selectedsites, elitesites, maxfunccounter, maxstag);
    // auto end = chrono::high_resolution_clock::now();
    // fout2 << chrono::duration_cast<chrono::nanoseconds>(end - begin).count() << ";";
    if (INTERSECTION_0) h->intsort(bm, beestotal, true);
    else h->sortswarm(bm, beestotal);
    h->sortmembees(bm);

    double result = h->runcycle(bm, maxiter, fout, fout2, fout3);
    cout << "best value = " << result << endl;
    int dim = bm.getDim();
    double distance = 0;
    for (int i = 0; i < dim; i++)
    {
        cout << "coordinate " << i << ": " << h->membees[0].x[i] << endl;
        distance += (h->membees[0].x[i] - bm.getGlobMinX()[i]) * (h->membees[0].x[i] - bm.getGlobMinX()[i]);
    }
    distance = sqrt(distance);
    cout << "distance: " << distance << endl;
    fout << distance << ";";
    delete h;
    return result;
}

bool testBench(const BM &bm, ofstream &fout, int param[8], ofstream &fout2, ofstream &fout3)
{
    cout << "*************Testing benchmark**********" << endl;
    cout << bm;
    fout << bm.getDesc() << ";";
    auto begin = chrono::high_resolution_clock::now();
    double v = findMin(bm, fout, param, fout2, fout3);
    cout << "the difference is " << v - bm.getGlobMinY() << endl;
    fout << v - bm.getGlobMinY() << ";";
    auto end = chrono::high_resolution_clock::now();
    cout << "time since the beginning: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms" << endl;
    fout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << ";";
    cout << "****************************************" << endl
         << endl;
}

main()
{
    omp_set_num_threads(6);
    ofstream fout("data.csv", ios::app);
    ofstream fout2("execution_time.csv", ios::app);
    int h = 19;
    for (int i = 0; i < 5; i++)
    {
        string n1 = "squadresults";
        string n2 = to_string(i);
        string n3 = ".csv";
        string name = n1 + n2 + n3;
        ofstream fout3(name, ios::app);
        int param[8] = {100, 50, 100, 5, 3, 3, 5000, 100};
        if (i > 20)
        {
            param[0] = 200;
        }
        if (i > 40)
        {
            param[0] = 50;
        }
        if (i > 60)
        {
            param[0] = 100;
            param[1] = 15;
            param[2] = 30;
            param[3] = 5;
            param[4] = 3;
            param[5] = 3;
        }
        if (i > 80)
        {
            param[0] = 200;
            param[1] = 15;
            param[2] = 30;
            param[3] = 5;
            param[4] = 3;
            param[5] = 3;
        }
        if (i > 100)
        {
            param[0] = 50;
            param[1] = 15;
            param[2] = 30;
            param[3] = 5;
            param[4] = 3;
            param[5] = 3;
        }
        h++;
        if (h == 20)
        {
            h = 0;
            fout << "scouts;" << "selectedbees;" << "elitebees;" << "selectedsquads;" << "elitesquads;";
            fout << "funccounter;" << "maxiter;" << "maxstag;" << "INTERSECTION_0;" << "INTERSECTION_SCOUTS;";
            fout << "INTERSECTION_SQUADS;" << "RANGESHRINK;" << "INITIALSHRINK;" << endl;
            fout << param[0] << ";" << param[1] << ";" << param[2] << ";" << param[3] << ";" << param[4] << ";";
            fout << param[5] << ";" << param[6] << ";" << param[7] << ";" << INTERSECTION_0 << ";" << INTERSECTION_SCOUTS << ";";
            fout << INTERSECTION_SQUADS << ";" << RANGESHRINK << ";" << INITIALSHRINK << ";" << endl;
        }

        auto begin = chrono::high_resolution_clock::now();
        Benchmarks<double> tests;
        int t = 0;
        // Alpine1Benchmark<double> ab(3);
        // BealeBenchmark<double> bb;
        // testBench(bb, fout, param, fout2);
        for (auto bm : tests)
        {
            testBench(*bm, fout, param, fout2, fout3);
            // fout2 << endl;
            // if (t < 10)
            // {
            //     testBench(*bm, fout, param, fout2);
            //     fout2 << endl;
            // }
            // else
            // {
            //     break;
            // }
            // t++;
        }
        auto end = chrono::high_resolution_clock::now();
        cout << "TIME TOTAL: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms" << endl;
        fout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << ";";
        fout << endl;
        // fout2 << endl;
    }
    fout << endl;
}
