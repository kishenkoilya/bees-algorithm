/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tutorial.cpp
 *
 * Created on November 20, 2017, 9:13 AM
 */
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "testfuncs/benchmarks.hpp"
#include <stdlib.h>
#include <time.h>

using BM = Benchmark<double>;

class doublebee
{
    public:
    int dim;
    double v = std::numeric_limits<double>::max();
    std::vector<double> x;
    double shrink = 1.0;
        doublebee(){}
        doublebee(const BM& bm)
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

        void flyto(const BM& bm, std::vector<double> leaderpos, double rangeshrink = 1.0)
        {
            dim = bm.getDim();
            shrink = rangeshrink;
            for (int i = 0; i < dim; i++)
            {
                double a = bm.getBounds()[i].first;
                double b = bm.getBounds()[i].second;
                double minval = leaderpos[i] - ((b-a)*shrink);
                if (minval < a) minval = a;
                double maxval = leaderpos[i] + ((b-a)*shrink);
                if (maxval > b) maxval = b;
                x[i] = minval + ((double)rand() / (RAND_MAX + 1.0)) * (maxval - minval);

            }
            v = bm.calcFunc(x);
        }
        ~doublebee(){
    x.clear();
    x.shrink_to_fit();
}
};

class membee
{
public:
    int dim;
    double v = std::numeric_limits<double>::max();
    std::vector<double> x;
    double rangeshrink;
    int funccounter;
    membee(){}
    membee(const BM& bm, doublebee bee, int func, int sites)
    {
        dim = bm.getDim();
        x.resize(dim); 
        x.reserve(dim + 1);
        x = bee.x;
        v = bee.v;
        funccounter = func;
        if(bee.shrink == 1.0)
        {
            rangeshrink = 1.0 / sites * 2.0;
        }
        else
        {
            rangeshrink = bee.shrink / 2.0;
        }
    }
    ~membee(){
    x.clear();
    x.shrink_to_fit();
}
};

class hive
{
    public:
        std::vector<doublebee> swarm;
        std::vector<membee> membees;
        int scoutbeecount;
        int selectedbeecount;
        int bestbeecount;
        int selsitescount;
        int bestsitescount;
        int beestotal;
        int maxfunccounter;
        int maxstag;
        double rangeshrinker;
        hive(const BM& bm, int scouts, int selbees, int bestbees, int selsites, int bestsites, double shrinker, int maxfunccount, int maxstagnation)
        {
            scoutbeecount = scouts;
            selectedbeecount = selbees;
            bestbeecount = bestbees;
            selsitescount = selsites;
            bestsitescount = bestsites;
            rangeshrinker = shrinker;
            maxfunccounter = maxfunccount;
            maxstag = maxstagnation;
            beestotal = selectedbeecount * selsitescount + bestbeecount * bestsitescount + scoutbeecount;
            swarm.resize(beestotal);
            swarm.reserve(beestotal + 1);
            membees.resize(selsitescount + bestsitescount);
            membees.reserve(selsitescount + bestsitescount + 1);
            srand(time(NULL));
            for (int i = 0; i < beestotal; i++)//hive creation
            {
                swarm[i] = doublebee(bm);
            }
            for (int i = 0; i < bestsitescount + selsitescount; i++)
            {
                membees[i] = membee(bm, swarm[i], 0, bestsitescount + selsitescount);
            }
        }

        void sortswarm(const BM& bm, int sortrange)
        {
            for (int i = 0; i < bestsitescount + selsitescount; i++)
            {
                membee bestbee = membees[i];
                doublebee bestswarmbee = swarm[i];
                int bestnum = i;
                bool ind = false;
                for (int j = i+1; j < sortrange; j++)
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
                    for (int k = 0; k < bestsitescount + selsitescount; k++)
                    {
                        if (bestswarmbee.v < membees[k].v) 
                        {
                            membee reserve = membees[k];
                            membees[k] = membee(bm, bestswarmbee, 0, bestsitescount + selsitescount);
                            for(int j = bestsitescount + selsitescount - 1; j > k + 1; j--)
                            {
                                membees[j] = membees[j - 1];
                            }
                            if (k < bestsitescount + selsitescount - 1) {
                                membees[k + 1] = reserve;
                            }
                            reserve.~membee();
                            break;
                        }
                    }
                }
                bestswarmbee.~doublebee();
                bestbee.~membee();
            }
        }

        void sortmembees()
        {
            for(int i = 0; i < bestsitescount + selsitescount; i++)
            {
                membee best = membees[i];
                int bestnum = i;
                bool ind = false;
                for(int j = i+1; j < bestsitescount + selsitescount; j++)
                {
                    if(membees[j].v < best.v)
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
        }

        void beesatwork(const BM& bm)
        {
            for(int i = 0; i < bestsitescount + selsitescount; i++)//squad search cycle
            {
                membee bestsbee = membees[i];
                doublebee bestsitebee = swarm[i];
                if(i < bestsitescount)//best spots squads
                {
                    bool ind = false;
                    for(int j = scoutbeecount + i * bestbeecount; j < scoutbeecount + (i + 1) * bestbeecount; j++)
                    {
                        swarm[j].flyto(bm, membees[i].x, membees[i].rangeshrink);
                        if(swarm[j].v < bestsbee.v && swarm[j].v < bestsitebee.v)
                        {
                            bestsitebee = swarm[j];
                            ind = true;
                        }
                    }
                    if (ind)
                    {
                        membees[i] = membee(bm, bestsitebee, membees[i].funccounter + 1, bestsitescount + selsitescount);
                    }
                    else
                    {
                        membees[i].funccounter++;
                        if (membees[i].funccounter >= maxfunccounter)
                        {
                            membees[i].funccounter = 0;
                            membees[i].rangeshrink *= rangeshrinker;
                        }
                    }
                }
                else//selected sites squads
                {
                    bool ind = false;
                    for(int j = scoutbeecount + bestsitescount * bestbeecount + (i - bestsitescount) * selectedbeecount; j < scoutbeecount + bestsitescount * bestbeecount + (i - bestsitescount + 1) * selectedbeecount; j++)
                    {
                        swarm[j].flyto(bm, membees[i].x, membees[i].rangeshrink);
                        if(swarm[j].v < bestsbee.v && swarm[j].v < bestsitebee.v)
                        {
                            bestsitebee = swarm[j];
                            ind = true;
                        }
                    }
                    if (ind)
                    {
                        membees[i] = membee(bm, bestsitebee, membees[i].funccounter + 1, bestsitescount + selsitescount);
                    }
                    else
                    {
                        membees[i].funccounter++;
                        if (membees[i].funccounter >= maxfunccounter)
                        {
                            membees[i].funccounter = 0;
                            membees[i].rangeshrink *= rangeshrinker;
                        }
                    }
                }
            }

            int dim = bm.getDim();
            std::vector<double> center(dim);//variables for scouts

            for(int i = 0; i < dim; i++)
            {   
                double a = bm.getBounds()[i].first;
                double b = bm.getBounds()[i].second;           
                center[i] = (a + b) / 2.0;
            }

            for(int i = 0; i < scoutbeecount; i++)
            {
                swarm[i].flyto(bm, center);
            }
        }

        double runcycle(const BM& bm, int maxiter)
        {
            double bestfit = std::numeric_limits<double>::max();;
            int stagnationcount = 0;
            int z;
            for (z = 0; z < maxiter; z++)
            {
                beesatwork(bm);
                sortswarm(bm, scoutbeecount);
                sortmembees();
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
                    std::cout << "algorithm reached stagnation on iteration N" << z << std::endl;
                    break;
                }
            }
            return bestfit;
        }
        ~hive()
        {
            for (int i = 0; i < beestotal; i++)
            {
                swarm[i].~doublebee();
            }
            for (int i = 0; i < bestsitescount + selsitescount; i++)
            {
                membees[i].~membee();
            }
            membees.clear();
            membees.shrink_to_fit();
            swarm.clear();
            swarm.shrink_to_fit();
        }
};

double findMin(const BM& bm, int cnt) {
    const int scoutbeecount = 300;//my parameters
    const int selectedbeecount = 10;
    const int bestbeecount = 30;
    const int selsitescount = 15;
    const int bestsitescount = 5;
    const int maxfunccounter = 1;
    const double rangeshrinker = 0.5;
    const int maxiter = cnt;
    const int maxstag = 30;
    const int beestotal = selectedbeecount * selsitescount + bestbeecount * bestsitescount + scoutbeecount;
    hive *h = new hive(bm, scoutbeecount, selectedbeecount, bestbeecount, selsitescount, bestsitescount, rangeshrinker, maxfunccounter, maxstag);
    h->sortswarm(bm, beestotal);
    double result = h->runcycle(bm, maxiter);
    std::cout << "best value = " << result << std::endl;
    int dim = bm.getDim();
    for (int i = 0; i < dim; i++)
    {
        std::cout << "coordinate " << i << ": " << h->membees[0].x[i] << std::endl;
    }
    delete h;
    return result;
}

bool testBench(const BM& bm) {
    const int ntries = 100;
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm, ntries);
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

main() {
    ZakharovBenchmark<double> zb(3);
    testBench(zb);
    Benchmarks<double> tests;
    for (auto bm : tests) {
        testBench(*bm);
    }
}
