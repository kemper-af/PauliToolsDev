#include<string>
#include<fstream>
#include<iostream>
#include<map>
#include<unordered_map>
#include<set>
#include<vector>
#include<algorithm>
#include<cassert>
#include<cstring>
#include<stdio.h>
#include<tuple>
#include<random>
#include<functional>
#include<complex>
  using namespace std;

#include "cartan/Pauli.h"
#include "cartan/algebra_types.h"
#include "cartan/dla_generator.h"
#include "cartan/models.h"
#include "cartan/involutions.h"

typedef long Pint;
typedef map<Pint, double> algebra_element;
typedef map<Pint, cdouble> c_algebra_element;



cdouble II(0,1);


#ifdef _OPENMP
    #include<omp.h>
#endif





#include "helper_functions.cc"
//#include "iterative_optimize.cc"

enum modelType{ TFXY, HEIS };



template<typename container>
void dump_algebra(container& algebra)
{    
    cout << "\t";
    for(const auto& term : algebra)
    {
        cout << term << " ";
    }
    cout << endl;
}





void optimize_KHK(const vector<string>& ham, const set<Pint>& algebra_k, const set<Pint>& algebra_h);



int main(int argc, char **argv)
{
#ifdef _OPENMP
    {
    cout << "\nOpenMP Information: ---------------------" << endl;
    #pragma omp parallel
    {
        //int my_th_id = omp_get_thread_num();
        //printf("Hello from thread %i\n", my_th_id);
        int numthreads = omp_get_num_threads();
        #pragma omp single
        {
            printf("There are %i threads\n", numthreads);
        }
    }
    cout << "-----------------------------------------\n" << endl;
    }
#endif

    // ====================== READ COMMAND LINE ================

    if(not (argc==3))
    {
        throw std::invalid_argument("Usage: iterative_cartan Nq MODEL");
    }

    int newNq = atoi(argv[1]);
    if(newNq < 2 or newNq > 24)
    {
        cout << "Nq: " << Nq << endl;
        throw std::invalid_argument("Wrong value for Nq.");
    }
    setNq(newNq);
    cout << endl;


    modelType model;
    {
    if(strcmp(argv[2], "TFXY")==0)
    {
        model = TFXY;
        cout << "Using TFXY model" << endl;
    }
    else if(strcmp(argv[2], "HEIS")==0)
    {
        model = HEIS;
        cout << "Using Heisenberg model" << endl;
        if(Nq > 8)
        {
            cout << "You set up Heisenberg for Nq > 8\nThis is not OK" << endl;
            throw std::invalid_argument("Problem size too large");
        }
    }
    else
        throw std::invalid_argument("Expected a model argument: TFXY or HEIS");
    }

    cout << "model: " << model << endl;

    // =======================================================

    /* Get Ham from function */
    vector<string> ham;
    switch(model)
    {
        case HEIS : ham = heisenberg(false); break;
        case TFXY : ham = tfxy(false); break;
    }

    cout << "Hamiltonian vector: " << endl;
    for(const auto& term : ham)
    {
        cout << term << ", ";
    }
    cout << "\n" << endl;

    Pauli p;


    // Fill the Hamiltonian algebra objects
    algebra algebra_g, algebra_h;

    for(const auto& hamterm : ham)
    {
        algebra_g.insert(Pauli(hamterm));
    }

    // See what's in it
    cout << "Hamiltonian algebra: " << endl;
    dump_algebra(algebra_g);


    // Produce an initial guess for algebra_h
    cout << "\nInitial guess for h:" << endl;
    if(model == TFXY)
        for(size_t i=0; i < Nq; i++)
        {
            algebra_h.insert(Pauli(ham[i]));
        }
    else if(model == HEIS)
        for(size_t i=0; i < ham.size(); i += 6)
        {
            algebra_h.insert(Pauli(ham[i]));
            algebra_h.insert(Pauli(ham[i+1]));
        }


    dump_algebra(algebra_h);
    cout << "\n" << endl;
    


    // Produce the full algebra
    auto&& [algebra_k, algebra_m ] = get_algebra_by_commuting(
            &algebra_g, &algebra_h, parityY);







    // Examine final h
    cout << "\nFinal results for h:" << endl;
    dump_algebra(algebra_h);
    cout << "\n" << endl;


    set<Pint> algebra_b;
    if(model == TFXY)
    {
        // For TFXY, b = h
        algebra_b = algebra_h;
    }
    else if(model == HEIS)
    {
        // For Heisenberg, use XX---- YY---- --XX-- --YY-- etc
        for(const auto& h: algebra_h)
        {
            
            /*
            dump_pauli(h);
            cout << endl;
            cout << countX(h, Nq) << " " <<  countY(h, Nq) << " " <<  countZ(h, Nq) << endl;
            */
            if(countZ(h, Nq) > 0)
                continue;
            
            if(countX(h, Nq) == 2 and countY(h, Nq) == 0)
                algebra_b.insert(h);

            if(countX(h, Nq) == 0 and countY(h, Nq) == 2)
                algebra_b.insert(h);
        }
    }


    cout << "b set: " << endl;
    dump_algebra(algebra_b);
    cout << endl;

    // Group the k by b ------------------------------
    map<Pint,set<Pint>> kgrouping;
    auto algebra_k_copy(algebra_k);
    for(const auto& b: algebra_b)
    {
        set<Pint> bk;
        for(const auto& k: algebra_k_copy)
        {
            Pint comm = symplectic_binary_comm(b, k);
            if(comm==0)
                continue;
            else
                bk.insert(k);
        }

        // Now remove the found ones
        for(const auto& k: bk)
        {
            algebra_k_copy.erase(k);
        }

        kgrouping[b] = bk;
    }

    for(const auto& b_bk: kgrouping)
    {
        cout << "b term: ";
        dump_pauli(b_bk.first);
        cout << endl << "\t";
        for(const auto& k: b_bk.second)
        {
            dump_pauli(k);
            cout << " ";
        }
        cout << endl;
    }
    // ------------------------------------------

    cout << "===================== Going into b optimization =====================" << endl;


    // Produce the H that will be our target.  As the optimization happens, this will
    // get updated to K_{n-1}+ ... K2+ K1+ H K1 K2 ... K_{n-1}

    // Make the Hamiltonian algebra element
    algebra_element KHK;
    for(const auto& hamterm : ham)
    {
        Pint term = string_to_binary_rep(hamterm);
        KHK[term] = 1.0;
    }
    

    int idx=0;
    map<Pint, vector<double>> thetas;
    for(const auto& b_bk: kgrouping)
    {
        cout << "\nWorking on the b string: ";
        dump_pauli(b_bk.first);
        cout << endl;

        if(b_bk.second.size() == 0)
        {
            cout << "\tReached the last element!" << endl;
        }
        else 
        {
            vector<Pint> Kgroup;
            for(const auto& k : b_bk.second)
                Kgroup.push_back(k);

            algebra_element b;
            b[b_bk.first] = 1.0;
            const auto& kset = b_bk.second;
            //KLUGE
            vector<double> result_thetas;// = optimize_KHK(KHK, Kgroup, algebra_h, b);
            thetas[b_bk.first] = result_thetas;


            // Now update KHK
            for(int i=Kgroup.size()-1; i >= 0; i--)
            {
                KHK = Ad(-result_thetas[i], Kgroup[i], KHK);
            }

            // Truncate small terms
            // This is awful code
            {
                vector<Pint> to_be_removed;
                for(const auto& pauli : KHK)
                {
                    if(abs(pauli.second) < 1e-7)
                        to_be_removed.push_back(pauli.first);
                }

                for(const auto& pauli: to_be_removed)
                    KHK.erase(pauli);
            }
        }


        // Check if KHK commutes with b1
        /*
        cout << "KHK: " << endl;
        for(const auto& pauli : KHK)
        {
            cout << pauli.second << "\t";
            dump_pauli(pauli.first);
            cout << " ";
            dump_pauli(b_bk.first);
            cout << " ";
            int comm = symplectic_binary(b_bk.first, pauli.first);
            cout << comm << endl;

        }
        cout << endl;
        */


        double threshold = 1e-4;
        double residual = 0.;
        for(const auto& term : KHK)
        {
            Pint pauli = term.first;
            double coef = term.second;

            if(algebra_h.count(pauli) == 0)
                residual += abs(term.second);
        }
        cout << "\tResidual: " << residual << endl;

        //if(fabs(residual) < 1e-8)
            //break;

        // KLUGE
        //break;
    }
   

    cout << "\nOptimization achieved." << endl;
    //cout << "Total number of calls to cost function: " << fcalls << endl;


    // Test the final overlap
    {
        algebra_element KHK;
        for(const auto& hamterm : ham)
        {
            Pint term = string_to_binary_rep(hamterm);
            KHK[term] = 1.0;
        }

        for(const auto& b_bk: kgrouping)
        {
            if(b_bk.second.size() == 0)
                continue;

            vector<Pint> Kgroup;
            for(const auto& k : b_bk.second)
                Kgroup.push_back(k);

            for(int i=Kgroup.size()-1; i >= 0; i--)
            {
                KHK = Ad(-thetas[b_bk.first][i], Kgroup[i], KHK);
            }

        }

        double threshold = 1e-4;
        double residual = 0.;
        for(const auto& term : KHK)
        {
            Pint pauli = term.first;
            double coef = term.second;

            if(algebra_h.count(pauli) == 0)
                residual += abs(term.second);

            if(abs(coef) < threshold)
                continue;
            //cout << "\t" << coef << " ";
            //dump_pauli(pauli);
            //cout << " + " << endl;

        }
        cout << "\tResidual: " << residual << endl;
    }


    return 0;
}






