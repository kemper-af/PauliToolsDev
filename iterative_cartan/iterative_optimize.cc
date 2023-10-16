
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
typedef dlib::matrix<double,0,1> column_vector;

static long fcalls = 0;

column_vector costfcn_grad_cvector(const column_vector& thetas, const vector<Pint>& Kgroup, const algebra_element& v, const algebra_element& H)
{
    // This function will obtain the derivatives of the cost function for each element of
    // the group

    column_vector derivs(Kgroup.size());

    // Optimization: don't recompute m every time, just keep track of
    // what modification one makes every iteration

    //for(int j=0; j < Kgroup.size(); j++)
    int j=0;
    algebra_element m1, m2;
    {

        // m1 = \pi_{i=0}^{i=j} Ad_ki  (v) 
        // m2 = \pi_{i<j}  Ad_ki* (H) 

        m1 = v;
        if(j > 0 and false) // Will be optimized away
        {
            for(int i=0; i < j; i++) // From i=0 to i=j-1
            {
                m1 = Ad(thetas(i), Kgroup[i], m1);
            }
        }

        m2 = H;
        for(int i=Kgroup.size()-1; i >= j; i--) // NB: is in reverse order from i=N-1 to i=j
        {
            m2 = Ad(-thetas(i), Kgroup[i], m2);
        }
        

        // Given m1, m2, need to compute
        // term1 =  i tr(kj m1 m2)
        // term2 = -i tr(kj m2  m1)


        // Determine kj m2
        c_algebra_element kj_m2;

        for(const auto& term : m2)
        {
            auto newPauli = symplectic_binary_product(Kgroup[j], term.first);
            kj_m2[ newPauli.first ]       = term.second * newPauli.second;
        }


        c_algebra_element kj_m1;

        for(const auto& term : m1)
        {
            auto newPauli = symplectic_binary_product(Kgroup[j], term.first);
            kj_m1[ newPauli.first ] = term.second * newPauli.second;

        }

        cdouble term1 = inner_product(kj_m1, m2);
        cdouble term2 = inner_product(kj_m2, m1);
        //cout << "d: " << II*(term1-term2) << endl;

        derivs(j) = -imag(term1 - term2);

    }

    // For j larger than 1, start modifying m1 and m2
    for(int j=1; j < Kgroup.size(); j++)
    {
        // Add an Ad for m1
        m1 = Ad(thetas(j-1), Kgroup[j-1], m1);

        // Remove an Ad for m2
        m2 = Ad(thetas(j-1), Kgroup[j-1], m2);

        // Determine kj m2
        c_algebra_element kj_m2;

        for(const auto& term : m2)
        {
            auto newPauli = symplectic_binary_product(Kgroup[j], term.first);
            kj_m2[ newPauli.first ]       = term.second * newPauli.second;
        }


        c_algebra_element kj_m1;

        for(const auto& term : m1)
        {
            auto newPauli = symplectic_binary_product(Kgroup[j], term.first);
            kj_m1[ newPauli.first ] = term.second * newPauli.second;

        }

        cdouble term1 = inner_product(kj_m1, m2);
        cdouble term2 = inner_product(kj_m2, m1);
        //cout << "d: " << II*(term1-term2) << endl;

        derivs(j) = -imag(term1 - term2);
    }


    return derivs;
}

double costfcn_abs2_gradient( const column_vector& thetas, const vector<Pint>& Kgroup, const algebra_element& v, const algebra_element& H)
{
    if(fcalls++ % 1000 == 0)
        cout << "\tCalls to cost function: " << fcalls-1 << endl;

    double dd = dlib::length(
            costfcn_grad_cvector(thetas, Kgroup, v, H)
            );
    return dd * dd;
}

double costfcn_cvector(const column_vector& thetas, const vector<Pint>& Kgroup, const algebra_element& v, const algebra_element& H)
{
    if(fcalls++ % 1000 == 0)
        cout << "\tCalls to cost function: " << fcalls-1 << endl;

    // Produce KvK
    algebra_element KvK = v;
    //dump_algebra_element(v);

    //for(int i=Kgroup.size()-1; i >=0 ; i--)
    for(int i=0; i < Kgroup.size(); i++)
    {
        KvK = Ad(thetas(i), Kgroup[i], KvK);
    }

    //cout << "KvK: " << endl;
    //dump_algebra_element(KvK);

    double dot = inner_product(H, KvK);
    return dot;
}


double costfcn(const vector<double>& thetas, const vector<Pint>& Kgroup, const algebra_element& v, const algebra_element& H)
{
    // Produce KvK
    algebra_element KvK = v;
    //dump_algebra_element(v);

    //for(int i=Kgroup.size()-1; i >=0 ; i--)
    for(int i=0; i < Kgroup.size(); i++)
    {
        KvK = Ad(thetas[i], Kgroup[i], KvK);
    }

    //cout << "KvK: " << endl;
    //dump_algebra_element(KvK);

    double dot = inner_product(H, KvK);
    return dot;
}


vector<double> optimize_KHK(const algebra_element& H, const vector<Pint>& Kgroup, const set<Pint>& algebra_h, const algebra_element& v)
{    
    //cout << "Hello from optimize_KHK" << endl;
    //cout << "Kgroup size: " << Kgroup.size() << endl;
    //cout << "thetas size: " << thetas.size() << endl;

    //=========== VERSION 2 USING VECTOR OF ANGLES AND Ks =======
    // Generate group elements for k
    vector<double> thetas;
    for(int i=0; i < Kgroup.size(); i++)
        thetas.push_back(0.1); // Random angle



    // Produce KvK
    algebra_element KvK = v;

    //for(int i=Kgroup.size()-1; i >= 0; i--)
    for(int i=0; i < Kgroup.size(); i++)
    {
        KvK = Ad(thetas[i], Kgroup[i], KvK);
    }

    //cout << "KvK: " << endl;
    //dump_algebra_element(KvK);

    double dot = inner_product(H, KvK);
    //cout << "\n\tInner product tr(H, KvK) : " << dot << endl;
    //cout << endl;



    // **** DLIB STUFF ***********************************
   
    
    dlib::matrix<double,0,1> thetas_column_vector;
    thetas_column_vector.set_size(thetas.size());
    for(int i=0; i < thetas.size(); i++)
        thetas_column_vector(i) = thetas[i];
   
    // Make a curried cost function
    using namespace::placeholders; //for _1, _2, _3
    using curried_costfcn = function<double(const column_vector&)>;
    using curried_grad_costfcn = function<column_vector(const column_vector&)>;
    curried_costfcn costfcn_single = bind(costfcn_cvector, _1, Kgroup, v, H);
    curried_costfcn costfcn_abs2_grad_single = bind(costfcn_abs2_gradient, _1, Kgroup, v, H);
    curried_grad_costfcn costfcn_grad_single = bind(costfcn_grad_cvector, _1, Kgroup, v, H);

    //cout << "Starting optimization" << endl;

    if(false)
    {
        /* Derivatives from function */
        auto lexd = costfcn_grad_single(thetas_column_vector);

        /* Manually compute derivatives */
        //size_t ktest[] = {0,7,Kgroup.size()-1};
        cout << "Derivatives: " << endl;
        auto dlibd = dlib::derivative(costfcn_single)(thetas_column_vector);
        for(int i=0; i < Kgroup.size(); i++)
        {
            cout << "d: " << dlibd(i) << " " << lexd(i) << endl;
        }
        cout << dlib::length(dlibd - lexd) << endl;

    }
    //return thetas;


    
    /*
    dlib::find_min_bobyqa(costfcn_single,
            thetas_column_vector,
            2*thetas.size() + 1,
            dlib::uniform_matrix<double>(thetas.size(), 1, -M_PI),
            dlib::uniform_matrix<double>(thetas.size(), 1, M_PI),
            0.1,
            1e-8,
            1000000
            );
            */
        
    double tols[] = {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6};
    if(true)
    {
        for(const auto tol : tols)
        {
            cout << "tol: " << tol << endl;
            /*
            dlib::find_max_box_constrained(
                    dlib::lbfgs_search_strategy(10),
                    //dlib::objective_delta_stop_strategy(1e-8),
                    dlib::gradient_norm_stop_strategy(tol), //.be_verbose(),
                    costfcn_single,
                    costfcn_grad_single,
                    //dlib::derivative(costfcn_single),
                    thetas_column_vector,
                    dlib::uniform_matrix<double>(thetas.size(), 1, -M_PI),
                    dlib::uniform_matrix<double>(thetas.size(), 1, M_PI)
                    );
                    */
            /*
            dlib::find_min_bobyqa(
                    costfcn_abs2_grad_single,
                    thetas_column_vector,
                    2*thetas.size() + 1,
                    dlib::uniform_matrix<double>(thetas.size(), 1, -M_PI),
                    dlib::uniform_matrix<double>(thetas.size(), 1, M_PI),
                    0.1,
                    1e-8,
                    1000000
            );
            */
            dlib::find_min_using_approximate_derivatives(
                    dlib::lbfgs_search_strategy(10),
                    dlib::objective_delta_stop_strategy(tol).be_verbose(),
                    costfcn_abs2_grad_single,
                    thetas_column_vector,
                    -1);

            cout << "Function calls: " << fcalls << endl;
        }
    }

    

            
    // Convert it back
    for(int i=0; i < thetas.size(); i++)
        thetas[i] = thetas_column_vector(i);


    // ****** DLIB DONE **********************************


    /*
    algebra_element KHK = H;
    //for(int i=0; i < Kgroup.size(); i++)
    for(int i=Kgroup.size()-1; i >= 0; i--)
    {
        ////if(abs(thetas[i]) < 1e-7)
            //continue;

        KHK = Ad(-thetas[i], Kgroup[i], KHK);
    }
    //dump_algebra_element(KHK, 1e-4);
    //dump_algebra(algebra_h);

    // Output KHK which must lie in h and compute
    // the residual
    double threshold = 1e-4;
    double residual = 0.;
    for(const auto& term : KHK)
    {
        int pauli = term.first;
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
    */

    return thetas;
}


