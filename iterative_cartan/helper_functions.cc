
map<string, double> ham_from_file(string filename)
{
    map<string, double> ham;
    ifstream f(filename);
    string pauli;
    double coef;
    while(f >> coef >> pauli)
    {
        //cout << coef << " " << pauli << endl;
        ham[pauli] = coef;
    }
    f.close();

    return ham;
}



pair<Pint,double> Ad_symplectic_binary_product(const Pint c1, const Pint c2);
algebra_element Ad(const double theta, const Pint k, const algebra_element& b);

/*
Pint symplectic_binary(const Pint c1, const Pint c2)
{
    const int N = Nq;
    Pint a1 = (c1 & a_mask);
    Pint b1 = (c1 & b_mask) >> N;
    Pint a2 = (c2 & a_mask) ;
    Pint b2 = (c2 & b_mask) >> N;


    if(__builtin_parity(a1 & b2) ^ __builtin_parity(a2 & b1))
        return c1 ^ c2;
    else
        return 0;
}
*/

void print_int_as_binary(Pint a, int maxbits)
{
    Pint x=1;
    for(int i=maxbits-1; i >= 0; i--)
        cout << bool(a & (x << i));
}

/*
inline bool involution_countY(Pint p, int N)
{
    pauli_int c = p.code;
    return __builtin_parity(  (c & a_mask) & ((c & b_mask) >> N) );
    return __builtin_popcount( (c & a_mask) & ((c & b_mask) >> N) );
}

inline bool involution_countZ(Pint c, int N)
{
    return __builtin_parity(  ~(c & a_mask) & ((c & b_mask) >> N) );
    //return __builtin_popcount(  ~(c & a_mask) & ((c & b_mask) >> N) );
}

inline bool involution_countX(Pint c, int N)
{
    return __builtin_parity(  (c & a_mask) & ~((c & b_mask) >> N) );
    //print_int_as_binary(c, 2*N);
    //cout << " ";
    //print_int_as_binary( (c & a_mask) , N);
    //cout << " ";
    //print_int_as_binary( ~((c & b_mask)>>N) , N);
    //cout << " ";
    //return __builtin_popcount( (c & a_mask) & ~((c & b_mask) >> N) );
}

inline int countY(Pauli p, int N)
{
    pauli_int c = p.code;
    //return __builtin_parity(  (c & a_mask) & ((c & b_mask) >> N) );
    return __builtin_popcount( (c & a_mask) & ((c & b_mask) >> N) );
}

inline int countZ(Pint c, int N)
{
    //return __builtin_parity(  ~(c & a_mask) & ((c & b_mask) >> N) );
    return __builtin_popcount(  ~(c & a_mask) & ((c & b_mask) >> N) );
}

inline int countX(Pint c, int N)
{
    //return __builtin_parity(  (c & a_mask) & ~((c & b_mask) >> N) );
    //print_int_as_binary(c, 2*N);
    //cout << " ";
    //print_int_as_binary( (c & a_mask) , N);
    //cout << " ";
    //print_int_as_binary( ~((c & b_mask)>>N) , N);
    //cout << " ";
    return __builtin_popcount( (c & a_mask) & ~((c & b_mask) >> N) );
}

*/



void dump_pauli(const Pint encoded_pauli, ostream& out)
{
    Pint x = 1;
    for(int i=0; i < Nq; i++)
    {
        bool a = encoded_pauli & (x << i);
        bool b = encoded_pauli & (x << (i+Nq));
 
        if(a and b)
            out << 'Y';
        else if(a)
            out << 'X';
        else if(b)
            out << 'Z';
        else
            out << '-';
    }
}

void dump_pauli(const Pint encoded_pauli)
{
    dump_pauli(encoded_pauli, cout);
}




template<typename T>
void dump_algebra_element(const map<Pint, T>& g, double threshold=1e-8)
{
    for(const auto& term : g)
    {
        Pint pauli = term.first;
        T coef = term.second;
        if(abs(coef) < threshold)
            continue;
        cout << "\t" << coef << " ";
        dump_pauli(pauli);
        cout << " + " << endl;
    }
}

pair<Pint,double> Ad_symplectic_binary_product(const pair<Pint,double> c1_pair, const pair<Pint,double> c2_pair)
{
    const int N = Nq;

    Pint c1 = c1_pair.first;
    Pint c2 = c2_pair.first;

    return Ad_symplectic_binary_product(c1, c2);

}

pair<Pint,cdouble> symplectic_binary_product(const Pint c1, const Pint c2)
{
    /*
     * Work out the product with the sign included
     * A B = C
     * The terms that have a negative sign are:
     *     Y X = -iZ 
     *     Z Y = -iX 
     *     X Z = -iY 
     * the X component is a, the Z component is b. Printed with the
     * convention above, that is written b|a. With that,
     * the three cases above are
     *     1|1  0|1  =  - 1|0
     *     1|0  1|1  =  - 0|1
     *     0|1  1|0  =  - 1|1
     * in logic gates,
     *   ( b1 & a1) & (~b2 & a2)
     *   ( b1 &~a1) & ( b2 & a2)
     *   (~b1 & a1) & ( b2 & ~a2)
     *
     * These result in either a -1 (per qubit). The total number of
     * negative signs is the parity of this operation.
     *
     * The total number of Is is the number of not-equal Paulis
     */



    const int N = Nq;


    Pint a1 = (c1 & a_mask);
    Pint b1 = (c1 & b_mask) >> N;
    Pint a2 = (c2 & a_mask) ;
    Pint b2 = (c2 & b_mask) >> N;

//#define verbose_ad1
    
    
#ifdef verbose_ad1
    cout << "Going into product." << endl;
    cout << "op1:\t";
    print_int_as_binary(a1, Nq);
    cout << "|";
    print_int_as_binary(b1, Nq);
    cout << endl;
    cout << "op2:\t";
    print_int_as_binary(a2, Nq);
    cout << "|";
    print_int_as_binary(b2, Nq);
    cout << endl;

    cout << "(a1 & ( b1) & a2 & ~b2) : ";
    print_int_as_binary((a1 & ( b1) & a2 & ~b2), Nq);
    cout << endl;
    cout << "(~a1 & b1 & (a2) & b2) : ";
    print_int_as_binary((~a1 & b1 & (a2) & b2), Nq);
    cout << endl;
    cout << "((a1) & ~b1 & ~a2 & (b2)) : ";
    print_int_as_binary(((a1) & b1 & ~a2 & (b2)), Nq);
    cout << endl;
    cout << "a1 & a2 & (~b1) & (~b2): "; // X equality
    print_int_as_binary(a1 & a2 & (~b1) & (~b2), Nq);
    cout << endl;
    cout << "a1 & a2 & b1 & b2: ";      // Y equality
    print_int_as_binary(a1 & a2 & b1 & b2, Nq);
    cout << endl;
    cout << "(~a1) & (~a2) & (b1 & b2) : "; // Z equality
    print_int_as_binary((~a1) & (~a2) & (b1 & b2), Nq);
    cout << endl;
#endif   
    
    

    Pint product = c1 ^ c2;


    Pint a3 = (product & a_mask);
    Pint b3 = (product & b_mask) >> N;

    
    //int p1 = __builtin_popcount(a1 | b1);
    //int p2 = __builtin_popcount(a2 | b2);
    //int p3 = __builtin_popcount(a3 | b3);

    // Count overlapping Paulis
    int p4 = __builtin_popcount((a1 | b1) & (a2 | b2));

    Pint sign_counter = (a1 & (b1) & a2 & ~b2) |    // Y X = -iZ
                       (~a1 & b1 & (a2) & b2) |     // Z Y = -iX
                       ((a1) & ~b1 & ~a2 & (b2)) ;  // X Z = -iY
    Pint equal_paulis = 
                       a1 & a2 & (~b1) & (~b2) | // X equality
                       (~a1) & (~a2) & b1 & b2 | // Y equality
                       a1 & a2 & b1 & b2;        // Z equality

#ifdef verbose_ad1
    cout << "op3:\t";
    print_int_as_binary(a3, Nq);
    cout << "|";
    print_int_as_binary(b3, Nq);
    cout << endl;
    
    cout << "Sign counter: ";
    print_int_as_binary(sign_counter,2*Nq);
    cout << " " << __builtin_parity(sign_counter) << endl;
                       
    
    cout << "Paulis in p1 : " << p1 << endl;
    cout << "Paulis in p2 : " << p2 << endl;
    cout << "Paulis in p3 : " << p3 << endl;
    cout << "Overlapping Paulis : " << p4 << endl;
    cout << "Negative signs: " << __builtin_popcount(sign_counter) << endl;
    cout << "Sign parity: " << __builtin_parity(sign_counter) << endl;
    cout << "Equal Paulis: " << __builtin_popcount(equal_paulis) << endl;
#endif   
    

    int icount = (p4 -  __builtin_popcount(equal_paulis)) % 4;
    while(icount < 0)
        icount += 4;
    while(icount >= 4)
        icount -= 4;

    const cdouble sign_arr[4] = {1, II, -1, -II};
    cdouble sign = sign_arr[icount];


                      
    if(__builtin_parity(sign_counter))
        sign *= -1;

    return pair(product, sign);
}

pair<Pint,double> Ad_symplectic_binary_product(const Pint c1, const Pint c2)
{
    /*
     * Work out the product with the sign included
     * A B = C
     * (-i) (iA) (iB) = (iC)
     * The terms that have a negative sign are:
     *     X Y = iZ -> iX iY = -iZ
     *     Y Z = iX -> iY iZ = -iX
     *     Z X = iY -> iZ iX = -iY
     * the X component is a, the Z component is b. Printed with the
     * convention above, that is written b|a. With that,
     * the three cases above are
     *     0|1  1|1  =  - 1|0
     *     1|1  1|0  =  - 0|1
     *     1|0  0|1  =  - 1|1
     * in logic gates,
     *   (~b1 & a1) & (b2 & a2)
     *   ( b1 & a1) & (b2 & ~a2)
     *   ( b1 &~a1) & (~b2 & a2)
     *
     * These result in either 1 or 0 (per qubit). The total number of
     * negative signs is the parity of this operation.
     */



    const int N = Nq;


    Pint a1 = (c1 & a_mask);
    Pint b1 = (c1 & b_mask) >> N;
    Pint a2 = (c2 & a_mask) ;
    Pint b2 = (c2 & b_mask) >> N;

//#define verbose_ad1
    
    
#ifdef verbose_ad1
    cout << "Going into commuting." << endl;
    cout << "op1:\t";
    print_int_as_binary(a1, Nq);
    cout << "|";
    print_int_as_binary(b1, Nq);
    cout << endl;
    cout << "op2:\t";
    print_int_as_binary(a2, Nq);
    cout << "|";
    print_int_as_binary(b2, Nq);
    cout << endl;

    cout << "(a1 & (~b1) & a2 & b2) : ";
    print_int_as_binary((a1 & (~b1) & a2 & b2), Nq);
    cout << endl;
    cout << "(a1 & b1 & (~a2) & b2) : ";
    print_int_as_binary((a1 & b1 & (~a2) & b2), Nq);
    cout << endl;
    cout << "((~a1) & b1 & a2 & (~b2)) : ";
    print_int_as_binary(((~a1) & b1 & a2 & (~b2)), Nq);
    cout << endl;
    cout << "a1 & a2 & (~b1) & (~b2): "; // X equality
    print_int_as_binary(a1 & a2 & (~b1) & (~b2), Nq);
    cout << endl;
    cout << "a1 & a2 & b1 & b2: ";      // Y equality
    print_int_as_binary(a1 & a2 & b1 & b2, Nq);
    cout << endl;
    cout << "(~a1) & (~a2) & (b1 & b2) : "; // Z equality
    print_int_as_binary((~a1) & (~a2) & (b1 & b2), Nq);
    cout << endl;
    
#endif   
    

    Pint product = c1 ^ c2;

    
    if(not(__builtin_parity(a1 & b2) ^ __builtin_parity(a2 & b1))) // They commute
    {
        return pair(product, 0.0);
    }
    

    Pint a3 = (product & a_mask);
    Pint b3 = (product & b_mask) >> N;

    
    int p1 = __builtin_popcount(a1 | b1);
    int p2 = __builtin_popcount(a2 | b2);
    int p3 = __builtin_popcount(a3 | b3);

    int sign = 1;
    Pint sign_counter = (a1 & (~b1) & a2 & b2) | 
                       (a1 & b1 & (~a2) & b2) | 
                       ((~a1) & b1 & a2 & (~b2)) |
                       a1 & a2 & (~b1) & (~b2) |
                       (~a1) & (~a2) & b1 & b2 |
                       a1 & a2 & b1 & b2;

#ifdef verbose_ad1
    cout << "op3:\t";
    print_int_as_binary(a3, Nq);
    cout << "|";
    print_int_as_binary(b3, Nq);
    cout << endl;
    
    cout << "Sign counter: ";
    print_int_as_binary(sign_counter,2*Nq);
    cout << " " << __builtin_parity(sign_counter) << endl;
                       
    
    cout << "p1 : " << p1 << endl;
    cout << "p2 : " << p2 << endl;
    cout << "p3 : " << p3 << endl;
    cout << "p1 + p2 - p3 + 1: " << p1 + p2 - p3 + 1 << endl;
    

#endif   
    

    if(((p1 + p2 - p3 - 1)/2) % 2)
        sign *= -1;
                      
    if(__builtin_parity(sign_counter))
        sign *= -1;

    return pair(product, sign);
}

algebra_element Ad(const pair<Pint, double>& k, const algebra_element& b)
{
    return Ad(k.second, k.first, b);

    /*
    // Implements Ad_(exp(i th k) (b)
    double cos2th = cos(2*k.second), sin2th = sin(2*k.second);

    // Make a copy
    algebra_element c;

    for(const auto& x : b)
    {
        if(abs(x.second) < 1e-8) // Coefficient was 0
            continue;



        auto y = Ad_symplectic_binary_product(k, x);

        //cout << "Product: ";
        //dump_pauli(y.first);
        //cout << '\t' << y.second << endl;
        if(abs(y.second) < 1e-8) // Either they commute or the coefficient was 0
        {



            // They commuted
            if(c.count(x.first) > 0)
                c[x.first] += x.second;
            else
                c[x.first] = x.second;

            continue;
        }



        if(c.count(x.first) > 0)
            c[x.first] += x.second * cos2th;
        else
            c[x.first] = x.second * cos2th;

        if(c.count(y.first) > 0)
            c[y.first] += y.second * sin2th;
        else
            c[y.first] = y.second * sin2th;
    }


    return c;
    */
}


algebra_element Ad(const double theta, const Pint k, const algebra_element& b)
{
    // Implements Ad_(exp(i th k) (b)
    double cos2th = cos(2*theta), sin2th = sin(2*theta);

    // Make a copy
    algebra_element c;

#ifdef verbose_ad
    cout << "Ad with ";
    dump_pauli(k);
    cout << endl;
#endif

    for(const auto& bx : b)
    {
#ifdef verbose_ad
        cout << "New element is now:\n";
        dump_algebra_element(c);
        cout << "element done.";
        cout << endl;
#endif

        if(abs(bx.second) < 1e-8) // Coefficient was 0
            continue;

        
#ifdef verbose_ad
        cout << "\nEvaluating element of b:";
        dump_pauli(bx.first);
        cout << ", " << bx.second;
        cout << endl;
#endif
        
        auto y = Ad_symplectic_binary_product(k, bx.first);

        //cout << "Product: ";
        //dump_pauli(y.first);
        //cout << '\t' << y.second << endl;
        if(abs(y.second) < 1e-8) // Either they commute or the coefficient was 0
        {

            
#ifdef verbose_ad
            cout << "Found a commuting one : " ;
            dump_pauli(bx.first);
            cout << " which gives : " ;
            dump_pauli(y.first);
            cout << endl;
#endif
            

            // They commuted
            if(c.count(bx.first) > 0)
                c[bx.first] += bx.second;
            else
                c[bx.first] = bx.second;

            continue;
        }

        
#ifdef verbose_ad
        cout << "Found a non-commuting one : ";
        dump_pauli(bx.first);
        cout << " which gives : " ;
        dump_pauli(y.first);
        cout << endl;
        cout << "The product has sign " << y.second << endl;
        cout << "Adding coefs " << bx.second * cos2th << " and " << bx.second * sin2th << endl;
#endif
        

        if(c.count(bx.first) > 0)
            c[bx.first] += bx.second * cos2th;
        else
            c[bx.first] = bx.second * cos2th;

        if(c.count(y.first) > 0)
            c[y.first] += bx.second * y.second * sin2th;
        else
            c[y.first] = bx.second * y.second * sin2th;
    }


    //cout << "I'm done, getting out!" << endl;
    return c;
}

double inner_product(const vector<double>& theta, const vector<Pint> b, const algebra_element& a)
{
    double dot=0;
    for(const auto& aterm : a)
    {
        Pint akey = aterm.first;
        vector<Pint>::const_iterator b_iter = find(b.begin(), b.end(), akey);
        if(b_iter != b.end()) // Found one!
        {
            const auto& b_idx = *b_iter;
            
            
            /*
            cout << "Matched " << aterm.second;
            dump_pauli(akey);
            cout << " and " << theta[b_idx];
            dump_pauli(b_idx);
            cout << endl;
            */
            

            dot += aterm.second * theta[b_idx];
        }
    }
    return dot;

}

double inner_product(const algebra_element& a, const algebra_element& b)
{
    // For now do the dumb thing
    double dot=0;
    for(const auto& aterm : a)
    {
        Pint akey = aterm.first;
        algebra_element::const_iterator b_iter = b.find(akey);
        if(b_iter != b.end()) // Found one!
        {
            const auto& bterm = *b_iter;
            
            
            
            /*
            cout << "Matched " << aterm.second;
            dump_pauli(akey);
            cout << " and " << bterm.second;
            dump_pauli(bterm.first);
            cout << endl;
            */

            dot += aterm.second * bterm.second;
            
        }
    }
    return dot;
}

cdouble inner_product(const c_algebra_element& a, const algebra_element& b)
{
    // For now do the dumb thing
    cdouble dot=0;
    for(const auto& aterm : a)
    {
        Pint akey = aterm.first;
        algebra_element::const_iterator b_iter = b.find(akey);
        if(b_iter != b.end()) // Found one!
        {
            const auto& bterm = *b_iter;
            
            if(false)
            if(abs(aterm.second) > 1e-8 and abs(bterm.second) > 1e-8)
            {
                dump_pauli(akey);
                cout << "\t" << aterm.second << " * " << bterm.second << endl;
            }
            
            
            /*
            cout << "Matched " << aterm.second;
            dump_pauli(akey);
            cout << " and " << bterm.second;
            dump_pauli(bterm.first);
            cout << endl;
            */

            dot += aterm.second * bterm.second;
            
        }
    }
    return dot;
}
