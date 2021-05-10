# myGebra
A simple Algebra Linear library that contains all Matrix main properties and the main Vector operations

# Test Code

```
#include <iostream>
#include <math.h>
//#include <boost/timer/timer.hpp>
#include "Matrix.h"
#include "LA_Vector.h"
#include <vector>
#include <boost/format.hpp>

using namespace L_Algebra;
using namespace std;

int main()
{

    vector<int> vec;
    vec.push_back(76);
    vec.push_back(76);
    vec.push_back(76);
    vec.push_back(76);
    vec.push_back(76);

    Vector vv(vec);

    int sd = 87, ds = 56;

    Fraction ffr = 10;


    Matrix b(3,4,3);
    Matrix c{5,5,3};

    Matrix a = {{-5, 5, -6, 1, 0}, {0, -5, 10, -3, 3}, {1, 11, 6, 1, 7}, {4, 5, -9, 9, -7}, {-5, 10, 0, -4, 4}};
    Matrix s = {{5, 5, -6, 1, 0}, {3, 4, 5, 7, 8}, {1, 11, 6, 1, 7}, {4, 5, -9, 9, -7}, {5, 10, 0, -4, 4}};
    Matrix s1 = {{5, 5, -6, 1, 0}, {3, 4, 5, 7, 8}, {1, 11, 6, 1, 7}, {4, 5, -9, 9, -7}, {5, 10, 0, -4, 4}};

    cout << a * 23;


    Matrix sw = {{-5}};

    Matrix d = {{1, 0, 2}, {2, 3, 7}};//, {-2, 2, 1, 7}, {-2, 3, 4, 1} };
    Matrix e = {{1, 1}, {0, 0} };
    Matrix g = {{0, 1}, {1, 0} };
    Matrix h = {{1, 0}, {0, 1} };
    Matrix i = {{1, 1}, {0, 1 } };

    // cout << turnMatricesIntoLinearCombination({e, g, h, i});

    try
    {
        //  cout << boost::format("%1% %3%") % 36 % 77 % 34;
    }
    catch (exception& e)
    {
        cout << e.what();
    }


    Matrix f = { {4, 0, 7, 6}, {1, 0, 7, 7}, {8, 0, 8, 8}};//, {-1, -4, -5, 0} };
    Matrix ff = { {4, 2, 7, 6, 5, 6}, {1, 7, 7, 7, 8, 0}, {8, 2, 8, 8, 9, 1}, {-1, -4, -5, 0, 1, 5} };

    Matrix mx1 = { {4, 1, 3, 1}, {3, 1, 3, 0}, {5, 1, 4, 1} };
    Matrix mx11 = { {1, 4, 8, 2}, {1, 4, 4, 9}, {1, 4, 4, 3}, {1, 4, 5, 5} };



   // cout << f << endl << endl;

   // vector<Vector> test = null_space(mx11);

    //cout << f.gaussJordanElimination();

//    for(auto e : test)
//        cout << e << endl;
//
//    cout << endl << nullity(f);


    b(0,2) = 4;
    b(1,2) = 5;
    b(1,3) = 2;
    b(2,0) = -8;
    b(2,3) = 9;
    b(0,0) = 1;
    b(0,1) = 2;



    //cout << mx11 << endl << endl;
    //vector<Vector> test3 = null_space(mx11);

//        for(auto e : test3)
//        cout << e << endl;


    //  cout << mx11.determinant();


    /*

     Vector lav1 = {1, 2, 1};
    Vector lav2 = {2, 9, 0};
    Vector lav = {3, 3, 4};

    Vector lav1 = {1, 5, 3};
    Vector lav2 = {-2, 6, 2};
    Vector lav = {3, -1, 1};

    Vector lav1 = {1, 2, -1};
    Vector lav2 = {6, 4, 2};
    Vector lav3 = {9, 2, 7};

    Vector lav1 = {3, 6, -9, -3};
    Vector lav2 = {6, -2, 5, 1};
    Vector lav3 = {-1, -2, 3, 1};
    Vector lav4 = {2, 3, 0, -2, 0};

    Vector lav3 = {3, 2, 1};
    */

    // cout << p.gaussJordanElimination();

    Matrix mx({ {3, 1, 1, 1}, {5, 2, 3, -2}});//,{-1, -2, 3, 1}});

    //  cout << mx.gaussJordanElimination();


    initializer_list<initializer_list<Fraction>> A = { {1, 3}, {1, -2} };
    initializer_list<Vector> B = { {3, 5}, {1, 2} };
    initializer_list<Vector> C = {{1, 0, 0, 0, }, {-2, 1, 0, 0, }, {5, 3, 0, 0}, {0, 0, 1, 0}, {3, 0, 0, 0} };
    //  Vector vec = {3, 2};

    Matrix gt(A);
    Matrix wz = { {0, 0, 0, 2, 9, 6}, {0, 0, 0, 4, 5, 8} };
    Matrix wzf = { {3, 2, 9, 2, 9, 6}, {6, 4, 5, 4, 5, 8} };
    Matrix z = { {1, 3, -2, 0, 2, 0}, {2, 6, -5, -2, 4, -3}, {0, 0, 5, 10, 0, 15}, {2, 6, 0, 8, 4, 18} };

//    cout << gt;

    Matrix dz = { {4, 1, 5, 1, 7, 8, 2}, {6, 3, 3, 5, 2, 3, 1}};//, {0, 0, 5, 10, 0, 15}, {2, 6, 0, 8, 4, 18} };

    Matrix fz = { {1, 3, 4, 4}, {2, 3, 5, 4}, {9, 1, 7, 2}};// {-1, -4, -5, 0} };
    Matrix tfz = { {1, 3, 4, 4, 1}, {2, 3, 5, 4, 5}, {9, 1, 7, 2, 3}};// {-1, -4, -5, 0} };

    Matrix khan = { {1, 1, 2, 3, 2}, {1, 1, 3, 1, 4} };
    Matrix kha = { {2, 0, 2}, {-1, 0, -1}, {-1, 0, -1} };

//    boost::timer::cpu_timer timer;
//    wz.gaussJordanElimination();
    //  timer.stop();


    //  cout << timer.format();

    Vector lav1 = {0, -2, 2};
    Vector lav2 = {1, 3, -1};
    Vector lav3 = {9, 0, 0};
    Vector lav4 = {4, 0, 2};
    Vector v = { 0, 0, 0};

    Matrix p = { {4, 0}, {-2, -2} };
    Matrix ph = { {1, -1}, {2, 3} };
    Matrix ph1 = { {0, 2}, {1, 4} };
    Matrix ph2 = { {-1, 5}, {7, 1} };
    Matrix ph21 = { {6, -8}, {-1, -8} };
    Matrix ph3 = { {6, 0}, {3, 8} };
    Matrix ph0 = { {0, 0}, {0, 0} };

    Fraction fr1(27, 17);
    Fraction fr2(43, 34);
    Fraction fr3(-29, 306);

    Matrix mcf(3, 3, {2, 3, 5, 6, 4, 5, 5, 8, 9});

    double db = 10.0 / 3;

    Fraction frt;



    // cout << frt;

    // cout << s << endl;


    try
    {
//        cout << s.main_diagonal() << endl;
//        cout << s.secondary_diagonal() << endl;

        //cout << coordinate_vector_relative_to_basis({ {0,1,0}, { {-4,5}, 0, {3,5}, }, { {3,5}, 0, {4,5} } }, {1,1,1});

        //cout << change_basis(vec, A, B);

        //cout << kha.gaussJordanElimination() << endl;

        //vector<Vector> v = null_space(kha);
        //  cout << coordinate_vector_relative_to_basis({ lav1, lav2,lav3}, lav4);

        // for(auto e : v)
        //    cout << e << endl;

        //  cout << endl << khan.rank();
    }
    catch(exception& e)
    {
        cout << e.what();
    }

//cout << lav2 * (lav ->* lav1);

}

```
