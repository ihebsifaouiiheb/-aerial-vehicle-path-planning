#include<iostream>
#include"ilcplex/ilocplex.h";

#include<chrono>
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <float.h>
#include<iomanip>

using namespace std;
typedef IloArray<IloNumVarArray> NumVar2D;
typedef IloArray<NumVar2D> NumVar3D;


/*float Uij(int i, int j) {
    return sqrt((v[i].first - v[j].first) * (v[i].first - v[j].first) + (v[i].second - v[j].second) * (v[i].second - v[j].second));
}*/


int main()
{
#pragma region Problem Data
    freopen("CMT1.txt", "r", stdin);
    int startingNodeIndex, endingNodeIndex, Dimension;
    int R = 2, G = 1000, T = 100;
    vector<pair<float, float>> v;
    //int nnbreR = 2;
    //int nG = 1;
    //int nT = 1;

    //int* nbreR = new int[R] {1, 2};
    //int* G = new int[nG] {1000};
    //int* T = new int[nT] {100};

    string s;
    while (s != "DIMENSION") {
        cin >> s;
    }
    cin >> s;
    cin >> Dimension;

    while (s != "NODE_COORD_SECTION") {
        cin >> s;
    }

    /*bool*** x;
    int **t = new int*[Dimension];

    x = new bool** [Dimension];
    for (int i = 0; i < Dimension; ++i) {
        x[i] = new bool* [Dimension];
        for (int y = 0; y < R; ++y) {
            x[i][y] = new bool[R];
            for (int z = 0; z < R; ++z) {
                x[i][y][z] = 0;
            }
        }
    }

    for (int i = 0; i < R; i++)
        t[i] = new int[R];*/


    int colWidth = 10;
    cout << setfill('*') << setw(3 * colWidth) << "*" << endl;
    cout << setfill(' ') << fixed;
    cout << setw(colWidth) << "Node" << setw(colWidth) << "Abscissa" << setw(colWidth) << "Ordinate" << endl;
    cout << setfill('*') << setw(3 * colWidth) << "*" << endl;
    cout << setfill(' ') << fixed;

    float n, z, y;
    for (int i = 0; i < Dimension; ++i) {
        cin >> n >> z >> y;
        cout << setprecision(0) << setw(colWidth) << n << setprecision(4) << setw(colWidth)
            << (int)z << setw(colWidth) << (int)y << endl;
        v.push_back({ z,y });
    }
    startingNodeIndex = 0;
    endingNodeIndex = v.size() - 1;
    //n = 5;
    //cout << "V[" << n << "] = " << v[n].first <<" "<< v[n].second << " " << Uij(n, n+1);

    float** U = new float* [Dimension];
    for (int i = 0; i < Dimension; i++) {
        U[i] = new float[Dimension];
        for (int j = 0; j < Dimension; ++j) {
            U[i][j] = sqrt((v[i].first - v[j].first) * (v[i].first - v[j].first) +
                (v[i].second - v[j].second) * (v[i].second - v[j].second));
        }
    }
#pragma endregion

    IloEnv env;
    IloModel Model(env);

#pragma region Decision variables
    NumVar2D t(env, Dimension);
    for (int i = 0; i < Dimension; ++i) {
        t[i] = IloNumVarArray(env, R, 0, T, ILOINT);
    }

    NumVar3D X(env, Dimension);
    for (int i = 0; i < Dimension; ++i) {
        X[i] = NumVar2D(env, Dimension);
        for (int j = 0; j < Dimension; ++j) {
            X[i][j] = IloNumVarArray(env, R, 0, 1, ILOBOOL);
            /*for (int k = 1; k <= R; ++k) {
                //X[i][j][k] = IloNumVar(env, 0, 1, ILOBOOL);
            }*/
        }
    }
#pragma endregion
    

#pragma region Object Function
    IloExpr exp0(env);

    for (int i = 1; i < Dimension - 1; ++i) {
        for (int k = 0; k < R; ++k) {
            exp0 += t[i][k];
        }
    }

    Model.add(IloMinimize(env, exp0));
#pragma endregion

#pragma region Constraints
    // constraint 2
    for (int i = 1; i < Dimension-1; ++i) {
        for (int k = 0; k < R; ++k) {
            IloExpr exp1(env), exp2(env);
            for (int j = 0; j < Dimension; ++j) {
                exp1 += X[j][i][k];
                exp2 += X[i][j][k];
            }
            Model.add(exp1 == exp2);
        }
    }

    // constraint 3
    for (int i = 1; i < Dimension - 1; ++i) {
        IloExpr exp3(env);
        for (int k = 0; k < R; ++k) {
            for (int j = 0; j < Dimension; ++j) {
                exp3 += X[i][j][k];
            }
            Model.add(exp3 == 1);
        }
    }

    //constraint 4
    for (int k = 0; k < R; ++k) {
        IloExpr exp4(env);
        for (int j = 0; j < Dimension; ++j) {
            exp4 += X[startingNodeIndex][j][k];
        }
        Model.add(exp4 == 1);
    }
    
    //constraint 5
    for (int k = 0; k < R; ++k) {
        IloExpr exp5(env);
        for (int j = 0; j < Dimension; ++j) {
            exp5 += X[j][endingNodeIndex][k];
        }
        Model.add(exp5 == 1);
    }

    //constraint 6
    for (int k = 0; k < R; ++k) {
        IloExpr exp6(env);
        for (int i = 0; i < Dimension; ++i) {
            for (int j = 0; j < Dimension; ++j) {
                exp6 += X[i][j][k] * U[i][j];
            }
        }
        Model.add(exp6 <= T);
    }

    //constraint 7
    for (int k = 0; k < R; ++k) {
        for (int j = 0; j < Dimension; ++j) {
            for (int i = 0; i < Dimension - 1; ++i) {
                IloExpr exp7(env);
                exp7 = t[i][k] + U[i][j] - (1 - X[i][j][k]) * G - t[j][k];
                Model.add(exp7 <= 0);
            }
        }
    }

    //constraint 8
    for (int k = 0; k < R; ++k) {
        for (int i = 0; i < Dimension; ++i) {
            IloExpr exp8(env);
            exp8 = t[i][k];
            Model.add(exp8 >= 0);
        }
    }
#pragma endregion

    IloCplex cplex(Model);
    cplex.setOut(env.getNullStream()); // pour eviter les commentaire el zeydin (ysimplifi w ye7seb)
    if (!cplex.solve()) {
        env.error() << "Failed for some reason" << endl;
        throw(-1);
    }

    double obj = cplex.getObjValue();

    cout << "\n\n\t The objective value is: " << obj << endl;
    
    for (int i = 0; i < Dimension; ++i) {
        for (int j = 0; j < Dimension; ++j) {
            for (int k = 0; k < R; ++k) {
                double Xval = cplex.getValue(X[i][j][k]);
                cout << "\t\t\t X[" << i + 1 << "][" << j + 1 << "][" << k + 1 << "] = " << Xval << endl;
            }
        }
    }

    for (int i = 0; i < Dimension; ++i) {
        for (int k = 0; k < R; ++k) {
            double tval = cplex.getValue(t[i][k]);
            cout << "\t\t\t t[" << i + 1 << "][" << k + 1 << "] = " << tval << endl;
        }
    }

    return 0;
}
