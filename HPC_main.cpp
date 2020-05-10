//include all the libraries and functions needed.
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <string>
#include <stdlib.h>
#include "cblas.h"
#include <iomanip>
#include <boost/program_options.hpp>
#include <mpi.h>
//#include<cassert>
     using namespace std;
     namespace po = boost::program_options;

     #define F77NAME(x) x##_
     extern "C" {
         //define the Lapack subroutines used in order to solve the matrix-vector operations required
         void F77NAME(dpbsv) (const char& UPLO, const int& N, const int& KD, const int& NRHS, double* AB, const int& LDAB, double* B, const int& LDB, int& INFO);
         void F77NAME(daxpy) (const int& N, const double& DA, double* DX, const int& INCX, double* DY, const int& INCY);
         void F77NAME(dscal) (int N, double DA, double* DX, int INCX);
         void F77NAME(dtbmv) (const char& UPLO, const char& TRANS, const char& DIAG, const int& N, const int& K, double* A, const int& LDA, double* x, const int INCX);
         void F77NAME(dsbmv) (const char& UPLO, const int& N, const int& K, const double& alpha, double* AB, const int& LDA, double* x, const int& INCX, const double& beta, double* Y, const int& INCY);
         void F77NAME(dgemv) (const char& TRANS, const int& M, const int& N, const double& alpha, double* A, const int& lda, double* X, const int& INCX, const double& beta, double* Y, const int& INCY);
	 }

     //function that allows matrices to be printed - it is used in order to check the values of the matrices computed
     void PrintMatrix(int ncols, int row, double PM[]) {
         for (int i = 0; i < row; i++) {
            cout << "[";
             for (int j = 0; j < ncols; j++) {
                 cout << PM[row*j+i] << " ";
             }
             cout << "]" << endl;
         }
         cout << "\n" << endl;
      }
     //function to compute and assembly the global stiffnes matrix for the number of elements selected
     void AssemblyStiff (int ncols, double A, double E, double l, double I, double K[]) {
       double m2[6] = {(2*(A * E))/l, 0, (24 * E * I)/(l*l*l), 0, 0, (8 * E * I)/l};                                          //define the submatrix with the non-zero values of the upper triangle in the matrix
       double m1[9] = {-(A * E)/l, 0, 0, 0, -(12 * E * I)/(l*l*l), -(6 * E * I)/(l*l), 0, (6 * E * I)/(l*l), (2 * E * I)/l};  //define the submatrix with the non-zero values of the lower triangle in the matrix

       int It1 = 9; int It2 = 0;  //initialize the iterators for both submatrices
       for (int i = 0; i < ncols; i++) {
            for (int j = 0; j < 6; j++) {
                if (j >= 5 - (i%3)) {
                   if (i%3 == 0) {It2 = 0;}
                   K[6*i + j] = m2[It2];  //assign the correspondent values of in the stiffness matrix
                   ++It2;
                }
                else if (j >= 2 - (i%3)) {
                        if (i%3 == 0 && It1 == 9) {It1 = 0;}
                        K[6*i + j] = m1[It1];  //assign the corresponent values in the stiffness matrix
                        ++It1;
                }
                else {
                 K[6*i + j] = 0;  //fill with zeros the empty gaps in the matrix
                }
            }
        }
     }
    //function to assembly the global mass matrix for the number of elements selected - this matrix is used for the dynamic problem
     void AssemblyMass (int ncols, double A, double rho, double l, double alpham, double dt, double M[]) {
       double mm2[3] = {rho*A*l*(1 / (dt*dt)), rho*A*l*(1 / (dt*dt)), alpham*2.0*(l*l*l*rho*A) * (1 / (dt*dt))};  //define the submatrix where the real values are stores
       int It1 = 0;
       for (int i = 0; i < ncols; i++) {
            if (i%3 == 0) {It1 = 0;}
                M[i*6 + 5] = mm2[It1];  //assign the correspondent values to the positions in the matrix (only diagonal non-zero)
                It1++;
        }
     }
     //function to assembly the global mass matrix for the number of elements selected - this matrix is used for the parallelisation
     void AssembleMassMatrix (double* M, const int dof, double l, double mult) {
       double alpha = 1.0 / 24.0;
       for (int i = 0; i < dof; ++i) {
           if (i % 3 < 2) {
              M[i] = mult;
           } 
           else if (i % 3 == 2) {
                   M[i] = 2 * mult * alpha * l * l;
           }
       }
     }
     //function to compute the vector of forces taking into account the distributed and concentrated loads - used in the dynamic problem
     void ForceMatrix (int nx, double qx, double qy, double Fy, double l, double F[]) {
         double matrixf[3] = {qx*l, qy*l, 0.0};  //define subvector with the 3-value diagonal (since it repeats every 3 positions)
         int itr1 = 0;
         for (int i = 0; i < (nx-1)*3; i++) {
             if (i%3 == 0) {itr1 = 0;}
                F[i] = matrixf[itr1];
                itr1++;
             if (i == ((nx-1)*3)/2) {
                 F[i] += Fy;  //add the concentrated load in case we are in the central node
             }
          }
     }
     //function to compute the vector of forces taking into account the distributed and concentrated loads - used in the static problem
     void ForceVector (int ncols, double l, double qy, double Fy, double* F) {
        double fm[3] = {0.0, qy * l, 0.0}; //define subvector with the 3-value diagonal
        int it = -3;  //initialize the iterator
        for (int i = 0; i < ncols; i++) {
            if (i%3 == 0) {it += 3;}
               F[i] = fm[i-it];
        }
        F[(ncols-1)/2] += Fy;  //add the concentrated load in the case of the central node
      }
     //function to compute the vector of forces taking into account the distributed and concentrated loads - used un the parallelisation
     void ForceMatrixPar (const int dofl, double qy, double Fy, double l, double F[], int position) {
        //double matrixf[3] = {qx*l, qy*l, 0.0};
        for (int i = 0; i < dofl; ++i) {
            if (i%3 == 1) {
				F[i] = qy * l;
			}
            else {
               F[i] = 0;              
            }
        }
		F[position] += Fy;
     }
	 
    //begin the main function that initialize the code 
    int main(int argc, const char* argv[]) {  
      //define the input parameters to be selected in order to solve the problem
      po::options_description desc("Solve the cantilever problem for the values and conditions selected. ");
      desc.add_options()
      ("length", po::value<double>(),       "Length of the cantilever beam.")  //input the length of the beam
      ("nr_elements", po::value<int>(),     "Number of elements selected.")  //input the number of elements for the Finite Element solver
      ("area", po::value<double>(),         "Area of the cross-section.")  //input the value of the cross-sectional area desired
      ("mInertia", po::value<double>(),     "Moment of inertia.")  //input the Moment of Inertia of the cross-section been used
      ("yMod", po::value<double>(),         "Young's modulus of the material.")  //input the Young's Modulud of the material that forms the beam
      ("static",                            "Type of problem to solve, static.")  //select in case you want to run the static solver (task 1)
      ("dynamic",                           "Type of problem to solve, dynamic.")  //select in case you want to run the dynamic solver (tasks 2 & 3)
      ("explicit",                          "Type of finite difference method to use, explicit scheme.")  //select in case you want to run the explicit -dynamic- solver (task 2)
      ("implicit",                          "Type of finite difference method to use, implicit scheme.")  //select in case you want to run the implicit -dynamic- solver (task 3)
      ("parallel",                          "Use Parallel in order to compute the solution.")  //select for computing the problem in parallel mode
      ("time", po::value<double>(),         "Maximum time.")  //select the maximum time for solving the dynamic problem
      ("timeSS", po::value<double>(),       "Time at which steady-state starts.") 
      ("Nt", po::value<double>(),           "Timestep for the solver.")  //select the number of time-elements (timestep in the discretisation) for the dynamic solver
      ("density", po::value<double>(),      "Density of the material of the beam.")  //input the density of the material of the beam
      ("use-blas",                          "Use the BLAS implementation.")  //allows the code to run BLAS in order to do matrix-array calculations
      ("help",                              "Print help message.");  //help message to facilitate program use

      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);

      //if statement - prints the different inputs to be selected and information about them in case the user needs them
      if (vm.count("help")) {
         cout << "Solve the cantilever problem for the values and conditions selected." << endl;
         cout << desc << endl;
         return 0;
      }
      //convert the arguments into variables to be used through the problem
      double    L     = vm.count("length")       ? vm["length"]     .as<double>()  : 10.0;
      const int nx    = vm.count("nr_elements")  ? vm["nr_elements"].as<int>()     : 24;
      double    A     = vm.count("area")         ? vm["area"]       .as<double>()  : 0.012;
      double    I     = vm.count("mInertia")     ? vm["mInertia"]   .as<double>()  : 1.44e-05;
      double    E     = vm.count("yMod")         ? vm["yMod"]       .as<double>()  : 210e+09;
      double    T     = vm.count("time")         ? vm["time"]       .as<double>()  : 10;
      double    Nt    = vm.count("Nt")           ? vm["Nt"]         .as<double>()  : 100000;
	    double    T1    = vm.count("timeSS")       ? vm["timeSS"]     .as<double>()  : 5.0;
      double    rho   = vm.count("density")      ? vm["density"]    .as<double>()  : 7850;      

      //if statement - check and validate the number of inputs given
      if (argc < 7){
         cout << "Error: Not enough or invalid arguments, please try again\n";  //print error message
         exit(0);
      }
      else {
          for(int i = 1; i < argc; i++) {
             //cout << argv[i] << endl;
             if (argv[i] <= 0){
                 cout << "Error: Invalid arguments, please select positive numbers\n";  //print error message
                 exit(0);
             }
      }
     // if statement - check that the number of elements selected allows a the concentrated load to be placed at a node
     if (nx % 2 == 0) {
        cout << "Values validated" << endl;  //print message that informs that the input phase has been completed successfuly

        double l = L / nx;  //length of a single element
        const int nrows = 6;  //number of rows (number of degrees of freedom per element)
        const int nnodes = nx+1;  //number of nodes (number of elements + 1)
        const int dim2 = 3*nnodes;  //number of degrees of freedom for the whole beam, including boundary conditions
        const int ncols = dim2 - 6; //number of DoFs for the whole beam excluding boundary conditions
      
      //define the array with the locations of the nodes along the beam
      double nodeloc[nnodes-2];
      for (int i = 0; i < nnodes-2; i++) {
        nodeloc[i] = (i+1)*l;  //compute the lengthwise distance along the beam for that particular node
      }

       //declare the matrices and submatrices to be used in the global stiffness matrix assembly
       double K[ncols*6] = {};
       double* m1 = new double[9];
       double* m2 = new double[6];
       //compute the final global stiffness matrix using the function shown above
       AssemblyStiff (ncols, A, E, l, I, K);
       //PrintMatrix (ncols, 6, K);  //print-matrix function used to check values
   
       //if statement - check what type of problems would the user like to solve (static or dynamic)
       if (vm.count("static")) {

          double* F = new double[ncols];  //declare the vector of forces
          double qy = -1000.0;  //define the value of the distributed load  [N/m]
          double Fy = -1000.0;  //define the value of the concentrated load [N]
          double qx = 0.0;      //define the value of the x-direction distributed load [0]
          //compute the complete vector of force s for the number inputs selected
          ForceVector (ncols, l, qy, Fy, F);

          int info;  //declare variable that allows to see if the subroutine is working properly  
          double* u = new double[ncols];  //declare the vector of displacements (solution to the problem)
          cblas_dcopy(ncols, F, 1, u, 1);  //make a copy of the Force vector to the displacements vector 

          //compute the solution of the equation Ax = B for x, which is the vector of displacements
          F77NAME(dpbsv) ('U', ncols, 5, 1, K, 6, u, ncols, info);
          //if statement - display and error message if the subroutine is having trouble solvinf the Ax = B problem
          if (info) {
             cout << "Error solving, ID: " << info << endl;
          }
          else {

          //output the displacement values to a text file in order to plot the solution
          ofstream displData("task1_displ.txt"); //generate the .txt file
          for (int k = 0; k < ncols/3; k++) {
              displData << fixed << setprecision(10) << u[3*k + 1] << " " << nodeloc[k] << "\r\n";  //append the results for each node to the .txt file in columns with a separation
          }
          //output the input values given at the start to solve the problem
          ofstream displInput("task1_info.txt");  //generate the .txt file
          displInput << L << "\r\n";  //place the length value
          displInput << nx << "\r\n";  //place the number of elements value
          displInput << A << "\r\n";  //place the area value
          displInput << E << "\r\n";  //place the Young's modulus value
          displInput << I << "\r\n";  //place the moment of inertia value
          displInput << qy << "\r\n";  //place the distributed load value
          displInput << Fy << "\r\n";  //place the concentrated load value

          cout << "Calculations finished successfuly." << endl;
          }
          return 0;
        }  
        //if statement - initialize the solution to the dynamic problem
        else if (vm.count("dynamic")){
             //declare the assembled mass matrix and its inverse
             double M[6*ncols] = {};
             double invM[6*ncols] = {};

             double alpham = 0.04166667;  //define constant alpha to be used in the mass matrix
             //define the values of the loads
             double qx = 0;  //x-direction distributed load
             double Fy = -1000.0;  //y-direction concetrated load
             double qy = -1000.0;  //-y-direction concentrated load
  
             double dt = T / Nt;

             AssemblyMass (ncols, A, rho, l, alpham, dt, M);

    if (vm.count("explicit")) {
      if (vm.count("parallel")) {
        // N.B in the current configuration - there is little-to-no point to parallelising the program. This however provides a good basis to model the beam as a half, (i.e. in earlier tasks), and then split the halves into parallel constituent parts, being careful of boundary conditions.
        int err = 0;
        // Extend for parallelisation over for more processes?
    
        // Initialise the MPI environment.
        err = MPI_Init(NULL, NULL);
    
        if (err != MPI_SUCCESS) {
            cout << "Failed to initialise MPI environment." << endl;
            //return -1;
        }
    
        // Number of processes.
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
        // Rank of local process.
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
        ofstream displData4;
        if (world_rank == 0) {
    
           displData4.open("task4_displ.txt");
        
           if (!displData4.good()) {
              cout << "Error: Unable to write to output file: task4_displ.txt" << endl;
              return -1;
           }
        
        }
    
    // Define number of degrees of freedom in each partition.
    const int dof = ((nx / world_size) + 1) * 3;
    
    double K[6 * dof] = {};
    double M[dof] = {};
    double F[dof] = {};
    
    double u[dof] = {0.0};
    double ul[dof] = {0.0};
    
    double mult = rho * A * l / (dt * dt);
    int index = dof - 5;
    
    AssemblyStiff(dof, A, E, l, I, K);
    //AssemblyMass(dof, A, rho, l, alpham, dt, M);
    AssembleMassMatrix(M, dof, l, mult);
	
    double subMatrix[9] = {-1 * A * E / l, 0, 0,
                            0, -12 * E * I / (l * l * l), -6 * E * I / (l * l), 0,
                            6 * E * I / (l * l), 2 * E * I / l};
    
    for (int i = 0; i < dof; ++i) {
        K[6 * i + 5] +=  -2.0 * M[i];
    }   

    
    double t = 0.0;
    double p;
    
    double receive [3];
    
    //int tag = 0;
    
    while (t <= T) {
        
        if (t <= T1) {
            p = t / T1;
            ForceMatrixPar (dof, qy * p, Fy * p, l, F, index);
        } else {
            ForceMatrixPar (dof, qy, Fy, l, F, index);
        }

        F77NAME(dsbmv) ('U', dof, 5, -1.0, K, 6, u, 1, 1.0, F, 1);      
        
        if (world_rank == 0) {
            
            MPI_Ssend(&u[dof - 6], 3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(receive, 3, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        } else if (world_rank == 1) {
            
            MPI_Recv(receive, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Ssend(&u[dof - 6], 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
        }
        
        F77NAME(dgemv) ('N', 3, 3, -1, subMatrix, 3, receive, 1, 1, &F[dof - 3], 1);
        
        F77NAME(dsbmv) ('U', dof, 0, -1, M, 1, ul, 1, 1, F, 1);
        //cblas_dcopy (dof, u, 1, ul, 1);
        
        for (int i = 0; i < dof; ++i) {
            u[i] = F[i] / M[i];
        }

        cblas_dcopy (dof, u, 1, ul, 1);
        
        if (world_rank == 0) {
			
			     displData4 << fixed << setprecision(10) << u[index] << " " << t << "\r\n";
            
            /*fOutput2.precision(5);
            fOutput2.width(6);
            fOutput2 << fixed << t << " ";
            
            fOutput2.precision(10);
            fOutput2.width(11);
            fOutput2 << u[index] << "\n";*/

        }
        
        t += dt;
    }
    
    if (world_rank == 0) {
        displData4.close();
    }
       MPI_Finalize();
}
else {
     if (vm.count("timeSS")) {
	             double M[6*ncols] = {};
                 double invM[6*ncols] = {};

                 double alpham = 0.04166667;
    
                 double qx = 0;
                 double Fy = -1000.0;
                 double qy = -1000.0;

                 double dt = T / Nt;

                 AssemblyMass (ncols, A, rho, l, alpham, dt, M);
 
                          for (int i = 0; i < 6*ncols; i++) {
                             if (M[i] != 0) {
                             invM[i] = 1.0 / M[i];
                             //cout << invM[i] << '\n' << endl;
                             }
                             else {M[i] = 0;}
                          }
  
                          double ul[ncols] = {0.0};
                          double u[ncols] = {0.0};
                          double uu[ncols] = {};
                          cblas_dcopy (ncols, u, 1, uu, 1);
    
                          double F[ncols] = {};
                          double varFy, varqy;

                          F77NAME (daxpy) (6*ncols, -2.0, M, 1, K, 1);			  
                          //PrintMatrix(6*ncols, 6, K);

                          ofstream displData2("task2_displ.txt");
                          if (!displData2.good()) {
                             cout << "Error: unable to write to output file: task2_displ.txt" << endl;
                             return -1;
                          }

                          double t = 0.0;
                          while (t <= T) {
                               if (t <= T1) {
                                  varFy = (Fy/T1) * t;
                                  varqy = (qy/T1) * t;
                               }
                               else {
                                  varFy = Fy;
                                  varqy = qy;
                               }

                               ForceMatrix (nx, qx, varqy, varFy, l, F);

                               F77NAME (dsbmv) ('U', ncols, 5, -1.0, K, 6, u, 1, 1.0, F, 1);
                               //cblas_dcopy (ncols, term1, 1, rhs, 1);
                               F77NAME (dsbmv) ('U', ncols, 5, -1.0, M, 6, ul, 1, 1.0, F, 1);
                               //cblas_dcopy (ncols, F, 1, term1, 1); 
      
                               for (int i = 0; i < ncols; i++) {
                                   uu[i] = invM[6*i + 5] * F[i];
                               }
      
                               displData2 << fixed << setprecision(10) << uu[ncols/2] << " " << t << "\r\n";
                               //delete[]ul;
                               if (t > T1) {
                                uStor[t] = uu[ncols/2];
                               }

                               cblas_dcopy (ncols, u, 1, ul, 1);
                               cblas_dcopy (ncols, uu, 1, u, 1);
                               t += dt;
                          }
                          displData2.close();
                        }
                        else {

                          double M[6*ncols] = {};
                          double invM[6*ncols] = {};

                          double alpham = 0.04166667;
    
                          double qx = 0;
                          double Fy = -1000.0;
                          double qy = -1000.0;

                          double dt = T / Nt;

                          AssemblyMass (ncols, A, rho, l, alpham, dt, M);
 
                          for (int i = 0; i < 6*ncols; i++) {
                             if (M[i] != 0) {
                             invM[i] = 1.0 / M[i];
                             //cout << invM[i] << '\n' << endl;
                             }
                             else {M[i] = 0;}
                          }
  
                          double ul[ncols] = {0.0};
                          double u[ncols] = {0.0};
                          double uu[ncols] = {};
                          cblas_dcopy (ncols, u, 1, uu, 1);
    
                          double F[ncols] = {};
                          double varFy, varqy;

                          F77NAME (daxpy) (6*ncols, -2.0, M, 1, K, 1);        
                          //PrintMatrix(6*ncols, 6, K);

                          ofstream displData2f("task2f.txt");
                          if (!displData2f.good()) {
                             cout << "Error: unable to write to output file: task2f.txt" << endl;
                             return -1;
                          }

                          for (int T1 = 1; T1 < 10; ++i){
                          double t = 0.0;
                          while (t <= T) {
                               if (t <= T1) {
                                  varFy = (Fy/T1) * t;
                                  varqy = (qy/T1) * t;
                               }
                               else {
                                  varFy = Fy;
                                  varqy = qy;
                               }

                               ForceMatrix (nx, qx, varqy, varFy, l, F);

                               F77NAME (dsbmv) ('U', ncols, 5, -1.0, K, 6, u, 1, 1.0, F, 1);
                               //cblas_dcopy (ncols, term1, 1, rhs, 1);
                               F77NAME (dsbmv) ('U', ncols, 5, -1.0, M, 6, ul, 1, 1.0, F, 1);
                               //cblas_dcopy (ncols, F, 1, term1, 1); 
      
                               for (int i = 0; i < ncols; i++) {
                                   uu[i] = invM[6*i + 5] * F[i];
                               }
      
                               displData2f << fixed << setprecision(10) << uu[ncols/2] << " " << t << "\r\n";
                               //delete[]ul;
                              

                               cblas_dcopy (ncols, u, 1, ul, 1);
                               cblas_dcopy (ncols, uu, 1, u, 1);
                               t += dt;
                          }
                          displData2.close();
                        }




                        }
                      }

     else if (vm.count("implicit")) {
         double gamma = 0.5;
         double beta = 0.25;

         double udot[ncols] = {0.0};
         double udotu[ncols] = {0.0};
         double uddot[ncols] = {0.0};
         double uddotu[ncols] = {0.0};
         double ul[ncols] = {0.0};
         double u[ncols] = {0.0};
         double uu[ncols] = {0.0};
    
         double F[ncols] = {0.0};
         double tb[ncols] = {0.0};
         double uddotf[ncols] = {0.0};

         double varFy, varqy;

         F77NAME (daxpy) (6*ncols, 4.0, M, 1, K, 1);
         //PrintMatrix (ncols, 6, K);
         double Kt[6*ncols] = {};
         double T1 = 5.0;
         
         ofstream displData3("task3_displ.txt");
         if (!displData3.good()) {
            cout << "Error: unable to write to output file: task3_displ.txt" << endl;
            return -1;
         }

         double bdt2 = 1.0/(beta*dt*dt);
         double bdt = 1.0/(beta*dt);
         double b2 = 1.0/(2*beta);

         double t = 0.0;
         while (t <= T) {
          if (t <= T1) {
           varFy = (Fy/T1) * t;
           varqy = (qy/T1) * t;
          }
          else {
          varFy = Fy;
          varqy = qy;
          }
         
         ForceMatrix (nx, qx, varqy, varFy, l, F);
         //PrintMatrix(ncols, 1, F);

         /*cblas_dscal(ncols, 1 / (beta*dt), udot, 1);
         cblas_dscal(ncols, ((1/(2*beta))-1), uddot, 1);
         cblas_dscal(ncols, (1/(beta*dt*dt)), u, 1);*/

         for (int i = 0; i < ncols; i++) {
          tb[i] = bdt2*u[i] + bdt*udot[i] + (b2-1)*uddot[i];    
          //cout << tb[i] << " ";    
         }

         F77NAME (dsbmv) ('U', ncols, 5, (dt*dt), M, 6, tb, 1, 1.0, F, 1);
         cblas_dcopy(ncols, F, 1, uu, 1);
         int info;
         cblas_dcopy(ncols*6, K, 1, Kt, 1);
         F77NAME(dpbsv) ('U', ncols, 5, 1, Kt, 6, uu, ncols, info);

         displData3 << fixed << setprecision(10) << uu[ncols/2] << " " << t << "\r\n";
         //cblas_dscal(ncols, (1/(beta*dt*dt)), uu, 1);

         for (int i = 0; i < ncols; i++) {
         uddotu[i] = bdt2*uu[i] - bdt2*u[i] - bdt*udot[i] - (b2-1)*uddot[i];
         }

         /*cblas_dscal(ncols, (dt*(1-gamma)), uddotf, 1);
         cblas_dscal(ncols, dt*gamma, uddotu, 1);*/

         for (int i = 0; i < ncols; i++) {
         udotu[i] = udot[i] + dt*(1-gamma)*uddot[i] + dt*gamma*uddotu[i];
         }
         
         //cblas_dcopy (ncols, u, 1, ul, 1);
         cblas_dcopy (ncols, uu, 1, u, 1);
         //cblas_dcopy (ncols, udot, 1, udotl, 1);
         cblas_dcopy (ncols, udotu, 1, udot, 1);
         //cblas_dcopy (ncols, uddot, 1, uddotl, 1);
         cblas_dcopy (ncols, uddotu, 1, uddot, 1);
         //} 
         t += dt;
         }
         displData3.close();
         }


     else {cout << "No valid solving method selected." << endl;}
   }
   else {cout << "Invalid argument for the type of problem to solve." << endl;}
    }
 else { cout << "Select an even number of elements." << endl;}
      return 0;
  }
}
