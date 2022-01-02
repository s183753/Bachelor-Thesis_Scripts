/*----------------------------------------------------------------------------------- 
    Victor Lorentzen, s183753
    BSc. Physics and Nanotechnology, DTU
    Bachelor's Thesis: Ray Tracing in Tokamak
-----------------------------------------------------------------------------------*/
#include <iostream>
#include <math.h>
#include <vector>

// Physical constants
double m_e = 9.1093837015e-31;                  // Mass of electron in [kg]
double q = 1.602176634e-19;                     // Charge of electron in [C]
double c = 299792458;                           // Speed of light in [m/s]
double epsilon_0 = 8.8541878128e-12;            // Vacuum permitivity [F/m]
double omega = 2*M_PI*140e9;                    // Angular frequency [1/s] or [Hz]

// Variables
int N = 1;                                      // Number of iterations
double delta_t = 0.1,k0 = 1/sqrt(2),            // Preset for parameters    
L_x =5, L_y =5, L_x_min = 2, L_y_min = -L_y;    
int runs = 1,dim = 3;
std::vector<int> mode = {1};
double theta = M_PI/2;                          // Angle between 
double accuracy = 0.01;                         // Numerical accuracy
double delta[3] = {accuracy,accuracy,accuracy}; // Step size for numerical approximations
double B_0 = 1.6;                                 // Magnetic field strength at start of tokamak
//double n_0 = 2.43e19;                           // Electron density constant [1/m^3]
double n_end = 1e20;                          // Electron density peak [1/m^3]
//double n_end = 5e20;
double n_0 = 2e20;                              // TESTING EBW
double T_end = 10*q;                            // Temperature peak [keV]
std::vector<std::vector<double>> r_init,k_init; // Vectors of initial conditions
std::vector<double> B_vector(dim);              // Magnetic Field Vector

//--------------------------------------------------------------------------------------------------
// GRIDS OF ELECTRON DENSITY, MAGNETIC FIELD AND TEMPERATURE
// Grids symmetric following cylindrical angle, grid coords will be R and Z

double mur = 5, muz = 5, sigmar = 2, sigmaz = 2;
double rs = 0,zs = 0;
double muPeakr = 8;
double muPeakz = 5;
double sigmaPeakr = 0.5, sigmaPeakz = 0.5;

// ELECTRON DENSITY GRID
const int N_nGrid = 100;
double nGrid[N_nGrid][N_nGrid];
void make_grid_n_e(){
    

    for (int r = 0; r < N_nGrid; r++){
        for (int z = 0; z < N_nGrid; z++){
            rs = (double)r / 10;
            zs = (double)z / 10;
            
            //nGrid[r][z] = 0.5*n_0 + 0.5*n_0*r/100; // For testing EBW
 
            nGrid[r][z] = n_end*exp((double)-1/2*((rs-mur)*(rs-mur)/(sigmar*sigmar) + (zs-muz)*(zs-muz)/(sigmaz*sigmaz))); // FOR GAUSSIAN

            //nGrid[r][z] += n_end*exp((double)-1/2*((rs-muPeakr)*(rs-muPeakr)/(sigmaPeakr*sigmaPeakr) + (zs-muPeakz)*(zs-muPeakz)/(sigmaPeakz*sigmaPeakz)));   // FOR GAUSSIAN WITH ADD PEAK
            if (nGrid[r][z] < 0) nGrid[r][z] = 0;
        }
    }
}

// MAGNETIC FIELD GRID
const int N_BGrid = 100;
double BGrid[N_BGrid][N_BGrid];
void make_grid_B(){
    double rad = 0;
    //rad = (N_BGrid)/(L_x - L_x_min);
    rad = (L_x - L_x_min)/N_BGrid;
    int spot = 50;
    B_0 = 5*(L_x*spot - spot*L_x_min + L_x_min*N_BGrid)/N_BGrid;
    for (int r = 0; r < N_BGrid; r++){
        for (int z = 0; z < N_BGrid; z++){
            
            BGrid[r][z] = B_0/(r*rad + L_x_min);                  // Magnetic field as 1/R 

            //BGrid[r][z] = 0.004*B_end*r;                  // LINEAR B-Field used for testing
            //BGrid[r][z] = 3;                              // CONSTANT B-field used for testing
            if (nGrid[r][z] < 0) nGrid[r][z] = 0;
        }
    }
}

// TEMPERATURE GRID
const int N_TGrid = 100;
double TGrid[N_TGrid][N_TGrid];
void make_grid_T(){
    double mur = 5, muz = 5, sigmar = 2, sigmaz = 2;
    double rs = 0,zs = 0;
    for (int r = 0; r < N_TGrid; r++){
        for (int z = 0; z < N_TGrid; z++){
            rs = (double)r / 10;
            zs = (double)z / 10;
            //TGrid[r][z] = T_end;                      // CONSTANT TEMPERATURE USED FOR TESTING
            //TGrid[r][z] = T_end + T_end*r;            // LINEAR TEMPERATURE USED FOR TESTING
            TGrid[r][z] = T_end*exp((double)-1/2*((rs-mur)*(rs-mur)/(sigmar*sigmar) + (zs-muz)*(zs-muz)/(sigmaz*sigmaz)));
            //TGrid[r][z] += T_end*exp((double)-1/2*((rs-muPeakr)*(rs-muPeakr)/(sigmaPeakr*sigmaPeakr) + (zs-muPeakz)*(zs-muPeakz)/(sigmaPeakz*sigmaPeakz)));
            if (TGrid[r][z] < 0) TGrid[r][z] = 0;
        }
    }
}

// FINDING GRID COORDINATES FROM REAL SPACE POSITION VECTOR
// USING LINEAR TRANSLATION
std::vector<int> grid_coords(std::vector<double> r, int gridID, double interpolation[]){
    double slopex = (N_nGrid)/(L_x - L_x_min);                  // slope of radial transformation
    double interceptx = (N_nGrid) - slopex*(L_x + 0.01*L_x);    // intercept of radial transformation

    double slopey = (N_nGrid)/(L_y - L_y_min);                  // slope of z transformation
    double intercepty = (N_nGrid) - slopey*(L_y + 0.01*L_y);    // intercept of z transformation
    std::vector<int> output(2);
    
    // radius from x and z components of r-vector
    output[0] = floor(slopex*sqrt(r[0]*r[0] + r[2]*r[2]) + interceptx);     
    interpolation[0] = slopex*sqrt(r[0]*r[0] + r[2]*r[2]) + interceptx - output[0];

    // if out of bounds of grids, coords become edge of grid
    if (output[0] < 0) {
        output[0] = 0;
        interpolation[0] = 0;
    }
    if (output[0] > N_nGrid - 2) {
        output[0] = N_nGrid - 2;
        interpolation[0] = 0;
    }
    
    // y of r-vector becomes Z
    output[1] = floor(slopey*r[1] + intercepty);
    interpolation[1] = slopey*r[1] + intercepty - output[1];

    // if out of bounds of grids, coords become edge of grid
    if (output[1] < 0) {
        output[1] = 0;
        interpolation[1] = 0;
    }
    if (output[1] > N_nGrid - 2) {
        output[1] = N_nGrid - 2;
        interpolation[1] = 0;
    }
    
    return output;
}

// Length of a vector, doenst have to be k
double k_value(std::vector<double> k_vector){
    double k = 0;
    for (int i = 0; i < int(k_vector.size()); i++){
        k += k_vector[i]*k_vector[i];
    }
    k = sqrt(k);
    return k;
}

// k*k vector dot product with self
double k_square(std::vector<double> k_vector){
    double k = 0;
    for (int i = 0; i < int(k_vector.size()); i++){
        k += k_vector[i]*k_vector[i];
    }
    return k;
}

// corresponding electron density from real space coordinates
double n_e(std::vector<double> r){
    double interpolation[2];
    std::vector<int> coords = grid_coords(r,1,interpolation);

    return nGrid[coords[0]][coords[1]] + interpolation[0]*(nGrid[coords[0]+1][coords[1]] - nGrid[coords[0]][coords[1]]) + interpolation[1]*(nGrid[coords[0]][coords[1]+1] - nGrid[coords[0]][coords[1]]);
}

// numerical gradient of electron density
double grad_n_e(std::vector<double> r, int id){

    std::vector<double> corrP = {r[0], r[1], r[2]};
    std::vector<double> corrM = {r[0], r[1], r[2]};
    corrP[id] += delta[id];
    corrM[id] -= delta[id];

    return (n_e(corrP) - n_e(corrM))/(2*delta[id]);      
}

// corresponding magnetic field from real space coordinates
double B(std::vector<double> r){
    double interpolation[2];
    std::vector<int> coords = grid_coords(r,0,interpolation);

    return BGrid[coords[0]][coords[1]] + interpolation[0]*(BGrid[coords[0] + 1][coords[1]] - BGrid[coords[0]][coords[1]]) + interpolation[1]*(BGrid[coords[0]][coords[1]+1] - BGrid[coords[0]][coords[1]]);
}

// numerical gradient of magnetic field 
double grad_B(std::vector<double> r, int id){
    
    std::vector<double> corrP = {r[0], r[1], r[2]};
    std::vector<double> corrM = {r[0], r[1], r[2]};
    corrP[id] += delta[id];
    corrM[id] -= delta[id];
    
    return (B(corrP) - B(corrM))/(2*delta[id]);      
}

// corresponding electron density from real space coordinates
double T_e(std::vector<double> r){
    double interpolation[2];
    std::vector<int> coords = grid_coords(r,0,interpolation);

    return TGrid[coords[0]][coords[1]] + interpolation[0]*(TGrid[coords[0] + 1][coords[1]] - TGrid[coords[0]][coords[1]]) + interpolation[1]*(TGrid[coords[0]][coords[1]+1] - TGrid[coords[0]][coords[1]]);
}

// numerical gradient of temperature
double grad_T_e(std::vector<double> r, int id){

    std::vector<double> corrP = {r[0], r[1], r[2]};
    std::vector<double> corrM = {r[0], r[1], r[2]};
    corrP[id] += delta[id];
    corrM[id] -= delta[id];

    return (T_e(corrP) - T_e(corrM))/(2*delta[id]);      
}

// Function returning 1, 0 or -1 from the sign of the inputed value
int sign(double value){
    if (value > 0) return 1;
    else if (value < 0) return -1;
    else return 0;
}

std::vector<double> crossProduct(std::vector<double> v_A, std::vector<double> v_B) {
    std::vector<double> output;
    output.push_back(v_A[1] * v_B[2] - v_A[2] * v_B[1]);
    output.push_back(v_A[3] * v_B[0] - v_A[0] * v_B[2]);
    output.push_back(v_A[0] * v_B[1] - v_A[1] * v_B[0]);
    return output;
}

double dotProduct(std::vector<double> v_A, std::vector<double> v_B){
    double product = 0;
    for (int i = 0; i < int(v_A.size()); i++){
        product += v_A[i]*v_B[i];
    }
    return product;
}

double omegaP(std::vector<double> r){           // Plasma frequency squared
    return (n_e(r)*q*q)/(m_e*epsilon_0);
}


// -------- NOTATION STEMS FROM ELECTROMAGNETIC WAVES FOR THERMONUCLEAR FUSION ------------
// Functions for easier notation throughout

double Gamma(double X, double Y, double angle){
    return sqrt(pow(Y,4)*pow(sin(angle),4) + 4*Y*Y*(1-X)*(1-X)*cos(angle)*cos(angle));
}

double A(double X){
    return 2*X*(1-X);
}

double B_m(double X, double Y, double angle, int m){
    return 2*(1-X) - Y*Y*sin(angle)*sin(angle) + m*Gamma(X,Y,angle);
}

double G(double X, double Y, double angle, int m){
    return A(X)/B_m(X,Y,angle,m);
}

double dD_dk(std::vector<double> k_vector){
    return -2*k_value(k_vector)*c*c/(omega*omega);
}

double dG_dX_more(double X, double Y, double angle, int m){
    double B_this = B_m(X,Y,angle,m);
    double A_this = A(X);
    return (2*(1-2*X)*B_this + 2*A_this + m*A_this*4*Y*Y*(1-X)*cos(angle)*cos(angle)/Gamma(X,Y,angle))/(B_this*B_this);
}

double dD_dX(double X, double Y, double angle, int m){
    return -dG_dX_more(X,Y,angle,m);
}

double dG_dYY_more(double X, double Y, double angle, int m){
    return A(X)*(sin(angle)*sin(angle) - m*(Y*Y*pow(sin(angle),4) + 2*(1-X)*(1-X)*cos(angle)*cos(angle))/Gamma(X,Y,angle))/(pow(B_m(X,Y,angle,m),2));
}

double dD_dYY(double X, double Y, double angle, int m){
    return -dG_dYY_more(X,Y,angle,m);
}

double dD_domega(std::vector<double> k_vector, double X, double Y, double angle, int m){
    return 2*k_square(k_vector)*c*c/pow(omega,3) + 2*X*dG_dX_more(X,Y,angle,m)/omega + 2*Y*Y*dG_dYY_more(X,Y,angle,m)/omega;
}

double dD_dtheta(double X, double Y, double angle, int m){
    return - G(X,Y,angle,m)*2*sin(angle)*cos(angle)*Y*Y/B_m(X,Y,angle,m) * (1-m*(Y*Y*pow(sin(angle),2)-2*pow(1-X,2)/Gamma(X,Y,angle)));
}

std::vector<double> grad_X(std::vector<double> r, double X){
    std::vector<double> output;
    for (int i = 0; i < dim; i++){
        output.push_back(X*grad_n_e(r,i)/n_e(r));
    }
    return output;
}

std::vector<double> grad_YY(std::vector<double> r, double Y){
    std::vector<double> output;
    for (int i = 0; i < dim; i++){
        output.push_back(2*Y*Y*grad_B(r,i)/B(r));
    }
    return output;
}

std::vector<double> grad_theta(std::vector<double> r, std::vector<double> k, double angle, std::vector<double> B_vector){
    std::vector<double> output;
    for (int i = 0; i < dim; i++){
        output.push_back(-1/(k_value(k)*sin(angle))*(k[i]*grad_B(r,i)/B(r)));
    }
    return output;
}

// -------------- DISPERSION RELATION FOR EBW ---------------------------------------------

double X(std::vector<double> r){
    return n_e(r)*q*q/(m_e*epsilon_0*omega*omega);
}

double Y(std::vector<double> r){
    return q*B(r)/(m_e*omega);
}

double D(double omega_p_2, double omega_c){
    return omega_c/omega * omega_p_2/(omega*omega - omega_c*omega_c);
}

double S(double omega_p_2, double omega_c){
    return 1 - omega_p_2/(omega*omega - omega_c*omega_c);
}

double P(double X){
    return 1 - X;
}

double l_Te_2(std::vector<double> r, double omega_p_2, double omega_c){
    return 3*omega_p_2*omega_c*omega_c/((4*omega_c*omega_c - omega*omega)*(omega*omega-omega_c*omega_c))*(m_e*T_e(r)/(q*q*B(r)*B(r)));
}

double dP_domega( std::vector<double> r){
    return 2*omegaP(r)/pow(omega,3);
}

double dd_domega_B (std::vector<double> r,double omega_p_2, double omega_c){
    return q*q*q*B(r)*n_e(r)*(q*q*B(r)*B(r) - 3*omega*omega*m_e*m_e)/(omega*omega*epsilon_0*pow(q*q*B(r)*B(r)-omega*omega*m_e*m_e,2));
}

double dD2_domega (std::vector<double> r){
    return -2*pow(q,6)*B(r)*B(r)*n_e(r)*n_e(r)*(q*q*B(r)*B(r)-3*omega*omega*m_e*m_e)/(pow(omega,3)*pow(epsilon_0,2)*pow(q*q*B(r)*B(r)-omega*omega*m_e*m_e,3));
}

double dS_domega_B (double omega_p_2,double omega_c){
    return 2*omega_p_2*omega/pow(omega*omega - omega_c*omega_c,2);
}

double dS2_domega (std::vector<double> r){
    
    return 2*S(omegaP(r),q*B(r)/m_e)*dS_domega_B(omegaP(r),q*B(r)/m_e);
}

double dl_domega (std::vector<double> r){
    return -30*T_e(r)*q*q*(q*q*B(r)*B(r) - 2*omega*omega*m_e*m_e/5)*m_e*m_e*m_e*m_e*omega*n_e(r)/((epsilon_0*pow(4*q*q*B(r)*B(r)-omega*omega*m_e*m_e,2))*pow(q*q*B(r)*B(r)-omega*omega*m_e*m_e,2));
}

double dD_domega_B(std::vector<double> r,std::vector<double> k,double omega_p_2,double omega_c){
    return dl_domega(r)*pow(sin(theta),4)*k_square(k)*k_square(k) + dS_domega_B(omegaP(r),q*B(r)/m_e)*pow(sin(theta),2)*k_square(k) + dP_domega(r)*pow(cos(theta),2)*k_square(k) - omega*omega/(c*c)*(dS2_domega(r) - dD2_domega(r)) - 2*omega/(c*c)*(S(omegaP(r),q*B(r)/m_e) - D(omegaP(r),q*B(r)/m_e));
}

// WITH ONLY DOUBLE NOT VECTOR NOTATION ---------------------------------------------------------------------------------

double D(std::vector<double> r){
    return Y(r)*X(r)/(1-pow(Y(r),2));
}

double S(std::vector<double> r){
    return 1 - X(r)/(1 - pow(Y(r),2));
}

double P(std::vector<double> r){
    return 1 - X(r);
}

double grad_X(std::vector<double> r, int id){
    return X(r)*grad_n_e(r,id)/n_e(r);
}

double grad_Y(std::vector<double> r, int id){
    return Y(r)*grad_B(r,id)/B(r);
}

double grad_YY(std::vector<double> r, int id){
    return 2*pow(Y(r),2)*grad_B(r,id)/B(r);
}

double grad_Y4(std::vector<double> r, int id){
    return 4*pow(Y(r),4)*grad_B(r,id)/B(r);
}

double grad_l(std::vector<double> r, int id){
    return 3/(omega*omega*m_e)*((grad_X(r,id)*T_e(r) + X(r)*grad_T_e(r,id))/((4*pow(Y(r),2) - 1)*(1-pow(Y(r),2))) - (5*grad_YY(r,id)-4*grad_Y4(r,id))*T_e(r)*X(r)/pow(5*pow(Y(r),2)-4*pow(Y(r),4)-1,2) );
}

double grad_P(std::vector<double> r, int id){
    return -grad_X(r,id);
}

double grad_S(std::vector<double> r, int id){
    return -(grad_X(r,id)/(1-pow(Y(r),2)) + X(r)*(grad_YY(r,id))/pow((1-pow(Y(r),2)),2));
}

double grad_S2(std::vector<double> r, int id){
    return 2*S(r)*grad_S(r,id);
}

double grad_D(std::vector<double> r, int id){
    return -grad_S(r,id)*Y(r)+X(r)/(1-pow(Y(r),2))*grad_Y(r,id);
}

double grad_D2(std::vector<double> r, int id){
    return 2*D(r)*grad_D(r,id);
}


//----------------------------------- INITIAL CONDITIONS FOR k -------------------------------------------------------------------------------------------

double k_BW(std::vector<double> sp){        // SOLUTION FOR WAVE NUMBER X-EBW USING "-"
    double wp2,wc,s,d,lte,p;
    wp2 = omegaP(sp);
    wc = q*B(sp)/m_e;
    s=S(wp2,wc);
    lte = l_Te_2(sp,wp2,wc);
    d = D(wp2,wc);
    p = P(wp2/(omega*omega));
    double a,b,co;

    a = lte*pow(sin(theta),4);
    b = s*pow(sin(theta),2) + p*pow(cos(theta),2);
    co = -omega*omega/(c*c)*(s*s - d*d);
    

    return sqrt(b*(-1 - sqrt(1 - 4*a*co/(b*b)))/(2*a));
}

double k_BX(std::vector<double> sp){        // SOLUTION FOR WAVE NUMBER X-EBW USING "+"
    double wp2,wc,s,d,lte,p;
    wp2 = omegaP(sp);
    wc = q*B(sp)/m_e;
    s=S(wp2,wc);
    lte = l_Te_2(sp,wp2,wc);
    d = D(wp2,wc);
    p = P(wp2/(omega*omega));
    double a,b,co;

    a = lte*pow(sin(theta),4);
    b = s*pow(sin(theta),2) + p*pow(cos(theta),2);
    co = -omega*omega/(c*c)*(s*s - d*d);
    
    return sqrt(b*(-1 + sqrt(1 - 4*a*co/(b*b)))/(2*a));
}

double k_O(std::vector<double> sp){         // SOLUTION FOR WAVE NUMBER O-mode 
    return sqrt((omega*omega - n_e(sp)*q*q/(m_e*epsilon_0))/(c*c));
}

double k_X(std::vector<double> sp){         //SOLUTION FOR WAVE NUMBER X-mode Appleton-Hartree
    double X = omegaP(sp)/(omega*omega);
    double Y = q*B(sp)/m_e/(omega);
    return sqrt(omega*omega/(c*c)*(1-2*X*(1-X)/(2*(1-X) - Y*Y*pow(sin(theta),2)-Gamma(X,Y,theta))));
}

// ----------------------------------------------------------------------------------------

void Print_grids(){
    for (int i = 0; i < N_nGrid; i++){
        for (int j = 0; j < N_nGrid; j++){
            std::cout << nGrid[i][j] << std::endl;
            std::cout << BGrid[i][j] << std::endl;
            std::cout << TGrid[i][j] << std::endl;
        }
    }
}

int main(){
    r_init.push_back({0,0});
    k_init.push_back({k0,k0});
    mode.pop_back();

    /* Input parameters from Python "ptrace":
        N,       delta_t,   L_x_min,L_x,        L_y_min,L_y,    runs,            dim,           theta
    Iterations,  Ray Steps, radius dimensions,  Z- dimensions   number of waves, dim = 3: 3D,   Angle to B-field
    */
    std::cin >> N >> delta_t >> L_x_min >> L_x >> L_y_min >> L_y >> runs >> dim >> theta;
    
    int tempss;
    // mode: 1 -> O-mode, 0 -> EBW, -1 -> X-mode
    for (int i = 0; i < runs; i++){
        std::cin >> tempss;
        mode.push_back(tempss);
    }

    
    // INITIAL VECTORS EMPTY
    r_init.pop_back();
    k_init.pop_back();

    // INITIALIZING GRIDS
    make_grid_B();
    make_grid_n_e();
    make_grid_T();

    // GRIDS TO PYTHON
    Print_grids();

    // DATA STORAGE VECTORS
    std::vector<std::vector<std::vector<double>>> r(runs);     // position vector 
    std::vector<std::vector<std::vector<double>>> k(runs);     // wave vector
    std::vector<std::vector<double>> omega_p_2(runs);          // Plasma frequency squared values stored for plotting
    std::vector<std::vector<double>> omega_c(runs);            // Cyclotron frequency for changing B-field


    double input[dim*runs];
    for (int i = 0; i < dim*runs; i++){
        std::cin >> input[i];
    }
    std::vector<double> temp;
    for (int i = 0; i < runs; i++){ 
        for (int j = 0; j < dim; j++){
            temp.push_back(input[j+i*dim]);
        }
        r_init.push_back(temp);
        temp.clear();
    }

    // CALCULATING INITIAL WAVE NUMBER FOR EACH MODE
    double k0[runs]; 
    for (int i = 0; i < runs; i++){
        if(mode[i] == 0){                                       // EBW
            k0[i] = k_BW(r_init[i]);
        }
        else if (mode[i] == 1){                                 // O-mode
            k0[i] = sqrt((omega*omega - n_e(r_init[i])*q*q/(m_e*epsilon_0))/(c*c));
        }
        else if (mode[i] == -1){                                // X-mode
            k0[i] = k_BX(r_init[i]);
        }
        std::cout << k0[i] << std::endl;
    }

    for (int i = 0; i < dim*runs; i++){
        std::cin >> input[i];
    }

    for (int i = 0; i < runs; i++){
        for (int j = 0; j < dim; j++){
            temp.push_back(input[j+i*dim]);
        }
        k_init.push_back(temp);
        temp.clear();
    }

    for (int i = 0; i < runs; i++){     
        r[i].push_back(r_init[i]);               // Initial of r-vector
        k[i].push_back(k_init[i]);               // Initial of k-vector
    }
    // END OF INITIALIZATION
    
    // START OF RAY TRACING
    double k_x,k_x2;                            // x-coordinate of wave vector and change in x-coordinate       
    double k_y,k_y2;                            // y-coordinate of wave vector and change in y-coordinate
    double k_z,k_z2 = 0;                        // z-coordinate  
    double r_x,r_x2;                            // x-coordinate of position vector and change in x-coordinate  
    double r_y,r_y2;                            // y-coordinate of position vector and change in y-coordinate
    double r_z,r_z2 = 0;                        // z-coordinate
    double len_r = 0;                           // Length of r-vector
    double dDdomega;
    double X, Y;
    std::vector<double> kcrossB(dim);
    double lenCross;
    std::vector<double> e_k(dim),e_theta(dim);
    int m;
    std::vector<double> temp1(dim), temp2(dim), temp3(dim);
    double denom;
    std::vector<std::vector<double>> k_theoX(runs);
    std::vector<std::vector<double>> k_theoB(runs);
    std::vector<std::vector<double>> k_theoO(runs);

    

    for (int i = 0; i < N; i++){            // Calculating using k^(n+1) = k^n + Delta_k*delta_t, same for r
        for (int j = 0; j < runs; j++){
            omega_p_2[j].push_back(omegaP(r[j].back()));        // Plasma frequency squared omega_p^2
            omega_c[j].push_back((q*B(r[j].back()))/m_e);       // Cyclotron frequency omega_c
            
            X = omega_p_2[j][i]/(omega*omega);
            Y = omega_c[j][i]/omega;
            
            m = mode[j];
            std::cout << m << std::endl;
            k_x  = (k[j][i][0]);                                // 1. coordinate of k^n
            k_y  = (k[j][i][1]);                                // 2. coordinate of k^n
            k_z  = (k[j][i][2]);
            r_x  = (r[j][i][0]);                                // 1. coordinate of r^n
            r_y  = (r[j][i][1]);                                // 2. coordinate of r^n
            r_z  = (r[j][i][2]);
            
            len_r = sqrt(r_x*r_x + r_z*r_z);
            
            B_vector[0] = -r_z/len_r*B(r[j][i]);
            B_vector[1] = 0;
            B_vector[2] = r_x/len_r*B(r[j][i]);

            theta = acos(dotProduct(B_vector,k[j].back())/(k_value(B_vector)*k_value(k[j].back())));

            //if (m == 0){
            if (m == 0 || m == -1){
                dDdomega = dD_domega_B(r[j].back(),k[j].back(),omega_p_2[j].back(),omega_c[j].back());
                
                r_x2 = l_Te_2(r[j].back(),omega_p_2[j].back(),omega_c[j].back())*pow(sin(theta),4)*4*k_square(k[j].back())*k_x + S(omega_p_2[j].back(),omega_c[j].back())*2*sin(theta)*sin(theta)*k_x + P(X)*cos(theta)*cos(theta)*2*k_x;
                r_y2 = l_Te_2(r[j].back(),omega_p_2[j].back(),omega_c[j].back())*pow(sin(theta),4)*4*k_square(k[j].back())*k_y + S(omega_p_2[j].back(),omega_c[j].back())*2*sin(theta)*sin(theta)*k_y + P(X)*cos(theta)*cos(theta)*2*k_y;
                r_z2 = l_Te_2(r[j].back(),omega_p_2[j].back(),omega_c[j].back())*pow(sin(theta),4)*4*k_square(k[j].back())*k_z + S(omega_p_2[j].back(),omega_c[j].back())*2*sin(theta)*sin(theta)*k_z + P(X)*cos(theta)*cos(theta)*2*k_z;
                temp1[0] = r_x2;
                temp1[1] = r_y2;
                temp1[2] = r_z2;
                denom = k_value(temp1);
                
                r_x2 = -sign(dDdomega)*r_x2/denom;
                r_y2 = -sign(dDdomega)*r_y2/denom;
                r_z2 = -sign(dDdomega)*r_z2/denom;

                k_x2 = sign(dDdomega)*(grad_l(r[j].back(),0)*pow(sin(theta),4)*pow(k_square(k[j].back()),2) + grad_S(r[j].back(),0)*k_square(k[j].back())*pow(sin(theta),2) + grad_P(r[j].back(),0)*pow(cos(theta),2)*k_square(k[j].back()) - pow(omega/c,2)*(grad_S2(r[j].back(),0) - grad_D2(r[j].back(),0)))/denom;
                k_y2 = sign(dDdomega)*(grad_l(r[j].back(),1)*pow(sin(theta),4)*pow(k_square(k[j].back()),2) + grad_S(r[j].back(),1)*k_square(k[j].back())*pow(sin(theta),2) + grad_P(r[j].back(),1)*pow(cos(theta),2)*k_square(k[j].back()) - pow(omega/c,2)*(grad_S2(r[j].back(),1) - grad_D2(r[j].back(),1)))/denom;
                k_z2 = sign(dDdomega)*(grad_l(r[j].back(),2)*pow(sin(theta),4)*pow(k_square(k[j].back()),2) + grad_S(r[j].back(),2)*k_square(k[j].back())*pow(sin(theta),2) + grad_P(r[j].back(),2)*pow(cos(theta),2)*k_square(k[j].back()) - pow(omega/c,2)*(grad_S2(r[j].back(),2) - grad_D2(r[j].back(),2)))/denom;
            }
            else {
                kcrossB = crossProduct(k[j].back(),crossProduct(k[j].back(),B_vector));
                lenCross = k_value(kcrossB);

                for (int a = 0; a < dim; a++){
                    e_k[a] = k[j][i][a]/k_value(k[j][i]);
                    e_theta[a] = kcrossB[a]/lenCross;
                }
                denom = sqrt(dD_dk(k[j][i])*dD_dk(k[j][i]) + pow(dD_dtheta(X,Y,theta,m)/k_value(k[j][i]),2));

                r_x2 = -sign(dD_domega(k[j][i],X,Y,theta,m))*(dD_dk(k[j][i])*e_k[0] + dD_dtheta(X,Y,theta,m)*e_theta[0]/k_value(k[j][i]))/denom;
                r_y2 = -sign(dD_domega(k[j][i],X,Y,theta,m))*(dD_dk(k[j][i])*e_k[1] + dD_dtheta(X,Y,theta,m)*e_theta[1]/k_value(k[j][i]))/denom;
                r_z2 = -sign(dD_domega(k[j][i],X,Y,theta,m))*(dD_dk(k[j][i])*e_k[2] + dD_dtheta(X,Y,theta,m)*e_theta[2]/k_value(k[j][i]))/denom;

                k_x2 = sign(dD_domega(k[j][i],X,Y,theta,m))*(dD_dX(X,Y,theta,m)*grad_X(r[j][i],X)[0] + dD_dYY(X,Y,theta,m)*grad_YY(r[j][i],Y)[0] + dD_dtheta(X,Y,theta,m)*grad_theta(r[j][i],k[j][i],theta,B_vector)[0])/denom;
                k_y2 = sign(dD_domega(k[j][i],X,Y,theta,m))*(dD_dX(X,Y,theta,m)*grad_X(r[j][i],X)[1] + dD_dYY(X,Y,theta,m)*grad_YY(r[j][i],Y)[1] + dD_dtheta(X,Y,theta,m)*grad_theta(r[j][i],k[j][i],theta,B_vector)[1])/denom;
                k_z2 = sign(dD_domega(k[j][i],X,Y,theta,m))*(dD_dX(X,Y,theta,m)*grad_X(r[j][i],X)[2] + dD_dYY(X,Y,theta,m)*grad_YY(r[j][i],Y)[2] + dD_dtheta(X,Y,theta,m)*grad_theta(r[j][i],k[j][i],theta,B_vector)[2])/denom;
                
                
                /*
                denom = sqrt(dD_dk(k[j][i])*dD_dk(k[j][i]) + pow(dD_dtheta(X(r[j].back()),Y(r[j].back()),theta,m)/k_value(k[j][i]),2));

                r_x2 = -sign(dD_domega(k[j][i],X(r[j].back()),Y(r[j].back()),theta,m))*(dD_dk(k[j][i])*e_k[0] + dD_dtheta(X(r[j].back()),Y(r[j].back()),theta,m)*e_theta[0]/k_value(k[j][i]))/denom;
                r_y2 = -sign(dD_domega(k[j][i],X(r[j].back()),Y(r[j].back()),theta,m))*(dD_dk(k[j][i])*e_k[1] + dD_dtheta(X(r[j].back()),Y(r[j].back()),theta,m)*e_theta[1]/k_value(k[j][i]))/denom;
                r_z2 = -sign(dD_domega(k[j][i],X(r[j].back()),Y(r[j].back()),theta,m))*(dD_dk(k[j][i])*e_k[2] + dD_dtheta(X(r[j].back()),Y(r[j].back()),theta,m)*e_theta[2]/k_value(k[j][i]))/denom;

                k_x2 = sign(dD_domega(k[j][i],X(r[j].back()),Y(r[j].back()),theta,m))*(dD_dX(X(r[j].back()),Y(r[j].back()),theta,m)*grad_X(r[j][i],X(r[j].back()))[0] + dD_dYY(X(r[j].back()),Y(r[j].back()),theta,m)*grad_YY(r[j][i],Y(r[j].back()))[0] + dD_dtheta(X(r[j].back()),Y(r[j].back()),theta,m)*grad_theta(r[j][i],k[j][i],theta,B_vector)[0])/denom;
                k_y2 = sign(dD_domega(k[j][i],X(r[j].back()),Y(r[j].back()),theta,m))*(dD_dX(X(r[j].back()),Y(r[j].back()),theta,m)*grad_X(r[j][i],X(r[j].back()))[1] + dD_dYY(X(r[j].back()),Y(r[j].back()),theta,m)*grad_YY(r[j][i],Y(r[j].back()))[1] + dD_dtheta(X(r[j].back()),Y(r[j].back()),theta,m)*grad_theta(r[j][i],k[j][i],theta,B_vector)[1])/denom;
                k_z2 = sign(dD_domega(k[j][i],X(r[j].back()),Y(r[j].back()),theta,m))*(dD_dX(X(r[j].back()),Y(r[j].back()),theta,m)*grad_X(r[j][i],X(r[j].back()))[2] + dD_dYY(X(r[j].back()),Y(r[j].back()),theta,m)*grad_YY(r[j][i],Y(r[j].back()))[2] + dD_dtheta(X(r[j].back()),Y(r[j].back()),theta,m)*grad_theta(r[j][i],k[j][i],theta,B_vector)[2])/denom;
                */
            }
            k_theoB[j].push_back(k_BW(r[j].back()));
            k_theoX[j].push_back(k_BX(r[j].back()));
            k_theoO[j].push_back(k_O(r[j].back()));
            /*std::cout << r[j].back()[0] << std::endl;
            std::cout << r[j].back()[1] << std::endl;
            std::cout << r[j].back()[2] << std::endl;
            std::cout << k[j].back()[0] << std::endl;
            std::cout << k[j].back()[1] << std::endl;
            std::cout << k[j].back()[2] << std::endl;*/
            std::cout << mode[j] << std::endl;
            std::cout << theta << std::endl;
            //std::cout << n_e(r[j].back()) << std::endl;

            k[j].push_back({k_x + k_x2*delta_t, k_y + k_y2*delta_t, k_z + k_z2*delta_t});      // Storing current wave vector (k) WITH Z-coordinate
            r[j].push_back({r_x + r_x2*delta_t, r_y + r_y2*delta_t, r_z + r_z2*delta_t});      // Storing current position vector (r)
        }
    }

    
    double closestECR = 100, closestUH = 100, closestLH = 100;
    int n = 1000;
    std::vector<double> RUH(n);
    std::vector<double> RECR(n);
    std::vector<double> RLH(n);
    double testingECR,testingUH,testingLH;
    std::vector<double> spacing = {L_x_min,L_y_min,0};
    for (int Z = 0; Z < n; Z++){
        closestECR = 100; closestUH = 100; closestLH = 100;
        for (int R = 0; R < n; R++){
            testingUH = fabs(omegaP(spacing)/(omega*omega) + q*q*pow(B(spacing),2)/(m_e*m_e*omega*omega) - 1);
            if (testingUH < closestUH && !(testingUH != testingUH)) {
                closestUH = testingUH;
                RUH[Z] = spacing[0];
            }
            testingECR = fabs(q*B(spacing)/(m_e*omega) - 1);
            if (testingECR < closestECR && !(testingECR != testingECR)) {
                closestECR = testingECR;
                RECR[Z] = spacing[0];
            }
            testingLH = fabs(omegaP(spacing)/(omega*omega) + q*B(spacing)/(m_e*omega) - 1);
            if (testingLH < closestLH && !(testingLH != testingLH)) {
                closestLH = testingLH;
                RLH[Z] = spacing[0];
            }

            spacing[0] += (L_x-L_x_min)/n;
        }
        spacing[1] += (L_y-L_y_min)/n;
        spacing[0] = L_x_min;
    }
    /*
    std::vector<double> BS(n), XS(n),OS(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){

        }
    }
    */

    // PRINING TO PYTHON

    for (int i = 0; i < n; i++){
        std::cout << RECR[i] << std::endl;
        std::cout << RUH[i] << std::endl;
        std::cout << RLH[i] << std::endl;
    }


    // Printing cyclotron- and plasma frequencies
    for (int i = 0; i < N; i++){
        for (int j = 0; j < runs; j++){
            //std::cout << "Frequencies:" << std::endl;
            std::cout << omega_c[j][i] << std::endl;
            std::cout.flush();
            std::cout << omega_p_2[j][i] << std::endl;
            std::cout.flush();
        }
    }

    // Printing k
    for (int i = 0; i < N; i++){
        for (int run = 0; run < runs; run++){
            for (int j = 0; j < dim; j++){
                //std::cout << run << " ";
                std::cout << k[run][i][j] << std::endl;
                std::cout.flush();
            }
        }
    }

    // Printing r
    for (int i = 0; i < N; i++){
        for (int run = 0; run < runs; run++){
            for (int j = 0; j < dim; j++){
                //std::cout << "r " << dim << " at " << run << " ";
                std::cout << r[run][i][j] << std::endl;
                std::cout.flush();
            }
        }
    }

    for (int i = 0; i < N; i++){
        for (int run = 0; run < runs; run++){
            //std::cout << run << " ";
            std::cout << k_theoX[run][i] << std::endl;
            std::cout << k_theoB[run][i] << std::endl;
            std::cout << k_theoO[run][i] << std::endl;
            std::cout.flush();

        }
    }

    return 0;
}