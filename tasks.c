
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define precision 0.000001 //how precise we want it to be
#define TRUE 0             //True or false
#define FALSE 1
#define GAMMA 1.4            //default gamma value
#define M_PI 3.141592653589  //PI
#define INIT 0               //initialise the values
#define EXCLUDE_FIRST_LINE 1 //ignore the first line of data
#define STEPS 10             //number of steps between points

//stores a given point for interpolation
typedef struct
{
    double x;
    double y;
} point_t;

typedef long unsigned int lu_t; //typedef for number of bytes

//safely allocate memory
void *safe_malloc(lu_t num_bytes)
{
    void *ptr = malloc(num_bytes);
    if (ptr == NULL)
    {
        printf("ERROR: malloc(%lu)\n", num_bytes);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

//safely open file
FILE *safe_fopen(const char *path, const char *mode)
{
    FILE *fp = fopen(path, mode);
    if (fp == NULL)
    {
        perror("file open error");
        exit(EXIT_FAILURE);
    }
    return fp;
}

//the function
double f(double x0, double M, double t, double gamma)
{
    double func;
    //converting to radians
    double b = x0 * (M_PI / 180);
    double theta = t * (M_PI / 180);
    double frac = ((M * M * (sin(b) * sin(b)) - 1) / (M * M * (gamma + cos(2 * b)) + 2));
    func = 2 * (1 / tan(b)) * frac - tan(theta);
    return func;
}

//an approximation to the derivatie
double f_derivative(double x0, double M, double t, double gamma)
{
    double x1, x2, y1, y2;
    //using our degree of precision
    //making an approximation using first principles
    x1 = x0 - precision;
    x2 = x0 + precision;
    y1 = f(x1, M, t, gamma);
    y2 = f(x2, M, t, gamma);

    return (y2 - y1) / (x2 - x1);
}
//runs the main newton-rhapson
double shockwave_value(double t, double M, double x0, int maxitr, int *theta_max, int *theta_max_found, double gamma)
{

    int success_flag = FALSE; //success flag for if method was successful
    // int what_i;               //which value of i does it converge for

    double x1; // x1 values

    for (int i = 1; i <= maxitr; i++)
    {
        x1 = x0 - f(x0, M, t, gamma) / f_derivative(x0, M, t, gamma); //x0 - delta
                                                                      // printf(" At Iteration no. %3d, x = %.6lf\n", i, x1);
        if ((fabs(f(x0, M, t, gamma) / f_derivative(x0, M, t, gamma)) < precision) && (success_flag == FALSE))
        {

            success_flag = TRUE;
            // what_i = i;
        }
        //set x0 to the next value of x
        x0 = x1;
    }

    //check to see if this is true, print what the root is after however many iterations
    if (success_flag == TRUE)
    {
        //printf("Root = %.6lf, after %d iterations\n", x1, what_i);
        return x1;
    }
    //if false then return then print the error statement
    else
    {
        //printf("Either not enough iterations to reach precision or does not converge\n");

        if (*theta_max_found == FALSE)
        {
            *theta_max = (int)t;
            *theta_max_found = TRUE;
        }
        // printf("%d\n", theta_max);

        return FALSE;
    }
}

void shockwave_calculation(double m, double beta_u, double beta_l, int question, double gamma)
{

    int theta;
    double *b_l; //lower
    double *b_u; //upper
    int theta_max = INIT;
    int theta_max_found = FALSE;
    //write to this file
    FILE *fp = safe_fopen("out_shock.csv", "a");

    //find the value of theta_max
    for (int i = 1; i <= 90; i++)
    {
        theta = (double)i;
        shockwave_value(theta, m, 89, 50, &theta_max, &theta_max_found, gamma);
    }
    //allocating memory

    b_u = (double *)safe_malloc((theta_max + 1) * sizeof(double));
    b_l = (double *)safe_malloc((theta_max + 1) * sizeof(double));

    //find the roots
    for (int i = 0; i < theta_max; i++)
    {
        //if intial value
        if (i == 0)
        {
            b_l[0] = asin(1 / m) * 180 / M_PI; //convert to degrees and the boundary conditions
            b_u[0] = 90.0;
        }
        else
        { //evaluate at these points for the upper and lower angles
            theta = (double)i;
            b_l[i] = shockwave_value(theta, m, (double)beta_l, 50, &theta_max, &theta_max_found, gamma);
            b_u[i] = shockwave_value(theta, m, (double)beta_u, 50, &theta_max, &theta_max_found, gamma);
        }
        //error checking has to realistically be in this bound
        if ((b_u[i] > 0) && (b_u[i] < 90))
        { //check which question it is
            if (question == 3)
            {
                fprintf(fp, "%.6lf,%d,%.6lf,%.6lf\n", m, i, b_l[i], b_u[i]);
            }
        }
    }
    fclose(fp);
    theta_max_found = FALSE;
    //free memory
    free(b_u);
    free(b_l);
}

void shockwave(const char *q2_file)
{
    int N_count = INIT; //number of lines of data
    int i = INIT;       //counting variables
    int j = INIT;
    char ch;           //character stream
    double *input_val; //list of values for M for part c

    double M, theta, beta_l, beta_u, gamma; //parameters for initial question
    double M2 = INIT;                       //paramaters for part B m-value

    FILE *fp = safe_fopen(q2_file, "r");

    //checking for the number of lines
    for (ch = getc(fp); ch != EOF; ch = getc(fp))
    {
        if (ch == '\n')
        { // Increment entries_count if this character is newline
            //checking for the number of lines in the data
            N_count++;
        }
    }

    input_val = (double *)safe_malloc((N_count) * sizeof(double));

    fclose(fp);
    fp = safe_fopen(q2_file, "r");
    //read in the values of the variables from the CSV

    for (ch = getc(fp); ch != EOF; ch = getc(fp))
    {
        if (ch == '\n')
        {
            if (i == 0)
            {
                fscanf(fp, "%lf,%lf,%lf,%lf,%lf", &M, &theta, &beta_l, &beta_u, &gamma);
                //j++;
            }
            //if line number is 2
            else if (i == 2)
            {
                fscanf(fp, "%lf", &M2);
            }
            //line number more than 4
            else if (i >= 4)
            {

                fscanf(fp, "%lf", &input_val[j]);
                j++;
            }
            //add a new line to number of lines
            i++;
        }
    }

    fclose(fp);

    //part a
    //modify the values of initial guesses
    beta_u = 89;
    beta_l = 1;
    shockwave_value(theta, M, beta_u, 50, 0, 0, gamma);
    shockwave_value(theta, M, beta_l, 50, 0, 0, gamma);

    //part b
    shockwave_calculation(M2, beta_u, beta_l, 2, gamma);
    //part c
    fp = safe_fopen("out_shock.csv", "w");
    fprintf(fp, "M,theta,beta_ lower,beta_upper\n");
    fclose(fp); //close file stream

    for (int count = 0; count < j - 1; count++)
    {
        //printf("%lf\n", input_val[count]);
        shockwave_calculation(input_val[count], beta_u, beta_l, 3, gamma);
    }

    free(input_val); //free memory
}

void thomas(double *ai, double *bi, double *ci, double *qi,
            double *a_star, double *q_star, double *xi, int N_count)
{

    for (int i = 0; i < N_count; i++)
    {
        //for the first value do this (formulas given/proved in report)
        if (i == 0)
        {
            a_star[0] = ai[0];
            q_star[0] = qi[0];
        }
        else
        {
            //for all other values do this
            a_star[i] = ai[i] - ((ci[i] * bi[i - 1]) / a_star[i - 1]);
            q_star[i] = qi[i] - ((ci[i] * q_star[i - 1]) / a_star[i - 1]);
        }
    }
    //perform back substution
    for (int i = N_count - 1; i > 0; i--)
    {
        if (i == (N_count - 1))
        {
            xi[i] = q_star[i] / a_star[i];
        }
        else
        {
            xi[i] = (q_star[i] - (bi[i] * xi[i + 1])) / a_star[i];
        }
    }
}

void linalgbsys(const char *q4_file)
{
    //declare arrays to be read in
    double *ai;
    double *bi;
    double *ci;
    double *qi;

    double *a_star;
    double *q_star;
    double *xi;

    int N_count = INIT; // the size of the A matrix will be this^2
    int i = INIT;       //counting variable
    char c;

    FILE *fp = safe_fopen(q4_file, "r");

    //checking for the number of lines
    for (c = getc(fp); c != EOF; c = getc(fp))
    {
        if (c == '\n')
        { // Increment entries_count if this character is newline
            //checking for the number of lines in the data
            N_count++;
        }
    }
    N_count = N_count - EXCLUDE_FIRST_LINE;
    fclose(fp);
    //allocate memory

    ai = (double *)safe_malloc(N_count * sizeof(double));
    bi = (double *)safe_malloc(N_count * sizeof(double));
    ci = (double *)safe_malloc(N_count * sizeof(double));
    qi = (double *)safe_malloc(N_count * sizeof(double));
    a_star = (double *)safe_malloc(N_count * sizeof(double));
    q_star = (double *)safe_malloc(N_count * sizeof(double));
    xi = (double *)safe_malloc(N_count * sizeof(double));

    //initialise the solution matrix for printing such that no conditional jumps
    for (int i = 0; i < N_count; i++)
    {
        xi[i] = INIT;
        ai[i] = INIT;
        bi[i] = INIT;
        ci[i] = INIT;
        qi[i] = INIT;
        a_star[i] = INIT;
        q_star[i] = INIT;
    }

    fp = safe_fopen(q4_file, "r");
    //read in the numbers from the file
    for (c = getc(fp); c != EOF; c = getc(fp))
    {
        if (c == '\n')
        {
            fscanf(fp, "%lf,%lf,%lf,%lf", &ai[i], &bi[i], &ci[i], &qi[i]);
            i++;
        }
    }
    fclose(fp);

    //execute function
    thomas(ai, bi, ci, qi,
           a_star, q_star, xi, N_count);

    fp = safe_fopen("out_linalsys.csv", "w");

    for (int i = 0; i < N_count; i++)
    {
        //please use the print function to assess if valgrind gives error and doesn't print to file
        //printf("%.6lf\n", xi[i]);
        fprintf(fp, "%.6lf\n", xi[i]);
    }

    fclose(fp);

    //free memory
    free(ai);
    free(bi);
    free(ci);
    free(qi);
    free(a_star);
    free(q_star);
    free(xi);
}

//calculate the value of S at a given point
double calculate_S(double ai, double bi, double ci, double di, double x, double xi)
{
    double s;
    s = ai + (bi * (x - xi)) + ci * ((x - xi) * (x - xi)) + di * ((x - xi) * (x - xi) * (x - xi));
    return s;
}

void interp(const char *q5_file, const double xo)
{
    //the structure to hold the point
    point_t *point;
    //interpolation coefficients
    double *h;
    double *a;
    double *b;
    double *c;
    double *d;

    //for solving tri_diagonal
    double *ci;
    double *bi;
    double *ai;
    double *qi;
    double *a_star;
    double *q_star;
    double *xi;
    int num_elem;            //number of elements
    double current_S = INIT; //hold the current value of S

    int N_count = INIT; //number of lines counted
    int i = INIT;       //counting variable
    char ch;            //file steam chracter

    double guess = xo; //the value to be found

    FILE *fp = safe_fopen(q5_file, "r");

    //checking for the number of lines
    for (ch = getc(fp); ch != EOF; ch = getc(fp))
    {
        if (ch == '\n')
        { // Increment entries_count if this character is newline
            //checking for the number of lines in the data
            N_count++;
        }
    }

    fclose(fp);
    //allocate memory

    point = (point_t *)safe_malloc(N_count * sizeof(point_t));
    h = (double *)safe_malloc(N_count * sizeof(double));
    a = (double *)safe_malloc(N_count * sizeof(double));
    b = (double *)safe_malloc(N_count * sizeof(double));
    c = (double *)safe_malloc(N_count * sizeof(double));
    d = (double *)safe_malloc(N_count * sizeof(double));

    ai = (double *)safe_malloc(N_count * sizeof(double));
    bi = (double *)safe_malloc(N_count * sizeof(double));
    ci = (double *)safe_malloc(N_count * sizeof(double));
    qi = (double *)safe_malloc(N_count * sizeof(double));
    a_star = (double *)safe_malloc(N_count * sizeof(double));
    q_star = (double *)safe_malloc(N_count * sizeof(double));
    xi = (double *)safe_malloc(N_count * sizeof(double));

    //initialise these values
    for (int i = 0; i < N_count; i++)
    {
        a[i] = INIT;
        b[i] = INIT;
        c[i] = INIT;
        d[i] = INIT;
        h[i] = INIT;
        point[i].x = INIT;
        point[i].y = INIT;

        ai[i] = INIT;
        bi[i] = INIT;
        ci[i] = INIT;
        qi[i] = INIT;
        a_star[i] = INIT;
        q_star[i] = INIT;
        xi[i] = INIT;
    }

    fp = safe_fopen("in_interp.csv", "r");

    //read in the values from the file

    for (ch = getc(fp); ch != EOF; ch = getc(fp))
    {
        if (ch == '\n')
        {
            fscanf(fp, "%lf,%lf", &point[i].x, &point[i].y);
            i++;
        }
    }
    fclose(fp);
    num_elem = N_count - 1;

    //cubic spline will have 4n conditions
    //calculate the h value
    for (i = 0; i < num_elem; i++)
    {

        a[i] = point[i].y;

        h[i] = point[i + 1].x - point[i].x;
    }
    //set last values of the matrix
    c[0] = 0.0;
    c[num_elem] = 0.0;

    for (int i = 0; i < num_elem; i++)
    {
        if (i == 0)
        {
            qi[0] = 0.0;
        }
        else
        {
            qi[i] = 3.0 * (a[i + 1] - a[i]) / h[i] - 3.0 * (a[i] - a[i - 1]) / h[i - 1];
        }
    }
    qi[num_elem] = 0;
    //making the big A matrix, elements of the tridiagonal to solve
    //initial
    ai[0] = (double)1;
    bi[0] = 0;
    ci[0] = 0;
    //middle values
    for (i = 1; i < num_elem; i++)
    {
        ci[i] = h[i - 1];
        bi[i] = h[i];
        ai[i] = 2.0 * (h[i - 1] + (h[i]));
    }
    //final values
    ai[num_elem] = (double)1;
    bi[num_elem] = h[num_elem];
    ci[num_elem] = h[num_elem];

    //THOMAS
    for (int i = 0; i < num_elem; i++)
    {
        //for the first value do this (formulas given/proved in report)
        if (i == 0)
        {
            a_star[0] = ai[0];
            q_star[0] = qi[0];
        }
        else
        {
            //for all other values do this
            a_star[i] = ai[i] - ((ci[i] * bi[i - 1]) / a_star[i - 1]);
            q_star[i] = qi[i] - ((ci[i] * q_star[i - 1]) / a_star[i - 1]);
        }
    }
    //perform back substution
    for (int i = num_elem - 1; i > 0; i--)
    {
        if (i == (num_elem - 1))
        {
            xi[i] = q_star[i] / a_star[i];
        }
        else
        {
            xi[i] = (q_star[i] - (bi[i] * xi[i + 1])) / a_star[i];
        }
    }
    //THOMAS END

    //assign values of solved column to c
    for (int i = 0; i < num_elem; i++)
    {
        c[i] = xi[i];
    }

    //find values of b and d
    for (int i = 0; i < num_elem; i++)
    {
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        b[i] = (1.0 / h[i]) * (a[i + 1] - a[i]) - (h[i] / 3.0) * (2.0 * c[i] + c[i + 1]);
    }
    //output file stream
    fp = safe_fopen("out_interp.csv", "w+");
    //print to file
    fprintf(fp, "xo,f(xo)\n");
    for (int i = 0; i < num_elem - 1; i++)
    {

        //check for the point between the first and next point
        if ((point[i].x > guess) && (point[i + 1].x < guess))
        {
            current_S = calculate_S(a[i], b[i], c[i], d[i], guess, point[i].x);
            //please use the print function to assess if valgrind gives error and doesn't print to file
            //printf("%.6lf,%.6lf\n", guess, current_S);
            fprintf(fp, "%.6lf,%.6lf\n", guess, current_S);
        }
        //check for the value between the first point and the one before it
        else if ((point[i + 1].x > guess) && (point[i].x < guess))
        {
            current_S = calculate_S(a[i], b[i], c[i], d[i], guess, point[i].x);
            //please use the print function to assess if valgrind gives error and doesn't print to file
            //printf("%.6lf,%.6lf\n", guess, current_S);
            fprintf(fp, "%.6lf,%.6lf\n", guess, current_S);
        }
    }

    fclose(fp);

    //Outputting values with stepsize for the graph with stepsize 10
    // printf("xo,f(xo)\n");
    // for (int i = 0; i < num_elem - 1; i++)
    // {
    //     double stepsize = (point[i + 1].x - point[i].x) / STEPS;

    //     for (int j = 0; j < STEPS; j++)
    //     {
    //         printf("%lf,%lf\n", point[i].x + j * stepsize, calculate_S(a[i], b[i], c[i], d[i], point[i].x + j * stepsize, point[i].x));
    //     }
    // }

    //REMEMBER TO FREE
    //free memory
    free(point);
    free(h);
    free(a);
    free(b);
    free(c);
    free(d);
    free(ai);
    free(bi);
    free(ci);
    free(qi);
    free(a_star);
    free(q_star);
    free(xi);
}

//first order upwiind scheme
double rhs_scheme(int i, double *f_n, double delta_x, int Nx, double c, char scheme)
{
    //upwinding scheme
    if (scheme == 'u')
    { //boundary value  at i == 0
        double partial_x;
        if (i == 0)
        {
            partial_x = (f_n[i + 1] - f_n[i]) / delta_x;
        }
        else
        {
            partial_x = (f_n[i] - f_n[i - 1]) / delta_x;
        }
        return -1.0 * c * partial_x;
    }
    //central scheme
    else if (scheme == 'c')
    {
        double partial_x;
        //value at boundary i == 0
        if (i == 0)
        {
            partial_x = (f_n[i + 1] - f_n[i]) / delta_x;
        }
        //value at boundary i == Nx
        else if (i == Nx)
        {
            partial_x = (f_n[i] - f_n[i - 1]) / delta_x;
        }
        else
        {

            partial_x = (f_n[i + 1] - f_n[i - 1]) / (2 * delta_x);
        }

        return -1.0 * c * partial_x;
    }
    else
    {
        return FALSE;
    }
}

//RK2 routine
void rk_2(double delta_x, double delta_t, double *f_n_1, double *f_n_half, double *f_n, double c, int Nx, char scheme, double time_point, int out_iter)
{
    double timestep_count = INIT; //number of time levels to go through
    int count_iter = INIT;
    while (count_iter < out_iter)
    {
        if (timestep_count <= time_point)
        {
            for (int i = 0; i <= Nx; i++)
            {
                f_n_half[i] = 0;
                f_n_1[i] = 0;
            }
            //assign values of f_n_half
            for (int i = 0; i <= Nx; i++)
            {
                f_n_half[i] = f_n[i] + (delta_t) * (rhs_scheme(i, f_n, delta_x, Nx, c, scheme));
            }
            //assign values of f_n_1
            for (int i = 0; i <= Nx; i++)
            {
                f_n_1[i] = f_n[i] + (0.5) * (delta_t) * ((rhs_scheme(i, f_n_half, delta_x, Nx, c, scheme) + (rhs_scheme(i, f_n, delta_x, Nx, c, scheme))));
            }
            //assign for the next point
            for (int i = 0; i <= Nx; i++)
            {
                f_n[i] = f_n_1[i];
            }
            //increment the count
            timestep_count = timestep_count + delta_t;
        }
        count_iter++;
    }
}

void waveeqn(const char *q6_file)
{
    //variables read in as parameters
    double c;
    int Nx;
    double CFL;
    int out_iter;

    char ch;        //character for file reading
    int i = INIT;   //counting array
    double delta_x; //delta_x
    double delta_t; //delta_t
    double *f_n;    //arrays for the function values
    double *f_n_half;
    double *f_n_1;
    double value_t = 0.2; //by default find 0.2

    char scheme; //type of scheme

    FILE *fp = safe_fopen(q6_file, "r");

    //read in the values from the file

    for (ch = getc(fp); ch != EOF; ch = getc(fp))
    {
        if (ch == '\n')
        {
            fscanf(fp, "%lf,%d,%lf,%d", &c, &Nx, &CFL, &out_iter);
            i++;
        }
    }
    fclose(fp);

    f_n = (double *)safe_malloc((Nx + 1) * sizeof(double));
    f_n_half = (double *)safe_malloc((Nx + 1) * sizeof(double));
    f_n_1 = (double *)safe_malloc((Nx + 1) * sizeof(double));

    delta_x = 1.0 / Nx; //assign delta x as 1/Nx

    //the value is equal to CFL which is 1 in this instance
    delta_t = (CFL * delta_x) / c;

    for (int i = 0; i <= Nx; i++)
    {
        f_n[i] = INIT;
        f_n_1[i] = INIT;
        f_n_half[i] = INIT;
    }
    //initial conditions
    for (int i = 0; i <= Nx; i++)
    {

        if (delta_x * (i) < 0.125)
        {
            f_n[i] = 0;
        }
        else if ((delta_x * (i) >= 0.125) && (delta_x * (i) <= 0.375))
        {
            f_n[i] = 0.5 * (1.0 - cos(8.0 * M_PI * (((i * delta_x)) - 0.125)));
        }
        else if ((delta_x * (i) < 0.375) && (delta_x * (i) <= 1.0))
        {
            f_n[i] = 0;
        }
    }

    fp = safe_fopen("out_waveeqn_1U.csv", "w");

    value_t = 0.2;
    //call upwinding scheme
    scheme = 'u';

    rk_2(delta_x, delta_t, f_n_1, f_n_half, f_n, c, Nx, scheme, value_t, out_iter);
    fprintf(fp, "x,f(x)\n");

    for (int i = 0; i <= Nx; i++)
    {
        fprintf(fp, "%lf,%lf\n", delta_x * i, f_n_1[i]);
    }

    fclose(fp);

    //reset all the matricies for the central scheme
    for (int i = 0; i <= Nx; i++)
    {
        f_n[i] = INIT;
        f_n_1[i] = INIT;
        f_n_half[i] = INIT;
    }
    for (int i = 0; i <= Nx; i++)
    {

        if (delta_x * (i) < 0.125)
        {
            f_n[i] = 0;
        }
        else if ((delta_x * (i) >= 0.125) && (delta_x * (i) <= 0.375))
        {
            f_n[i] = 0.5 * (1.0 - cos(8.0 * M_PI * (((i * delta_x)) - 0.125)));
        }
        else if ((delta_x * (i) < 0.375) && (delta_x * (i) <= 1.0))
        {
            f_n[i] = 0;
        }
    }
    fp = safe_fopen("out_waveeqn_2C.csv", "w");
    //call upwinding scheme
    scheme = 'c';

    rk_2(delta_x, delta_t, f_n_1, f_n_half, f_n, c, Nx, scheme, value_t, out_iter);
    fprintf(fp, "x,f(x)\n");

    for (int i = 0; i <= Nx; i++)
    {
        fprintf(fp, "%lf,%lf\n", delta_x * i, f_n_1[i]);
    }
    //reset and clear arrays
    for (int i = 0; i <= Nx; i++)
    {
        f_n[i] = INIT;
        f_n_1[i] = INIT;
        f_n_half[i] = INIT;
    }
    //exact solution plotter performs translation by c*value_t
    for (int i = 0; i <= Nx; i++)
    {
        if (delta_x * (i) < (0.125 + (c * value_t)))
        {
            f_n[i] = 0;
        }
        else if ((delta_x * (i) >= (0.125 + (c * value_t))) && (delta_x * (i) <= (0.375 + (c * value_t))))
        {
            f_n[i] = 0.5 * (1.0 - cos(8.0 * M_PI * ((((i * delta_x) - (c * value_t))) - 0.125)));
        }
        else if ((delta_x * (i) < (0.375 + (c * value_t))) && (delta_x * (i) <= (1.0)))
        {
            f_n[i] = 0;
        }
    }

    // printf("x,f(x)\n");

    // for (int i = 0; i <= Nx; i++)
    // {
    //     printf("%lf,%lf\n", delta_x * i, f_n[i]);
    // }

    fclose(fp);
    //free memory
    free(f_n);
    free(f_n_half);
    free(f_n_1);
}