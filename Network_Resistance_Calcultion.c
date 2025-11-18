#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

// limiting number for number of edges per single particle node (memory  allocation purpose)
#define Z 50

int main()
{
    double Resistance, Resistivity ; // network Resistance and Resistivity
    double Resistance_probe; // probe - particle resistance

    double V_supply = 1.0; // arbitrary Voltage supply value
    double I_probe_Right, I_probe_Left; // Currents flowing trough Left and Right probes
    int N = 0; // Total Number of Nodes
    int M = 0; // Total Number of Edges
    int M1,M2; // start Indexes of left probe-particles edges (M1) and  right probe-particles edges (M2)
    int kZ = 20;  // probe nodes can have kZ times more edges (probe contacts) than  particle nodes

    double SizeX = 0.0; // variables to determine margins of sample box
    double SizeY = 0.0;
    double SizeZ = 0.0;
    double O_X = 10000.0;
    double O_Y = 10000.0;
    double O_Z = 10000.0;

    double R_edge = 1.0; // particle-particle contact resistance

    int i, j, k, l;

    char s[20], s2[30] , s3[30];
    char line[256];
    int  source, target;
    double x_read, y_read, z_read, D_read;
    int node_read, ContactArea_read;

    char s1[20] = "exp_36_V3"; // filename label of graph-like structure derived from FIB-SEM experiment data

    strcpy(s2, s1);
    strcat(s2, ".txt");
    printf("%s\n", s2);

    FILE *f1;

    // reading graph structure from the file

    f1 = fopen(s2,"r");
    fgets(line, sizeof(line), f1);

    node_read = 0;

	double *x = NULL; // pointer to array of X coordinates of particles centers
    double *y = NULL; // pointer to array of  Y coordinates of particles centers
    double *z = NULL; //  pointer to array of  Z coordinates of particles centers
    double *radius = NULL;   //  pointer to array of  particles radii
    double *V = NULL; //   pointer to array of node potentials
    double *In = NULL; //   pointer to array of  net node currents, must approach to zero

    int *numn = NULL; // pointer to array of numbers of  nearest neighbors for each node
    int **nn = NULL;  //pointer to 2D array  of indexes of nearest neighbors for each node
    int **edges = NULL; //pointer to 2D array of indexes of connected edges for each node

    double *R = NULL; // pointer to array of edge resistances
    int *E_origin = NULL; //  pointer to array of edge ORIGIN node indexes
    int *E_end = NULL; // pointer to array of edge END node indexes
    double *Ie = NULL; //  pointer to array of edge currents

    double *t1 = NULL; // temporary pointers
    int    *t2 = NULL;
    int   **t3 = NULL;
    double *tR = NULL;
    double *tIe = NULL;

    int newN;

    // reading "Nodes" part of file

	while (node_read == 0)
        {
            fscanf(f1,"%lf\n",&x_read);
            fscanf(f1,"%lf\n",&y_read);
            fscanf(f1,"%lf\n",&z_read);
            fscanf(f1,"%lf\n",&D_read);
            fscanf(f1,"%s\n",&s);

            node_read = strcmp(s,"node");
            if (node_read != 0) break;

            newN = N + 1;

            t1 = realloc(x, newN * sizeof(double)); if (!t1) exit(1); x = t1;
            t1 = realloc(y, newN * sizeof(double)); if (!t1) exit(1); y = t1;
            t1 = realloc(z, newN * sizeof(double)); if (!t1) exit(1); z = t1;
            t1 = realloc(radius, newN * sizeof(double)); if (!t1) exit(1); radius = t1;
            t1 = realloc(V, newN * sizeof(double)); if (!t1) exit(1); V = t1;
            t1 = realloc(In, newN * sizeof(double)); if (!t1) exit(1); In = t1;

            t2 = realloc(numn, newN * sizeof(int)); if (!t2) exit(1); numn = t2;

            t3 = realloc(nn, newN * sizeof(*nn)); if (!t3) exit(1); nn = t3;
            nn[N] = malloc(Z * sizeof(int));
            if (!nn[N]) exit(1);

            t3 = realloc(edges, newN * sizeof(*edges)); if (!t3) exit(1); edges = t3;
            edges[N] = malloc(Z * sizeof(int));
            if (!edges[N]) exit(1);

            x[N] = x_read * 1e-1;
            y[N] = y_read * 1e-1;
            z[N] = z_read * 1e-1;
            radius[N] = 0.5 * D_read * 1e-1;

            if (x[N] > SizeX) SizeX = x[N];
            if (y[N] > SizeY) SizeY = y[N];
            if (z[N] > SizeZ) SizeZ = z[N];

            if (x[N] < O_X) O_X = x[N];
            if (y[N] < O_Y) O_Y = y[N];
            if (z[N] < O_Z) O_Z = z[N];

            numn[N] = 0;
            for (l = 0; l < Z; l++) nn[N][l] = 0;
            for (l = 0; l < Z; l++) edges[N][l] = 0;

            V[N] = 0.0;
            In[N] = 0.0;

            N = newN;
        }

        printf("Nodes N = %d\n", N);

        // reading "Edges" part of file

         while    ( fgets(line, sizeof(line), f1) != NULL)
        {
            sscanf(line,"%d\n", &source);
            fscanf(f1,"%d\n", &target);
            fscanf(f1,"%d\n",&ContactArea_read);
            fgets(line, sizeof(line), f1);

            i = source-1; j = target-1;

            nn[i][numn[i]] = j;
            edges[i][numn[i]] = M;
            E_origin = realloc(E_origin, (M + 1) * sizeof(int));
            if (!E_origin) exit(1);
            E_origin[M] = i;
            numn[i]++;

            nn[j][numn[j]] = i;
            edges[j][numn[j]] = M;
            E_end = realloc(E_end, (M + 1) * sizeof(int));
            if (!E_end) exit(1);
            E_end[M] = j;
            numn[j]++;

            tR = realloc(R, (M + 1) * sizeof(double));
            if (!tR) exit(1);
            R = tR;

            tIe = realloc(Ie, (M + 1) * sizeof(double));
            if (!tIe) exit(1);
            Ie = tIe;

            R[M] = R_edge;
            Ie[M] = 0.0;

            M++;
        }

    fclose(f1);

    printf("Edges M = %d\n", M);

    Resistance_probe = 0.00001; // probe - particle contact resistance

 // Probe node 1 on the Left  (X axis) --------------------------------------------------------------

    newN = N + 1;

    t1 = realloc(x, newN * sizeof(double)); if (!t1) exit(1); x = t1;
    t1 = realloc(y, newN * sizeof(double)); if (!t1) exit(1); y = t1;
    t1 = realloc(z, newN * sizeof(double)); if (!t1) exit(1); z = t1;
    t1 = realloc(radius, newN * sizeof(double)); if (!t1) exit(1); radius = t1;
    t1 = realloc(V, newN * sizeof(double)); if (!t1) exit(1); V = t1;
    t1 = realloc(In, newN * sizeof(double)); if (!t1) exit(1); In = t1;

    t2 = realloc(numn, newN * sizeof(int)); if (!t2) exit(1); numn = t2;

    t3 = realloc(nn, newN * sizeof(*nn)); if (!t3) exit(1); nn = t3;
    nn[N] = malloc( Z*kZ* sizeof(int));
    if (!nn[N]) exit(1);

    t3 = realloc(edges, newN * sizeof(*edges)); if (!t3) exit(1); edges = t3;
    edges[N] = malloc( Z*kZ*sizeof(int));
    if (!edges[N]) exit(1);

    x[N] = -30.0;
    y[N] = SizeY/2.;
    z[N] = SizeZ/2.;

    numn[N] = 0;
    for(l=0; l<Z*kZ; l++) nn[N][l] = 0;
    for(l=0; l<Z*kZ; l++) edges[N][l] = 0;

    V[N] = V_supply;
    In[N] = 0.0;

    N = newN;

 // --------------------------------------------

    i = N-1;

    M1 = M;

   // Probe node 1 connecting particles on the left face of sample
    for(k=0; k<N-1; k++)
        if ( x[k]-2.*radius[k]<0.0) // vicinity criteria
		    {
                j = k;

                nn[i][numn[i]] = j;
                edges[i][numn[i]] = M;
                E_origin = realloc(E_origin, (M + 1) * sizeof(int));
                if (!E_origin) exit(1);
                E_origin[M] = j;
                numn[i]++;

                nn[j][numn[j]] = i;
                edges[j][numn[j]] = M;
                E_end = realloc(E_end, (M + 1) * sizeof(int));
                if (!E_end) exit(1);
                E_end[M] = i;
                numn[j]++;

                tR = realloc(R, (M + 1) * sizeof(double));
                if (!tR) exit(1);
                R = tR;

                tIe = realloc(Ie, (M + 1) * sizeof(double));
                if (!tIe) exit(1);
                Ie = tIe;

				R[M] = Resistance_probe;
				Ie[M] = 0.0;

				M++;
		    }

	// Probe node 2 on the Right (X axis)

    newN = N + 1;

    t1 = realloc(x, newN * sizeof(double)); if (!t1) exit(1); x = t1;
    t1 = realloc(y, newN * sizeof(double)); if (!t1) exit(1); y = t1;
    t1 = realloc(z, newN * sizeof(double)); if (!t1) exit(1); z = t1;
    t1 = realloc(radius, newN * sizeof(double)); if (!t1) exit(1); radius = t1;
    t1 = realloc(V, newN * sizeof(double)); if (!t1) exit(1); V = t1;
    t1 = realloc(In, newN * sizeof(double)); if (!t1) exit(1); In = t1;

    t2 = realloc(numn, newN * sizeof(int)); if (!t2) exit(1); numn = t2;

    t3 = realloc(nn, newN * sizeof(*nn)); if (!t3) exit(1); nn = t3;
    nn[N] = malloc( Z*kZ* sizeof(int));
    if (!nn[N]) exit(1);

    t3 = realloc(edges, newN * sizeof(*edges)); if (!t3) exit(1); edges = t3;
    edges[N] = malloc( Z*kZ*sizeof(int));
    if (!edges[N]) exit(1);

    x[N] = SizeX+30.0;
    y[N] = SizeY/2.;
    z[N] = SizeZ/2.;

    numn[N] = 0;
    for(l=0; l<Z*kZ; l++) nn[N][l] = 0;
    for(l=0; l<Z*kZ; l++) edges[N][l] = 0;

    V[N] = 0.0;
    In[N] = 0.0;

    N = newN;

    M2 = M;
    i = N-1;

    // Probe node 2 connecting particles on the right face of sample
    for(k=0; k<N-2; k++)
        if ( x[k]+2.*radius[k]>SizeX) // vicinity criteria
		    {
			    j = k;

                nn[i][numn[i]] = j;
                edges[i][numn[i]] = M;
                E_origin = realloc(E_origin, (M + 1) * sizeof(int));
                if (!E_origin) exit(1);
                E_origin[M] = j;
                numn[i]++;

                nn[j][numn[j]] = i;
                edges[j][numn[j]] = M;
                E_end = realloc(E_end, (M + 1) * sizeof(int));
                if (!E_end) exit(1);
                E_end[M] = i;
                numn[j]++;

                tR = realloc(R, (M + 1) * sizeof(double));
                if (!tR) exit(1);
                R = tR;

                tIe = realloc(Ie, (M + 1) * sizeof(double));
                if (!tIe) exit(1);
                Ie = tIe;

				R[M] = Resistance_probe;
				Ie[M] = 0.0;

				M++;
		    }


    //---- Main Calculations --------------------------


    for(i=0; i<N-2; i++)
        {
           V[i] = 0.0; //V_supply*(1.-x[i]/SizeX);
        }

    int S = 0;
    int Si = 0;
    int Ss = 10;

    double step = 0.0;
    double dstep = 1.0;

    strcpy(s3, "output_"); // output file for visualization in OVITO (Optional feature)
    strcat(s3, s1);
    strcat(s3, ".xyz");

    FILE *f;
	f = fopen(s3,"w");

	FILE *f_iterations;
	f_iterations = fopen("R_iterations.txt","w");

	double Dipole_magn, Dipole_orientX, Dipole_orientY, Dipole_orientZ  ; // variables for visualization of edges in OVITO

	double S_V_R, S_1_R, S_Ie_N;
	double S_Ie,  S_Ie_av, dI_I;

    S_Ie = 1.0; // stopping criteria 1 (in paper)
    dI_I = 1.0; // stopping criteria 2 (test)

    //  while (dI_I > 0.000001) // stopping criteria 2 (test)
    while (S_Ie > 0.000001) // stopping criteria 1 (in paper)
	{
         if ( Si == Ss )
         {
            S_Ie = 0.0;
            for(i=0; i<N-2; i++)
		    {
			   S_Ie_N = 0.0;

			   for (j = 0; j < numn[i]; j++)
               {
                 S_Ie_N += (V[nn[i][j]]-V[i])/R[edges[i][j]]; // deviations from the current balance per node
               }

               S_Ie +=fabs(S_Ie_N); // sum of absolute values of deviations from the current balance for all nodes
		    }
            S_Ie_av = S_Ie/(N-2);

            I_probe_Left = 0.0; I_probe_Right = 0.0;
            for(i=M1; i<M2; i++) I_probe_Left += Ie[i];
            for(i=M2; i<M; i++) I_probe_Right += Ie[i];

            Resistance = (V[N-2]-V[N-1])/I_probe_Right;

            dI_I = (I_probe_Left-I_probe_Right)/I_probe_Left;

            printf(" Step %d Resistance = %lf I_Left= %lf I_Right = %lf dI/I = %lf Net nodes current  = %lf \n",S,Resistance,I_probe_Left, I_probe_Right, (I_probe_Left-I_probe_Right)/I_probe_Left,  S_Ie);

            fprintf(f_iterations,"%d %lf %lf\n",S*Ss,Resistance,S_Ie);

			S++;
			Si = 0;
          }

        for(i=0; i<N-2; i++)
		{
			S_V_R = 0.0; S_1_R = 0.0;

			for (j = 0; j < numn[i]; j++)
			{
              S_V_R += V[nn[i][j]]/R[edges[i][j]];
              S_1_R += 1.0/R[edges[i][j]];
			}

			V[i] = S_V_R/S_1_R; // main equation for Node potential calculation
		}

        for(i=0; i<M; i++)
			{
				k = E_origin[i];
				l = E_end[i];

                Ie[i] = fabs((V[l]-V[k])/R[i]); // absolute value of current flowing trough Edge
            }

    	step += dstep;
		Si++;
	}

	// saving data to OVITO-readable format

    fprintf(f,"%d\n",M1);

    fprintf(f,"Lattice=\"%d %d %d %d %d %d %d %d %d\"",1.0,0,0,0,1.0,0,0,0,1.0);
    fprintf(f," Properties=pos:I:3:Dipole:R:1:Dipole Orientation:R:2:Charge:R:1 Time=%lf\n",step);

    for(i=0; i<M1; i++)
    {
        if (V[E_origin[i]]>V[E_end[i]])
        {
            k = E_origin[i];
            l = E_end[i];
        }
        else
        {
            l = E_origin[i];
            k = E_end[i];
        }

        Dipole_magn = sqrt((x[l]-x[k])*(x[l]-x[k]) + (y[l]-y[k])*(y[l]-y[k]) + (z[l]-z[k])*(z[l]-z[k]));
        Dipole_orientX =  x[l]-x[k];
        Dipole_orientY =  y[l]-y[k];
        Dipole_orientZ =  z[l]-z[k];

        fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %10.8lf\n",x[k],y[k],z[k], Dipole_magn, Dipole_orientX, Dipole_orientY, Dipole_orientZ , Ie[i]);
    }

	fclose(f);
    fclose(f_iterations);

	I_probe_Right = 0.0;
    for(i=M2; i<M; i++) I_probe_Right += Ie[i];

    Resistance = (V[N-2]-V[N-1])/I_probe_Right; // final value of Net Resistance

    printf("Number of probe contacts on the Left  = %d\n",  M2 - M1);
    printf("Number of probe contacts on the Right = %d\n", M - M2);

    // final value of Net Resistivity taking into account dimensions of sample
	Resistivity = Resistance * ((SizeY-O_Y)*1.0e-6)  * ((SizeZ-O_Z)*1.0e-6) / ((SizeX-O_X)*1.0e-6);

    // saving results to the file
    FILE *f_results;
	f_results = fopen("Results.txt","a");

    fprintf(f_results,"File: %s Nodes = %d  Edges = %d Left contacts = %d Right contacts = %d NN = %lf Error  = %lf R_net = %lf Resistivity = %lE \n", s2,  N-2, M1, M2-M1, M-M2, 2.*M1/(N-2.), S_Ie, Resistance, Resistivity  );

	fclose(f_results);

	printf(" O_x = %lf Size_x = %lf\n", O_X, SizeX);
	printf(" O_y = %lf Size_y = %lf\n", O_Y, SizeY);
	printf(" O_z = %lf Size_z = %lf\n", O_Z, SizeZ);

	printf(" Resistance = %lf\n",Resistance );
	printf(" Resistivity  = %lE\n",Resistivity  );

    free(x);
	free(y);
	free(z);
	free(radius);
	free(V);
	free(In);
	free(numn);
	for (int i = 0; i < N; i++)
    {
        free(nn[i]);     // free each row
        free(edges[i]);  // free each row
    }
	free(nn);
	free(edges);
    free(R);
    free(E_origin);
    free(E_end);
	free(Ie);

    return 0;
}
