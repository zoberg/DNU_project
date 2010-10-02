#include <R.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "nrutil.h"

#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
#define imin(a,b) (iminarg1=(a),iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

extern double *dvector();
extern void free_dvector();
extern double fabs(), exp(), log();

double ComputePrimal(int mtrg, int input_dim, double gamma, int *support, double *alpha, double *xi)
{
   int i, j;
   double temp = 0;
   for (i=1; i <= mtrg; i++) 
      temp +=  log(1.0+exp(xi[i]));
   for (j=1; j <= input_dim; j++)
      if (support[j] == 1)
         temp += gamma*fabs(alpha[j]);
   return(temp);
}

void ComputeFAndDelta(int j, int mtrg,  double *xi, int *training_index, int *target, double **input, double *F, double *Delta)
{
   int i;
   double temp1, temp2 = 0, temp = 0;

   if (j == 0) {
      for (i=1; i <= mtrg; i++) {
          temp1 = exp(xi[i]);
          temp += temp1*target[training_index[i]]/(1+temp1);
          temp2 += temp1/((1+temp1)*(1+temp1));
      }
   } else {
      for (i=1; i <= mtrg; i++) {
          temp1 = exp(xi[i]);
          temp += temp1*target[training_index[i]]*input[training_index[i]][j]/(1+temp1);
          temp2 += temp1*input[training_index[i]][j]*input[training_index[i]][j]/((1+temp1)*(1+temp1));
      }
   }
   *F = temp;
   *Delta = temp2;
}

void UpdateXi(int j, int mtrg,  double *xi, int *training_index, int *target, double var_old, double var_new, double **input)
{
   int i;

   if (j == 0) {
       for (i=1; i <= mtrg; i++) 
           xi[i] += target[training_index[i]]*(var_new-var_old);

   } else {
       for (i=1; i <= mtrg; i++) 
           xi[i] += target[training_index[i]]*(var_old-var_new)*input[training_index[i]][j];
  }
}

void InitializeVariables(int mtrg, int input_dim, int *support,  double *alpha, double *xi,  double *bias)
{	
	int i,j;

        for (i=1; i <= mtrg; i++) 
	      xi[i]  = 0; 
        for (j=1; j <= input_dim; j++)
              support[j] = alpha[j] = 0;

        *bias = 0;
}

int FindMaxViolator(int k, int mtrg, int input_dim, double *xi, int *training_index, int *target, int *support, double **input, double *alpha, double gamma, double *F, double *Delta, double tol)
{

   int i, j, v;
   double MaxViol, viol_j, F_t, Delta_t;
   double dmaxarg1, dmaxarg2, dminarg1, dminarg2;


   if (k) {
       v = 0;
       ComputeFAndDelta(0, mtrg, xi, training_index, target, input,  &F_t, &Delta_t);
       MaxViol = fabs(F_t);

       for (j=1; j <= input_dim; j++) {
            if (support[j] == 1) {
                  ComputeFAndDelta(j, mtrg, xi, training_index, target, input,  &F_t, &Delta_t);
                  if (alpha[j] > 0)
                      viol_j =  fabs(gamma-F_t);
                  else
                      viol_j = fabs(gamma+F_t);

                  if (viol_j > MaxViol) {
                       MaxViol = viol_j; v = j;
                       *F = F_t; *Delta = Delta_t;
                  } 
            }
       }
    } else {
      v = -1;
      MaxViol = tol;
      for (j=1; j <= input_dim; j++) {
            if (support[j] == 0) {
                  ComputeFAndDelta(j, mtrg, xi, training_index, target, input,  &F_t, &Delta_t);
                    viol_j = DMAX(F_t-gamma, -F_t-gamma);
                    if (viol_j < 0)
                         viol_j = 0;
                    if (viol_j > MaxViol) {
                       MaxViol = viol_j; v = j;
                       *F = F_t; *Delta = Delta_t;
                    }
           }
       }
    }

    if (MaxViol > tol)
        return (v);
    else 
        return(-1);
}
     

void OptimizeForAlpha(int v, int mtrg, int *support, int *training_index, int *target, double *alpha, double *bias, double *xi, double **input, double *F, double *Delta,  double gamma, double tol, int input_dim)
{
    int i, flag;
    double temp, alpha_cur, alpha_new, slope_v, slope_0;
    double  Fv_at_0, Deltav_at_0, L, H;
    double *temp_xi;

    L = H = flag = 0;

    if (v == 0) {

       L = -DBL_MAX; H = DBL_MAX;

       slope_v = *F;

       alpha_cur = *bias;

    } else if (support[v] != 0) {
        
        temp_xi = dvector(1, mtrg);

        for (i=1; i <= mtrg; i++)
            temp_xi[i] = xi[i];

        alpha_cur = alpha[v]; alpha_new = 0;

        UpdateXi(v, mtrg, temp_xi, training_index, target, alpha_cur, alpha_new, input);
        ComputeFAndDelta(v, mtrg, temp_xi, training_index, target, input, &Fv_at_0, &Deltav_at_0);

        if ((alpha_cur > 0 && (gamma-*F) > 0 && (gamma-Fv_at_0) > 0 ) || (alpha_cur < 0 && (-gamma-*F) < 0 && (-gamma-Fv_at_0) < 0 )) {
            
            flag = 1;
            for (i=1; i <= mtrg; i++)
                xi[i] = temp_xi[i];
 
            *F = Fv_at_0; *Delta = Deltav_at_0;
            support[v] = alpha[v]  = 0;

            if (alpha_cur > 0) {
                slope_v = -gamma-Fv_at_0;
                L = -DBL_MAX; H = 0;
            } else {
                slope_v = gamma-Fv_at_0;
                L = 0; H = DBL_MAX;
            }

            alpha_cur = alpha_new = 0;
            if (-gamma-*F < 0 && gamma-*F > 0)
                slope_v = 0;
        } else {

            if (alpha[v] > 0) {

                slope_v = gamma-*F;

                if  (gamma-*F > 0 && gamma-Fv_at_0 < 0) {

                    L = 0; H = alpha[v];

                } else if (gamma-*F < 0) {

                    L = alpha[v]; H = DBL_MAX;
                }

            } else {

                slope_v = -gamma-*F;

               if (-gamma-*F < 0 && -gamma-Fv_at_0 > 0) {

                    L = alpha[v]; H = 0;

               } else if (-gamma-*F > 0) {

                    L = -DBL_MAX; H = alpha[v];

               }
            }
        }

        free_dvector(temp_xi, 1, mtrg);        

    } else {
       
         alpha_cur = alpha[v];
   
        if (gamma-*F < 0) {

             L = 0; H = DBL_MAX;

             slope_v = gamma-*F;
     
        } else if (-gamma-*F > 0) {

             L = -DBL_MAX; H = 0;
             
             slope_v = -gamma-*F;

        }

    }

    while (fabs(slope_v) > .1*tol) {

        alpha_new = alpha_cur - slope_v/(*Delta);

        if (alpha_new <= L || alpha_new >= H)
         
             alpha_new = (L+H)/2.0;
        
        UpdateXi(v, mtrg,  xi, training_index, target, alpha_cur, alpha_new, input);
        ComputeFAndDelta(v, mtrg,  xi, training_index, target, input, F, Delta);

        if (v != 0) {
        
           if (alpha_new > 0)

              slope_v = gamma-*F;
           else
   
             slope_v = -gamma-*F;

        } else

           slope_v = *F;
   
        if (slope_v > .1*tol)
 
             H = alpha_new;

        else if (slope_v < -.1*tol)

             L = alpha_new;

        alpha_cur = alpha_new;

     };


  

        
    if (v > 0) {
       alpha[v] = alpha_new;
       if (fabs(alpha[v]) > 0) 
           support[v] = 1;
    } else
       *bias = alpha_new;

}


void SparseLOGREGTrain(double tol, int mtrg, int input_dim, double **input, int *target, int *training_index, double *alpha, int *support,  double gamma_old, double gamma,  double *bias,  double *xi)
{
	int count, flag, i, j, k, id, p_maxviol;
	double primal, temp, MaxViol, Delta_v;
         double dmaxarg1, dmaxarg2, dminarg1, dminarg2;
        double F, Delta;

        if (gamma_old < 1e-6) 
            InitializeVariables(mtrg, input_dim, support, alpha, xi,  bias);
	
        id = FindMaxViolator(0, mtrg, input_dim, xi, training_index, target, support, input, alpha, gamma, &F, &Delta, tol);

        if (id == -1)
           id = FindMaxViolator(1, mtrg, input_dim, xi, training_index, target,support,  input, alpha, gamma, &F, &Delta, tol);

        while ( id >= 0) {
           do {
              OptimizeForAlpha(id, mtrg, support, training_index, target, alpha, bias, xi, input, &F, &Delta,  gamma, tol, input_dim);
              id = FindMaxViolator(1, mtrg, input_dim, xi, training_index, target, support, input, alpha, gamma, &F, &Delta, tol);
           } while (id >= 0);
          id = FindMaxViolator(0, mtrg, input_dim, xi, training_index, target, support, input, alpha, gamma, &F, &Delta, tol);
        };
}




double LogregValidate(int m, int mtrg, int input_dim, double **input, int *target, int *training_index, int *val_index, double *alpha, int *support, double *xi, double *b, double gamma,  double *val_output)
{
	int i, j;
	static int k = 0;
	double temp, cost;
	
	cost = 0;
	for (i=1; i <= m-mtrg; i++) {
		k++;
		temp = 0;
		for (j=1; j <= input_dim; j++)
			temp += alpha[j]*input[val_index[i]][j];
                temp -= *b;

                val_output[i] = temp;
		
                if (target[val_index[i]] > 0) {
                    if (temp < 0)
                        cost += 1.0; 
                } else {
                    if (temp > 0)
                        cost += 1.0;
                }
	}
	return(cost);
}


void my_fun(double *zi, double *ansi, int *ni)
{
	int i, j;
	double ans = ansi[0];
	int m = ni[0];
	int n = ni[1];
	double **z = dmatrix(1, m, 1, n);
	//printf("zi=\n");
	//for(i=1; i<=n*n; i++){
	//	printf(" %f", zi[i-1]);
	//}
	printf("Ans=%f\n", ans);
	printf("Length=%d%d\n", n, m);
	/*for (i=1; i <= m; i++){
		for (j=1; j <= n; j++){
			z[i][j] = zi[i + (j-1)*m - 1];
			//printf("input z[%d]=%f\n", j + (i-1)*n - 1, zi[j + (i-1)*n - 1] );
			//printf("z[%d][%d]=%f\n", i, j, z[i][j]);
		}
	}*/
	my_transform(zi, z, m, n);
	for (i=1; i <= m; i++){
		for(j=1; j <= n; j++){
			ans = ans + z[i][j];
			//printf("Current ans[%d][%d]=%f\n", i, j, ans);
		}
	}
	ansi[0] = ans;
	printf("Summ=%f\n", ans);
}


void my_transform(double *input, double **output, int m, int n)
{
	int i, j;
	for (i=1; i <= m; i++){
		for (j=1; j <= n; j++){
			output[i][j] = input[i + (j-1)*m - 1];
		}
	};
}

void my_transform_int(int *input, int **output, int m, int n)
{
	int i, j;
	for (i=1; i <= m; i++){
		for (j=1; j <= n; j++){
			output[i][j] = input[i + (j-1)*m - 1];
		}
	};
}

void my_transform_int_back(int **input, int *output, int m, int n)
{
	int i, j;
	for (i=1; i <= m; i++){
		for (j=1; j <= n; j++){
			output[i + (j-1)*m - 1] = input[i][j];
		}
	};
}

void my_v_transform(int *input, int *output, int n)
{
	int i;
	for (i=1; i <= n; i++){
		output[i] = input[i - 1];
		}
}

void my_v_transform_double(double *input, double *output, int n)
{
	int i;
	for (i=1; i <= n; i++){
		output[i] = input[i - 1];
		}
}

void my_v_transform_double_back(double *input, double *output, int n)
{
	int i;
	for (i=1; i <= n; i++){
		output[i-1] = input[i];
		}
}

void my_indexx(double *idxp, double *idxn, int *mplusi, int *mnegi)
{
		
	long id=-2;	
	int i, mplus=mplusi[0], mneg=mnegi[0];
	float *ranpos, *ranneg;
	unsigned long *temp_idxp, *temp_idxn;

	temp_idxp = lvector(1, mplus);
        temp_idxn = lvector(1, mneg);
	ranpos = vector(1,mplus);
        ranneg = vector(1,mneg);

	for (i=1; i <= mplus; i++)
        	ranpos[i] = ran1(&id);
        for (i=1; i <= mneg; i++)
        	ranneg[i] = ran1(&id);
	
	indexx(mplus, ranpos, temp_idxp);
	indexx(mneg, ranneg, temp_idxn);

	for (i=1; i <= mplus; i++)
        	idxp[i-1] = temp_idxp[i];
        for (i=1; i <= mneg; i++)
        	idxn[i-1] = temp_idxn[i];	

	free_lvector(temp_idxp, 1, mplus);
	free_lvector(temp_idxn, 1, mneg);
	free_vector(ranpos, 1,mplus);
        free_vector(ranneg, 1,mneg);
}


void my_train(int *Intkfoldi, int *mplusi, int *mnegi, int *val_indexi, int *training_indexi, int *pos_indexi, int *neg_indexi, int *idxpi, int *idxni)
{
	int posex_per_fold, negex_per_fold, m_train=62;
	int  posex_in_last_fold, negex_in_last_fold;
	int i, k, k1, k2, id1, id2, Intkfold=Intkfoldi[0], mplus=mplusi[0], mneg=mnegi[0];
	int **training_index, **val_index, *neg_index, *pos_index, *idxn, *idxp;

	training_index = imatrix(1,Intkfold,1,m_train);
	val_index = imatrix(1,Intkfold,1,m_train);
	neg_index = ivector(1, mneg);
	pos_index = ivector(1, mplus);
	idxn = ivector(1, mneg);
	idxp = ivector(1, mplus);

	my_transform_int(training_indexi, training_index, Intkfold, m_train);
	my_transform_int(val_indexi, val_index, Intkfold, m_train);

	my_v_transform(neg_indexi, neg_index, mneg);
	my_v_transform(pos_indexi, pos_index, mplus);
	my_v_transform(idxni, idxn, mneg);
	my_v_transform(idxpi, idxp, mplus);

	posex_per_fold =  floor(mplus / (double)Intkfold);
	negex_per_fold =  floor(mneg / (double) Intkfold);
	posex_in_last_fold = mplus - (Intkfold-1)*posex_per_fold;
	negex_in_last_fold = mneg - (Intkfold-1)*negex_per_fold;

	for (k=1; k <= Intkfold; k++) {
                   id1 = id2 = 1;
                   k1 = (k-1)*posex_per_fold+1;
                   k2 = (k-1)*negex_per_fold+1;
                   if (k < Intkfold) {
                        for (i=1; i <= mplus; i++) {
                            if (i >= k1 && i <= k1 + posex_per_fold - 1) {
                                 val_index[k][id2] = pos_index[idxp[i]]; id2++;
                             } else {
                                 training_index[k][id1] = pos_index[idxp[i]];id1++;
                             }
                        }
                        for (i=1; i <= mneg; i++) {
                             if (i >= k2 && i <= k2 + negex_per_fold - 1) {
                                  val_index[k][id2] = neg_index[idxn[i]]; id2++;
                             } else {
                                  training_index[k][id1] = neg_index[idxn[i]];id1++;
                             }
                         }
                   } else {
                         for (i=1; i <= mplus; i++) {
                             if (i >= k1 && i <= k1 + posex_in_last_fold - 1) {
                                val_index[k][id2] = pos_index[idxp[i]]; id2++;
                             } else {
                                training_index[k][id1] = pos_index[idxp[i]];id1++;
                             }
                         }
                         for (i=1; i <= mneg; i++) {
                              if (i >= k2 && i <= k2 + negex_in_last_fold - 1) {
                                  val_index[k][id2] = neg_index[idxn[i]]; id2++;
                              } else {
                                  training_index[k][id1] = neg_index[idxn[i]];id1++;
                              }
                         }
                    }
	}
	my_transform_int_back(training_index, training_indexi, Intkfold, m_train);
	my_transform_int_back(val_index, val_indexi, Intkfold, m_train);

}

void my_SLR(double *tol, int *input_dim, double *input, int *total_y, int *training_index, int *trgex_per_fold, double *gamma, double *gamma_old, double *alpha, double *bias)
{
	int m = trgex_per_fold[0], j;
	double **total_input = dmatrix(1, m, 1, input_dim[0]);
	int *support = ivector(1,input_dim[0]);
	my_transform(input, total_input, m, input_dim[0]);
	double *xi = dvector(1, m);
	for (j=1; j <= input_dim[0]; j++) {
		support[j]=0;
	}

	for (j=1; j <= m; j++) {
		xi[j] = 0;
	}

	SparseLOGREGTrain(tol[0], trgex_per_fold[0], input_dim[0], total_input, total_y, training_index, alpha, support, gamma_old[0], gamma[0], &bias[0], xi);

	free_dmatrix(total_input, 1, m, 1, input_dim[0]);
	free_dvector(xi, 1, m);
	free_ivector(support, 1, input_dim[0]);	
}
