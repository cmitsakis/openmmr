// Code copied from https://www.csie.ntu.edu.tw/~cjlin/papers/online_ranking/
// Licensed under the terms of the BSD-3-Clause license.
// You can find the full text of the license in the file "weng_lin_license".
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rank.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

#define max(x,y) (((x)>(y))?(x):(y))
#define PREDICTIONS_SIZE 30
#define COMPETITIVE_TEST_NUM 1688
#define RARE_CASE 5
#define M_PI 3.14159265358979323846264338327950288

int method;
int rare_players;
int competitive_test = 0;
int competitive_test_game = 0;

/*void exit_with_help()
{
	printf(
	"Usage: rank method data_file\n"
	"method 0: BT model (full)\n"
	"method 1: BT model (partial)\n"
	"method 2: PL model\n"
	"method 3: TM model (full)\n"
	"method 4: TM model (partial)\n"
	"method 5: Glicko (two-player games only)\n"
	);
	exit(1);
}*/


double phi(double x)
{
	return exp(-0.5*x*x)/sqrt(2*M_PI);
}

double Phi(double x)
{
	double z = fabs(x/sqrt(2.0));
	double t = 1/(1+0.5*z);
	double res  = t * exp (-z * z - 1.26551223 + t * (1.00002368 + t *
(0.37409196 + t * (0.09678418 + t * (-0.18628806 + t *
(0.27886807 + t * (-1.13520398 + t * (1.48851587 + t *
				      (-0.82215223 + t * 0.17087277))))))))) ;
	if (x >=0) 
		res = 2.0-res;
	return res/2.0;
}

// the calculation of V, W, Vt, Wt follows TrueSkill

double min_value = 2.222758749e-162; // 1.0e-50;

double V(double x, double t)
{
	double xt = x-t;
	double denom = Phi(x-t);
	// 2.222758749e-162)
	if (denom < min_value)
		//	if (xt < -10)
		return -xt;
	else
	{
		double tmp= phi(xt)/denom;
		if (isnan((float) tmp))
			fprintf(stderr, "nan %lf %lf %lf %lf\n", xt, phi(xt), Phi(xt), tmp);
		return tmp;
	}
}

double W(double x, double t)
{
	double xt = x-t;
        double denom = Phi(xt);
	if (denom < min_value)
	//	if (xt < -10)
	{
		if (x < 0)
			return(1);
		else
			return(0);
	}
	else
		return V(x, t) * (V(x, t) + xt);
}

double Vt(double x, double t)
{
	double xx = fabs(x);
	double b = (Phi(t-xx)-Phi(-t-xx));

	if (b < 1e-5)
	{
		if (x < 0)
			return -x-t;
		else
			return -x+t;
	}
	else
	{
		double a = phi(-t-xx) - phi(t-xx);
		if (x < 0)
			return -a/b;
		else
			return a/b;
	}
}

double Wt(double x, double t)
{
	double xx = fabs(x);
	double b = (Phi(t-xx)-Phi(-t-xx));

	if (b < min_value)
		return 1.0;
	else
		return ((t-xx)*phi(t-xx)+(t+xx)*phi(-t-xx))/b + Vt(x,t)*Vt(x,t);
}

void evaluate(int num_teams, int *rank,
	      double *team_mu, double *team_sigmasq, 
	      int *predictions) 
{
	int i, j;
	for (i=0; i<num_teams; i++)
	{
		for (j=i+1; j<num_teams; j++)
		{
			if (rank[i] != rank[j])
			{
				int hard = 0;
				if (abs(rank[i]-rank[j]) <= 2)
				{
					hard = 1;
					predictions[5]++;
				}
				if (rare_players == 1)
					predictions[10]++;
				if (competitive_test_game == 1)
					predictions[15]++;

				predictions[0]++;
				if ((rank[i]-rank[j])*(team_mu[i] - team_mu[j]) < 0)
				{
					predictions[1]++;
					predictions[6]+=hard;
					predictions[11]+=rare_players;
					predictions[16]+=competitive_test_game;
				}
				else if ((rank[i]-rank[j])*(team_mu[i] - team_mu[j]) > 0)
				{
					predictions[2]++;
					predictions[7]+=hard;
					predictions[12]+=rare_players;
					predictions[17]+=competitive_test_game;

				}
				double team_i = team_mu[i]-3*sqrt(team_sigmasq[i]);
				double team_j = team_mu[j]-3*sqrt(team_sigmasq[j]);
				if ((rank[i]-rank[j])*(team_i - team_j) < 0)
				{
					predictions[3]++;
					predictions[8]+=hard;
					predictions[13]+=rare_players;
					predictions[18]+=competitive_test_game;
				}
				else if ((rank[i]-rank[j])*(team_i - team_j) > 0)
				{
					predictions[4]++;
					predictions[9]+=hard;
					predictions[14]+=rare_players;
					predictions[19]+=competitive_test_game;
				}
			}
		}
	}
}


void find_team_mu_sigma(int num_teams, int *team_num_players, 
			double *team_mu, double *team_sigmasq, 
			double **mu, double **sigma)
{
	int i, j;
	for (i=0; i<num_teams; i++)
	{
		team_mu[i] = mu[i][0];
		team_sigmasq[i] = sigma[i][0]*sigma[i][0];
		for (j=1; j<team_num_players[i]; j++)
		{
			team_mu[i] += mu[i][j];
			team_sigmasq[i] += sigma[i][j]*sigma[i][j];
		}
		//		printf("mu %3.5f %3.5f ", team_mu[i], sqrt(team_sigmasq[i]));
	}
	//	printf("\n");
}

// rank starts with 0. If scores are [3 1 4]
// then rank is [1 2 0]

void update_BT_full(double beta, 
		    int num_teams, int *rank, 
		    int *team_num_players, 
		    double **mu,
		    double **sigma, int *predictions)
{
	int i, j, q;
	double *team_mu = Malloc(double, num_teams);
	double *team_sigmasq = Malloc(double, num_teams);

	//printf("num_teams=%d team_num_players[0]=%d\nmu[0][0]=%3.5f mu[0][1]=%3.5f mu[1][0]=%3.5f\nsigma[0][0]=%3.5f sigma[0][1]=%3.5f\n", num_teams, team_num_players[0], mu[0][0], mu[0][1], mu[1][0], sigma[0][0], sigma[0][1]); // for debugging
	find_team_mu_sigma(num_teams, team_num_players, team_mu, team_sigmasq, mu, sigma);

	evaluate(num_teams, rank, team_mu, team_sigmasq, predictions); 

	// update

	for (i=0; i<num_teams; i++)
	{
		double Omega = 0;
		double Delta = 0;

		for (q=0; q<num_teams; q++)
		{
			if (q == i) continue;

			double ciq = sqrt(team_sigmasq[i] +
					  team_sigmasq[q] + 2*beta*beta);
			double piq = 1/(1+exp((team_mu[q]-team_mu[i])/ciq));
			double sigsq_to_ciq = team_sigmasq[i]/ciq;
			double s = 0;
			if (rank[q] > rank[i])
				s = 1;
			else if (rank[q] == rank[i])
				s = 0.5;

			Omega += sigsq_to_ciq* (s-piq);
			double gamma = sqrt(team_sigmasq[i])/ciq;
			// double gamma = 1;
			// double gamma = 1.0/num_teams;
			Delta += gamma*sigsq_to_ciq/ciq*piq*(1-piq);
		}
		for (j=0; j<team_num_players[i]; j++)
		{
			double sigmaijsq = sigma[i][j]*sigma[i][j];
			mu[i][j] += sigmaijsq/team_sigmasq[i]*Omega;
			sigma[i][j] *= sqrt(max(1-sigmaijsq/team_sigmasq[i]*Delta, 0.0001));
		}
		//		printf("%3.5f %3.5f ", mu[i][0], sigma[i][0]);
	}
	//	printf("\n");

	free(team_mu);
	free(team_sigmasq);
}


void update_BT_partial(double beta, 
		       int num_teams, int *rank, 
		       int *rank_idx,
		       int *team_num_players, 
		       double **mu, double **sigma, 
		       int *predictions)
{
	int i, j, q;
	int q12[2];
	int a;
	double *team_mu = Malloc(double, num_teams);
	double *team_sigmasq = Malloc(double, num_teams);

	find_team_mu_sigma(num_teams, team_num_players, team_mu, team_sigmasq, mu, sigma);

	evaluate(num_teams, rank, team_mu, team_sigmasq, predictions); 

	// update

	for (a=0; a<num_teams; a++)
	{
		i = rank_idx[a];

		double Omega = 0;
		double Delta = 0;

		int q12_size = 0;
		if (a >= 1)
		{
			q12[q12_size] = rank_idx[a-1];
			q12_size++;
		}
		if (a <= num_teams-2)
		{
			q12[q12_size] = rank_idx[a+1];
			q12_size++;
		}

		for (j=0; j<q12_size; j++)
		{
			q = q12[j];

			if (q == i) continue;

			double ciq = sqrt(team_sigmasq[i] +
					  team_sigmasq[q] + 2*beta*beta);
			double piq = 1/(1+exp((team_mu[q]-team_mu[i])/ciq));
			double sigsq_to_ciq = team_sigmasq[i]/ciq;
			double s = 0;
			if (rank[q] > rank[i])
				s = 1;
			else if (rank[q] == rank[i])
				s = 0.5;

			Omega += sigsq_to_ciq*(s-piq);
			double gamma = sqrt(team_sigmasq[i])/ciq;
			// double gamma = 1;
			//			double gamma = 1.0/num_teams;
			Delta += gamma*sigsq_to_ciq/ciq*piq*(1-piq);

		}
		for (j=0; j<team_num_players[i]; j++)
		{
			double sigmaijsq = sigma[i][j]*sigma[i][j];
			mu[i][j] += sigmaijsq/team_sigmasq[i]*Omega;
			sigma[i][j] *= sqrt(max(1-sigmaijsq/team_sigmasq[i]*Delta, 0.0001));
		}
		//		printf("%3.5f %3.5f ", mu[i][0], sigma[i][0]);
	}
	//	printf("\n");

	free(team_mu);
	free(team_sigmasq);
}

void update_PL(double beta, 
	       int num_teams, int *rank, 
	       int *team_num_players, 
	       double **mu,
	       double **sigma, int *predictions)
{
	int i, j, q;
	double *team_mu = Malloc(double, num_teams);
	double *team_sigmasq = Malloc(double, num_teams);

	find_team_mu_sigma(num_teams, team_num_players, team_mu, team_sigmasq, mu, sigma);

	evaluate(num_teams, rank, team_mu, team_sigmasq, predictions); 

	// update

	double c = 0;
	for (i=0; i<num_teams; i++)
		c += team_sigmasq[i] + beta*beta;
	c = sqrt(c);
	int *A = Malloc(int, num_teams);
	double *sum_q = Malloc(double, num_teams);
	for (i=0; i<num_teams; i++)
	{
		A[i] = 0;
		sum_q[i] = 0;
	}
	for (i=0; i<num_teams; i++)
	{
		double tmp = exp(team_mu[i]/c);
		for (q=0; q<num_teams; q++)
		{
			if (rank[i] >= rank[q])
				sum_q[q] += tmp;
			if (rank[i] == rank[q])
				A[q]++;
		}
	}

	for (i=0; i<num_teams; i++)
	{
		double Omega = 0;
		double Delta = 0;
		double tmp1 = exp(team_mu[i]/c);
		for (q=0; q<num_teams; q++)
		{
			double tmp = tmp1/sum_q[q];
			if (rank[q] <= rank[i])
			{
				Delta += tmp*(1-tmp)/(double) A[q];
				if (q==i)
					Omega += (1-tmp)/(double) A[q];
				else
					Omega -= tmp/(double) A[q];
			}
		}
		Omega *= team_sigmasq[i]/c;
		Delta *= team_sigmasq[i]/c/c;

		double gamma = sqrt(team_sigmasq[i])/c;
		// using smaller Delta
		Delta *= gamma;

		for (j=0; j<team_num_players[i]; j++)
		{
			double sigmaijsq = sigma[i][j]*sigma[i][j];
			mu[i][j] += sigmaijsq/team_sigmasq[i]*Omega;
			sigma[i][j] *= sqrt(max(1-sigmaijsq/team_sigmasq[i]*Delta, 0.0001));
		}
		//		printf("%3.5f %3.5f ", mu[i][0], sigma[i][0]);
	}
	//	printf("\n");

	free(team_mu);
	free(team_sigmasq);
	free(sum_q);
	free(A);
}


void update_TM_full(double beta, 
		    int num_teams, int *rank, 
		    int *team_num_players, 
		    double **mu,
		    double **sigma, int *predictions)
{
	double epsilon = 0.1;

	int i, j, q;
	double *team_mu = Malloc(double, num_teams);
	double *team_sigmasq = Malloc(double, num_teams);

	find_team_mu_sigma(num_teams, team_num_players, team_mu, team_sigmasq, mu, sigma);

	evaluate(num_teams, rank, team_mu, team_sigmasq, predictions); 

	// update

	for (i=0; i<num_teams; i++)
	{
		double Omega = 0;
		double Delta = 0;

		for (q=0; q<num_teams; q++)
		{
			if (q == i) continue;

			double ciq = sqrt(team_sigmasq[i] +
					  team_sigmasq[q] + 2*beta*beta);
			double tmp = (team_mu[i]-team_mu[q])/ciq;
			double sigsq_to_ciq = team_sigmasq[i]/ciq;
			double gamma = sqrt(team_sigmasq[i])/ciq;
			if (rank[q] > rank[i])
			{
				Omega += sigsq_to_ciq*V(tmp, epsilon/ciq);
				Delta += gamma*sigsq_to_ciq/ciq*W(tmp, epsilon/ciq);
			}
			else if (rank[q] < rank[i])
			{
				Omega += -sigsq_to_ciq*V(-tmp, epsilon/ciq);
				Delta += gamma*sigsq_to_ciq/ciq*W(-tmp, epsilon/ciq);
			}
			else 
			{
				Omega += sigsq_to_ciq*Vt(tmp, epsilon/ciq);
				Delta += gamma*sigsq_to_ciq/ciq*Wt(tmp, epsilon/ciq);
			}
		}
		for (j=0; j<team_num_players[i]; j++)
		{
			double sigmaijsq = sigma[i][j]*sigma[i][j];
			mu[i][j] += sigmaijsq/team_sigmasq[i]*Omega;
			sigma[i][j] *= sqrt(max(1-sigmaijsq/team_sigmasq[i]*Delta, 0.0001));
		}
		//		printf("%3.5f %3.5f ", mu[i][0], sigma[i][0]);
	}
	//	printf("\n");

	free(team_mu);
	free(team_sigmasq);
}


void update_TM_partial(double beta, 
			  int num_teams, int *rank, 
			  int *rank_idx,
			  int *team_num_players, 
			  double **mu,
			  double **sigma, int *predictions)
{
	double epsilon = 0.1;

	int i, j, q;
	int q12[2];
	int a;
	double *team_mu = Malloc(double, num_teams);
	double *team_sigmasq = Malloc(double, num_teams);

	find_team_mu_sigma(num_teams, team_num_players, team_mu, team_sigmasq, mu, sigma);

	evaluate(num_teams, rank, team_mu, team_sigmasq, predictions); 

	// update

	for (a=0; a<num_teams; a++)
	{
		i = rank_idx[a];

		double Omega = 0;
		double Delta = 0;

		int q12_size = 0;
		if (a >= 1)
		{
			q12[q12_size] = rank_idx[a-1];
			q12_size++;
		}
		if (a <= num_teams-2)
		{
			q12[q12_size] = rank_idx[a+1];
			q12_size++;
		}

		for (j=0; j<q12_size; j++)
		{
			q = q12[j];

			if (q == i) continue;

			double ciq = 2*sqrt(team_sigmasq[i] +
					  team_sigmasq[q] + 2*beta*beta);
			double tmp = (team_mu[i]-team_mu[q])/ciq;
			double sigsq_to_ciq = team_sigmasq[i]/ciq;
			if (rank[q] > rank[i])
			{
				Omega += sigsq_to_ciq*V(tmp, epsilon/ciq);
				Delta += sigsq_to_ciq/ciq*W(tmp, epsilon/ciq);
			}
			else if (rank[q] < rank[i])
			{
				Omega += -sigsq_to_ciq*V(-tmp, epsilon/ciq);
				Delta += sigsq_to_ciq/ciq*W(-tmp, epsilon/ciq);
			}
			else 
			{
				Omega += sigsq_to_ciq*Vt(tmp, epsilon/ciq);
				Delta += sigsq_to_ciq/ciq*Wt(tmp, epsilon/ciq);
			}
			// Delta *= sqrt(team_sigmasq[i])/ciq;
			//			Delta *= team_sigmasq[i]/ciq/ciq;
		}
		for (j=0; j<team_num_players[i]; j++)
		{
			double sigmaijsq = sigma[i][j]*sigma[i][j];
			mu[i][j] += sigmaijsq/team_sigmasq[i]*Omega;
			sigma[i][j] *= sqrt(max(1-sigmaijsq/team_sigmasq[i]*Delta, 0.0001));
		}
		//		printf("%3.5f %3.5f ", mu[i][0], sigma[i][0]);
	}
	//	printf("\n");

	free(team_mu);
	free(team_sigmasq);
}

void update_Glicko(double beta, 
		   int num_teams, int *rank, 
		   int *team_num_players, 
		   double **mu,
		   double **sigma, int *predictions)
{
	int i, j;

	if (num_teams != 2)
		printf("Glicko supports only two-player games\n");
	for (i=0; i<1; i++)
		if (team_num_players[i] != 1)
			printf("Each team can have one player\n");

	double *team_mu = Malloc(double, num_teams);
	double *team_sigmasq = Malloc(double, num_teams);

	find_team_mu_sigma(num_teams, team_num_players, team_mu, team_sigmasq, mu, sigma);

	evaluate(num_teams, rank, team_mu, team_sigmasq, predictions); 

	const double q = log(10)/400;
	double g[2], E[2], delta_sq[2];

	for (i=0; i<2; i++)
		g[i] = 1/sqrt(1+3*q*q*team_sigmasq[i]/M_PI/M_PI);

	printf("g[0] %g g[1] %g\n", g[0], g[1]);

	// update

	for (i=0; i<2; i++)
	{
		if (i==0)
			j = 1;
		else
			j = 0;
		
		E[i] = 1/(1+pow(10.0, -g[j]*(team_mu[i]-team_mu[j])/400));
		delta_sq[i] = 1/(q*q*g[j]*g[j]*E[i]*(1-E[i]));
	}

	printf("rank %d %d\n", rank[0], rank[1]);
	printf("old mu,sigma %3.5f %3.5f %3.5f %3.5f\n", mu[0][0], mu[1][0], sigma[0][0], sigma[1][0]);

	for (i=0; i<2; i++)
	{
		if (i==0)
			j = 1;
		else
			j = 0;
		
		E[i] = 1/(1+pow(10.0, -g[j]*(team_mu[i]-team_mu[j])/400));
		delta_sq[i] = 1/(q*q*g[j]*g[j]*E[i]*(1-E[i]));
	}

	for (i=0; i<2; i++)
	{
		if (i==0)
			j = 1;
		else
			j = 0;
		double s = 0;
		if (rank[i] < rank[j])
			s = 1;
		else if (rank[i] == rank[j])
			s = 0.5;
		mu[i][0] += q/(1/(team_sigmasq[i]+1/delta_sq[i]))*g[j]*g[j]*(s-E[i]);
		sigma[i][0] = sqrt(1/(1/team_sigmasq[i]+1/delta_sq[i]));
	}
	printf("new mu,sigma %3.5f %3.5f %3.5f %3.5f\n", mu[0][0], mu[1][0], sigma[0][0], sigma[1][0]);

	free(team_mu);
	free(team_sigmasq);
}

/*int main(int argc,char **argv)
{
	int i, j, k;
	FILE *fp;

	double init_mu = 25;
	double init_sigma = init_mu/3;
	double beta = init_sigma*0.5;

	if (argc != 3)
		exit_with_help();
	method = atoi(argv[1]);
	fp=fopen(argv[2],"r");
	
	if(fp==NULL)
	{
		fprintf(stderr,"can't open file %s\n", argv[2]);
		exit(1);
	}

	// check if we will do competitive tests
	FILE *ct_infile = NULL;
	FILE *ct_outfile = NULL;
	int competitive_test_cases = 0;
	if (strcmp(argv[2], "../data/HeadToHead.data") == 0)
	{
		competitive_test = 1;
		ct_infile = fopen("HeadToHead.ctcases_ts", "r");
		ct_outfile = fopen("HeadToHead.ctcases_ours", "w");
	}

	int num_games, num_players;
	fscanf(fp,"%d\n",&num_games);
	fscanf(fp,"%d\n",&num_players);

	double *players_mu = Malloc(double, num_players);
	double *players_sigma = Malloc(double, num_players);
	int *players_occurances = Malloc(int, num_players);

	int *predictions = Malloc(int, PREDICTIONS_SIZE);

	for (i=0; i<num_players; i++)
	{
		players_mu[i] = init_mu;
		players_sigma[i] = init_sigma;
		players_occurances[i] = 0;
	}

	for (i=0; i<num_games; i++)
	{
		int num_teams, player_id;
		int num_players_ateam;
		int num_players_agame = 0;
		fscanf(fp,"%d\n",&num_teams);
		double **team_players_mu = Malloc(double *, num_teams);
		double **team_players_sigma = Malloc(double *, num_teams);
		int **team_players_id = Malloc(int *, num_teams);
		int *team_scores = Malloc(int, num_teams);
		int *rank = Malloc(int, num_teams);
		int *rank_idx = Malloc(int, num_teams);
		int *team_num_players = Malloc(int, num_teams);
		double avg_occurances = 0;

		for (j=0; j<num_teams; j++)
		{
			fscanf(fp,"%d\n",&(team_scores[j]));
			fscanf(fp,"%d\n",&num_players_ateam);
			num_players_agame += num_players_ateam;
			team_num_players[j] = num_players_ateam;
			team_players_mu[j] = Malloc(double, num_players_ateam);
			team_players_sigma[j] = Malloc(double, num_players_ateam);
			team_players_id[j] = Malloc(int, num_players_ateam);

			for (k=0; k<num_players_ateam; k++)
			{
				fscanf(fp,"%d\n",&player_id);
				avg_occurances += (double) players_occurances[player_id];
				players_occurances[player_id]++;
				team_players_mu[j][k] = players_mu[player_id];
				team_players_sigma[j][k] = players_sigma[player_id];
				team_players_id[j][k] = player_id;
			}
		}

		// check occurances
		avg_occurances /= num_players_agame;
		if (avg_occurances <= RARE_CASE)
			rare_players = 1;
		else
			rare_players = 0;

		int old_errors = predictions[2];
		if (competitive_test == 1)
		{
			fscanf(ct_infile, "%d\n", &competitive_test_game);
		}

		// convert scores to rank; now scores are decreasing
		int s = 0;
		for (j=0; j<num_teams; j++)
		{
			if (j > 0)
				if (team_scores[j-1] > team_scores[j])
					s = j;
			rank[j]=s;
			rank_idx[j] = j;
		}

		// start evaluation at the 2nd game
		if (i==1)
			for (j=0; j<PREDICTIONS_SIZE; j++)
				predictions[j] = 0;
		
		switch (method)
		{
		case 0:
			update_BT_full(beta, num_teams, rank, 
				       team_num_players, team_players_mu,
				       team_players_sigma, predictions);
			break;
		case 1:
			update_BT_partial(beta, num_teams, rank, rank_idx,
					  team_num_players, team_players_mu,
					  team_players_sigma, predictions);
			break;
		case 2:
			update_PL(beta, num_teams, rank, 
				       team_num_players, team_players_mu,
				       team_players_sigma, predictions);
			break;
		case 3:
			update_TM_full(beta, num_teams, rank, 
					     team_num_players, team_players_mu,
					     team_players_sigma, predictions);
			break;
		case 4:
			update_TM_partial(beta, num_teams, rank, rank_idx,
					     team_num_players, team_players_mu,
					     team_players_sigma, predictions);
			break;
		case 5:
			update_Glicko(beta, num_teams, rank, 
				      team_num_players, team_players_mu,
				      team_players_sigma, predictions);
			break;
		default:      
			exit_with_help();
		}

		if (competitive_test == 1)
		{
			if (i>= 1 && predictions[2] > old_errors && competitive_test_cases < COMPETITIVE_TEST_NUM)
			{
				fprintf(ct_outfile, "1\n");
				competitive_test_cases ++;
			}
			else 
				fprintf(ct_outfile, "0\n"); 
		}

		for (j=0; j<num_teams; j++)
		{
			for (k=0; k<team_num_players[j]; k++)
			{
				player_id = team_players_id[j][k];

				players_mu[player_id] = team_players_mu[j][k];
				players_sigma[player_id] = team_players_sigma[j][k];
			}
			free(team_players_mu[j]);
			free(team_players_sigma[j]);
			free(team_players_id[j]);
		}
		free(team_players_mu);
		free(team_players_sigma);
		free(team_players_id);
		free(team_scores);
		free(rank);
		free(rank_idx);
	}

	printf("Using mu: correct %d\n"
	       "total %d error %g\n",
	       predictions[1], predictions[0],
	       1.0-(double) predictions[1]/ (double) predictions[0]);
	printf("Using mu-3sigma: correct %d\n"
	       "total %d error %g\n",
	       predictions[3], predictions[0],
	       1.0-(double) predictions[3]/ (double) predictions[0]);
	printf("Using mu (hard cases): correct %d\n"
	       "total %d error %g\n",
	       predictions[6], predictions[5],
	       1.0-(double) predictions[6]/ (double) predictions[5]);
	printf("Using mu-3sigma (hard cases): correct %d\n"
	       "total %d error %g\n",
	       predictions[8], predictions[5],
	       1.0 - (double) predictions[8]/ (double) predictions[5]);
	printf("Using mu (rare cases): correct %d\n"
	       "total %d error %g\n",
	       predictions[11], predictions[10],
	       1.0-(double) predictions[11]/ (double) predictions[10]);
	printf("Using mu-3sigma (rare cases): correct %d\n"
	       "total %d error %g\n",
	       predictions[13], predictions[10],
	       1.0 - (double) predictions[13]/ (double) predictions[10]);
	if (competitive_test == 1) 
	{
		printf("Using mu (competitive test): correct %d\n"
		       "total %d error %g\n",
		       predictions[16], predictions[15],
		       1.0-(double) predictions[16]/ (double) predictions[15]);
		printf("Using mu-3sigma (competitive test): correct %d\n"
		       "total %d error %g\n",
		       predictions[18], predictions[15],
		       1.0 - (double) predictions[18]/ (double) predictions[15]);
	}
	free(players_mu);
	free(players_sigma);
	free(players_occurances);
	free(predictions);
	fclose(fp);
}
*/
