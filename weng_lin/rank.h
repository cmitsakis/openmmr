void update_BT_full(double beta, 
		    int num_teams, int *rank, 
		    int *team_num_players, 
		    double **mu,
		    double **sigma, int *predictions);

void update_BT_partial(double beta, 
		       int num_teams, int *rank, 
		       int *rank_idx,
		       int *team_num_players, 
		       double **mu, double **sigma, 
		       int *predictions);

void update_PL(double beta, 
	       int num_teams, int *rank, 
	       int *team_num_players, 
	       double **mu,
	       double **sigma, int *predictions);

void update_TM_full(double beta, 
		    int num_teams, int *rank, 
		    int *team_num_players, 
		    double **mu,
		    double **sigma, int *predictions);
