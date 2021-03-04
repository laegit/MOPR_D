
/*
	This code was provided by PISA and adapted by Felipe Fernandes (2018) 
*/


/*===========================================================================*
 * hyp_ind.c: implements the unary hypervolume indicator as proposed in
 *            Zitzler, E., and Thiele, L. (1998): Multiobjective Optimization
 *            Using Evolutionary Algorithms - A Comparative Case Study.
 *            Parallel Problem Solving from Nature (PPSN-V), 292-301; a more
 *            detailed discussion can be found in
 *            Zitzler, E., Thiele, L., Laumanns, M., Fonseca, C., and
 *            Grunert da Fonseca, V (2003): Performance Assessment of
 *            Multiobjective Optimizers: An Analysis and Review. IEEE
 *            Transactions on Evolutionary Computation, 7(2), 117-132.
 *
 * Usage:
 *   hyp_ind [<param_file>] <data_file> <reference_set> <output_file>
 *
 *   <param_file> specifies the name of the parameter file for eps_ind; the
 *     file has the following format:
 *
 *       dim <integer>
 *       obj <+|-> <+|-> ...
 *       method <0|1>
 *       nadir <real> <real> ...
 *
 *     The first line defines the number of objectives, the second for each
 *     objective whether it is minimized (-) or maximized, the third
 *     line determines whether the hypervolume is calculated relative to
 *     the reference set (1) or not (0), and the last line gives the worst
 *     value for each objective (reference point).
 *     If the parameter file is omitted, the number of objectives is determined
 *     from the data file and it is assumed that all objectives are to be
 *     minimized, that the nadir point is (2.1, 2.1, ..., 2.1), and that a
 *     a reference set is given (method=1).
 *
 *   <data_file> specifies a file that contains the output of one or
 *     several runs of a selector/variator pair; the format corresponds to
 *     the one defined in the specification of the PISA monitor program.
 *
 *   <reference_set> is the name of a file that contains the reference set
 *     according to which the indicator values are calculated; the file
 *     format is the same as for the data file.
 *
 *   <output_file> defines the name of the file to which the computed
 *     indicator values are written to.
 *
 * IMPORTANT: In order to make the output of this tool consistent with
 *   the other indicator tools, for method 0 (no reference set) the
 *   negative hypervolume is outputted as indicator value. Thus,
 *   independently of which type of problem (minimization,
 *   maximization, mixed minimization/maximization) and of which type
 *   of method (with or without reference set) one considers, a lower
 *   indicator value corresponds to a better approximation set.
 *
 * Author:
 *   Eckart Zitzler, February 3, 2005 / last update August 9, 2005 */

#ifndef HYP_IND_CPP
#define HYP_IND_CPP

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define error(X,Y)  if (X) fprintf(stderr, Y "\n"), exit(1)

#define MAX_LINE_LENGTH  2048 /* maximal length of lines in the files */
#define MAX_STR_LENGTH  256 /* maximal length of strings in the files */


/*DEVE SER MINIMIZADO*/
class Hypervolume{

private:
	int  dim;  /* number of objectives */
	int  obj[NUMOBJETIVOS];  /* obj[i] = 0 means objective i is to be minimized */
	int  method;  /* 0 = no reference set, 1 = with respect to reference set */
	double nadir[NUMOBJETIVOS];;  /* reference point for hypervolume calculation */
	double *curr_run;

	int  dominates(double  *point1, double  *point2, int  no_objectives)
	    /* returns true if 'point1' dominates 'points2' with respect to the
	       to the first 'no_objectives' objectives
	    */
	{
	    int  i;
	    int  better_in_any_objective, worse_in_any_objective;
	    
	    better_in_any_objective = 0;
	    worse_in_any_objective = 0;
	    for (i = 0; i < no_objectives && !worse_in_any_objective; i++)
	      if (point1[i] > point2[i])
		better_in_any_objective = 1;
	      else if (point1[i] < point2[i])
		worse_in_any_objective = 1;
	    return (!worse_in_any_objective && better_in_any_objective);
	}

	int  weakly_dominates(double  *point1, double  *point2, int  no_objectives)
	    /* returns true if 'point1' weakly dominates 'points2' with respect to the
	       to the first 'no_objectives' objectives
	    */
	{
	    int  i;
	    int  worse_in_any_objective;
	    
	    worse_in_any_objective = 0;    
	    for (i = 0; i < no_objectives &&  !worse_in_any_objective; i++)
	      if (point1[i] < point2[i])
		worse_in_any_objective = 1;
	    return (!worse_in_any_objective);
	} 

	void  swap(double  *front, int  i, int  j)
	{
	    int  k;
	    double  temp;

	    for (k = 0; k < dim; k++) {
		temp = front[i * dim + k];
		front[i * dim + k] = front[j * dim + k];
		front[j * dim + k] = temp;
	    }
	}

	int  filter_nondominated_set(double  *front, int  no_points, int  no_objectives)
	    /* all nondominated points regarding the first 'no_objectives' dimensions
	       are collected; the points 0..no_points-1 in 'front' are
	       considered; the points in 'front' are resorted, such that points
	       [0..n-1] represent the nondominated points; n is returned
	    */
	{
	    int  i, j;
	    int  n;
	    
	    n = no_points;
	    i = 0;
	    while (i < n) {
		j = i + 1;
		while (j < n) {
		    if (dominates(&(front[i * dim]), &(front[j * dim]),
				  no_objectives)) {
			/* remove point 'j' */
			n--;
			swap(front, j, n);
		    }
		    else if (dominates(&(front[j * dim]), &(front[i * dim]),
				       no_objectives)) {
			/* remove point 'i'; ensure that the point copied to index 'i'
			   is considered in the next outer loop (thus, decrement i) */
			n--;
			swap(front, i, n);
			i--;
			break;
		    }
		    else
			j++;
		}
		i++;
	    }
	    return n;
	}

	double  surface_unchanged_to(double  *front, int  no_points, int  objective)
	     /* calculate next value regarding dimension 'objective'; consider
		points 0..no_points-1 in 'front'
	     */
	{
	  int     i;
	  double  min, value;

	  error(no_points < 1, "run-time error");
	  min = front[objective];
	  for (i = 1; i < no_points; i++) {
	    value = front[i * dim + objective];
	    if (value < min)  min = value;
	  }
	  
	  return min;
	} 

	int  reduce_nondominated_set(double  *front, int  no_points, int  objective,
				     double  threshold)
	    /* remove all points which have a value <= 'threshold' regarding the
	       dimension 'objective'; the points [0..no_points-1] in 'front' are
	       considered; 'front' is resorted, such that points [0..n-1] represent
	       the remaining points; 'n' is returned
	    */
	{
	    int  n;
	    int  i;
	    
	    n = no_points;
	    for (i = 0; i < n; i++)
		if (front[i * dim + objective] <= threshold) {
		    n--;
		    swap(front, i, n);
		}
	    
	    return n;
	} 

	double  calc_hypervolume(double  *front, int  no_points, int  no_objectives)
	{
	    int     n;
	    double  volume, distance;

	    volume = 0;
	    distance = 0;
	    n = no_points;
	    while (n > 0) {
		int     no_nondominated_points;
		double  temp_vol, temp_dist;

		no_nondominated_points = filter_nondominated_set(front, n,
								 no_objectives - 1);
		temp_vol = 0;
		if (no_objectives < 3) {
		    error(no_nondominated_points < 1, "run-time error");
		    temp_vol = front[0];
		}
		else
		    temp_vol = calc_hypervolume(front, no_nondominated_points,
						no_objectives - 1);
		temp_dist = surface_unchanged_to(front, n, no_objectives - 1);
		volume += temp_vol * (temp_dist - distance);
		distance = temp_dist;
		n = reduce_nondominated_set(front, n, no_objectives - 1, distance);
	    }
	    
	    return volume;
	}

	double  calc_ind_value(double  *a, int  size_a)
	{
	    int  i, k;
	    double  temp;
	    
	    /* re-calculate objective values relative to reference point */
	    for (i = 0; i < size_a; i++) {
	        for (k = 0; k < dim; k++) {
		    switch (obj[k]) {
		    case 0:
		        temp = nadir[k] - a[i * dim + k];
			error(temp < 0, "error in data or reference set file 4");
			a[i * dim + k] = temp;
			break;
		    default:
		        temp = a[i * dim + k] - nadir[k];
			error(temp < 0, "error in data or reference set file 3");
			a[i * dim + k] = temp;
			break;
		    }
		}
	    }
	    /* calculate indicator values */
	    return calc_hypervolume(a, size_a, dim);
	}

public:


	Hypervolume(){
		dim = NUMOBJETIVOS; 
		curr_run = (double *) malloc(dim * MAXARCSIZE * sizeof(double)); 
		for (int i=0; i<NUMOBJETIVOS; i++){
			obj[i] = 0;
			nadir[i] = 2.1;
		} 
		method = 0;  /* 0 = no reference set, 1 = with respect to reference set */
		
	}

	/*DEVE SER MINIMIZADO*/
	double getHypervolume(ParetoSet *current_position, std::vector<double> lbound, std::vector<double>  ubound){
		
		list<SolucaoEdgeSet *> sol = current_position->getElementos();
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		int cont = 0, k=0;
		while (i!=sol.end()){
			curr_run[k++] = 1.0+((double)(*i)->getObj(0)-lbound[0])/(ubound[0]-lbound[0]);
			curr_run[k++] = 1.0+((double)(*i)->getObj(1)-lbound[1])/(ubound[1]-lbound[1]);
			curr_run[k++] = 1.0+((double)(*i)->getObj(2)-lbound[2])/(ubound[2]-lbound[2]);
			// fprintf(stdout, "%.9e %.9e\n", 1.0+((double)(*i)->getObj(0)-lbound[0])/(ubound[0]-lbound[0]), 1.0+((double)(*i)->getObj(1)-lbound[1])/(ubound[1]-lbound[1]));
			cont++;
			i++;
		}
		double ind_value = calc_ind_value(curr_run, cont);
		return (-ind_value); 
		// fprintf(stdout, "Hypervolume = %.9e\n\n", -ind_value);
	}
};


#endif