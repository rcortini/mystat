#include <string.h>
#include <unistd.h>
#include "mylib.h"

const char *program_name = "mystat";

/* if columns were not specified, set defaults */
void init_cols (const int cflag, const size_t ncols, size_t **cols) {
  if (!cflag) {
    unsigned int i;
    *cols = (size_t *) malloc (ncols*sizeof (unsigned int));
    for (i=0; i<ncols; i++)
      (*cols) [i] = i;
  }
}


void print_usage () {
  printf ("Usage: %s [OPTIONS] <filename> <command> ...\n", program_name);
  printf ("Commands:\n");
  printf ("\taverage\n");
  printf ("\tvariance\n");
  printf ("\tdevst\n");
  printf ("\taverage_variance\n");
  printf ("\taverage_devst\n");
  printf ("\thistogram <bin_min> <bin_max> <n_bins>\n");
  printf ("\tlinear_fit\n");
  printf ("\tweighted_linear_fit\n");
  printf ("\tpolynomial_fit\n");
  printf ("\tblock_average <nblocks>\n");
  printf ("\tfft\n");
  printf ("\tpsd\n");
  printf ("\tcorrelation\n");
  printf ("Options:\n");
  printf ("\t-h : Print this help and exit\n");
  printf ("\t-v : Verbose\n");
  printf ("\t-c <n> : Column number of file\n");
}

int main (int argc, char *argv[]) {
  int c, cflag = 0, vflag = 0;
  size_t ncols, *cols;
  unsigned int retcode;
  char *filename, *command;
  FILE *f_in;
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
  (void) old_handler;

  /* parse command line options */
  while ((c = getopt (argc, argv, "c:h::v::")) != -1) {
    switch (c) {
      case 'h' :
	print_usage (program_name);
	exit (0);
	break;
      case 'v' :
	vflag = 1;
	break;
      case 'c' :
	cflag = 1;
	ncols = parse_ranges (optarg, &cols);
	break;
      default :
	err_message ("Incorrect usage\n");
	print_usage (program_name);
	exit (EXIT_FAILURE);
    }
  };

  /* check for proper invocation */
  if (argc-optind < 2) {
    err_message ("Incorrect usage\n");
    print_usage (program_name);
    exit (EXIT_FAILURE);
  }

  /* the first non-option argument is the filename to read */
  filename = argv [optind];
  if (strcmp (filename, "-")==0)
    f_in = stdin;
  else
    f_in = safe_fopen (filename, "r");

  /* the second non-option argument is the command name to execute */
  command = argv [optind+1];

  /* now check what is the command and act */
  if (strcmp (command, "average_variance")==0) {
    double **data;
    unsigned int i, N;

    init_cols (cflag, 1, &cols);
    if (!cflag) ncols = 1;
    N = read_data (f_in, ncols, cols, &data);
    for (i=0; i<ncols; i++) {
      double av, var;
      retcode = average_variance (N, data [i], &av, &var);
      printf ("%.8e %.8e ", av, var);
      if (retcode)
	break;
      free (data [i]);
    }
    printf ("\n");
    free (data);
  }
  else if (strcmp (command, "average_devst")==0) {
    double av, ds;
    double **data;
    unsigned int i, N;

    init_cols (cflag, 1, &cols);
    if (!cflag) ncols = 1;
    N = read_data (f_in, ncols, cols, &data);
    for (i=0; i<ncols; i++) {
      retcode = average_devst (N, data [i], &av, &ds);
      printf ("%.8e %.8e ", av, ds);
      if (retcode)
	break;
      free (data [i]);
    }
    printf ("\n");
    free (data);
  }
  else if (strcmp (command, "average")==0) {
    double av;
    double **data;
    unsigned int i, N;

    init_cols (cflag, 1, &cols);
    if (!cflag) ncols = 1;
    N = read_data (f_in, ncols, cols, &data);
    for (i=0; i<ncols; i++) {
      retcode = average (N, data [i], &av);
      printf ("%.8e ", av);
      if (retcode)
	break;
      free (data [i]);
    }
    printf ("\n");
    free (data);
  }
  else if (strcmp (command, "devst")==0) {
    double ds;
    double **data;
    unsigned int i, N;

    init_cols (cflag, 1, &cols);
    if (!cflag) ncols = 1;
    N = read_data (f_in, ncols, cols, &data);
    for (i=0; i<ncols; i++) {
      retcode = devst (N, data [i], &ds);
      printf ("%.8e ", ds);
      if (retcode)
	break;
      free (data [i]);
    }
    printf ("\n");
    free (data);
  }
  else if (strcmp (command, "variance")==0) {
    double var;
    double **data;
    unsigned int i, N;

    init_cols (cflag, 1, &cols);
    if (!cflag) ncols = 1;
    N = read_data (f_in, ncols, cols, &data);
    for (i=0; i<ncols; i++) {
      retcode = variance (N, data [i], &var);
      printf ("%.8e ", var);
      if (retcode)
	break;
      free (data [i]);
    }
    printf ("\n");
    free (data);
  }
  else if (strcmp (command, "histogram")==0) {
    unsigned int i, N, argstart=optind+2, nbins;
    double bin_min, bin_max, sum;
    double **data;
    gsl_histogram *hist;

    /* get input parameters */
    if (argc-optind < 5) {
      err_message ("Insufficient arguments for histogram\n");
      print_usage ();
      fclose (f_in);
      exit (EXIT_FAILURE);
    }
    bin_min = atof (argv [argstart]);
    bin_max = atof (argv [argstart+1]);
    nbins = atoi (argv [argstart+2]);

    /* read data */
    init_cols (cflag, 1, &cols);
    N = read_data (f_in, 1, cols, &data);

    /* calculate histogram */
    retcode = histogram (N, *data, nbins, bin_min, bin_max, &hist);

    /* get total number of values in the histogram */
    sum = gsl_histogram_sum (hist);

    /* output */
    for (i=0; i<nbins; i++) {
      double lower, upper, mean;

      /* get the range of the i-th bin in the histogram, and calculate mean */
      gsl_histogram_get_range (hist, i, &lower, &upper);
      mean = (upper+lower)/2.;

      /* write line in file */
      printf ("%f %f\n", mean, gsl_histogram_get (hist, i)/sum);
    }

    /* free memory */
    gsl_histogram_free (hist);
    free (*data);
    free (data);
  }
  else if (strcmp (command, "weighted_linear_fit")==0) {
    double **data;
    unsigned int i, N;
    double *w, *x, *y;
    linear_fit_results fit_results;

    /* check that we have a sufficient number of columns */
    if (ncols < 3) {
      err_message ("Insufficient number of columns: 3 needed, %u given\n", ncols);
      exit (EXIT_FAILURE);
    }

    init_cols (cflag, 3, &cols);
    N = read_data (f_in, 3, cols, &data);
    w = (double *) malloc (N*sizeof (double));

    /* check the values of the input sigmas  and assign weights*/
    for (i=0; i<N; i++) {
      double s = data [2][i];
      if (s==0) {
	err_message ("Zero variance point encountered in weighted linear fit!\n");
	return MYLIB_FAIL;
      }
      w [i] = 1./(s*s);
    }
    x = data[0];
    y = data[1];

    /* do the fit */
    weighted_linear_fit (N, x, y, w, &fit_results);
    retcode = fit_results.retcode;

    /* print result */
    print_linear_fit_results (&fit_results, vflag);

    /* success */
    free (data[0]);
    free (data[1]);
    free (data[2]);
    free (w);
    free (data);
  }
  else if (strcmp (command, "linear_fit")==0) {
    double **data;
    unsigned int N;
    double *x, *y;
    linear_fit_results fit_results;

    /* check that we have a sufficient number of columns */
    if (cflag)
      if (ncols < 2) {
	err_message ("Insufficient number of columns: 2 needed, %u given\n", ncols);
	exit (EXIT_FAILURE);
      }

    init_cols (cflag, 2, &cols);
    N = read_data (f_in, 2, cols, &data);
    x = data[0];
    y = data[1];

    /* do the fit */
    linear_fit (N, x, y, &fit_results);
    retcode = fit_results.retcode;

    /* print result */
    print_linear_fit_results (&fit_results, vflag);

    /* success */
    free (data[0]);
    free (data[1]);
    free (data);
  }
  else if (strcmp (command, "polynomial_fit")==0) {
    unsigned int N, degree, argstart=optind+2;
    double **data;
    double *x, *y;
    double x0;
    multifit_results *fit_results;

    /* get input parameters */
    if (argc-optind < 4) {
      err_message ("Insufficient arguments for polynomial fit\n");
      print_usage ();
      exit (EXIT_FAILURE);
    }
    x0 = atof (argv [argstart]);
    degree = atoi (argv [argstart+1]);
    fit_results = multifit_results_alloc (degree);

    /* check that we have a sufficient number of columns */
    init_cols (cflag, 2, &cols);
    if (cflag)
      if (ncols < 2) {
	err_message ("Insufficient number of columns: 2 needed, %u given\n", ncols);
	exit (EXIT_FAILURE);
      }

    /* read data from input */
    N = read_data (f_in, 2, cols, &data);
    x = data[0];
    y = data[1];

    /* do the fit */
    polynomial_fit (N, x, y, degree, x0, fit_results);
    retcode = fit_results->retcode;

    /* print result */
    print_multifit_results (fit_results, vflag);

    /* success */
    free (data[0]);
    free (data[1]);
    free (data);
    multifit_results_free (fit_results);
  }
  else if (strcmp (command, "block_average")==0) {
    unsigned int argstart=optind+2;
    double **data;
    unsigned int N, nblocks;
    block_average_results *ba_results;

    /* get input parameters */
    if (argc-optind < 3) {
      err_message ("Insufficient arguments for block averaging\n");
      print_usage ();
      exit (EXIT_FAILURE);
    }
    nblocks = atoi (argv [argstart]);

    /* allocate memory */
    ba_results = block_average_results_alloc (nblocks);

    /* read data */
    init_cols (cflag, 1, &cols);
    N = read_data (f_in, 1, cols, &data);

    /* do the block averaging */
    retcode = block_average (N, *data, nblocks, ba_results);

    /* print results */
    print_block_average_results (ba_results, vflag);

    /* success */
    free (*data);
    free (data);
    block_average_results_free (ba_results);
  }
  else if (strcmp (command, "fft")==0) {
    double **data, *fft_results;
    unsigned int N;

    /* read data */
    init_cols (cflag, 1, &cols);
    N = read_data (f_in, 1, cols, &data);
    fft_results = (double *) malloc (N*sizeof (double));

    /* perform the fft */
    retcode = fft (N, *data, fft_results);

    /* print results */
    print_fft_results (N, fft_results, vflag);

    /* free memory */
    free (*data);
    free (data);
    free (fft_results);
  }
  else if (strcmp (command, "psd")==0) {
    double **data, *psd_results;
    unsigned int i, N;

    /* read data */
    init_cols (cflag, 1, &cols);
    N = read_data (f_in, 1, cols, &data);
    psd_results = (double *) malloc (N*sizeof (double));

    retcode = psd (N, *data, psd_results);
    for (i=0; i<N; i++)
      printf ("%f\n", psd_results [i]);
    free (psd_results);
    free (*data);
    free (data);
  }
  else if (strcmp (command, "correlation")==0) {
    double **data, *corr;
    unsigned int i, N;

    /* read data */
    init_cols (cflag, 1, &cols);
    N = read_data (f_in, 1, cols, &data);
    corr = (double *) malloc (2*N*sizeof (double));
    retcode = correlation (N, *data, corr);
    for (i=0; i<N/2; i++)
      printf ("%f\n", corr [2*i]);
    free (corr);
    free (*data);
    free (data);
  }
  else {
    err_message ("Invalid command '%s'\n", command);
    print_usage ();
    exit (EXIT_FAILURE);
  }

  /* close input file */
  fclose (f_in);
  free (cols);

  return retcode;
}
