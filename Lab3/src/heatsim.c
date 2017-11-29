/*
 * heatsim.c
 *
 *  Created on: 2011-11-17
 *      Author: francis
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>

#include "config.h"
#include "part.h"
#include "grid.h"
#include "cart.h"
#include "image.h"
#include "heat.h"
#include "memory.h"
#include "util.h"

#define PROGNAME "heatsim"
#define DEFAULT_OUTPUT_PPM "heatsim.png"
#define DEFAULT_DIMX 1
#define DEFAULT_DIMY 1
#define DEFAULT_ITER 100
#define MAX_TEMP 1000.0
#define DIM_2D 2

typedef struct ctx {
	cart2d_t *cart;
	grid_t *global_grid;
	grid_t *curr_grid;
	grid_t *next_grid;
	grid_t *heat_grid;
	int numprocs;
	int rank;
	MPI_Comm comm2d;
	FILE *log;
	int verbose;
	int dims[DIM_2D];
	int isperiodic[DIM_2D];
	int coords[DIM_2D];
	int reorder;
	int north_peer;
	int south_peer;
	int east_peer;
	int west_peer;
	MPI_Datatype vector;
} ctx_t;

typedef struct command_opts {
	int dimx;
	int dimy;
	int iter;
	char *input;
	char *output;
	int verbose;
} opts_t;

static opts_t *global_opts = NULL;

__attribute__((noreturn))
static void usage(void) {
	fprintf(stderr, PROGNAME " " VERSION " " PACKAGE_NAME "\n");
	fprintf(stderr, "Usage: " PROGNAME " [OPTIONS] [COMMAND]\n");
	fprintf(stderr, "\nOptions:\n");
	fprintf(stderr, "  --help	this help\n");
	fprintf(stderr, "  --iter	number of iterations to perform\n");
	fprintf(stderr, "  --dimx	2d decomposition in x dimension\n");
	fprintf(stderr, "  --dimy	2d decomposition in y dimension\n");
	fprintf(stderr, "  --input  png input file\n");
	fprintf(stderr, "  --output ppm output file\n");
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

static void dump_opts(struct command_opts *opts) {
	printf("%10s %s\n", "option", "value");
	printf("%10s %d\n", "dimx", opts->dimx);
	printf("%10s %d\n", "dimy", opts->dimy);
	printf("%10s %d\n", "iter", opts->iter);
	printf("%10s %s\n", "input", opts->input);
	printf("%10s %s\n", "output", opts->output);
	printf("%10s %d\n", "verbose", opts->verbose);
}

void default_int_value(int *val, int def) {
	if (*val == 0)
		*val = def;
}

static int parse_opts(int argc, char **argv, struct command_opts *opts) {
	int idx;
	int opt;
	int ret = 0;

	struct option options[] = { { "help", 0, 0, 'h' },
			{ "iter", 1, 0, 'r' }, { "dimx", 1, 0, 'x' }, { "dimy",
					1, 0, 'y' }, { "input", 1, 0, 'i' }, {
					"output", 1, 0, 'o' }, { "verbose", 0,
					0, 'v' }, { 0, 0, 0, 0 } };

	memset(opts, 0, sizeof(struct command_opts));

	while ((opt = getopt_long(argc, argv, "hvx:y:l:", options, &idx)) != -1) {
		switch (opt) {
		case 'r':
			opts->iter = atoi(optarg);
			break;
		case 'y':
			opts->dimy = atoi(optarg);
			break;
		case 'x':
			opts->dimx = atoi(optarg);
			break;
		case 'i':
			if (asprintf(&opts->input, "%s", optarg) < 0)
				goto err;
			break;
		case 'o':
			if (asprintf(&opts->output, "%s", optarg) < 0)
				goto err;
			break;
		case 'h':
			usage();
			break;
		case 'v':
			opts->verbose = 1;
			break;
		default:
			printf("unknown option %c\n", opt);
			ret = -1;
			break;
		}
	}

	/* default values*/
	default_int_value(&opts->iter, DEFAULT_ITER);
	default_int_value(&opts->dimx, DEFAULT_DIMX);
	default_int_value(&opts->dimy, DEFAULT_DIMY);
	if (opts->output == NULL)
		if (asprintf(&opts->output, "%s", DEFAULT_OUTPUT_PPM) < 0)
			goto err;
	if (opts->input == NULL) {
		fprintf(stderr, "missing input file");
		goto err;
	}

	if (opts->dimx == 0 || opts->dimy == 0) {
		fprintf(stderr,
				"argument error: dimx and dimy must be greater than 0\n");
		ret = -1;
	}

	if (opts->verbose)
		dump_opts(opts);
	global_opts = opts;
	return ret;
	err:
	FREE(opts->input);
	FREE(opts->output);
	return -1;
}

FILE *open_logfile(int rank) {
	char str[255];
	sprintf(str, "out-%d", rank);
	FILE *f = fopen(str, "w+");
	return f;
}

ctx_t *make_ctx() {
	ctx_t *ctx = (ctx_t *) calloc(1, sizeof(ctx_t));
	return ctx;
}

void free_ctx(ctx_t *ctx) {
	if (ctx == NULL)
		return;
	free_grid(ctx->global_grid);
	free_grid(ctx->curr_grid);
	free_grid(ctx->next_grid);
	free_grid(ctx->heat_grid);
	free_cart2d(ctx->cart);
	if (ctx->log != NULL) {
		fflush(ctx->log);
		fclose(ctx->log);
	}
	FREE(ctx);
}

int init_ctx(ctx_t *ctx, opts_t *opts) {
	TODO("lab3");
	MPI_Comm_size(MPI_COMM_WORLD, &ctx->numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &ctx->rank);

	if (opts->dimx * opts->dimy != ctx->numprocs) {
		fprintf(stderr,
				"2D decomposition blocks must equal number of process\n");
		goto err;
	}
	ctx->log = open_logfile(ctx->rank);
	ctx->verbose = opts->verbose;
	ctx->dims[0] = opts->dimx;
	ctx->dims[1] = opts->dimy;
	ctx->isperiodic[0] = 1;
	ctx->isperiodic[1] = 1;
	ctx->reorder = 0;
	grid_t *new_grid = NULL;

	/* FIXME: create 2D cartesian communicator */
	MPI_Cart_create(MPI_COMM_WORLD, DIM_2D, ctx->dims, ctx->isperiodic, ctx->reorder, &ctx->comm2d);
	
	MPI_Cart_shift(ctx->comm2d, 1, 1, &ctx->north_peer, &ctx->south_peer);
	MPI_Cart_shift(ctx->comm2d, 0, 1, &ctx->west_peer, &ctx->east_peer);
	
	MPI_Cart_coords(ctx->comm2d, ctx->rank, DIM_2D, ctx->coords);

	/*
	 * FIXME: le processus rank=0 charge l'image du disque
	 * et transfert chaque section aux autres processus
	 */
	if( ctx->rank == 0 )
	{
		/* load input image */
		image_t *image = load_png(opts->input);
		if (image == NULL)
			goto err;

		/* select the red channel as the heat source */
		ctx->global_grid = grid_from_image(image, CHAN_RED);

		/* grid is normalized to one, multiply by MAX_TEMP */
		grid_multiply(ctx->global_grid, MAX_TEMP);

		/* 2D decomposition */
		ctx->cart = make_cart2d(ctx->global_grid->width,
				ctx->global_grid->height, opts->dimx, opts->dimy);
		cart2d_grid_split(ctx->cart, ctx->global_grid);

		/*
		* FIXME: send grid dimensions and data
		* Comment traiter le cas de rank=0 ?
		*/
		const n_procs = ctx->numprocs -1;  // number of processess but that of rank 0 (current)
		MPI_Request req[4 * n_procs];
		MPI_Status status[4 * n_procs];
		for (int rank=1; rank<ctx->numprocs; ++rank)
		{
			//
			int coordinates[DIM_2D];
			MPI_Cart_coords(ctx->comm2d, rank, DIM_2D, coordinates);
			grid_t *buf = cart2d_get_grid(ctx->cart, coordinates[0], coordinates[1]);

			// send grid dimensions
			int shift = rank - 1;
			MPI_Send(&buf->width, 1, MPI_INTEGER, rank, 0, ctx->comm2d, &req[shift+0]);
			MPI_Send(&buf->height, 1, MPI_INTEGER, rank, 1, ctx->comm2d, &req[shift+1]);
			MPI_Send(&buf->padding, 1, MPI_INTEGER, rank, 2, ctx->comm2d, &req[shift+2]);

			// send grid data
			MPI_Send(buf->dbl, buf->height * buf->width, MPI_INTEGER, rank, 3, ctx->comm2d, &req[shift+3]);
		}
		free(req);
		free(status);

		/*
		* FIXME: receive dimensions of the grid
		* store into new_grid
		*/
		int coordinates[DIM_2D];
		MPI_Cart_coords(ctx->comm2d, 0, DIM_2D, coordinates);
		new_grid = cart2d_get_grid(ctx->cart, coordinates[0], coordinates[1]);

	}else{

		MPI_Request req[4];
		MPI_Status status[4];

		int width, height, padding;
		MPI_Recv(&width, 1, MPI_INTEGER, rank, 0, ctx->comm2d, &req[0]);
		MPI_Recv(&height, 1, MPI_INTEGER, rank, 1, ctx->comm2d, &req[1]);
		MPI_Recv(&padding, 1, MPI_INTEGER, rank, 2, ctx->comm2d, &req[2]);

		//
		new_grid = make_grid(width, height, padding);
		MPI_Recv(new_grid->dbl, height * width, MPI_INTEGER, rank, 3, ctx->comm2d, &req[3]);

	}

	/* Utilisation temporaire de global_grid */
	//new_grid = ctx->global_grid;

	if (new_grid == NULL)
		goto err;
	/* set padding required for Runge-Kutta */
	ctx->curr_grid = grid_padding(new_grid, 1);
	ctx->next_grid = grid_padding(new_grid, 1);
	ctx->heat_grid = grid_padding(new_grid, 1);
	free_grid(new_grid);

	/* FIXME: create type vector to exchange columns */
	MPI_Type_vector(ctx->curr_grid->ph, ctx->curr_grid->pw, ctx->curr_grid->padding, MPI_DOUBLE, ctx->vector);
	MPI_Type_commit(ctx->vector);

	return 0;
	err: return -1;
}

void dump_ctx(ctx_t *ctx) {
	fprintf(ctx->log, "*** CONTEXT ***\n");
	fprintf(ctx->log, "rank=%d\n", ctx->rank);
	fprintf(ctx->log, "north=%d south=%d west=%d east=%d \n",
			ctx->north_peer, ctx->south_peer,
			ctx->east_peer, ctx->west_peer);
	fprintf(ctx->log, "***************\n");
}

void exchng2d(ctx_t *ctx) {
	/*
	 *  FIXME: Echanger les bordures avec les voisins
	 * 4 echanges doivent etre effectues
	 */
	
	grid_t *grid = ctx->next_grid;
	int width = grid->pw;
	int height = grid->ph;
	int *data = grid->dbl;
	MPI_Comm comm = ctx->comm2d;

	//
	MPI_Request req[8];
	MPI_Status status[8];

	// receive mpi message
	int *offset_recv_north = data;
	int *offset_recv_south = data + (height -1) * width;
	int *offset_recv_east = data + (width -1);
	int *offset_recv_west = data;
	
	MPI_Irecv(offset_recv_north, width, MPI_DOUBLE, ctx->north_peer, 0, comm1d, &req[0]);
	MPI_Irecv(offset_recv_south, width, MPI_DOUBLE, ctx->south_peer, 1, comm1d, &req[1]);
	MPI_Irecv(offset_recv_east, 1, ctx->vector, ctx->east_peer, 2, comm1d, &req[2]);
	MPI_Irecv(offset_recv_west, 1, ctx->vector, ctx->west_peer, 3, comm1d, &req[3]);

	// send mpi message
	int *offset_send_north = data + width;
	int *offset_send_south = data + (height -2) * width;
	int *offset_send_east = data + (width -2);
	int *offset_send_west = data + 1;
	
	MPI_Irecv(offset_send_north, width, MPI_DOUBLE, ctx->north_peer, 0, comm1d, &req[4]);
	MPI_Irecv(offset_send_south, width, MPI_DOUBLE, ctx->south_peer, 1, comm1d, &req[5]);
	MPI_Irecv(offset_send_east, 1, ctx->vector, ctx->east_peer, 2, comm1d, &req[6]);
	MPI_Irecv(offset_send_west, 1, ctx->vector, ctx->west_peer, 3, comm1d, &req[7]);

	MPI_Waitall(8, req, status);
	
}

int gather_result(ctx_t *ctx, opts_t *opts) {
	
	int ret = 0;
	grid_t *local_grid = grid_padding(ctx->next_grid, 0);
	if (local_grid == NULL)
		goto err;

	MPI_Comm comm = ctx->comm2d;
	
	/*
	 * FIXME: transfer simulation results from all process to rank=0
	 * use grid for this purpose
	 */
	if(ctx->rank != 0) {
		int size_grid = local_grid->height*local_grid->width;
		MPI_Send(local_grid->dbl, size_grid, MPI_DOUBLE, 0, 0, comm);
	} else {
		int rank_process;
        for(rank_process = 0; rank_process < ctx->numprocs; ++rank_process) {
			int coords_process[DIM_2D];
            MPI_Cart_coords(comm, rank_process, DIM_2D, coords_process);

			grid_t *tmp_grid = cart2d_get_grid(ctx->cart, coords_process[0], coords_process[1]);
            if(rank_process != 0) 
            {
                int size_grid = local_grid->height*local_grid->width;
                MPI_Recv(tmp_grid->dbl, size_grid, MPI_DOUBLE, rank_process, 0, comm, MPI_STATUS_IGNORE);
            }
            else
            {
                grid_copy(local_grid, tmp_grid);
            }
		}
	} 

	/* now we can merge all data blocks, reuse global_grid */
	if (ctx->rank == 0)
        cart2d_grid_merge(ctx->cart, ctx->global_grid);

	done: free_grid(local_grid);
	return ret;
	err: ret = -1;
	goto done;
}

int main(int argc, char **argv) {
	ctx_t *ctx = NULL;
	int rep, ret;
	opts_t opts;

	if (parse_opts(argc, argv, &opts) < 0) {
		printf("Error while parsing arguments\n");
		usage();
	}
	if (opts.verbose)
		dump_opts(&opts);

	MPI_Init(&argc, &argv);

	ctx = make_ctx();
	if (init_ctx(ctx, &opts) < 0)
		goto err;
	if (opts.verbose)
		dump_ctx(ctx);

	if (ctx->verbose) {
		fprintf(ctx->log, "heat grid\n");
		fdump_grid(ctx->heat_grid, ctx->log);
	}

	for (rep = 0; rep < opts.iter; rep++) {
		if (ctx->verbose) {
			fprintf(ctx->log, "iter %d\n", rep);
			fprintf(ctx->log, "start\n");
			fdump_grid(ctx->curr_grid, ctx->log);
		}

		grid_set_min(ctx->heat_grid, ctx->curr_grid);
		if (ctx->verbose) {
			fprintf(ctx->log, "grid_set_min\n");
			fdump_grid(ctx->curr_grid, ctx->log);
		}

		exchng2d(ctx);
		if (ctx->verbose) {
			fprintf(ctx->log, "exchng2d\n");
			fdump_grid(ctx->curr_grid, ctx->log);
		}

		heat_diffuse(ctx->curr_grid, ctx->next_grid);
		if (ctx->verbose) {
			fprintf(ctx->log, "heat_diffuse\n");
			fdump_grid(ctx->next_grid, ctx->log);
		}
		SWAP(ctx->curr_grid, ctx->next_grid);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (gather_result(ctx, &opts) < 0)
		goto err;

	if (ctx->rank == 0) {
		printf("saving...\n");
		if (save_grid_png(ctx->global_grid, opts.output) < 0) {
			printf("saving failed\n");
			goto err;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	ret = EXIT_SUCCESS;
done:
	free_ctx(ctx);
	MPI_Finalize();
	FREE(opts.input);
	FREE(opts.output);
	return ret;
err:
	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	ret = EXIT_FAILURE;
	goto done;
}

