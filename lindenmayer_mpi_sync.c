#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lindenmayer.h"
#include "lindenmayer_dp.h"
#include "pixmap.h"

#define INITIAL_EXPANDS 3

// Decomment to not write the image
// #define DONT_WRITE_IMAGE

#pragma pack(1)
typedef struct {
	pixel_t color;
	int x, y;
} mpi_pixel_t;
#pragma pack()

void expand_and_send_path(pixmap_t *p_pixmap, lindenmayer_system *p_lsystem, char *path,
               double start_x, double start_y, double start_angle, int scale,
               int previous_length, int total_length, coloring_f coloring_f)
{
	double x = start_x;
	double y = start_y;
	double angle = start_angle;
	mpi_pixel_t mpi_pixel;
	double a = x - (int)x;
	double b = y - (int)y;
	mpi_pixel.x = a <= 0.5 ? (int)x : (int)x + 1;
	mpi_pixel.y = b <= 0.5 ? (int)y : (int)y + 1;
	mpi_pixel.color = coloring_f(previous_length, total_length);
	MPI_Send(&mpi_pixel, sizeof(mpi_pixel_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
	for (int i = 0; path[i] != '\0'; ++i) {
		if (p_lsystem->is_forward[(int)path[i]]) {
			for (int j = 0; j < scale; ++j) {
				x += cos(angle);
				y += sin(angle);
				double a = x - (int)x;
				double b = y - (int)y;
				mpi_pixel.x = a <= 0.5 ? (int)x : (int)x + 1;
				mpi_pixel.y = b <= 0.5 ? (int)y : (int)y + 1;
				mpi_pixel.color = coloring_f(previous_length + i, total_length);
				MPI_Send(&mpi_pixel, sizeof(mpi_pixel_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
				// color_point(p_pixmap, x, y, pixel, blend_lighten);
			}
		} else if (path[i] == '+') {
				angle += p_lsystem->angle;
		} else if (path[i] == '-') {
				angle -= p_lsystem->angle;
		}
	}
	// Signal end
	mpi_pixel.x = -1;
	mpi_pixel.y = -1;
	MPI_Send(&mpi_pixel, sizeof(mpi_pixel_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int world_rank, world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if (world_rank == 0 && argc != 5) {
		fprintf(stderr, "Usage: %s curve_type iterations scaling coloring_type\n", argv[0]);
		fprintf(stderr, "Curve type is:\n");
		fprintf(stderr, "   0 = Dragon Curve\n");
		fprintf(stderr, "   1 = Koch Curve\n");
		fprintf(stderr, "   2 = Sierpinsky Triangle\n");
		fprintf(stderr, "   3 = Quadratic Gosper\n");
		fprintf(stderr, "   4 = Levy Curve\n");
		fprintf(stderr, "   5 = Pentaplexity\n");
		fprintf(stderr, "Coloring type is:\n");
		fprintf(stderr, "   0 = HSV coloring\n");
		fprintf(stderr, "   1 = Christmas coloring\n");
		return -1;
	}
	pixmap_t img;
	lindenmayer_system lsystem;
	coloring_f *p_coloring;

	// Chose the curve type
	switch (atoi(argv[1])) {
		case 0:
			initialize_dragon_curve(&lsystem);
			break;
		case 1:
			initialize_koch_curve(&lsystem);
			break;
		case 2:
			initialize_sierpinsky_triangle(&lsystem);
			break;
		case 3:
			initialize_quadratic_gosper(&lsystem);
			break;
		case 4:
			initialize_levy_curve(&lsystem);
			break;
		case 5:
			initialize_pentaplexity(&lsystem);
			break;
		default:
			initialize_dragon_curve(&lsystem);
	}
	// Choose the coloring type
	if (atoi(argv[4]) == 1) p_coloring = christmas_coloring;
	else p_coloring = hsv_coloring;

	// Find informations about the fractal using dynampic programming
	int n_threads = world_size - 1;
	int n_iterations = atoi(argv[2]);
	int scale = atoi(argv[3]);
	lindenmayer_dp_entry **dp = create_lindenmayer_dp_table(&lsystem, n_iterations);
	char *initially_expanded_path = expand_lsystem(&lsystem, INITIAL_EXPANDS);
	int initially_expanded_path_len = strlen(initially_expanded_path);
	int *starting = malloc(n_threads * sizeof(int));
	lindenmayer_dp_entry *entries = malloc(n_threads * sizeof(lindenmayer_dp_entry));
	for (int i = 0; i < n_threads; ++i) {
		starting[i] = i * initially_expanded_path_len / n_threads;
	}
	lindenmayer_dp_entry info = scan_rule(&lsystem, initially_expanded_path,
		dp[n_iterations - INITIAL_EXPANDS], compute_no_of_variables(&lsystem), 1,
		entries, starting, n_threads);

	// Initialize pixmap_t
	int height = (info.max_x - info.min_x + 10) * scale;
	int width = (info.max_y - info.min_y + 10) * scale;

	if (world_rank == 0) {
		initialize_pixmap(&img, width, height);
	}

	// Draw fractal
	if (world_rank == 0) {
		int have_finished = 0;
		MPI_Status status;
		while (have_finished < world_size - 1) {
			mpi_pixel_t mpi_pixel;
			MPI_Recv(&mpi_pixel, sizeof(mpi_pixel_t), MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			if (mpi_pixel.x == -1 && mpi_pixel.y == -1) ++have_finished;
			else color_point(&img, mpi_pixel.x, mpi_pixel.y, mpi_pixel.color, blend_lighten);
		}
	} else {
		// Expand the string
		int end = world_rank * initially_expanded_path_len / (world_size - 1);
		int len = end - starting[world_rank - 1];
		char *path = malloc((len + 1) * sizeof(char));
		memcpy(path, initially_expanded_path + starting[world_rank - 1], len);
		path[len] = '\0';
		for (int j = 0; j < n_iterations - INITIAL_EXPANDS; ++j) {
			char *tmp = path;
			path = expand_path(&lsystem, tmp);
			free(tmp);
		}
		expand_and_send_path(&img, &lsystem, path,
		  (-info.min_x + entries[world_rank - 1].x + 5) * scale,
		  (-info.min_y + entries[world_rank - 1].y + 5) * scale,
		  entries[world_rank - 1].angle, scale,
		  (world_rank - 1) * strlen(path), strlen(path) * n_threads, p_coloring);
		free(path);
	}

#ifndef DONT_WRITE_IMAGE
	if (world_rank == 0) {
		write_pixmap(&img, stdout);
	}
#endif
	// Free the used memory
	for (int i = 0; i <= n_iterations; ++i) free(dp[i]);
	free(dp);
	free(starting);
	free(initially_expanded_path);
	free(entries);
	clear_lsystem(&lsystem);
	if (world_rank == 0) {
		clear_pixmap(&img);
	}
	MPI_Finalize();
	return 0;
}
