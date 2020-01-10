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


typedef struct {
	mpi_pixel_t *data;
	int size, capacity;
} mpi_pixel_vector_t;

void initialize_pixel_vector(mpi_pixel_vector_t *v)
{
	v->data = malloc(2 * sizeof(mpi_pixel_t));
	v->size = 0;
	v->capacity = 2;
}

void pixel_vector_push_back(mpi_pixel_vector_t *v, double x, double y, pixel_t color)
{
	if (v->size == v->capacity) {
		v->capacity *= 2;
		v->data = realloc(v->data, 2 * v->capacity * sizeof(mpi_pixel_t));
	}
	int size = v->size;
	double a = x - (int)x;
	double b = y - (int)y;
	v->data[size].x = a <= 0.5 ? (int)x : (int)x + 1;
	v->data[size].y = b <= 0.5 ? (int)y : (int)y + 1;
	v->data[size].color = color;
	++v->size;
}

mpi_pixel_vector_t expand_and_send_path(pixmap_t *p_pixmap, lindenmayer_system *p_lsystem,
               char *path,
               double start_x, double start_y, double start_angle, int scale,
               int previous_length, int total_length, coloring_f coloring_f)
{
	double x = start_x;
	double y = start_y;
	double angle = start_angle;
	mpi_pixel_vector_t v;
	initialize_pixel_vector(&v);
	pixel_vector_push_back(&v, x, y, coloring_f(previous_length, total_length));
	for (int i = 0; path[i] != '\0'; ++i) {
		if (p_lsystem->is_forward[(int)path[i]]) {
			for (int j = 0; j < scale; ++j) {
				x += cos(angle);
				y += sin(angle);
				pixel_vector_push_back(&v, x, y, coloring_f(previous_length + i, total_length));
			}
		} else if (path[i] == '+') {
				angle += p_lsystem->angle;
		} else if (path[i] == '-') {
				angle -= p_lsystem->angle;
		}
	}
	return v;
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
	int n_iterations = atoi(argv[2]);
	int scale = atoi(argv[3]);
	lindenmayer_dp_entry **dp = create_lindenmayer_dp_table(&lsystem, n_iterations);
	char *initially_expanded_path = expand_lsystem(&lsystem, INITIAL_EXPANDS);
	int initially_expanded_path_len = strlen(initially_expanded_path);
	int *starting = malloc(world_size * sizeof(int));
	lindenmayer_dp_entry *entries = malloc(world_size * sizeof(lindenmayer_dp_entry));
	for (int i = 0; i < world_size; ++i) {
		starting[i] = i * initially_expanded_path_len / world_size;
	}
	lindenmayer_dp_entry info = scan_rule(&lsystem, initially_expanded_path,
		dp[n_iterations - INITIAL_EXPANDS], compute_no_of_variables(&lsystem), 1,
		entries, starting, world_size);

	// Initialize pixmap_t
	int height = (info.max_x - info.min_x + 10) * scale;
	int width = (info.max_y - info.min_y + 10) * scale;

	if (world_rank == 0) {
		initialize_pixmap(&img, width, height);
	}

	// Expand the string
	int end = (world_rank + 1) * initially_expanded_path_len / world_size;
	int len = end - starting[world_rank];
	char *path = malloc((len + 1) * sizeof(char));
	memcpy(path, initially_expanded_path + starting[world_rank], len);
	path[len] = '\0';
	for (int j = 0; j < n_iterations - INITIAL_EXPANDS; ++j) {
		char *tmp = path;
		path = expand_path(&lsystem, tmp);
		free(tmp);
	}
	mpi_pixel_vector_t v = expand_and_send_path(&img, &lsystem, path,
	  (-info.min_x + entries[world_rank].x + 5) * scale,
	  (-info.min_y + entries[world_rank].y + 5) * scale,
	  entries[world_rank].angle, scale,
	  world_rank * strlen(path), strlen(path) * world_size, p_coloring);
	free(path);
	if (world_rank == 0) {
		for (int i = 0; i < v.size; ++i) {
			color_point(&img, v.data[i].x, v.data[i].y, v.data[i].color, blend_lighten);
		}
		MPI_Status status;
		int buff_size = 2 * width * height / world_size * sizeof(mpi_pixel_t);
		int count;
		mpi_pixel_t *w = malloc(buff_size);
		for (int k = 1; k < world_size; ++k) {
			mpi_pixel_t mpi_pixel;
			MPI_Recv(w, buff_size, MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_BYTE, &count);
			count /= sizeof(mpi_pixel_t);
			for (int i = 0; i < count; ++i) {
				color_point(&img, w[i].x, w[i].y, w[i].color, blend_lighten);
			}
		}
		free(w);
		write_pixmap(&img, stdout);
	} else {
		MPI_Send(v.data, v.size * sizeof(mpi_pixel_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
	}
	free(v.data);

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
