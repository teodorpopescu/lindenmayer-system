#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lindenmayer.h"
#include "lindenmayer_dp.h"
#include "pixmap.h"

#define N_THREADS 4
#define INITIAL_EXPANDS 3

// Decomment to not write the image
// #define DONT_WRITE_IMAGE

void draw_path(pixmap_t *p_pixmap, lindenmayer_system *p_lsystem, char *path,
               double start_x, double start_y, double start_angle, int scale,
               int previous_length, int total_length, coloring_f coloring_f)
{
	double x = start_x;
	double y = start_y;
	double angle = start_angle;

	color_point(p_pixmap, x, y, coloring_f(previous_length, total_length), blend_lighten);
	for (int i = 0; path[i] != '\0'; ++i) {
		if (p_lsystem->is_forward[(int)path[i]]) {
			for (int j = 0; j < scale; ++j) {
				x += cos(angle);
				y += sin(angle);
				pixel_t pixel = coloring_f(previous_length + i, total_length);
				color_point(p_pixmap, x, y, pixel, blend_lighten);
			}
		} else if (path[i] == '+') {
				angle += p_lsystem->angle;
		} else if (path[i] == '-') {
				angle -= p_lsystem->angle;
		}
	}
}

typedef struct {
	int starting, ending, i;
	int n_iterations, scale;
	char *initially_expanded_path;
	lindenmayer_dp_entry *p_info, *p_entry;
	lindenmayer_system *p_lsystem;
	pixmap_t *p_pixmap;
	coloring_f *p_coloring;
} thread_info_t;

void *thread_function(void *p_info) {
	thread_info_t *p = (thread_info_t *)p_info;
	int end = p->ending;
	int len = end - p->starting;
	char *path = malloc((len + 1) * sizeof(char));
	memcpy(path, p->initially_expanded_path + p->starting, len);
	path[len] = '\0';
	for (int j = 0; j < p->n_iterations - INITIAL_EXPANDS; ++j) {
		char *tmp = path;
		path = expand_path(p->p_lsystem, tmp);
		free(tmp);
	}
	// Draw the lines
	draw_path(p->p_pixmap, p->p_lsystem, path,
		(-p->p_info->min_x + p->p_entry->x + 5) * p->scale,
		(-p->p_info->min_y + p->p_entry->y + 5) * p->scale, p->p_entry->angle, p->scale,
		p->i * strlen(path), strlen(path) * N_THREADS, p->p_coloring);
	free(path);

	return NULL;
}

int main(int argc, char *argv[])
{
	if (argc != 5) {
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
	int n_threads = N_THREADS;
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
	initialize_pixmap(&img, width, height);

	// Draw fractal
	pthread_t threads[N_THREADS];
	thread_info_t infos[N_THREADS];
	for (int i = 0; i < n_threads; ++i) {
		// Expand the string
		infos[i].starting = starting[i];
		infos[i].ending = (i + 1) * initially_expanded_path_len / n_threads;
		infos[i].i = i;
		infos[i].n_iterations = n_iterations;
		infos[i].scale = scale;
		infos[i].initially_expanded_path = initially_expanded_path;
		infos[i].p_info = &info;
		infos[i].p_entry = &entries[i];
		infos[i].p_lsystem = &lsystem;
		infos[i].p_pixmap = &img;
		infos[i].p_coloring = p_coloring;
		pthread_create(&threads[i], NULL, thread_function, &infos[i]);
	}
	for (int i = 0; i < n_threads; ++i) {
		pthread_join(threads[i], NULL);
	}
#ifndef DONT_WRITE_IMAGE
	write_pixmap(&img, stdout);
#endif

	// Free the used memory
	for (int i = 0; i <= n_iterations; ++i) free(dp[i]);
	free(dp);
	free(starting);
	free(initially_expanded_path);
	free(entries);
	clear_lsystem(&lsystem);
	clear_pixmap(&img);
	return 0;
}
