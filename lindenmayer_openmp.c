#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lindenmayer.h"
#include "lindenmayer_dp.h"
#include "pixmap.h"

#define INITIAL_EXPANDS 3
#define NUM_THREADS 4

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
	int n_iterations = atoi(argv[2]);
	int scale = atoi(argv[3]);
	lindenmayer_dp_entry **dp = create_lindenmayer_dp_table(&lsystem, n_iterations);
	char *initially_expanded_path = expand_lsystem(&lsystem, INITIAL_EXPANDS);
	int initially_expanded_path_len = strlen(initially_expanded_path);
	int *starting = malloc(NUM_THREADS * sizeof(int));
	lindenmayer_dp_entry *entries = malloc(NUM_THREADS * sizeof(lindenmayer_dp_entry));
	for (int i = 0; i < NUM_THREADS; ++i) {
		starting[i] = i * initially_expanded_path_len / NUM_THREADS;
	}
	lindenmayer_dp_entry info = scan_rule(&lsystem, initially_expanded_path,
		dp[n_iterations - INITIAL_EXPANDS], compute_no_of_variables(&lsystem), 1,
		entries, starting, NUM_THREADS);

	// Initialize pixmap_t
	int height = (info.max_x - info.min_x + 10) * scale;
	int width = (info.max_y - info.min_y + 10) * scale;
	initialize_pixmap(&img, width, height);

	// Draw fractal
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		int i = omp_get_thread_num();
		// Expand the string
		int end = (i + 1) * initially_expanded_path_len / NUM_THREADS;
		int len = end - starting[i];
		char *path = malloc((len + 1) * sizeof(char));
		memcpy(path, initially_expanded_path + starting[i], len);
		path[len] = '\0';
		for (int j = 0; j < n_iterations - INITIAL_EXPANDS; ++j) {
			char *tmp = path;
			path = expand_path(&lsystem, tmp);
			free(tmp);
		}
		// Draw the lines
		draw_path(&img, &lsystem, path, (-info.min_x + entries[i].x + 5) * scale,
			(-info.min_y + entries[i].y + 5) * scale, entries[i].angle, scale,
			i * strlen(path), strlen(path) * NUM_THREADS, p_coloring);
		free(path);
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
