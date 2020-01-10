#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lindenmayer.h"
#include "lindenmayer_dp.h"
#include "pixmap.h"

// Decomment to not write the image
// #define DONT_WRITE_IMAGE

void draw_path(pixmap_t *p_pixmap, lindenmayer_system *p_lsystem, char *path,
               double start_x, double start_y, double start_angle, int scale,
               coloring_f coloring_f)
{
	double x = start_x;
	double y = start_y;
	double angle = start_angle;
	int path_len = strlen(path);

	color_point(p_pixmap, x, y, coloring_f(0, path_len), blend_lighten);
	for (int i = 0; path[i] != '\0'; ++i) {
		if (p_lsystem->is_forward[(int)path[i]]) {
			for (int j = 0; j < scale; ++j) {
				x += cos(angle);
				y += sin(angle);
				pixel_t pixel = coloring_f(i, path_len);
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
	lindenmayer_dp_entry info = scan_rule(&lsystem, lsystem.start, dp[n_iterations],
		compute_no_of_variables(&lsystem), 1, NULL, NULL, 0);

	// Draw the fractal
	char *path = expand_lsystem(&lsystem, n_iterations);
	int height = (info.max_x - info.min_x + 10) * scale;
	int width = (info.max_y - info.min_y + 10) * scale;
	initialize_pixmap(&img, width, height);
	draw_path(&img, &lsystem, path, (-info.min_x + 5) * scale, (-info.min_y + 5) * scale, 0, scale, p_coloring);

#ifndef DONT_WRITE_IMAGE
	write_pixmap(&img, stdout);
#endif

	// Free the used memory
	for (int i = 0; i <= n_iterations; ++i) free(dp[i]);
	free(dp);
	clear_lsystem(&lsystem);
	free(path);
	clear_pixmap(&img);
	return 0;
}
