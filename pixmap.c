#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pixmap.h"

pixel_t hsv_coloring(double x, double period)
{
	pixel_t ans = {0, 0, 0};
	double h = x / period * 360;
	double tmp = 1.0 - fabs(fmod(h / 60, 2) - 1.0);
	if (h >= 0 && h < 60.0) {
		ans.r = 255; ans.g = tmp * 255;
	} else if (h < 120.0) {
		ans.r = tmp * 255; ans.g = 255;
	} else if (h < 160.0) {
		ans.g = 255; ans.b = tmp * 255;
	} else if (h  < 240.0) {
		ans.g = tmp * 255; ans.b = 255;
	} else if (h < 300.0) {
		ans.r = tmp * 255; ans.b = 255;
	} else if (h <= 360.0) {
		ans.r = 255; ans.b = tmp * 255;
	}
	return ans;
}

pixel_t christmas_coloring(double x, double period)
{
	pixel_t ans = {255, 255, 255};
	double value = fmod(x, period / 3) * 3 / period;
	if (value <= 0.5) {
		ans.g = (1.0 - 2 * value) * 255;
		ans.b = (1.0 - 2 * value) * 255;
	} else {
		ans.g = (2 * value - 1.0) * 255;
		ans.b = (2 * value - 1.0) * 255;
	}
	return ans;
}

pixel_t blend_overlay(pixel_t a, pixel_t b) {
	pixel_t ans;
	ans.r = a.r < 128 ? 2 * a.r * b.r / 255 : 255 - 2 * (255 - a.r) * (255 - b.r) / 255;
	ans.g = a.g < 128 ? 2 * a.g * b.g / 255 : 255 - 2 * (255 - a.g) * (255 - b.g) / 255;
	ans.b = a.b < 128 ? 2 * a.b * b.b / 255 : 255 - 2 * (255 - a.b) * (255 - b.b) / 255;
	return ans;
}

pixel_t blend_lighten(pixel_t a, pixel_t b) {
	pixel_t ans;
	ans.r = a.r < b.r ? b.r : a.r;
	ans.g = a.g < b.g ? b.g : a.g;
	ans.b = a.b < b.b ? b.b : a.b;
	return ans;
}

int initialize_pixmap(pixmap_t *p_pixmap, int width, int height)
{
	if (width <= 0 || height <= 0) {
		fprintf(stderr, "ERROR: Invalid height or width for pixmap.\n");
		return PIXMAP_ERROR;
	}

	p_pixmap->width = width;
	p_pixmap->height = height;

	/* Allocate the pixel matrix, catching the allocation errors */
	p_pixmap->pixels = malloc(height * sizeof(pixel_t *));
	if (p_pixmap->pixels == NULL) {
		fprintf(stderr, "ERROR: Not enough memory to allocate pixmap.\n");
		return PIXMAP_ERROR;
	}
	for (int i = 0; i < height; ++i) {
		p_pixmap->pixels[i] = malloc(width * sizeof(pixel_t));
		if (p_pixmap->pixels[i] == NULL) {
			// Free already allocated memory on error (preventing leaks)
			fprintf(stderr, "ERROR: Not enough memory to allocate pixmap.\n");
			for (int j = 0; j < i; ++j) free(p_pixmap->pixels[j]);
			free(p_pixmap->pixels);
			p_pixmap->pixels = NULL;
			return PIXMAP_ERROR;
		}
		memset(p_pixmap->pixels[i], 0, width * sizeof(pixel_t));
	}

	return PIXMAP_SUCCESS;
}

int clear_pixmap(pixmap_t *p_pixmap)
{
	if (p_pixmap->pixels == NULL) {
		fprintf(stderr, "ERROR: Deallocating unallocated pixmap.\n");
		return PIXMAP_ERROR;
	}

	// Deallocate the matrix
	for (int i = 0; i < p_pixmap->height; ++i) free(p_pixmap->pixels[i]);
	free(p_pixmap->pixels);
	p_pixmap->pixels = NULL;

	return PIXMAP_SUCCESS;
}

int write_pixmap(pixmap_t *p_pixmap, FILE *p_file)
{
	int e;
	if (p_pixmap->pixels == NULL) {
		fprintf(stderr, "ERROR: Writing unallocated pixmap.\n");
		return PIXMAP_ERROR;
	}

	/* Write the header */
	e = fprintf(p_file, "P6\n%d %d\n255\n", p_pixmap->width, p_pixmap->height);
	if (e < 0) {
		fprintf(stderr, "ERROR: While writing the header of the pixmap.\n");
		return PIXMAP_ERROR;
	}

	// Flush because fprintf and fwrite might print differently
	fflush(p_file);
	for (int i = 0; i < p_pixmap->height; ++i) {
		e = fwrite(p_pixmap->pixels[i], sizeof(pixel_t), p_pixmap->width, p_file);
		if (e != p_pixmap->width) {
			fprintf(stderr, "ERROR: While writing line %d.\n", i);
			return PIXMAP_ERROR;
		}
	}
	fflush(p_file);

	return PIXMAP_SUCCESS;
}

void color_point(pixmap_t *p_pixmap, double x, double y, pixel_t pixel, blend_f *f)
{
	// Calculate the distance to the nearest pixels
	double a = x - (int)x;
	double b = y - (int)y;
	// Approximate the given point with the closest pixel
	// Other option would be linear interpolation
	int discret_x = a <= 0.5 ? (int)x : (int)x + 1;
	int discret_y = b <= 0.5 ? (int)y : (int)y + 1;
	p_pixmap->pixels[discret_x][discret_y] =
		f(p_pixmap->pixels[discret_x][discret_y], pixel);
}
