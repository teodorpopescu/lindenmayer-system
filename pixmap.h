#ifndef PIXMAP_H
#define PIXMAP_H

#include <stdint.h>
#include <stdio.h>

#define PIXMAP_ERROR -1
#define PIXMAP_SUCCESS 0

#pragma pack(1)
typedef struct {
	uint8_t r;
	uint8_t g;
	uint8_t b;
} pixel_t;
#pragma pack()

typedef struct {
	int width, height;
	pixel_t **pixels;
} pixmap_t;

typedef pixel_t blend_f(pixel_t, pixel_t);
typedef pixel_t coloring_f(double x, double period);

/**
 *    Color in the hsv spectrum.
 */
pixel_t hsv_coloring(double x, double period);

/**
 *    Color the curve in alternativ white and red colors (alternating 3 times).
 */
pixel_t christmas_coloring(double x, double period);

/**
 *    Combine the pixels using the overlay blending mode.
 */
pixel_t blend_overlay(pixel_t a, pixel_t b);

/**
 *    Combine the pixels using the lighten blending mode.
 */
pixel_t blend_lighten(pixel_t a, pixel_t b);


/**
 *    Initialize the given pixmap by allocating the necessary memory given the
 * width and height. All the pixels will be black.
 *    @return PIXMAP_SUCCESS if successful or PIXMAP_ERROR otherwise
 */
int initialize_pixmap(pixmap_t *p_pixmap, int width, int height);

/**
 *    Free the memory used by the given pixmap. After this operation the pixmap
 * should no longer be used.
 *    @return PIXMAP_SUCCESS if successful or PIXMAP_ERROR otherwise
 */
int clear_pixmap(pixmap_t *p_pixmap);

/**
 *    Write the given pixmap as a PBM (portable bitmap format).
 *    @return PIXMAP_SUCCESS if successful or PIXMAP_ERROR otherwise
 */
int write_pixmap(pixmap_t *p_pixmap, FILE *p_file);

/**
 *    Color the given point (the size of a pixel). This will result in coloring
 * the neighbouring pixels in different ammounts. The blending mode indicated
 * by the given blending function will be used to blend the pixels with the
 * background.
 */
void color_point(pixmap_t *p_pixmap, double x, double y, pixel_t pixel, blend_f *f);

#endif
