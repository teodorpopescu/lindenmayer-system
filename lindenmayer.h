#ifndef LINDENMAYER_H
#define LINDENMAYER_H

#include <stdint.h>

typedef struct {
	char *rules[256];
	char *start;
	uint8_t is_forward[256];
	double angle;
} lindenmayer_system;

void initialize_dragon_curve(lindenmayer_system *p_lsystem);

void initialize_koch_curve(lindenmayer_system *p_lsystem);

void initialize_sierpinsky_triangle(lindenmayer_system *p_lsystem);

void initialize_quadratic_gosper(lindenmayer_system *p_lsystem);

void initialize_levy_curve(lindenmayer_system *p_lsystem);

void initialize_pentaplexity(lindenmayer_system *p_lsystem);

/**
 *    Deallocate the memory used by the given lindenmayer system. The result
 * might be undefined so it should no longer be used without initializing it
 * first
 */
void clear_lsystem(lindenmayer_system *p_lsystem);

/**
 *    Expand the given path using the given lindenmayer system. The returned
 * string will be allocated and should be deallocated by the user of this
 * function. The initial path will not be changed.
 */
char *expand_path(lindenmayer_system *p_lsystem, char *path);

/**
 *    Expand the given lindemayer system for n times.
 */
char *expand_lsystem(lindenmayer_system *p_lsystem, int n);
#endif
