#ifndef LINDENMAYER_DP_H
#define LINDENMAYER_DP_H

#include "lindenmayer.h"

typedef struct {
	char variable;
	double x, y, angle;
	double min_x, min_y;
	double max_x, max_y;
} lindenmayer_dp_entry;

int compute_no_of_variables(lindenmayer_system *p_lsystem);

/**
 *    Expand the rule based on the given lindenmayer system in the sense of
 * obtaining the boundaries and the final coordinate. This will be based on the
 * previous entries computed. do_expand signifies if the rules should expand
 * when possible or not.
 *    @polls should be NULL and n_polls should be 0 if you don't want to store
 * additional values. Otherwise, at indices that are in starting a poll will
 * be made.
 */
lindenmayer_dp_entry scan_rule(lindenmayer_system *p_lsystem, char *rule,
	lindenmayer_dp_entry *previous_entries, int n_variables, int do_expand,
	lindenmayer_dp_entry *polls, int *starting, int n_polls);

/**
 *    Return a matrix (n lines and no_of_variables columns) that contains
 * lindenmayer_dp_entries. Each entry represents, if we start to draw at (0, 0)
 * with starting angle = 0, where we will stop drawing and what angle will be
 * facing. Also, determine a tight window (but not the tightes possible) in
 * which we can draw the fractal defined by the lindenmayer system.
 */
lindenmayer_dp_entry **create_lindenmayer_dp_table(
	lindenmayer_system *p_lsystem, int n);

#endif
