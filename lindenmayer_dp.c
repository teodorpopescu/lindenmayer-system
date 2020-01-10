#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lindenmayer_dp.h"

#define PI 3.14159265359

double min(double x, double y)
{
	return x < y ? x : y;
}

double max(double x, double y)
{
	return x < y ? y : x;
}

int compute_no_of_variables(lindenmayer_system *p_lsystem)
{
	int ans = 0;
	for (int i = 0; i < 256; ++i) {
		if (p_lsystem->rules[i] != NULL) ++ans;
	}
	return ans;
}

void attach_variable_names(lindenmayer_system *p_lsystem, lindenmayer_dp_entry *v)
{
	for (int j = 0, i = 0; i < 256; ++i) {
		if (p_lsystem->rules[i] != NULL) v[j++].variable = (char)i;
	}
}

lindenmayer_dp_entry scan_rule(lindenmayer_system *p_lsystem, char *rule,
	lindenmayer_dp_entry *previous_entries, int n_variables, int do_expand,
	lindenmayer_dp_entry *polls, int *starting, int n_polls)
{
	lindenmayer_dp_entry ans;
	ans.x = ans.y = ans.angle = 0;
	ans.min_x = ans.max_x = 0;
	ans.min_y = ans.max_y = 0;
	for (int i_poll = 0, j = 0; rule[j] != '\0'; ++j) {
		if (polls != NULL && i_poll < n_polls && starting[i_poll] == j) {
			ans.angle = fmod(ans.angle, 2 * PI);
			polls[i_poll++] = ans;
		}
		if (p_lsystem->rules[(int)rule[j]] != NULL && do_expand) {
			int k;
			for (k = 0; k < n_variables; ++k) {
				if (previous_entries[k].variable == rule[j]) break;
			}
			double cos_tmp = cos(-ans.angle);
			double sin_tmp = sin(-ans.angle);
			double x1 = cos_tmp * previous_entries[k].max_x + sin_tmp * previous_entries[k].min_y;
			double x2 = cos_tmp * previous_entries[k].min_x + sin_tmp * previous_entries[k].min_y;
			double x3 = cos_tmp * previous_entries[k].max_x + sin_tmp * previous_entries[k].max_y;
			double x4 = cos_tmp * previous_entries[k].min_x + sin_tmp * previous_entries[k].max_y;
			double y1 = -sin_tmp * previous_entries[k].min_x + cos_tmp * previous_entries[k].min_y;
			double y2 = -sin_tmp * previous_entries[k].min_x + cos_tmp * previous_entries[k].max_y;
			double y3 = -sin_tmp * previous_entries[k].max_x + cos_tmp * previous_entries[k].min_y;
			double y4 = -sin_tmp * previous_entries[k].max_x + cos_tmp * previous_entries[k].max_y;
			double min_x = ans.x + min(min(x1, x2), min(x3, x4));
			double max_x = ans.x + max(max(x1, x2), max(x3, x4));
			double min_y = ans.y + min(min(y1, y2), min(y3, y4));
			double max_y = ans.y + max(max(y1, y2), max(y3, y4));
			ans.x += cos_tmp * previous_entries[k].x + sin_tmp * previous_entries[k].y;
			ans.y += -sin_tmp * previous_entries[k].x + cos_tmp * previous_entries[k].y;
			ans.min_x = min(ans.min_x, min_x);
			ans.max_x = max(ans.max_x, max_x);
			ans.min_y = min(ans.min_y, min_y);
			ans.max_y = max(ans.max_y, max_y);
			ans.angle += previous_entries[k].angle;
		} else if (p_lsystem->is_forward[(int)rule[j]]) {
			ans.x += cos(ans.angle);
			ans.y += sin(ans.angle);
			ans.min_x = min(ans.min_x, ans.x);
			ans.max_x = max(ans.max_x, ans.x);
			ans.min_y = min(ans.min_y, ans.y);
			ans.max_y = max(ans.max_y, ans.y);
		} else if (rule[j] == '+') {
			ans.angle += p_lsystem->angle;
		} else if (rule[j] == '-') {
			ans.angle -= p_lsystem->angle;
		}
	}
	ans.angle = fmod(ans.angle, 2 * PI);
	return ans;
}

lindenmayer_dp_entry **create_lindenmayer_dp_table(
	lindenmayer_system *p_lsystem, int n)
{
	lindenmayer_dp_entry **ans = malloc((n + 1) * sizeof(lindenmayer_dp_entry *));
	int n_variables = compute_no_of_variables(p_lsystem);
	ans[n] = malloc(n_variables * sizeof(lindenmayer_dp_entry));
	attach_variable_names(p_lsystem, ans[n]);

	// Base case
	if (n == 0) {
		for (int i = 0; i < n_variables; ++i) {
			ans[0][i].angle = ans[0][i].x = ans[0][i].y = 0;
			ans[0][i].min_x = ans[0][i].max_x = 0;
			ans[0][i].min_y = ans[0][i].max_y = 0;
		}
		return ans;
	}

	// Compute the current iteration based on the previous n
	lindenmayer_dp_entry **tmp = create_lindenmayer_dp_table(p_lsystem, n - 1);
	for (int i = 0; i < n; ++i) ans[i] = tmp[i];
	free(tmp);
	for (int i = 0; i < n_variables; ++i) {
		char tmp = ans[n][i].variable;
		ans[n][i] = scan_rule(p_lsystem, p_lsystem->rules[(int)ans[n][i].variable],
		                      ans[n - 1], n_variables, n > 1 ? 1 : 0, NULL, NULL, 0);
		ans[n][i].variable = tmp;
	}
	return ans;
}
