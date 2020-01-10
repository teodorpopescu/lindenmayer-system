#include "lindenmayer.h"

#include <stdlib.h>
#include <string.h>

#define PI 3.14159265359

void initialize_dragon_curve(lindenmayer_system *p_lsystem)
{
	for (int i = 0; i < 256; ++i) {
		p_lsystem->rules[i] = NULL;
		p_lsystem->is_forward[i] = 0;
	}
	p_lsystem->rules[(int)'X'] = malloc(6 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'X'], "X+YF+");
	p_lsystem->rules[(int)'Y'] = malloc(6 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'Y'], "-FX-Y");
	p_lsystem->start = malloc(3 * sizeof(char));
	strcpy(p_lsystem->start, "FX");
	p_lsystem->is_forward[(int)'F'] = 1;
	p_lsystem->angle = PI / 2;
}

void initialize_koch_curve(lindenmayer_system *p_lsystem)
{
	for (int i = 0; i < 256; ++i) {
		p_lsystem->rules[i] = NULL;
		p_lsystem->is_forward[i] = 0;
	}
	p_lsystem->rules[(int)'F'] = malloc(10 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'F'], "F+F-F-F+F");
	p_lsystem->start = malloc(2 * sizeof(char));
	strcpy(p_lsystem->start, "F");
	p_lsystem->is_forward[(int)'F'] = 1;
	p_lsystem->angle = PI / 2;
}

void initialize_sierpinsky_triangle(lindenmayer_system *p_lsystem)
{
	for (int i = 0; i < 256; ++i) {
		p_lsystem->rules[i] = NULL;
		p_lsystem->is_forward[i] = 0;
	}
	p_lsystem->rules[(int)'F'] = malloc(10 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'F'], "F-G+F+G-F");
	p_lsystem->rules[(int)'G'] = malloc(3 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'G'], "GG");
	p_lsystem->start = malloc(6 * sizeof(char));
	strcpy(p_lsystem->start, "F-G-G");
	p_lsystem->is_forward[(int)'F'] = 1;
	p_lsystem->is_forward[(int)'G'] = 1;
	p_lsystem->angle = PI * 2 / 3;
}

void initialize_quadratic_gosper(lindenmayer_system *p_lsystem)
{
	for (int i = 0; i < 256; ++i) {
		p_lsystem->rules[i] = NULL;
		p_lsystem->is_forward[i] = 0;
	}
	p_lsystem->rules[(int)'X'] = malloc(68 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'X'], "XFX-YF-YF+FX+FX-YF-YFFX+YF+FXFXYF-FX+YF+FXFX+YF-FXYF-YF-FX+FX+YFYF-");
	p_lsystem->rules[(int)'Y'] = malloc(68 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'Y'], "+FXFX-YF-YF+FX+FXYF+FX-YFYF-FX-YF+FXYFYF-FX-YFFX+FX+YF-YF-FX+FX+YFY");
	p_lsystem->start = malloc(4 * sizeof(char));
	strcpy(p_lsystem->start, "-YF");
	p_lsystem->is_forward[(int)'F'] = 1;
	p_lsystem->angle = PI / 2;
}

void initialize_levy_curve(lindenmayer_system *p_lsystem)
{
	for (int i = 0; i < 256; ++i) {
		p_lsystem->rules[i] = NULL;
		p_lsystem->is_forward[i] = 0;
	}
	p_lsystem->rules[(int)'F'] = malloc(7 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'F'], "-F++F-");
	p_lsystem->start = malloc(11 * sizeof(char));
	strcpy(p_lsystem->start, "F++F++F++F");
	p_lsystem->is_forward[(int)'F'] = 1;
	p_lsystem->angle = PI / 4;
}

void initialize_pentaplexity(lindenmayer_system *p_lsystem)
{
	for (int i = 0; i < 256; ++i) {
		p_lsystem->rules[i] = NULL;
		p_lsystem->is_forward[i] = 0;
	}
	p_lsystem->rules[(int)'F'] = malloc(19 * sizeof(char));
	strcpy(p_lsystem->rules[(int)'F'], "F++F++F+++++F-F++F");
	p_lsystem->start = malloc(14 * sizeof(char));
	strcpy(p_lsystem->start, "F++F++F++F++F");
	p_lsystem->is_forward[(int)'F'] = 1;
	p_lsystem->angle = PI / 5;
}


void clear_lsystem(lindenmayer_system *p_lsystem)
{
	for (int i = 0; i < 256; ++i) {
		if (p_lsystem->rules[i] != NULL) {
			free(p_lsystem->rules[i]);
			p_lsystem->rules[i] = NULL;
		}
	}
	free(p_lsystem->start);
}

char *expand_path(lindenmayer_system *p_lsystem, char *path)
{
	int new_size = 1; // 1 is for string terminator
	for (int i = 0; path[i] != '\0'; ++i) {
		if (p_lsystem->rules[(int)path[i]] == NULL) ++new_size;
		else new_size += strlen(p_lsystem->rules[(int)path[i]]);
	}
	char *new_path = malloc(new_size * sizeof(char));
	for (int j = 0, i = 0; path[i] != '\0'; ++i) {
		if (p_lsystem->rules[(int)path[i]] == NULL) {
			new_path[j++] = path[i];
		} else {
			for (int k = 0; p_lsystem->rules[(int)path[i]][k] != '\0'; ++k) {
				new_path[j++] = p_lsystem->rules[(int)path[i]][k];
			}
		}
	}
	new_path[new_size - 1] = '\0';
	return new_path;
}

char *expand_lsystem(lindenmayer_system *p_lsystem, int n)
{
	char *ans = expand_path(p_lsystem, p_lsystem->start);
	for (int i = 1; i < n; ++i) {
		char *tmp = ans;
		ans = expand_path(p_lsystem, tmp);
		free(tmp);
	}
	return ans;
}
