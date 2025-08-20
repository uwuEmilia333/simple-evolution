#include <stdlib.h>

float rand_float(float min, float max) {
    return min + (max - min) * (float)rand() / (float)RAND_MAX;
}

bool chancef(float chance) {
    return rand_float(0.0f, 1.0f) < chance;
}

bool chance(int numer, int denom) {
    return chancef((float)numer / (float)denom);
}
