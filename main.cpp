#include <stdio.h>

#define N 100

struct Organism {
    // TODO
};

enum class CellTag {
    Nothing,
    Food,
    Organism
};

struct Cell {
    CellTag tag;
    union {
        int food;
        Organism organism;
    } data;
};

static Cell grid[N][N];

int main() {
    puts("The current idea is to simulate an N*N grid where each cell either holds nothing, food, or an organism.");
    return 0;
}
