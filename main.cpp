#include <stdio.h>
#include <time.h>
#include <math.h>

#include <raylib.h>

#include "config.cpp"
#include "rand.cpp"

struct Organism {
    float genes[AMOUNT_GENES];
    Color color;
    int fullness;
};

enum class ActionRequest {
    Wait,
    MoveLeft,
    MoveRight,
    MoveUp,
    MoveDown,
    Eat,
    Reproduce
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

static Cell grid[GRID_SIZE][GRID_SIZE];
static long int iterations = 0;

Organism generate_organism() {
    Organism organism;
    for (size_t i = 0; i < AMOUNT_GENES; ++i)
        organism.genes[i] = GENE_BASE + VARIANCE * rand_float(-1.0, 1.0f);
    organism.color = ColorFromHSV(rand_float(0.0f, 360.0f), 1.0f, rand_float(0.6f, 0.8f));
    organism.fullness = ORGANISM_STOMACH_CAP;
    return organism;
}

Organism reproduce(const Organism *parent) {
    Organism child;
    for (size_t i = 0; i < AMOUNT_GENES; ++i) {
        child.genes[i] = parent->genes[i];
        if (chancef(MUTATION_CHANCE)) child.genes[i] += VARIANCE * rand_float(-1.0f, 1.0f);
    }
    child.color = parent->color;
    child.fullness = ORGANISM_STOMACH_CAP;
    return child;
}

ActionRequest act(const Organism *organism, int x, int y) {
    const size_t G = AMOUNT_GENES;
    auto g = [&](size_t i)->float { return fabsf(organism->genes[i % G]); };

    // allocate raw scores for: Wait, Left, Right, Up, Down, Eat, Reproduce
    float raw[7];
    for (int i = 0; i < 7; ++i) raw[i] = 0.01f; // small floor
    // mix genes into raw scores (cyclic)
    for (size_t i = 0; i < G; ++i) raw[i % 7] += g(i);

    // if eat impossible (no adjacent food) zero out eat weight
    bool food_adj = (x>0 && grid[x-1][y].tag==CellTag::Food) || (x+1<GRID_SIZE && grid[x+1][y].tag==CellTag::Food) ||
                    (y>0 && grid[x][y-1].tag==CellTag::Food) || (y+1<GRID_SIZE && grid[x][y+1].tag==CellTag::Food);
    if (!food_adj) raw[5] = 0.0f;

    // if reproduce impossible, zero out
    if (organism->fullness <= REPRODUCTION_COST) raw[6] = 0.0f;

    // invalidate directional moves that go into occupied/out-of-bounds
    if (x == 0 || grid[x-1][y].tag != CellTag::Nothing) raw[1] = 0.0f;
    if (x + 1 >= GRID_SIZE || grid[x+1][y].tag != CellTag::Nothing) raw[2] = 0.0f;
    if (y == 0 || grid[x][y-1].tag != CellTag::Nothing) raw[3] = 0.0f;
    if (y + 1 >= GRID_SIZE || grid[x][y+1].tag != CellTag::Nothing) raw[4] = 0.0f;

    // normalize and sample
    float sum = 0.0f;
    for (int i = 0; i < 7; ++i) sum += raw[i];
    if (sum <= 0.0f) return ActionRequest::Wait;

    float pick = rand_float(0.0f, sum);
    float acc = 0.0f;
    for (int i = 0; i < 7; ++i) {
        acc += raw[i];
        if (pick <= acc) {
            switch (i) {
                case 0: return ActionRequest::Wait;
                case 1: return ActionRequest::MoveLeft;
                case 2: return ActionRequest::MoveRight;
                case 3: return ActionRequest::MoveUp;
                case 4: return ActionRequest::MoveDown;
                case 5: return ActionRequest::Eat;
                default: return ActionRequest::Reproduce;
            }
        }
    }
    return ActionRequest::Wait;
}

void init_grid() {
    for (size_t y = 0; y < GRID_SIZE; ++y) {
        for (size_t x = 0; x < GRID_SIZE; ++x) {
            int r = rand() % 100;
            if (r < 90) {
                grid[x][y].tag = CellTag::Nothing;
            }
            else if (r < 99) {
                grid[x][y].tag = CellTag::Food;
                grid[x][y].data.food = rand() % 10;
            }
            else {
                grid[x][y].tag = CellTag::Organism;
                grid[x][y].data.organism = generate_organism();
            }
        }
    }
}

int count_cells(int (*value)(const Cell *)) {
    int total = 0;
    for (size_t y = 0; y < GRID_SIZE; ++y) {
        for (size_t x = 0; x < GRID_SIZE; ++x) {
            total += value(&grid[x][y]);
        }
    }
    return total;
}

void draw_board() {
    for (size_t y = 0; y < GRID_SIZE; ++y) {
        for (size_t x = 0; x < GRID_SIZE; ++x) {
            switch (grid[x][y].tag) {
            case CellTag::Nothing:
                DrawRectangle(x * SCALE, y * SCALE, SCALE, SCALE, DARKGRAY);
                break;
            case CellTag::Food:
                DrawRectangle(x * SCALE, y * SCALE, SCALE, SCALE, GOLD);
                break;
            case CellTag::Organism:
                DrawRectangle(x * SCALE, y * SCALE, SCALE, SCALE, grid[x][y].data.organism.color);
                break;
            }
        }
    }
}

void update_board() {
    iterations += 1;
    for (size_t y = 0; y < GRID_SIZE; ++y) {
        for (size_t x = 0; x < GRID_SIZE; ++x) {
            if (grid[x][y].tag != CellTag::Organism) continue;

            // work on a local copy of the organism
            Organism org = grid[x][y].data.organism;
            ActionRequest action = act(&org, (int)x, (int)y);

            auto write_back_here = [&]() {
                // put the possibly-updated organism back into its original cell
                grid[x][y].tag = CellTag::Organism;
                grid[x][y].data.organism = org;
            };

            switch (action) {
            case ActionRequest::Wait:
                org.fullness -= WAIT_COST;
                write_back_here();
                break;

            case ActionRequest::MoveLeft:
                org.fullness -= MOVE_COST;
                if (x > 0 && grid[x-1][y].tag == CellTag::Nothing) {
                    // place organism into left cell, clear original
                    grid[x-1][y].tag = CellTag::Organism;
                    grid[x-1][y].data.organism = org;
                    grid[x][y].tag = CellTag::Nothing;
                }
                else {
                    // failed attempt: stay but still pay cost
                    write_back_here();
                }
                break;

            case ActionRequest::MoveRight:
                org.fullness -= MOVE_COST;
                if (x + 1 < GRID_SIZE && grid[x+1][y].tag == CellTag::Nothing) {
                    grid[x+1][y].tag = CellTag::Organism;
                    grid[x+1][y].data.organism = org;
                    grid[x][y].tag = CellTag::Nothing;
                }
                else {
                    write_back_here();
                }
                break;

            case ActionRequest::MoveUp:
                org.fullness -= MOVE_COST;
                if (y > 0 && grid[x][y-1].tag == CellTag::Nothing) {
                    grid[x][y-1].tag = CellTag::Organism;
                    grid[x][y-1].data.organism = org;
                    grid[x][y].tag = CellTag::Nothing;
                }
                else {
                    write_back_here();
                }
                break;

            case ActionRequest::MoveDown:
                org.fullness -= MOVE_COST;
                if (y + 1 < GRID_SIZE && grid[x][y+1].tag == CellTag::Nothing) {
                    grid[x][y+1].tag = CellTag::Organism;
                    grid[x][y+1].data.organism = org;
                    grid[x][y].tag = CellTag::Nothing;
                }
                else {
                    write_back_here();
                }
                break;

            case ActionRequest::Eat: {
                // Attempt to eat in an order. If none available, pay a penalty (treat as WAIT).
                bool ate = false;

                // try left
                if (x > 0 && grid[x-1][y].tag == CellTag::Food) {
                    grid[x-1][y].data.food -= 1;
                    if (grid[x-1][y].data.food <= 0) grid[x-1][y].tag = CellTag::Nothing;
                    org.fullness = ORGANISM_STOMACH_CAP;
                    ate = true;
                }
                // try up
                else if (y > 0 && grid[x][y-1].tag == CellTag::Food) {
                    grid[x][y-1].data.food -= 1;
                    if (grid[x][y-1].data.food <= 0) grid[x][y-1].tag = CellTag::Nothing;
                    org.fullness = ORGANISM_STOMACH_CAP;
                    ate = true;
                }
                // try right
                else if (x + 1 < GRID_SIZE && grid[x+1][y].tag == CellTag::Food) {
                    grid[x+1][y].data.food -= 1;
                    if (grid[x+1][y].data.food <= 0) grid[x+1][y].tag = CellTag::Nothing;
                    org.fullness = ORGANISM_STOMACH_CAP;
                    ate = true;
                }
                // try down
                else if (y + 1 < GRID_SIZE && grid[x][y+1].tag == CellTag::Food) {
                    grid[x][y+1].data.food -= 1;
                    if (grid[x][y+1].data.food <= 0) grid[x][y+1].tag = CellTag::Nothing;
                    org.fullness = ORGANISM_STOMACH_CAP;
                    ate = true;
                }
                if (!ate) {
                    // failed attempt: organism still pays the (attempt) cost (use WAIT_COST)
                    org.fullness -= WAIT_COST;
                }

                // after eating (or failing), organism stays where it was (no movement)
                write_back_here();
                break;
            }

            case ActionRequest::Reproduce: {
                // Attempt reproduction: selecting only valid empty neighbours (bounds checked).
                // Important: only charge full REPRODUCTION_COST if reproduction actually places a child.
                // If reproduction fails, charge WAIT_COST instead.

                // collect valid empty neighbours
                int candidates_x[4];
                int candidates_y[4];
                int n = 0;
                if (x > 0 && grid[x-1][y].tag == CellTag::Nothing)      { candidates_x[n] = x - 1; candidates_y[n] = y; ++n; }
                if (y > 0 && grid[x][y-1].tag == CellTag::Nothing)      { candidates_x[n] = x;     candidates_y[n] = y - 1; ++n; }
                if (x + 1 < GRID_SIZE && grid[x+1][y].tag == CellTag::Nothing) { candidates_x[n] = x + 1; candidates_y[n] = y; ++n; }
                if (y + 1 < GRID_SIZE && grid[x][y+1].tag == CellTag::Nothing) { candidates_x[n] = x;     candidates_y[n] = y + 1; ++n; }

                if (n > 0 && org.fullness > REPRODUCTION_COST) {
                    // place child in a random empty neighbour
                    int pick = rand() % n;
                    int wx = candidates_x[pick];
                    int wy = candidates_y[pick];

                    Organism child = reproduce(&org);
                    grid[wx][wy].tag = CellTag::Organism;
                    grid[wx][wy].data.organism = child;

                    // parent pays full reproduction cost only on success
                    org.fullness -= REPRODUCTION_COST;
                }
                else {
                    // reproduction failed (no empty neighbour or not enough fullness): pay WAIT_COST instead
                    org.fullness -= WAIT_COST;
                }

                // parent remains where it was (if still alive)
                if (org.fullness > 0) write_back_here();
                else grid[x][y].tag = CellTag::Nothing;
                break;
            }

            default:
                // catch-all: put the organism back
                write_back_here();
                break;
            } // end switch

            // If the organism starved, clear the cell.
            if (grid[x][y].tag == CellTag::Organism && grid[x][y].data.organism.fullness <= 0) {
                grid[x][y].tag = CellTag::Nothing;
            }
        }
    }

    // --- spawn a little food each frame ---
    // Try up to GRID_SIZE*2 random attempts to find a cell that's Nothing or Food.
    // If found: if Nothing -> make a Food cell with 1 unit; if Food -> increment food by 1.
    // If grid is completely full of organisms, we'll just skip adding food.
    int attempts = GRID_SIZE * 2;
    while (attempts-- > 0) {
        int fx = rand() % GRID_SIZE;
        int fy = rand() % GRID_SIZE;
        if (grid[fx][fy].tag == CellTag::Nothing) {
            grid[fx][fy].tag = CellTag::Food;
            grid[fx][fy].data.food = 1;
            break;
        }
        else if (grid[fx][fy].tag == CellTag::Food) {
            grid[fx][fy].data.food += 1;
            break;
        }
        else {
            // It's an organism; try again
        }
    }
}

int main() {
    srand(time(NULL));
    init_grid();
    InitWindow(GRID_SIZE * SCALE, GRID_SIZE * SCALE, "Simple Organisms");
    SetTargetFPS(160);
    while (!WindowShouldClose()) {
        int organisms = count_cells([](const Cell *c) { return c->tag == CellTag::Organism ? 1 : 0; });
        if (organisms == 0) break;
        int food_total = count_cells([](const Cell *c) { return c->tag == CellTag::Food ? c->data.food : 0; });
        int food_cells = count_cells([](const Cell *c) { return c->tag == CellTag::Food ? 1 : 0; });
        BeginDrawing();
            draw_board();
        EndDrawing();
        update_board();
    }
    printf("all life has gone extinct after %d iterations\n", iterations);
    CloseWindow();
    return 0;
}
