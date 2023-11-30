#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* object){
	obj = object;
	origin = pos;
	cellsize = dims/cells;
	size = cells+1;
	nodes_length = size.product();
	nodes = new GridNode[nodes_length];
	node_area = cellsize.product();
}

// Copy constructor
Grid::Grid(const Grid& orig){}

Grid::~Grid(){
	delete[] nodes;
}

// Maps mass and velocity to the grid
void Grid::initializeMass() {
    // Reset the grid
    memset(nodes, 0, sizeof(GridNode) * nodes_length);

    // Map particle data to grid
    for (int i = 0; i < obj->size; i++) {
        Particle& p = obj->particles[i];
        // Convert particle position to grid coordinates
        p.grid_position = (p.position - origin) / cellsize;
        float ox = p.grid_position[0], oy = p.grid_position[1];

        // Shape function gives a blending radius of two;
        // so we do computations within a 2x2 square for each particle
        for (int idx = 0, y = oy - 1, y_end = y + 3; y <= y_end; y++) {
            if (y < 0 || y >= size[1]) continue;

            // Y-dimension interpolation
            float y_pos = oy - y;
            float wy = Grid::bspline(y_pos);
            float dy = Grid::bsplineSlope(y_pos);

            for (int x = ox - 1, x_end = x + 3; x <= x_end; x++, idx++) {
                if (x < 0 || x >= size[0]) continue;

                // X-dimension interpolation
                float x_pos = ox - x;
                float wx = Grid::bspline(x_pos);
                float dx = Grid::bsplineSlope(x_pos);

                // Final weight is the product of weights in each dimension
                float weight = wx * wy;
                p.weights[idx] = weight;

                // Weight gradient is a vector of partial derivatives
                p.weight_gradient[idx].setData(dx * wy, wx * dy);
                p.weight_gradient[idx] /= cellsize;

                // Interpolate mass
                nodes[(int)(y * size[0] + x)].mass += weight * p.mass;
            }
        }
    }
}


void Grid::initializeVelocities() {
    // Interpolate velocity after mass to conserve momentum
    for (int i = 0; i < obj->size; i++) {
        Particle& p = obj->particles[i];

        int ox = static_cast<int>(p.grid_position[0]);
        int oy = static_cast<int>(p.grid_position[1]);

        for (int idx = 0, y = oy - 1, y_end = y + 3; y <= y_end; y++) {
            for (int x = ox - 1, x_end = x + 3; x <= x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    // Interpolate velocity
                    int n = y * size[0] + x;
                    nodes[n].velocity += p.velocity * w * p.mass;
                    nodes[n].active = true;
                }
            }
        }
    }

    // Normalize velocity by mass for active nodes
    for (int i = 0; i < nodes_length; i++) {
        GridNode& node = nodes[i];
        if (node.active && node.mass > 0) {
            node.velocity /= node.mass;
        }
    }

    // Handle collisions
    collisionGrid();
}



// Maps volume from the grid to particles
// This should only be called once, at the beginning of the simulation
void Grid::calculateVolumes() const{
	//Estimate each particles volume (for force calculations)
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];

		int ox = p.grid_position[0],
			oy = p.grid_position[1];

		//First compute particle density
		p.density = 0;
		for (int idx=0, y=oy-1, y_end=y+3; y<=y_end; y++){
			for (int x=ox-1, x_end=x+3; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Node density is trivial
					p.density += w * nodes[(int) (y*size[0]+x)].mass;
				}
			}
		}

		p.density /= node_area;

		//Volume for each particle can be found from density
		p.volume = p.mass / p.density;
	}
}


void Grid::explicitVelocities(const Vector2f& gravity) {
    // Compute the forces and store them in velocity_new
    for (int i = 0; i < obj->size; i++) {
        Particle& p = obj->particles[i];

        // Solve for grid internal forces
        Matrix2f energy = p.energyDerivative();

        int ox = static_cast<int>(p.grid_position[0]);
        int oy = static_cast<int>(p.grid_position[1]);

        for (int idx = 0, y = oy - 1, y_end = y + 3; y <= y_end; y++) {
            for (int x = ox - 1, x_end = x + 3; x <= x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    // Weight the force onto nodes
                    int n = y * size[0] + x;
                    nodes[n].velocity_new += energy * p.weight_gradient[idx];
                }
            }
        }
    }

    // Compute velocities using Euler integration
    for (int i = 0; i < nodes_length; i++) {
        GridNode& node = nodes[i];
        if (node.active && node.mass > 0) {
            node.velocity_new = node.velocity + TIMESTEP * (gravity - node.velocity_new / node.mass);
        }
    }

    // Handle collisions
    collisionGrid();
}



#if ENABLE_IMPLICIT
void Grid::implicitVelocities() {
    // Initialize linear solve
    for (int idx = 0; idx < nodes_length; idx++) {
        GridNode& n = nodes[idx];
        n.imp_active = n.active;
        if (n.imp_active) {
            n.r.setData(n.velocity_new);
            n.err.setData(1);
        }
    }

    recomputeImplicitForces();
    for (int idx = 0; idx < nodes_length; idx++) {
        GridNode& n = nodes[idx];
        if (n.imp_active) {
            n.r = n.velocity_new - n.Er;
            n.p = n.r;
            n.rEr = n.r.dot(n.Er);
        }
    }

    recomputeImplicitForces();
    for (int idx = 0; idx < nodes_length; idx++) {
        GridNode& n = nodes[idx];
        if (n.imp_active) {
            n.Ep = n.Er;
        }
    }

    // Linear solve
    for (int i = 0; i < MAX_IMPLICIT_ITERS; i++) {
        bool done = true;
        for (int idx = 0; idx < nodes_length; idx++) {
            GridNode& n = nodes[idx];
            if (n.imp_active) {
                float div = n.Ep.dot(n.Ep);
                float alpha = n.rEr / div;
                n.err = alpha * n.p;
                float err = n.err.length();
                if (err < MAX_IMPLICIT_ERR || err > MIN_IMPLICIT_ERR || isnan(err)) {
                    n.imp_active = false;
                    continue;
                } else {
                    done = false;
                }
                n.velocity_new += n.err;
                n.r -= alpha * n.Ep;
            }
        }
        if (done) break;

        recomputeImplicitForces();
        for (int idx = 0; idx < nodes_length; idx++) {
            GridNode& n = nodes[idx];
            if (n.imp_active) {
                float temp = n.r.dot(n.Er);
                float beta = temp / n.rEr;
                n.rEr = temp;
                n.p *= beta;
                n.p += n.r;
                n.Ep *= beta;
                n.Ep += n.Er;
            }
        }
    }
}

void Grid::recomputeImplicitForces() {
    // Calculate delta force for each particle
    for (int i = 0; i < obj->size; i++) {
        Particle& p = obj->particles[i];
        int ox = p.grid_position[0], oy = p.grid_position[1];
        for (int idx = 0, y = oy - 1, y_end = y + 3; y <= y_end; y++) {
            for (int x = ox - 1, x_end = x + 3; x <= x_end; x++, idx++) {
                GridNode& n = nodes[(int)(y * size[0] + x)];
                if (n.imp_active) {
                    n.force += p.deltaForce(n.r, p.weight_gradient[idx]);
                }
            }
        }
    }

    // Compute Er for each node
    for (int idx = 0; idx < nodes_length; idx++) {
        GridNode& n = nodes[idx];
        if (n.imp_active) {
            n.Er = n.r - IMPLICIT_RATIO * TIMESTEP / n.mass * n.force;
        }
    }
}

#endif


void Grid::updateVelocities() const {
    for (int i = 0; i < obj->size; i++) {
        Particle& p = obj->particles[i];
        Vector2f pic, flip = p.velocity;
        Matrix2f& grad = p.velocity_gradient;
        grad.setData(0.0);
        p.density = 0;

        int ox = p.grid_position[0], oy = p.grid_position[1];
        for (int idx = 0, y = oy - 1, y_end = y + 3; y <= y_end; y++) {
            for (int x = ox - 1, x_end = x + 3; x <= x_end; x++, idx++) {
                float w = p.weights[idx];
                if (w > BSPLINE_EPSILON) {
                    GridNode &node = nodes[(int)(y * size[0] + x)];
                    pic += node.velocity_new * w;
                    flip += (node.velocity_new - node.velocity) * w;
                    grad += node.velocity_new.outer_product(p.weight_gradient[idx]);
                    p.density += w * node.mass;
                }
            }
        }
        p.velocity = flip * FLIP_PERCENT + pic * (1 - FLIP_PERCENT);
        p.density /= node_area;
    }
    collisionParticles();
}

void Grid::collisionGrid() {
    Vector2f delta_scale = Vector2f(TIMESTEP) / cellsize;

    for (int y = 0, idx = 0; y < size[1]; y++) {
        for (int x = 0; x < size[0]; x++, idx++) {
            GridNode &node = nodes[idx];
            if (node.active) {
                Vector2f new_pos = node.velocity_new * delta_scale + Vector2f(x, y);
                if (new_pos[0] < BSPLINE_RADIUS || new_pos[0] > size[0] - BSPLINE_RADIUS - 1) {
                    node.velocity_new[0] = 0;
                    node.velocity_new[1] *= STICKY;
                }
                if (new_pos[1] < BSPLINE_RADIUS || new_pos[1] > size[1] - BSPLINE_RADIUS - 1) {
                    node.velocity_new[0] *= STICKY;
                    node.velocity_new[1] = 0;
                }
            }
        }
    }
}


void Grid::collisionParticles() const {
    for (int i = 0; i < obj->size; i++) {
        Particle& p = obj->particles[i];
        Vector2f new_pos = p.grid_position + TIMESTEP * p.velocity / cellsize;
        // Collision with left and right borders
        if (new_pos[0] < BSPLINE_RADIUS - 1 || new_pos[0] > size[0] - BSPLINE_RADIUS)
            p.velocity[0] = -STICKY * p.velocity[0];
        // Collision with bottom and top borders
        if (new_pos[1] < BSPLINE_RADIUS - 1 || new_pos[1] > size[1] - BSPLINE_RADIUS)
            p.velocity[1] = -STICKY * p.velocity[1];
    }
}

void Grid::draw() {
    if (SUPPORTS_POINT_SMOOTH) glDisable(GL_POINT_SMOOTH);

    // Drawing grid nodes
    glPointSize(1);
    glColor3f(0.2, 0.2, 0.2);

    glBegin(GL_POINTS);
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++)
            glVertex2fv((origin + cellsize * Vector2f(i, j)).data);
    }
    glEnd();
}

