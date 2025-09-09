//intend to make dist3test.cpp a stable version, revert some changes to the float mechanics values distribution
//make it a website which can plot the heatmap and can also simulate and live plot the trajectory
//of a single serve (using the params in valid_serve.csv)

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <random>
#include <numeric> 
#include <algorithm> 

class PerlinNoise {
    std::vector<int> p;
public:
    PerlinNoise(unsigned int seed) {
        p.resize(256);
        std::iota(p.begin(), p.end(), 0);
        std::default_random_engine engine(seed);
        std::shuffle(p.begin(), p.end(), engine);
        p.insert(p.end(), p.begin(), p.end());
    }

    double noise(double x, double y, double z) const {
        int X = (int)floor(x) & 255;
        int Y = (int)floor(y) & 255;
        int Z = (int)floor(z) & 255;
        x -= floor(x);
        y -= floor(y);
        z -= floor(z);
        double u = fade(x);
        double v = fade(y);
        double w = fade(z);
        int A = p[X] + Y, AA = p[A] + Z, AB = p[A + 1] + Z;
        int B = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;

        return lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),
                                     grad(p[BA], x - 1, y, z)),
                             lerp(u, grad(p[AB], x, y - 1, z),
                                     grad(p[BB], x - 1, y - 1, z))),
                     lerp(v, lerp(u, grad(p[AA + 1], x, y, z - 1),
                                     grad(p[BA + 1], x - 1, y, z - 1)),
                             lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                                     grad(p[BB + 1], x - 1, y - 1, z - 1))));
    }

private:
    static double fade(double t) { return t * t * t * (t * (t * 6 - 15) + 10); }
    static double lerp(double t, double a, double b) { return a + t * (b - a); }
    static double grad(int hash, double x, double y, double z) {
        int h = hash & 15;
        double u = h < 8 ? x : y;
        double v = h < 4 ? y : h == 12 || h == 14 ? x : z;
        return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }
};


const double PI = 3.1415926535897932;

struct Vector3D {  
    double x = 0.0, y = 0.0, z = 0.0;
    Vector3D operator+(const Vector3D& other) const { return {x + other.x, y + other.y, z + other.z}; }
    Vector3D operator-(const Vector3D& other) const { return {x - other.x, y - other.y, z - other.z}; }
    Vector3D operator*(double scalar) const { return {x * scalar, y * scalar, z * scalar}; }
    Vector3D operator/(double scalar) const { return {x / scalar, y / scalar, z / scalar}; }
    double magnitude() const { return std::sqrt(x*x + y*y + z*z); }
};
Vector3D cross(const Vector3D& a, const Vector3D& b) {  
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

namespace constants {
    const double MASS = 0.27; //OG: 0.27
    const double DIAMETER = 0.207; //OG: 0.21
    const double RADIUS = DIAMETER / 2.0;
    const double AREA = PI * RADIUS * RADIUS;
    const double GRAVITY = 9.80665; //OG: 9.80665
    const double AIR_DENSITY_RHO = 1.1644; //OG: 1.1644
    const double AIR_VISCOSITY_MU = 1.872e-5; //OG: 1.872e-5
    const double DT = 0.001; //timestep, OG: 0.001, optimal : 0.001 to 0.0001, precise calc: 0.00001
    const double NET_X_POSITION = 9.0; //OG: 9.0
    const double NET_HEIGHT = 2.43; //OG:2.43
    const double COURT_LENGTH = 18.0; //OG:18.0
    const double COURT_WIDTH = 9.0;
    const double COR_NET_NORMAL = 0.15; // OG: 0.15, For velocity component perpendicular to the net (x)
    const double COR_NET_TANGENT = 0.5; // OG: 0.5, For velocity components parallel to the net (y, z)
    const double SPIN_TO_LIFT_FACTOR = 0.001;     // OG: 0.005, How effectively topspin is converted to upward velocity.
    const double SPIN_DAMPENING_FACTOR = 0.4;     // How much the ball's spin is reduced by the collision.

    const double FLOAT_SERVE_RPS_THRESHOLD = 0.75; // rps
    // scales the overall strength of the random force
    const double FLOAT_FORCE_MAGNITUDE = 0.2; // newtons!
    
    const double FLOAT_FORCE_FREQUENCY = 15.0; //  higher values = faster oscillation.
}

struct InitialConditions { Vector3D start_pos; double speed_mps, angle_deg, topspin_rps, slice_rps; };
struct SimulationResult { bool isValid; Vector3D finalPosition; };

double calculate_drag_coefficient(double speed) {  
    if (speed < 1e-6) return 0.0;
    const double re = (constants::AIR_DENSITY_RHO * speed * constants::DIAMETER) / constants::AIR_VISCOSITY_MU;
    const double x = log10(re);
    double cd = 0.0;
    if (x >= 3.000 && x < 4.699) {
        double m = x - 3.0; cd = 0.4800 - 0.0143 * m + 0.0084 * m * m - 0.0028 * m * m * m;
    } else if (x >= 4.699 && x < 5.176) {
        double m = x - 4.699; cd = 0.4500 - 0.0007 * m - 0.00001 * m * m - 0.0076 * m * m * m;
    } else if (x >= 5.176 && x < 5.398) {
        double m = x - 5.176; cd = 0.2800 - 0.4682 * m - 0.0152 * m * m + 0.7062 * m * m * m;
    } else if (x >= 5.398 && x < 6.000) {
        double m = x - 5.398; cd = 0.1800 + 0.0102 * m + 0.4543 * m * m - 0.3205 * m * m * m;
    } else if (x >= 6.000 && x <= 6.699) {
        double m = x - 6.0; cd = 0.2000 + 0.0001 * m * m * m;
    } else { cd = 0.20; }
    return cd;
}

Vector3D calculate_acceleration(const Vector3D& velocity, const Vector3D& angular_velocity, double time, const PerlinNoise& pnoise) {
    Vector3D force_gravity = {0.0, -constants::MASS * constants::GRAVITY, 0.0};
    double speed = velocity.magnitude();
    double cd = calculate_drag_coefficient(speed);
    Vector3D force_drag = velocity * (-0.5 * constants::AIR_DENSITY_RHO * constants::AREA * cd * speed);
    Vector3D force_magnus = cross(angular_velocity, velocity) * (0.5 * constants::AIR_DENSITY_RHO * constants::AREA * constants::RADIUS);
    
    Vector3D net_force = force_gravity + force_drag + force_magnus;

    double total_spin_rps = angular_velocity.magnitude() / (2.0 * PI);
    if (total_spin_rps <= constants::FLOAT_SERVE_RPS_THRESHOLD) {
        double noise_val1 = pnoise.noise(time * constants::FLOAT_FORCE_FREQUENCY, 0, 0);
        double noise_val2 = pnoise.noise(0, time * constants::FLOAT_FORCE_FREQUENCY, 0);
        
        Vector3D up = {0.0, 1.0, 0.0};
        Vector3D side = cross(velocity, up);
        if (side.magnitude() > 1e-6) {
            side = side / side.magnitude();
        }

        // combine random noise components ->create a force vector
        // ts perpendicular to the velocity.
        Vector3D random_direction = (up * noise_val1 + side * noise_val2);
        if (random_direction.magnitude() > 1e-6) {
            random_direction = random_direction / random_direction.magnitude();
        }

        // scale the force by the magnitude constant and the square of the speed
        // (aerodynamic forces scale with v^2)
        Vector3D force_float = random_direction * constants::FLOAT_FORCE_MAGNITUDE * (speed * speed / (25*25)); // Normalize to a 25m/s serve
        
        net_force = net_force + force_float;
    }

    return net_force / constants::MASS;
}

void handle_net_collision(Vector3D& velocity, Vector3D& angular_velocity) {
    double topspin_rad_s = -angular_velocity.z; 
    double vertical_kick = 0.0;
    if (topspin_rad_s > 0) { vertical_kick = topspin_rad_s * velocity.x * constants::SPIN_TO_LIFT_FACTOR; }
    velocity.x *= -constants::COR_NET_NORMAL;
    velocity.y = std::abs(velocity.y) * constants::COR_NET_TANGENT + vertical_kick;
    velocity.z *= constants::COR_NET_TANGENT;
    angular_velocity = angular_velocity * constants::SPIN_DAMPENING_FACTOR;
}

// create perlin noise
SimulationResult run_simulation(const InitialConditions& initial, const PerlinNoise& pnoise) {
    Vector3D position = initial.start_pos;
    double launch_angle_rad = initial.angle_deg * PI / 180.0;
    Vector3D velocity = {
        initial.speed_mps * std::cos(launch_angle_rad),
        initial.speed_mps * std::sin(launch_angle_rad), 0.0
    };
    Vector3D angular_velocity = {0.0, initial.slice_rps * 2.0 * PI, -initial.topspin_rps * 2.0 * PI };

    bool passed_net_successfully = false;
    bool stayed_in_bounds_laterally = true;
    bool landed_in_court = false;
    Vector3D last_position = position;
    bool net_interaction_handled = false;
    double time = 0.0;

    while (position.y > 0) {
        last_position = position;
        Vector3D acceleration = calculate_acceleration(velocity, angular_velocity, time, pnoise);
        velocity = velocity + acceleration * constants::DT;
        position = position + velocity * constants::DT;
        time += constants::DT;

        if (std::abs(position.z) > constants::COURT_WIDTH / 2.0) {
            stayed_in_bounds_laterally = false;
            break;
        }

        if (!net_interaction_handled && position.x >= constants::NET_X_POSITION) {
            net_interaction_handled = true;
            if ((position.y - constants::RADIUS) >= constants::NET_HEIGHT) {
                passed_net_successfully = true;
            } else if ((position.y + constants::RADIUS) > 0) {
                handle_net_collision(velocity, angular_velocity);
                if (velocity.y > 0) { passed_net_successfully = true; } 
                else { break; }
            } else { break; }
        }
        
        if (time > 5.0) break;
    }

    if (position.x > constants::NET_X_POSITION && position.x <= constants::COURT_LENGTH &&
        std::abs(position.z) <= constants::COURT_WIDTH / 2.0) {
        landed_in_court = true;
    }

    SimulationResult result;
    result.isValid = passed_net_successfully && landed_in_court && stayed_in_bounds_laterally;
    result.finalPosition = position;
    return result;
}

int main(){
    const int NUM_SIMULATIONS = 50000;
    static long long int x = pow(10, floor(log10(NUM_SIMULATIONS)/1.5));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> posZ_dist(-constants::COURT_WIDTH / 2.0, constants::COURT_WIDTH / 2.0);
    std::uniform_real_distribution<> speed_dist(16.0, 35.0);
    std::uniform_real_distribution<> angle_dist(-1.0, 20.0);
    
    std::uniform_real_distribution<> topspin_dist(0.0, 40.0); // adjust topspin -> more floats
    std::uniform_real_distribution<> slice_dist(-10.0, 10.0);

    std::ofstream outputFile;
    std::string filename = "valid_serves.csv"; // writing to file
    outputFile.open(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << "\n";
        return 1;
    }
    outputFile << "start_z,speed_mps,angle_deg,topspin_rps,slice_rps,final_x,final_y,final_z\n";
    int valid_serve_count = 0;
    std::cout << "Running " << NUM_SIMULATIONS << " simulations (with float serve physics)..." << "\n";

    for (int i = 0; i < NUM_SIMULATIONS; ++i) {
        PerlinNoise pnoise(rd());

        InitialConditions current_serve;
        current_serve.start_pos = {0.0, 2.8, posZ_dist(gen)};
        current_serve.speed_mps = speed_dist(gen);
        current_serve.angle_deg = angle_dist(gen);
        current_serve.topspin_rps = topspin_dist(gen);
        current_serve.slice_rps = slice_dist(gen);

        SimulationResult result = run_simulation(current_serve, pnoise);
        
        if (result.isValid) {
            valid_serve_count++;
            outputFile << std::fixed << std::setprecision(4)
                       << current_serve.start_pos.z << ","
                       << current_serve.speed_mps << ","
                       << current_serve.angle_deg << ","
                       << current_serve.topspin_rps << ","
                       << current_serve.slice_rps << ","
                       << result.finalPosition.x << ","
                       << result.finalPosition.y << ","
                       << result.finalPosition.z << "\n";
        }
        if ((i + 1) % x == 0) {  // prints progress every x, uses -O3 op
            std::cout << "Completed " << (i + 1) << "/" << NUM_SIMULATIONS << " simulations. Found " << valid_serve_count << " valid serves so far." << std::endl;
        }
    }
    outputFile.close();
    std::cout << "\n--- Simulation Finished ---" << "\n";
    std::cout << "Total simulations run: " << NUM_SIMULATIONS << "\n";
    std::cout << "Total valid serves found: " << valid_serve_count << "\n";
    std::cout << "Results saved to " << filename << "\n";
    return 0;
}