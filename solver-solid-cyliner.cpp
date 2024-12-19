#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <precice/precice.hpp>

using Vector = std::vector<double>;

class CircularDomain {
private:
    double radius;
    Vector nodes;
    Vector element_areas;

public:
    CircularDomain(double radius, int num_nodes) : radius(radius) {
        for (int i = 0; i < num_nodes; i++) {
            double theta = 2 * M_PI * i / num_nodes;
            nodes.push_back(radius * cos(theta));
            nodes.push_back(radius * sin(theta));
            element_areas.push_back(M_PI * radius * radius / num_nodes);
        }
    }

    const Vector& getNodes() const {
        return nodes;
    }

    const Vector& getElementAreas() const {
        return element_areas;
    }

    void printNodes() const {
        for (int i = 0; i < nodes.size(); i += 2) {
            std::cout << "Node " << i / 2 + 1 << ": (" << nodes[i] << ", " << nodes[i + 1] << ")" << std::endl;
        }
    }
};

struct DataContainer {
  void save_old_state(const Vector &vertices, const double &theta, const double &theta_dot, const double &time) {
    old_vertices = vertices;
    old_theta = theta;
    old_theta_dot = theta_dot;
    old_time = time;
  }

  void reload_old_state(Vector &vertices, double &theta, double &theta_dot, double &time) const {
    vertices = old_vertices;
    theta = old_theta;
    theta_dot = old_theta_dot;
    time = old_time;
  }

  Vector old_vertices;
  double old_theta;
  double old_theta_dot;
  double old_time;
};

class Solver {
public:
  Solver(const double moment_of_inertia) : moment_of_inertia(moment_of_inertia) {}

  void solve(const Vector &forces, const Vector &initial_vertices, Vector &vertices, double &theta, double &theta_dot, const double spring_constant, const double delta_t) const {
    double moment = 0;
    for (unsigned int i = 0; i < forces.size() / 2; ++i)
      moment += vertices[2 * i] * forces[2 * i + 1] - vertices[2 * i + 1] * forces[2 * i];

    const double theta_old = theta;

    theta = (1. / (1 - (spring_constant / moment_of_inertia) * std::pow(delta_t, 2))) *
            (std::pow(delta_t, 2) * moment / moment_of_inertia + delta_t * theta_dot + theta);

    theta_dot = (theta - theta_old) / delta_t;

    for (uint i = 0; i < vertices.size() / 2; ++i) {
      const double x_coord = initial_vertices[2 * i];
      vertices[2 * i] = x_coord * std::cos(theta) + initial_vertices[2 * i + 1] * std::sin(theta);
      vertices[2 * i + 1] = -x_coord * std::sin(theta) + initial_vertices[2 * i + 1] * std::cos(theta);
    }
    std::cout << "Theta: " << theta << " Theta dot: " << theta_dot << " Moment: " << moment << " 弹簧力: " << spring_constant * theta << std::endl;
  }

private:
  const double moment_of_inertia;
};

int main() {
  std::cout << "圆形刚体: 启动... \n";

  const std::string config_file_name("../precice-config.xml");
  const std::string solver_name("Solid");
  const std::string mesh_name("Solid-Mesh");
  const std::string data_write_name("Displacement");
  const std::string data_read_name("Force");

  constexpr double radius = 1.0;
  constexpr int num_nodes = 100;

  constexpr double density = 10000;
  constexpr double spring_constant = -25;

  CircularDomain circular_domain(radius, num_nodes);
  circular_domain.printNodes();

  constexpr double mass = M_PI * radius * radius * density;
  constexpr double inertia_moment = (1. / 2) * mass * radius * radius;

  precice::Participant precice(solver_name, config_file_name, 0, 1);

  const int dim = precice.getMeshDimensions(mesh_name);

  Vector forces(dim * num_nodes);
  Vector vertices = circular_domain.getNodes();
  Vector displacement(dim * num_nodes);
  std::vector<int> vertex_ids(num_nodes);
  double theta_dot = 0.0;
  double theta = 0.0;

  const Vector initial_vertices = vertices;

  precice.setMeshVertices(mesh_name, vertices, vertex_ids);

  for (uint i = 0; i < displacement.size(); ++i)
    displacement[i] = vertices[i] - initial_vertices[i];

  if (precice.requiresInitialData()) {
    precice.writeData(mesh_name, data_write_name, vertex_ids, displacement);
  }

  precice.initialize();

  DataContainer data_container;
  const Solver solver(inertia_moment);

  double time = 0;
  while (precice.isCouplingOngoing()) {

    if (precice.requiresWritingCheckpoint())
      data_container.save_old_state(vertices, theta, theta_dot, time);

    double dt = precice.getMaxTimeStepSize();
    time += dt;
    std::cout << "圆形刚体: t = " << time << "s \n";

    std::cout << "圆形刚体: 读取耦合数据 \n";
    precice.readData(mesh_name, data_read_name, vertex_ids, 0, forces);

    const double current_spring = spring_constant;
    solver.solve(forces, initial_vertices, vertices, theta, theta_dot, current_spring, dt);

    for (uint i = 0; i < displacement.size(); ++i)
      displacement[i] = vertices[i] - initial_vertices[i];

    std::cout << "圆形刚体: 写入耦合数据 \n";
    precice.writeData(mesh_name, data_write_name, vertex_ids, displacement);

    std::cout << "圆形刚体: 在时间上推进\n";
    precice.advance(dt);

    if (precice.requiresReadingCheckpoint())
      data_container.reload_old_state(vertices, theta, theta_dot, time);

    if (precice.isTimeWindowComplete()) {
    }
  }

  std::cout << "圆形刚体: 关闭...\n";

  return 0;
}
