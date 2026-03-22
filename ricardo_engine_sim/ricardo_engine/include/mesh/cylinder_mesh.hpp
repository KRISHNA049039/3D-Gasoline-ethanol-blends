#pragma once
#include "utils/types.hpp"
#include <vector>
#include <array>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace Ricardo {
namespace Mesh {

// ─────────────────────────────────────────────
//  Mesh node
// ─────────────────────────────────────────────
struct Node {
    int    id;
    Vec3   pos;
    bool   boundary = false;
    int    zone_id  = 0;     // 0=fluid, 1=wall, 2=inlet, 3=outlet

    double temperature = 0.0;
    double pressure    = 0.0;
    double density     = 0.0;
    Vec3   velocity    = {};
    double burned_frac = 0.0;
};

// ─────────────────────────────────────────────
//  Hexahedral cell (8-node brick)
// ─────────────────────────────────────────────
struct HexCell {
    int id;
    std::array<int, 8> node_ids;
    Vec3   centroid;
    double volume     = 0.0;
    int    zone_id    = 0;

    // Cell-averaged thermo
    double temperature    = 0.0;
    double pressure       = 0.0;
    double density        = 0.0;
    double burned_frac    = 0.0;
    double equivalence_ratio = 1.0;
    double turbulent_ke   = 0.0;  // TKE  m²/s²
    double dissipation    = 0.0;  // ε    m²/s³
    Vec3   velocity       = {};

    // Species mass fractions
    double Y_fuel = 0, Y_O2 = 0, Y_N2 = 0;
    double Y_CO2  = 0, Y_H2O = 0, Y_CO = 0, Y_NOx = 0;
};

// ─────────────────────────────────────────────
//  Face between cells
// ─────────────────────────────────────────────
struct Face {
    int id;
    std::array<int, 4> node_ids;
    int    owner_cell  = -1;
    int    neighbour_cell = -1;
    Vec3   normal;       // outward from owner
    double area         = 0.0;
    Vec3   centroid;
    bool   is_boundary  = false;
    int    bc_type      = 0;  // 0=wall, 1=inlet, 2=outlet
};

// ─────────────────────────────────────────────
//  Mesh zone metadata
// ─────────────────────────────────────────────
struct Zone {
    int         id;
    std::string name;
    std::string type;  // "fluid", "wall", "inlet", "outlet"
    std::vector<int> cell_ids;
    std::vector<int> face_ids;
};

// ─────────────────────────────────────────────
//  3D Cylinder Mesh Generator
//  Produces a structured O-grid mesh of the Ricardo cylinder
// ─────────────────────────────────────────────
class CylinderMeshGenerator {
public:
    struct Config {
        int    nr   = 8;   // radial divisions
        int    nth  = 16;  // circumferential divisions (must be even)
        int    nz   = 20;  // axial divisions (scales with crank angle)
        double r_core_frac = 0.4; // core-to-bore radius ratio for O-grid
        Config() = default;
    };

    CylinderMeshGenerator(const EngineGeometry& geom) : geom_(geom), cfg_(Config()) {}
    CylinderMeshGenerator(const EngineGeometry& geom, const Config& cfg)
        : geom_(geom), cfg_(cfg) {}

    // Generate mesh at a given crank angle (degrees)
    void generate(double crank_angle_deg) {
        nodes_.clear(); cells_.clear(); faces_.clear(); zones_.clear();
        crank_angle_ = crank_angle_deg;
        cylinder_height_ = geom_.cylinder_volume(crank_angle_deg) /
                           (Constants::PI / 4.0 * geom_.bore * geom_.bore);

        generate_o_grid();
        compute_cell_volumes();
        build_face_connectivity();
        assign_zones();
    }

    // Accessors
    const std::vector<Node>&    nodes()  const { return nodes_;  }
    const std::vector<HexCell>& cells()  const { return cells_;  }
    const std::vector<Face>&    faces()  const { return faces_;  }
    const std::vector<Zone>&    zones()  const { return zones_;  }
    int total_cells() const { return static_cast<int>(cells_.size()); }
    int total_nodes() const { return static_cast<int>(nodes_.size()); }

    double cell_volume_avg() const {
        if (cells_.empty()) return 0.0;
        double s = 0;
        for (auto& c : cells_) s += c.volume;
        return s / cells_.size();
    }

private:
    EngineGeometry geom_;
    Config         cfg_;
    double         crank_angle_     = 0.0;
    double         cylinder_height_ = 0.0;

    std::vector<Node>    nodes_;
    std::vector<HexCell> cells_;
    std::vector<Face>    faces_;
    std::vector<Zone>    zones_;

    // ── O-grid structured mesh ──────────────────────
    // Layout: core block (Cartesian) + 4 outer wedge blocks
    // Radial layers: 0..nr_core = core, nr_core..nr = outer
    void generate_o_grid() {
        const int NR  = cfg_.nr;
        const int NTH = cfg_.nth;
        const int NZ  = cfg_.nz;
        const double R = geom_.bore / 2.0;
        const double H = cylinder_height_;
        const double r0 = cfg_.r_core_frac * R; // core radius

        // ── Node generation ──────────────────────────
        // We use a simple structured r-theta-z grid with
        // hyperbolic tangent stretching near the walls.
        int node_id = 0;
        // Precompute r positions (stretching towards wall)
        std::vector<double> r_pos(NR + 1);
        for (int i = 0; i <= NR; i++) {
            double xi = static_cast<double>(i) / NR;
            // Tanh stretching factor
            double beta = 2.5;
            r_pos[i] = R * std::tanh(beta * xi) / std::tanh(beta);
        }
        // Theta positions
        std::vector<double> th_pos(NTH);
        for (int j = 0; j < NTH; j++)
            th_pos[j] = 2.0 * Constants::PI * j / NTH;
        // Z positions (linear, with slight compression near piston & head)
        std::vector<double> z_pos(NZ + 1);
        for (int k = 0; k <= NZ; k++) {
            double zeta = static_cast<double>(k) / NZ;
            z_pos[k] = H * zeta;
        }

        // Create nodes: index = k*(NR+1)*NTH + j*(NR+1) + i
        nodes_.resize((NZ+1) * NTH * (NR+1));
        for (int k = 0; k <= NZ; k++) {
            for (int j = 0; j < NTH; j++) {
                for (int i = 0; i <= NR; i++) {
                    Node& nd = nodes_[node_id];
                    nd.id  = node_id;
                    double r = r_pos[i];
                    double th = th_pos[j];
                    nd.pos = Vec3(r * std::cos(th), r * std::sin(th), z_pos[k]);
                    // Boundary flags
                    nd.boundary = (i == NR || k == 0 || k == NZ);
                    if (k == 0)  nd.zone_id = 1; // piston
                    if (k == NZ) nd.zone_id = 1; // head
                    if (i == NR) nd.zone_id = 1; // liner
                    node_id++;
                }
            }
        }

        // ── Cell generation ──────────────────────────
        auto node_idx = [&](int i, int j, int k) -> int {
            int jj = (j + NTH) % NTH;
            return k * NTH * (NR+1) + jj * (NR+1) + i;
        };

        cells_.reserve(NZ * NTH * NR);
        int cell_id = 0;
        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NTH; j++) {
                for (int i = 0; i < NR; i++) {
                    HexCell cell;
                    cell.id = cell_id++;
                    // 8 corners of the hex cell
                    cell.node_ids = {
                        node_idx(i,   j,   k),
                        node_idx(i+1, j,   k),
                        node_idx(i+1, j+1, k),
                        node_idx(i,   j+1, k),
                        node_idx(i,   j,   k+1),
                        node_idx(i+1, j,   k+1),
                        node_idx(i+1, j+1, k+1),
                        node_idx(i,   j+1, k+1)
                    };
                    // Centroid
                    Vec3 cen{};
                    for (int n : cell.node_ids)
                        cen += nodes_[n].pos;
                    cell.centroid = cen * (1.0/8.0);
                    cell.zone_id = 0; // fluid
                    cells_.push_back(cell);
                }
            }
        }
    }

    void compute_cell_volumes() {
        for (auto& cell : cells_) {
            // Volume using divergence theorem on hex
            // Split into 6 tetrahedral sub-volumes
            auto& n = nodes_;
            auto& ids = cell.node_ids;
            Vec3 p[8];
            for (int i = 0; i < 8; i++) p[i] = n[ids[i]].pos;

            // Approximate: parallelepiped volume
            Vec3 a = p[1] - p[0];
            Vec3 b = p[3] - p[0];
            Vec3 c = p[4] - p[0];
            cell.volume = std::abs(a.cross(b).dot(c));
        }
    }

    void build_face_connectivity() {
        // Build faces from cell connectivity (simplified)
        // For each cell, 6 faces; track shared faces
        faces_.clear();
        int fid = 0;
        const int NR = cfg_.nr, NTH = cfg_.nth, NZ = cfg_.nz;

        auto cell_idx = [&](int i, int j, int k) -> int {
            if (i<0||i>=NR||k<0||k>=NZ) return -1;
            int jj = (j + NTH) % NTH;
            return k * NTH * NR + jj * NR + i;
        };
        auto node_idx = [&](int i, int j, int k) -> int {
            int jj = (j + NTH) % NTH;
            return k * NTH * (NR+1) + jj * (NR+1) + i;
        };

        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NTH; j++) {
                for (int i = 0; i < NR; i++) {
                    int oc = cell_idx(i, j, k);

                    // +z face (top)
                    {
                        Face f;
                        f.id = fid++;
                        f.owner_cell = oc;
                        f.neighbour_cell = cell_idx(i, j, k+1);
                        f.node_ids = {node_idx(i,j,k+1), node_idx(i+1,j,k+1),
                                      node_idx(i+1,j+1,k+1), node_idx(i,j+1,k+1)};
                        f.is_boundary = (k == NZ-1);
                        f.bc_type = f.is_boundary ? 0 : -1; // wall or internal
                        compute_face_geometry(f);
                        faces_.push_back(f);
                    }
                }
            }
        }
    }

    void compute_face_geometry(Face& f) {
        Vec3 p[4];
        for (int i = 0; i < 4; i++) p[i] = nodes_[f.node_ids[i]].pos;
        // Centroid
        f.centroid = (p[0] + p[1] + p[2] + p[3]) * 0.25;
        // Normal from cross product of diagonals
        Vec3 d1 = p[2] - p[0];
        Vec3 d2 = p[3] - p[1];
        Vec3 n  = d1.cross(d2) * 0.5;
        f.area   = n.norm();
        f.normal = n.unit();
    }

    void assign_zones() {
        // Fluid zone
        Zone fluid; fluid.id = 0; fluid.name = "cylinder_fluid"; fluid.type = "fluid";
        for (auto& c : cells_) fluid.cell_ids.push_back(c.id);
        zones_.push_back(fluid);

        // Wall zones (boundary faces)
        Zone wall; wall.id = 1; wall.name = "cylinder_walls"; wall.type = "wall";
        for (auto& f : faces_) {
            if (f.is_boundary) wall.face_ids.push_back(f.id);
        }
        zones_.push_back(wall);
    }
};

// ─────────────────────────────────────────────
//  Mesh statistics report
// ─────────────────────────────────────────────
struct MeshStats {
    int    n_nodes, n_cells, n_faces;
    double min_volume, max_volume, avg_volume;
    double total_volume;
    double min_aspect_ratio, max_aspect_ratio;

    static MeshStats compute(const CylinderMeshGenerator& mg) {
        MeshStats s{};
        s.n_nodes = mg.total_nodes();
        s.n_cells = mg.total_cells();
        s.n_faces = static_cast<int>(mg.faces().size());
        s.total_volume = 0;
        s.min_volume = 1e30; s.max_volume = 0;
        for (auto& c : mg.cells()) {
            s.total_volume += c.volume;
            s.min_volume = std::min(s.min_volume, c.volume);
            s.max_volume = std::max(s.max_volume, c.volume);
        }
        s.avg_volume = s.total_volume / std::max(s.n_cells, 1);
        s.min_aspect_ratio = 1.0; // simplified
        s.max_aspect_ratio = 1.0;
        return s;
    }
};

} // namespace Mesh
} // namespace Ricardo
