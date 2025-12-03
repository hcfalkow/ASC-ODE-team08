#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  Vec<D> m_gravity=0.0;
public:

  // ADDED: distance constraints:
  struct DistanceConstraint { // between two masses
        size_t i, j; // indices of two masses
        double length;
    };
    std::vector<DistanceConstraint> m_constraints; // store distance constraints

    // Add a distance constraint between two masses
    void addDistanceConstraint(size_t i, size_t j, double length) {
        m_constraints.push_back({i,j,length}); // Add the constraint to the list
    }

    auto & constraints() { return m_constraints; }
  // END ADDED

  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;
  
  // ADDED: print distance constraints
  ost << "constraints:" << std::endl;
    for (auto c : mss.constraints())
        ost << "i = " << c.i << ", j = " << c.j << ", length = " << c.length << std::endl;
  // END ADDED

  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D*mss.masses().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;

    auto xmat = x.asMatrix(mss.masses().size(), D);
    auto fmat = f.asMatrix(mss.masses().size(), D);

    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    for (auto spring : mss.springs())
      {
        auto [c1,c2] = spring.connectors;
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        double force = spring.stiffness * (norm(p1-p2)-spring.length);
        Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);
        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force*dir12;
      }

    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) *= 1.0/mss.masses()[i].mass;


    // ADDED: distance constraints forces
    const size_t C = mss.constraints().size();
    if (C == 0) return; // nothing to do

    size_t N = mss.masses().size();
    size_t dim = D * N;

    // --- Build Mass Matrix (diagonal) ---
    Matrix<double> M(dim, dim);
    M = 0.0;
    for (size_t i = 0; i < N; ++i)
        for (int d = 0; d < D; ++d)
            M(i*D+d, i*D+d) = mss.masses()[i].mass;

    // --- Build Constraint Jacobian J (C x dimX) ---
    Matrix<double> J(C, dim);
    J = 0.0;
    for (size_t c = 0; c < C; ++c)
    {
        auto & con = mss.constraints()[c];
        Vec<D> diff = xmat.row(con.i) - xmat.row(con.j);
        for (int d = 0; d < D; ++d)
        {
            J(c, con.i*D + d) = 2.0 * diff(d);
            J(c, con.j*D + d) = -2.0 * diff(d);
        }
    }

    // --- Build RHS vector F (unscaled forces) ---
    Vector<> F(dim);
    for (size_t i = 0; i < N; ++i)
        for (int d = 0; d < D; ++d)
            F(i*D + d) = fmat.row(i)(d) * mss.masses()[i].mass;

    // --- Build Augmented System [M J^T; J 0] ---
    Matrix<double> A(dim + C, dim + C);
    A = 0.0;
    Vector<> rhs(dim + C);
    rhs = 0.0;

    // Top-left M
    for (size_t i = 0; i < dim; ++i)
        for (size_t j = 0; j < dim; ++j)
            A(i,j) = M(i,j);

    // Top-right J^T
    for (size_t i = 0; i < dim; ++i)
        for (size_t j = 0; j < C; ++j)
            A(i, dim + j) = J(j,i);

    // Bottom-left J
    for (size_t i = 0; i < C; ++i)
        for (size_t j = 0; j < dim; ++j)
            A(dim + i, j) = J(i,j);

    // RHS
    for (size_t i = 0; i < dim; ++i) rhs(i) = F(i);
    for (size_t i = 0; i < C; ++i) rhs(dim + i) = 0.0;

      // --- Solve the system (dense LU, simple version) ---
    Vector<> sol(dim + C);
    sol = rhs; // initialize solution

    // Forward elimination
    for (size_t k = 0; k < dim + C; ++k)
    {
        if (std::abs(A(k,k)) < 1e-12)
            throw std::runtime_error("Zero pivot in LU solve");

        for (size_t i = k+1; i < dim + C; ++i)
        {
            double factor = A(i,k)/A(k,k);
            A(i,k) = factor; // store L
            for (size_t j = k+1; j < dim + C; ++j)
                A(i,j) -= factor*A(k,j);
            sol(i) -= factor*sol(k);
        }
    }

    // Back substitution
    for (int i = int(dim + C) - 1; i >= 0; --i)
    {
        for (size_t j = i+1; j < dim + C; ++j)
            sol(i) -= A(i,j)*sol(j);
        sol(i) /= A(i,i);
    }
    

    // --- Write accelerations back to fmat ---
    for (size_t i = 0; i < N; ++i)
        for (int d = 0; d < D; ++d)
            fmat.row(i)(d) = sol(i*D + d);
  }
  // End ADDED
  
  // old numerical differentiation:
  // virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  // {
  //   // TODO: exact differentiation
  //   double eps = 1e-8;
  //   Vector<> xl(dimX()), xr(dimX()), fl(dimF()), fr(dimF());
  //   for (size_t i = 0; i < dimX(); i++)
  //     {
  //       xl = x;
  //       xl(i) -= eps;
  //       xr = x;
  //       xr(i) += eps;
  //       evaluate (xl, fl);
  //       evaluate (xr, fr);
  //       df.col(i) = 1/(2*eps) * (fr-fl);
  //     }
  // }
  virtual void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override // exact differentiation
  {
      df = 0.0; // Initialize the derivative matrix to zero

      auto xmat = x.asMatrix(mss.masses().size(), D); // Reshape x into a matrix of size (N x D)
      const size_t N = mss.masses().size(); // Number of masses

      for (auto & spring : mss.springs()) // Loop over each spring
      {
          auto [c1, c2] = spring.connectors; // Get the connectors of the spring

          // Only springs between two masses contribute
          if (c1.type != Connector::MASS || c2.type != Connector::MASS) // Check if both connectors are masses
              continue;

          size_t i = c1.nr; // Index of the first mass
          size_t j = c2.nr; // Index of the second mass

          Vec<D> p1 = xmat.row(i); // Position of the first mass
          Vec<D> p2 = xmat.row(j); // Position of the second mass

          Vec<D> r = p2 - p1; // Vector from mass i to mass j
          double l = norm(r); // Length of the spring
          if (l < 1e-12) continue; // Avoid division by zero

          Vec<D> d = (1.0 / l) * r; // Unit direction vector from mass i to mass j

          double k = spring.stiffness; // Spring stiffness
          double L = spring.length; // Rest length of the spring
          double a = k * (1.0 - L / l); // Coefficient for the identity part
          double b = k * (L / l); // Coefficient for the outer product part

          // Build the D×D stiffness matrix explicitly
          double K[D][D]; // Stiffness matrix
          for (int r1 = 0; r1 < D; ++r1)
              for (int c1 = 0; c1 < D; ++c1)
                  K[r1][c1] = a * (r1 == c1 ? 1.0 : 0.0) + b * d(r1) * d(c1); // K = a*I + b*(d⊗d)

          // divide by masses
          double inv_mi = 1.0 / mss.masses()[i].mass;
          double inv_mj = 1.0 / mss.masses()[j].mass;

          // Fill df manually using the stiffness matrix K
          for (int r1 = 0; r1 < D; ++r1)
          {
              for (int c1 = 0; c1 < D; ++c1)
              {
                  size_t ri = i*D + r1;
                  size_t rj = j*D + r1;
                  size_t ci = i*D + c1;
                  size_t cj = j*D + c1;

                  df(ri, ci) -= K[r1][c1] * inv_mi; // ∂ai/∂pi
                  df(ri, cj) += K[r1][c1] * inv_mi; // ∂ai/∂pj
                  df(rj, ci) += K[r1][c1] * inv_mj; // ∂aj/∂pi
                  df(rj, cj) -= K[r1][c1] * inv_mj; // ∂aj/∂pj
              }
          }
      }
  }


  
};

#endif
