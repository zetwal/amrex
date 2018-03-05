#include "AMReX_SPACE.H"
#include "AMReX_RealVect.H"
#include "AMReX_Utility.H"
using std::ostream;
using std::istream;
using std::ws;

namespace amrex
{

  RealVect tm;
  const RealVect RealVect::Unit(AMREX_D_DECL6(1.0,1.0,1.0,1.0,1.0,1.0));

  const RealVect RealVect::Zero(AMREX_D_DECL6(0.0,0.0,0.0,0.0,0.0,0.0));

  const Real*
  RealVect::dataPtr() const
  {
    return vect;
  }

  Real*
  RealVect::dataPtr()
  {
    return vect;
  }

  RealVect::RealVect (AMREX_D_DECL6(Real i, Real j, Real k, Real m, Real n, Real p))
  {
    AMREX_D_EXPR6(vect[0] = i, vect[1] = j, vect[2] = k, vect[3] = m, vect[4] = n, vect[5] = p);
  }

  RealVect::RealVect (const std::vector<Real>& vr )
  {
    AMREX_D_EXPR6(vect[0]=vr[0], vect[1]=vr[1], vect[2] = vr[2], vect[3] = vr[3], vect[4] = vr[4], vect[5] = vr[5]);
  }

  RealVect::RealVect ()
  {
    AMREX_D_EXPR6(vect[0]=0.0, vect[1]=0.0, vect[2] = 0.0, vect[3] = 0.0, vect[4] = 0.0, vect[5] = 0.0);
  }

  RealVect&
  RealVect::operator= (const RealVect &iv)
  {
    AMREX_D_EXPR6(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2], vect[3]=iv.vect[3], vect[4]=iv.vect[4], vect[5]=iv.vect[5]);
    return *this;
  }

  Real RealVect::dotProduct(const RealVect& a_rhs) const
  {
    return AMREX_D_TERM6(vect[0]*a_rhs.vect[0], +
                         vect[1]*a_rhs.vect[1], +
                         vect[2]*a_rhs.vect[2], +
                         vect[3]*a_rhs.vect[3], +
                         vect[4]*a_rhs.vect[4], +
                         vect[5]*a_rhs.vect[5]);
  }

  bool
  RealVect::operator== (const RealVect& p) const
  {
    return AMREX_D_TERM6(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2], && vect[3] == p[3], && vect[4] == p[4], && vect[5] == p[5]);
  }

  bool
  RealVect::operator!= (const RealVect& p) const
  {
    return AMREX_D_TERM6(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2], || vect[3] != p[3], || vect[4] != p[4], || vect[5] != p[5]);
  }

  RealVect&
  RealVect::operator+= (Real s)
  {
    AMREX_D_EXPR6(vect[0] += s, vect[1] += s, vect[2] += s, vect[3] += s, vect[4] += s, vect[5] += s);
    return *this;
  }

  RealVect&
  RealVect::operator+= (const RealVect& p)
  {
    AMREX_D_EXPR6(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2], vect[3] += p[3], vect[4] += p[4], vect[5] += p[5]);
    return *this;
  }

  RealVect&
  RealVect::operator*= (Real s)
  {
    AMREX_D_EXPR6(vect[0] *= s, vect[1] *= s, vect[2] *= s, vect[3] *= s, vect[4] *= s, vect[5] *= s);
    return *this;
  }

  RealVect&
  RealVect::operator*= (const RealVect &p)
  {
    AMREX_D_EXPR6(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2], vect[3] *= p[3], vect[4] *= p[4], vect[5] *= p[5]);
    return *this;
  }

  RealVect
  RealVect::operator* (Real s) const
  {
    RealVect v(AMREX_D_DECL6(vect[0]*s, vect[1]*s, vect[2]*s, vect[3]*s, vect[4]*s, vect[5]*s));
    return v;
  }

  RealVect
  RealVect::operator- (Real s) const
  {
    RealVect v(AMREX_D_DECL6(vect[0]-s, vect[1]-s, vect[2]-s, vect[3]-s, vect[4]-s, vect[5]-s));
    return v;
  }

  RealVect
  RealVect::operator+ (Real s) const
  {
    RealVect v(AMREX_D_DECL6(vect[0]+s, vect[1]+s, vect[2]+s, vect[3]+s, vect[4]+s, vect[5]+s));
    return v;
  }

  RealVect&
  RealVect::operator/= (Real s)
  {
    AMREX_D_EXPR6(vect[0] /= s, vect[1] /= s, vect[2] /= s, vect[3] /= s, vect[4] /= s, vect[5] /= s);
    return *this;
  }

  RealVect&
  RealVect::operator/= (const RealVect& p)
  {
    AMREX_D_EXPR6(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2], vect[3] /= p[3], vect[4] /= p[4], vect[5] /= p[5]);
    return *this;
  }

  RealVect
  RealVect::operator/ (Real s) const
  {
    RealVect result( AMREX_D_DECL6( vect[0] / s, vect[1] / s, vect[2] / s, vect[3] / s, vect[4] / s, vect[5] / s));
    return result ;
  }

  int
  RealVect::minDir(const bool& a_doAbs) const
  {
    int mDir = 0;
    for (int idir=0; idir<SpaceDim; idir++)
      {
        if (a_doAbs)
          {
            if (std::abs(vect[idir]) < std::abs(vect[mDir]))
              {
                mDir = idir;
              }
          }
        else
          {
            if (vect[idir] < vect[mDir])
              {
                mDir = idir;
              }
          }
      }
    return mDir;
  }

  int
  RealVect::maxDir(const bool& a_doAbs) const
  {
    int mDir = 0;
    for (int idir=0; idir<SpaceDim; idir++)
      {
        if (a_doAbs)
          {
            if (std::abs(vect[idir]) > std::abs(vect[mDir]))
              {
                mDir = idir;
              }
          }
        else
          {
            if (vect[idir] > vect[mDir])
              {
                mDir = idir;
              }
          }
      }
    return mDir;
  }

  RealVect
  BASISREALV (int dir)
  {
    assert(dir >= 0 && dir < SpaceDim);
    RealVect tmp = RealVect::Zero ;
    tmp.vect[dir] = 1;
    return tmp;
  }

  RealVect
  operator/ (Real            s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(s/p[0], s/p[1], s/p[2], s/p[3], s/p[4], s/p[5]));
  }
  RealVect
  operator+ (Real            s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(p[0] + s, p[1] + s, p[2] + s, p[3] + s, p[4] + s, p[5] + s));
  }

  RealVect
  operator- (Real            s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(s - p[0], s - p[1], s - p[2], s - p[3], s - p[4], s - p[5]));
  }

  RealVect
  operator* (Real            s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(s * p[0], s * p[1], s * p[2], s * p[3], s * p[4], s * p[5]));
  }

  RealVect
  operator/ (const RealVect& s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(s[0] / p[0], s[1] /p[1], s[2] / p[2], s[3] / p[3], s[4] / p[4], s[5] / p[5]));
  }

  RealVect
  operator+ (const RealVect& s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(p[0] + s[0], p[1] +s[1], p[2] + s[2], p[3] + s[3], p[4] + s[4], p[5] + s[5]));
  }

  RealVect
  operator- (const RealVect& s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(s[0] - p[0], s[1] - p[1], s[2] - p[2], s[3] - p[3], s[4] - p[4], s[5] - p[5]));
  }

  RealVect
  operator* (const RealVect& s,
             const RealVect& p)
  {
    return RealVect(AMREX_D_DECL6(p[0] * s[0], p[1] *s[1], p[2] * s[2], p[3] * s[3], p[4] * s[4], p[5] * s[5]));
  }

  std::ostream&
  operator<< (std::ostream& ostr, const RealVect& p)
  {
    ostr << "(" << AMREX_D_TERM6 ( p[0] ,<< "," << p[1], << "," << p[2]) , << "," << p[3]) , << "," << p[4]) , << "," << p[5]) << ")" ;
    return ostr;
  }

  std::istream&
  operator>> (std::istream& is,
              RealVect&      iv)
  {
    is >> std::ws;
    char c;
    is >> c;

    if (c == '(')
    {
        AMREX_D_EXPR6(is >> iv[0],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[1],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[2],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[3],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[4],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[5]);
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        amrex::Error("operator>>(istream&,IntVect&): expected \'(\'");
    }

    if (is.fail())
        amrex::Error("operator>>(istream&,IntVect&) failed");

    return is;
}
} //namespace amrex
