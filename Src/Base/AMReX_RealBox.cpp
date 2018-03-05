
#include <iostream>
#include <string>

#include <AMReX_RealBox.H>

namespace amrex {

RealBox::RealBox (const Box&  bx,
                  const Real* dx,
                  const Real* base)
{
    const int* blo = bx.loVect();
    const int* bhi = bx.hiVect();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = base[i] + dx[i]*blo[i];
        int shft = (bx.type(i) == IndexType::CELL ? 1 : 0);
        xhi[i] = base[i] + dx[i]*(bhi[i]+ shft);
    }   
}

RealBox::RealBox ()
{
    AMREX_D_TERM6(xlo[0] , = xlo[1] , = xlo[2] , = xlo[3] , = xlo[4] , = xlo[5] ) = 0.;
    AMREX_D_TERM6(xhi[0] , = xhi[1] , = xhi[2] , = xhi[3] , = xhi[4] , = xhi[5] ) = -1.;
}

RealBox::RealBox (const Real* a_lo,
                  const Real* a_hi)
{
    AMREX_D_EXPR6( xlo[0] = a_lo[0] , xlo[1] = a_lo[1] , xlo[2] = a_lo[2],
                   xlo[3] = a_lo[3] , xlo[4] = a_lo[4] , xlo[5] = a_lo[5]);
    AMREX_D_EXPR6(xhi[0] = a_hi[0] , xhi[1] = a_hi[1] , xhi[2] = a_hi[2],
                  xhi[3] = a_hi[3] , xhi[4] = a_hi[4] , xhi[5] = a_hi[5]);
}

RealBox::RealBox (const std::array<Real,AMREX_SPACEDIM>& a_lo,
                  const std::array<Real,AMREX_SPACEDIM>& a_hi)
{
    AMREX_D_EXPR6(xlo[0] = a_lo[0] , xlo[1] = a_lo[1] , xlo[2] = a_lo[2],
                  xlo[3] = a_lo[3] , xlo[4] = a_lo[4] , xlo[5] = a_lo[5]);
    AMREX_D_EXPR6(xhi[0] = a_hi[0] , xhi[1] = a_hi[1] , xhi[2] = a_hi[2],
                  xhi[3] = a_hi[3] , xhi[4] = a_hi[4] , xhi[5] = a_hi[5]);
}

RealBox::RealBox (AMREX_D_DECL6(Real x0, Real y0, Real z0),
                  AMREX_D_DECL6(Real x1, Real y1, Real z1))
{
    AMREX_D_EXPR6(xlo[0] = x0 , xlo[1] = y0 , xlo[2] = z0),
                  xlo[3] = x0 , xlo[4] = y0 , xlo[5] = z0);
    AMREX_D_EXPR6(xhi[0] = x1 , xhi[1] = y1 , xhi[2] = z1,
                  xhi[3] = x1 , xhi[4] = y1 , xhi[5] = z1);
}

bool
RealBox::contains (const RealBox& rb, Real eps) const
{
    return contains(rb.xlo, eps) && contains(rb.xhi, eps);
}

bool
RealBox::ok () const
{
    return (length(0) >= 0.0)
#if (AMREX_SPACEDIM > 1)
        && (length(1) >= 0.0)
#endif   
#if (AMREX_SPACEDIM > 2)
        && (length(2) >= 0.0)
#endif
#if (AMREX_SPACEDIM > 3)
        && (length(3) >= 0.0)
#endif
#if (AMREX_SPACEDIM > 4)
        && (length(4) >= 0.0)
#endif
#if (AMREX_SPACEDIM > 5)
        && (length(5) >= 0.0)
#endif
   ;
}

Real
RealBox::volume () const
{
    if (ok()) return AMREX_D_TERM6(length(0), *length(1), *length(2), *length(3), *length(4), *length(5));
    return 0.0;
}

bool
RealBox::contains (const Real* point, Real eps) const
{
    return  AMREX_D_TERM6((xlo[0]-eps < point[0]) && (point[0] < xhi[0]+eps),
                       && (xlo[1]-eps < point[1]) && (point[1] < xhi[1]+eps),
                       && (xlo[2]-eps < point[2]) && (point[2] < xhi[2]+eps),
                       && (xlo[3]-eps < point[3]) && (point[3] < xhi[3]+eps),
                       && (xlo[4]-eps < point[4]) && (point[4] < xhi[4]+eps),
                       && (xlo[5]-eps < point[5]) && (point[5] < xhi[5]+eps));
}

std::ostream&
operator << (std::ostream &os, const RealBox& b)
{
    os << "(RealBox ";
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        os << b.lo(i) << ' ' << b.hi(i) << ' ';
    os << ')';
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator >> (std::istream &is, RealBox& b)
{
    is.ignore(BL_IGNORE_MAX,'(');

    std::string s;

    is >> s;

    if (s != "RealBox")
    {
        std::cerr << "unexpected token in RealBox: " << s << '\n';
        amrex::Abort();
    }

    Real lo[AMREX_SPACEDIM];
    Real hi[AMREX_SPACEDIM];
#ifdef BL_USE_FLOAT
    double dlotemp, dhitemp;
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        is >> dlotemp >> dhitemp;
        lo[i] = dlotemp;
        hi[i] = dhitemp;
    }
#else
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is >> lo[i] >> hi[i];
#endif

    is.ignore(BL_IGNORE_MAX, ')');

    b = RealBox(lo,hi);

    return is;
}

}
