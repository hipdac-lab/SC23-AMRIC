#ifndef AMREX_EB2_IF_BOX_H_
#define AMREX_EB2_IF_BOX_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <algorithm>
#include <limits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

class BoxIF
    : GPUable
{
public:

    //! inside: is the fluid inside the box?
    BoxIF (const RealArray& a_lo, const RealArray& a_hi, bool a_inside)
        : m_lo(makeXDim3(a_lo)),
          m_hi(makeXDim3(a_hi)),
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        Real r = std::numeric_limits<Real>::lowest();
        AMREX_D_TERM(r = amrex::max(r, x-m_hi.x, -(x-m_lo.x));,
                     r = amrex::max(r, y-m_hi.y, -(y-m_lo.y));,
                     r = amrex::max(r, z-m_hi.z, -(z-m_lo.z)););
        return r*m_sign;
    }

    inline Real operator() (const RealArray& p) const noexcept
    {
        return this->operator() (AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:

    XDim3     m_lo;
    XDim3     m_hi;
    //
    Real      m_sign;
};

}

#endif
