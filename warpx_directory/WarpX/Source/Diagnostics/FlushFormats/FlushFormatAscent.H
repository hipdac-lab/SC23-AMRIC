#ifndef WARPX_FLUSHFORMATASCENT_H_
#define WARPX_FLUSHFORMATASCENT_H_

#include "FlushFormat.H"

#include "Diagnostics/ParticleDiag/ParticleDiag_fwd.H"

#ifdef AMREX_USE_ASCENT
#   include <AMReX_Conduit_Blueprint.H>
#endif
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#ifdef AMREX_USE_ASCENT
#   include <ascent.hpp>
#endif

#include <string>

/**
 * \brief This class aims at dumping performing in-situ diagnostics with ASCENT.
 * In particular, function WriteToFile takes fields and particles as input arguments,
 * and calls amrex functions to do the in-situ visualization.
 */
class FlushFormatAscent : public FlushFormat
{
public:
    /** Do in-situ visualization for field and particle data */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
        const std::string prefix, int file_min_digits,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        const bool use_pinned_pc = false,
        bool isBTD = false, int snapshotID = -1,
        int bufferID = 1, int numBuffers = 1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false,
        const amrex::Vector<int>& totalParticlesFlushedAlready = amrex::Vector<int>() ) const override;

#ifdef AMREX_USE_ASCENT
    /** \brief Do in-situ visualization for particle data.
     * \param[in] particle_diags Each element of this vector handles output of 1 species.
     * \param[out] a_bp_mesh blueprint mesh generated from the container
     * Only compile if AMREX_USE_ASCENT because we need to pass a conduit class
     */
    void WriteParticles(const amrex::Vector<ParticleDiag>& particle_diags, conduit::Node& a_bp_mesh) const;
#endif

    ~FlushFormatAscent() {}
};

#endif // WARPX_FLUSHFORMATASCENT_H_
