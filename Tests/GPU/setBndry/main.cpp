#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        IntVect n_cell;
        IntVect max_grid_size;
        int Ncomp = 1;
        int Nghost = 1;
        {
            ParmParse pp;

            pp.get("n_cell",n_cell);
            pp.get("max_grid_size",max_grid_size);
            pp.query("nghost",Nghost);
            pp.query("ncomp",Ncomp);
        }

        BoxArray ba;
        { 
            IntVect dom_lo(AMREX_D_DECL(          0,           0,           0));
            IntVect dom_hi(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1));
            Box domain(dom_lo, dom_hi);

            ba.define(domain);
            ba.maxSize(max_grid_size);

            amrex::Print() << " Using domain: " << domain
                           << " max_grid_size: " << max_grid_size
                           << " and Nghost: " << Nghost << std::endl;
        }
        DistributionMapping dm(ba);

        MultiFab ori;
        MultiFab check;
        MultiFab oriFused;
        MultiFab checkFused;
        MultiFab control;

        ori.define(ba, dm, Ncomp, Nghost);
        check.define(ba, dm, Ncomp, Nghost);
        oriFused.define(ba, dm, Ncomp, Nghost);
        checkFused.define(ba, dm, Ncomp, Nghost);
        control.define(ba, dm, Ncomp, Nghost);

        ori.setVal(0.0);
        check.setVal(0.0);
        oriFused.setVal(0.0);
        checkFused.setVal(0.0);
        control.setVal(-1.0);

        {
            BL_PROFILE_REGION("First run");
            ori.setBndry(1.0, 0, Ncomp);
            check.setBndryTestCell(1.0, 0, Ncomp);
            oriFused.setBndryFusedDiff(1.0, 0, Ncomp);
            checkFused.setBndryFusedTestCell(1.0, 0, Ncomp);

            // Fused without the launch to time setup overhead.
            control.setBndryFusedDiffEmpty(1.0, 0, Ncomp);
            control.setBndryFusedTestCellEmpty(1.0, 0, Ncomp);
        }

        {
            BL_PROFILE_REGION("Test");
            control.setVal(0.0);

            ori.setBndry(1.0, 0, Ncomp);
            check.setBndryTestCell(1.0, 0, Ncomp);
            oriFused.setBndryFusedDiff(1.0, 0, Ncomp);
            checkFused.setBndryFusedTestCell(1.0, 0, Ncomp);

            // Fused without the launch to time setup overhead.
            control.setBndryFusedDiffEmpty(1.0, 0, Ncomp);
            control.setBndryFusedTestCellEmpty(1.0, 0, Ncomp);
        }

        // Sum of all cells (Checks that boundary is set).
        for (MFIter mfi(ori); mfi.isValid(); ++mfi)
        {
            int idx = mfi.tileIndex();
            int ori_sum = 0.0;
            int check_sum = 0.0;
            int orif_sum = 0.0;
            int checkf_sum = 0.0;
 
            Box gbx = mfi.fabbox();
            FArrayBox& fab_ori = ori[mfi];
            FArrayBox& fab_check = check[mfi];
            FArrayBox& fab_orif = oriFused[mfi];
            FArrayBox& fab_checkf = checkFused[mfi];

            ori_sum = fab_ori.sum(gbx, 0, Ncomp);
            check_sum = fab_check.sum(gbx, 0, Ncomp);
            orif_sum = fab_orif.sum(gbx, 0, Ncomp);
            checkf_sum = fab_checkf.sum(gbx, 0, Ncomp);

            amrex::Print() << "Box " << idx << ":  " << ori_sum << " boundary cells set." << std::endl; 

            if (ori_sum != check_sum)
            {
                amrex::Print() << "Error in cell check on box " << idx << ": "
                               << check_sum << " != " << ori_sum << std::endl;
            }
            
            if (ori_sum != orif_sum)
            {
                amrex::Print() << "Error in original fused on box " << idx << ": "
                               << orif_sum << " != " << ori_sum << std::endl;
            }

            if (ori_sum != checkf_sum)
            {
                amrex::Print() << "Error in cell check fused on box " << idx << ": "
                               << checkf_sum << " != " << ori_sum << std::endl;
            }
        }

        // Sum of internal cells (Checks that the interior wasn't set.) 
        // Combined, these should confirm the correct result.
        {
            int ori_sum = ori.sum();
            int check_sum = check.sum();
            int orif_sum = oriFused.sum();
            int checkf_sum = checkFused.sum();

            amrex::Print() << "**Sum " << ":  " << ori_sum << " internal cells set." << std::endl; 

            if (ori_sum != check_sum)
            {
                amrex::Print() << "Error in cell check sum: "
                               << check_sum << " != " << ori_sum << std::endl;
            }
            
            if (ori_sum != orif_sum)
            {
                amrex::Print() << "Error in original fused sum: "
                               << orif_sum << " != " << ori_sum << std::endl;
            }

            if (ori_sum != checkf_sum)
            {
                amrex::Print() << "Error in cell check fused sum: "
                               << checkf_sum << " != " << ori_sum << std::endl;
            }
        }
    }
    amrex::Finalize();

    return 0;
}
