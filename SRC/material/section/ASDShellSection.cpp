/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Original implementation: Massimo Petracca (ASDEA)
//
// See ASDShellSection.h for the design notes, and
// OpenSees-Testing/new-asd-elements/ASDShellSection/ for the design
// document and the numpy oracle that gates the section algebra.

#include <ASDShellSection.h>
#include <PlateRebarMaterial.h>
#include <PlateFiberMaterial.h>
#include <UniaxialMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>
#include <cstdlib>

namespace {

    // Gauss-Lobatto rules on [-1, 1], n = 3..10. Values validated by the
    // red-gate oracle (exact up to degree 2n-3). Symmetric halves.
    constexpr int LOBATTO_MIN = 3;
    constexpr int LOBATTO_MAX = 10;

    const double* lobattoNodes(int n)
    {
        static const double x3[] = { -1.0, 0.0, 1.0 };
        static const double x4[] = { -1.0, -0.4472135954999579, 0.4472135954999579, 1.0 };
        static const double x5[] = { -1.0, -0.6546536707079771, 0.0, 0.6546536707079771, 1.0 };
        static const double x6[] = { -1.0, -0.7650553239294647, -0.2852315164806451,
                                      0.2852315164806451, 0.7650553239294647, 1.0 };
        static const double x7[] = { -1.0, -0.8302238962785670, -0.4688487934707142, 0.0,
                                      0.4688487934707142, 0.8302238962785670, 1.0 };
        static const double x8[] = { -1.0, -0.8717401485096066, -0.5917001814331423, -0.2092992179024789,
                                      0.2092992179024789, 0.5917001814331423, 0.8717401485096066, 1.0 };
        static const double x9[] = { -1.0, -0.8997579954114602, -0.6771862795107377, -0.3631174638261782, 0.0,
                                      0.3631174638261782, 0.6771862795107377, 0.8997579954114602, 1.0 };
        static const double x10[] = { -1.0, -0.9195339081664589, -0.7387738651055050, -0.4779249498104445,
                                      -0.1652789576663870, 0.1652789576663870, 0.4779249498104445,
                                       0.7387738651055050, 0.9195339081664589, 1.0 };
        static const double* xx[] = { x3, x4, x5, x6, x7, x8, x9, x10 };
        return xx[n - LOBATTO_MIN];
    }

    const double* lobattoWeights(int n)
    {
        static const double w3[] = { 1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0 };
        static const double w4[] = { 1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0 };
        static const double w5[] = { 0.1, 0.5444444444444444, 0.7111111111111111,
                                     0.5444444444444444, 0.1 };
        static const double w6[] = { 1.0 / 15.0, 0.3784749562978470, 0.5548583770354864,
                                     0.5548583770354864, 0.3784749562978470, 1.0 / 15.0 };
        static const double w7[] = { 1.0 / 21.0, 0.2768260473615659, 0.4317453812098627, 0.4876190476190476,
                                     0.4317453812098627, 0.2768260473615659, 1.0 / 21.0 };
        static const double w8[] = { 1.0 / 28.0, 0.2107042271435061, 0.3411226924835044, 0.4124587946587038,
                                     0.4124587946587038, 0.3411226924835044, 0.2107042271435061, 1.0 / 28.0 };
        static const double w9[] = { 1.0 / 36.0, 0.1654953615608055, 0.2745387125001617, 0.3464285109730463,
                                     0.3715192743764172, 0.3464285109730463, 0.2745387125001617,
                                     0.1654953615608055, 1.0 / 36.0 };
        static const double w10[] = { 1.0 / 45.0, 0.1333059908510701, 0.2248893420631264, 0.2920426836796838,
                                      0.3275397611838976, 0.3275397611838976, 0.2920426836796838,
                                      0.2248893420631264, 0.1333059908510701, 1.0 / 45.0 };
        static const double* ww[] = { w3, w4, w5, w6, w7, w8, w9, w10 };
        return ww[n - LOBATTO_MIN];
    }

    // reads the next input arg as a double only if it parses as one;
    // otherwise puts it back and returns false
    bool tryReadDouble(double& value)
    {
        if (OPS_GetNumRemainingInputArgs() < 1)
            return false;
        const char* s = OPS_GetString();
        char* end = nullptr;
        double v = strtod(s, &end);
        if (end == s || *end != '\0') {
            OPS_ResetCurrentInputArg(-1);
            return false;
        }
        value = v;
        return true;
    }

    bool tryReadInt(int& value)
    {
        double v;
        if (!tryReadDouble(v))
            return false;
        value = static_cast<int>(v);
        if (static_cast<double>(value) != v) {
            OPS_ResetCurrentInputArg(-1);
            return false;
        }
        return true;
    }

}

void* OPS_ASDShellSection()
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDShellSection - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr =
        "Want: section ASDShellSection $tag \\\n"
        "    -ply $ndMatTag $thickness <$nip> ... (>= 1 ply, bottom-to-top; nip in 3..10, default 5)\n"
        "    <-gap $thickness>                    (spacer between plies)\n"
        "    <-rebar $uniMatTag $dia $spacing $angle $z>\n"
        "    <-rebarArea $uniMatTag $areaPerWidth $angle $z>\n"
        "    <-rebarTop $uniMatTag $dia $spacing $angle $cover>\n"
        "    <-rebarBottom $uniMatTag $dia $spacing $angle $cover>\n"
        "                                         (z from the stack midplane; cover to the LAYER\n"
        "                                          CENTER from the stack face; angle in degrees)\n"
        "    <-offset $z0>                        (stack midplane position in the element frame)\n"
        "    <-shearCorrection $k>                (default 5/6; LayeredShell applies none: k = 1 reproduces it)\n"
        "    <-elasticShear <$S1 $S2>>            (freeze the transverse shear channels: no args = from\n"
        "                                          the initial ply tangents, or the two section stiffnesses)\n"
        "    <-rho $extraMassPerArea>\n";

    if (OPS_GetNumRemainingInputArgs() < 4) {
        opserr << "ASDShellSection ERROR: few arguments.\n" << descr;
        return 0;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
        opserr << "ASDShellSection ERROR: invalid tag.\n" << descr;
        return 0;
    }

    std::vector<ASDShellSection::StackItem> stack;
    std::vector<ASDShellSection::RebarItem> rebars;
    // -rebarTop/-rebarBottom are resolved after the stack is known
    struct PendingRebar { ASDShellSection::RebarItem item; double cover; bool top; };
    std::vector<PendingRebar> pending;
    double offset = 0.0;
    double k_shear = 5.0 / 6.0;
    int shear_mode = ASDShellSection::Shear_Integrated;
    double S1 = 0.0, S2 = 0.0;
    double rho_extra = 0.0;

    auto readRebarCommon = [&](const char* kw, ASDShellSection::RebarItem& item, bool by_area) -> bool {
        int matTag;
        int numData = 1;
        if (OPS_GetIntInput(&numData, &matTag) < 0) {
            opserr << "ASDShellSection ERROR: invalid material tag after " << kw << ".\n" << descr;
            return false;
        }
        item.mat = OPS_getUniaxialMaterial(matTag);
        if (item.mat == nullptr) {
            opserr << "ASDShellSection ERROR: uniaxial material " << matTag << " (" << kw << ") does not exist.\n";
            return false;
        }
        if (by_area) {
            double apw;
            if (OPS_GetDoubleInput(&numData, &apw) < 0 || apw <= 0.0) {
                opserr << "ASDShellSection ERROR: invalid area-per-unit-width after " << kw << ".\n" << descr;
                return false;
            }
            item.teq = apw;
        }
        else {
            double dia, spacing;
            if (OPS_GetDoubleInput(&numData, &dia) < 0 || dia <= 0.0) {
                opserr << "ASDShellSection ERROR: invalid diameter after " << kw << ".\n" << descr;
                return false;
            }
            if (OPS_GetDoubleInput(&numData, &spacing) < 0 || spacing <= 0.0) {
                opserr << "ASDShellSection ERROR: invalid spacing after " << kw << ".\n" << descr;
                return false;
            }
            item.teq = 0.25 * 3.141592653589793238462643 * dia * dia / spacing;
        }
        if (OPS_GetDoubleInput(&numData, &item.angle) < 0) {
            opserr << "ASDShellSection ERROR: invalid angle after " << kw << ".\n" << descr;
            return false;
        }
        return true;
    };

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* what = OPS_GetString();
        numData = 1;
        if (strcmp(what, "-ply") == 0) {
            ASDShellSection::StackItem item;
            int matTag;
            if (OPS_GetIntInput(&numData, &matTag) < 0) {
                opserr << "ASDShellSection ERROR: invalid material tag after -ply.\n" << descr;
                return 0;
            }
            item.mat = OPS_getNDMaterial(matTag);
            if (item.mat == nullptr) {
                opserr << "ASDShellSection ERROR: nD material " << matTag << " does not exist.\n";
                return 0;
            }
            if (OPS_GetDoubleInput(&numData, &item.t) < 0 || item.t <= 0.0) {
                opserr << "ASDShellSection ERROR: invalid thickness after -ply.\n" << descr;
                return 0;
            }
            item.nip = 5;
            int nip;
            if (tryReadInt(nip)) {
                if (nip < LOBATTO_MIN || nip > LOBATTO_MAX) {
                    opserr << "ASDShellSection ERROR: nip must be in " << LOBATTO_MIN << ".." << LOBATTO_MAX
                        << " (Lobatto n = 2 is the trapezoid and misintegrates z^2), got " << nip << ".\n" << descr;
                    return 0;
                }
                item.nip = nip;
            }
            stack.push_back(item);
        }
        else if (strcmp(what, "-gap") == 0) {
            ASDShellSection::StackItem item;
            if (OPS_GetDoubleInput(&numData, &item.t) < 0 || item.t <= 0.0) {
                opserr << "ASDShellSection ERROR: invalid thickness after -gap.\n" << descr;
                return 0;
            }
            item.nip = -1;
            stack.push_back(item);
        }
        else if (strcmp(what, "-rebar") == 0 || strcmp(what, "-rebarArea") == 0) {
            ASDShellSection::RebarItem item;
            if (!readRebarCommon(what, item, strcmp(what, "-rebarArea") == 0))
                return 0;
            if (OPS_GetDoubleInput(&numData, &item.z) < 0) {
                opserr << "ASDShellSection ERROR: invalid z after " << what << ".\n" << descr;
                return 0;
            }
            rebars.push_back(item);
        }
        else if (strcmp(what, "-rebarTop") == 0 || strcmp(what, "-rebarBottom") == 0) {
            PendingRebar p;
            p.top = (strcmp(what, "-rebarTop") == 0);
            if (!readRebarCommon(what, p.item, false))
                return 0;
            if (OPS_GetDoubleInput(&numData, &p.cover) < 0 || p.cover < 0.0) {
                opserr << "ASDShellSection ERROR: invalid cover after " << what << ".\n" << descr;
                return 0;
            }
            pending.push_back(p);
        }
        else if (strcmp(what, "-offset") == 0) {
            if (OPS_GetDoubleInput(&numData, &offset) < 0) {
                opserr << "ASDShellSection ERROR: invalid value after -offset.\n" << descr;
                return 0;
            }
        }
        else if (strcmp(what, "-shearCorrection") == 0) {
            if (OPS_GetDoubleInput(&numData, &k_shear) < 0 || k_shear <= 0.0) {
                opserr << "ASDShellSection ERROR: invalid value after -shearCorrection.\n" << descr;
                return 0;
            }
        }
        else if (strcmp(what, "-elasticShear") == 0) {
            if (tryReadDouble(S1)) {
                if (OPS_GetDoubleInput(&numData, &S2) < 0) {
                    opserr << "ASDShellSection ERROR: -elasticShear wants zero or two values.\n" << descr;
                    return 0;
                }
                if (S1 <= 0.0 || S2 <= 0.0) {
                    opserr << "ASDShellSection ERROR: -elasticShear stiffnesses must be positive.\n" << descr;
                    return 0;
                }
                shear_mode = ASDShellSection::Shear_ElasticUser;
            }
            else {
                shear_mode = ASDShellSection::Shear_ElasticAuto;
            }
        }
        else if (strcmp(what, "-rho") == 0) {
            if (OPS_GetDoubleInput(&numData, &rho_extra) < 0 || rho_extra < 0.0) {
                opserr << "ASDShellSection ERROR: invalid value after -rho.\n" << descr;
                return 0;
            }
        }
        else {
            opserr << "ASDShellSection ERROR: unknown keyword \"" << what << "\".\n" << descr;
            return 0;
        }
    }

    // checks on the stack
    int nply = 0;
    double H = 0.0;
    for (const auto& it : stack) {
        H += it.t;
        if (it.nip > 0)
            ++nply;
    }
    if (nply < 1) {
        opserr << "ASDShellSection ERROR: at least one -ply is required (a rebar-only "
            << "section has a singular transverse shear block).\n" << descr;
        return 0;
    }
    if (!stack.empty() && (stack.front().nip < 0 || stack.back().nip < 0)) {
        opserr << "ASDShellSection ERROR: the stack cannot start or end with a -gap.\n" << descr;
        return 0;
    }

    // resolve the -rebarTop/-rebarBottom covers (to the layer center)
    for (auto& p : pending) {
        p.item.z = p.top ? (0.5 * H - p.cover) : (-0.5 * H + p.cover);
        rebars.push_back(p.item);
    }

    // a rebar far outside the stack envelope is probably a typo, but it can
    // be intentional (external strengthening): warn, do not fail
    for (const auto& r : rebars) {
        if (std::abs(r.z) > 0.5 * H) {
            opserr << "ASDShellSection WARNING: section " << tag << " - a rebar layer sits at z = "
                << r.z << ", outside the stack envelope [" << -0.5 * H << ", " << 0.5 * H << "].\n";
        }
    }

    if (shear_mode == ASDShellSection::Shear_Integrated || shear_mode == ASDShellSection::Shear_ElasticAuto) {
        // nothing more to resolve here: the auto values need the wrapped
        // fibers and are computed in the constructor
    }

    return new ASDShellSection(tag, stack, rebars, offset, k_shear, shear_mode, S1, S2, rho_extra);
}

// static members
Vector ASDShellSection::stressResultant(8);
Matrix ASDShellSection::tangent(8, 8);
ID ASDShellSection::array(8);

ASDShellSection::ASDShellSection()
    : SectionForceDeformation(0, SEC_TAG_ASDShellSection)
    , strainResultant(8)
{
}

ASDShellSection::ASDShellSection(int tag,
    const std::vector<StackItem>& stack,
    const std::vector<RebarItem>& rebars,
    double offset,
    double k_shear,
    int shear_mode,
    double S1, double S2,
    double rho_extra)
    : SectionForceDeformation(tag, SEC_TAG_ASDShellSection)
    , m_offset(offset)
    , m_k(k_shear)
    , m_shear_mode(shear_mode)
    , m_S1(S1)
    , m_S2(S2)
    , m_rho_extra(rho_extra)
    , strainResultant(8)
{
    buildFibers(stack, rebars);
}

void ASDShellSection::buildFibers(const std::vector<StackItem>& stack,
    const std::vector<RebarItem>& rebars)
{
    // stack description
    m_h = 0.0;
    m_stack_t.clear();
    m_stack_nip.clear();
    for (const auto& it : stack) {
        m_stack_t.push_back(it.t);
        m_stack_nip.push_back(it.mat ? it.nip : -1);
        m_h += it.t;
    }
    // ply fibers, bottom-to-top
    m_fib_z.clear();
    m_fib_w.clear();
    m_fib_mat.clear();
    m_fib_rebar.clear();
    double zb = -0.5 * m_h;
    for (const auto& it : stack) {
        if (it.mat == nullptr) {
            zb += it.t;
            continue;
        }
        const double* x = lobattoNodes(it.nip);
        const double* w = lobattoWeights(it.nip);
        double zm = zb + 0.5 * it.t;
        for (int i = 0; i < it.nip; ++i) {
            NDMaterial* fib = it.mat->getCopy("PlateFiber");
            if (fib == nullptr && strcmp(it.mat->getType(), "ThreeDimensional") == 0) {
                // a plain 3D material: wrap it in the sigma_zz = 0 condenser
                fib = new PlateFiberMaterial(0, *it.mat);
            }
            if (fib == nullptr) {
                opserr << "ASDShellSection ERROR: section " << getTag() << " - nD material "
                    << it.mat->getTag() << " cannot return a PlateFiber material.\n";
                exit(-1);
            }
            m_fib_z.push_back(zm + 0.5 * it.t * x[i]);
            m_fib_w.push_back(0.5 * it.t * w[i]);
            m_fib_mat.push_back(fib);
            m_fib_rebar.push_back(0);
        }
        zb += it.t;
    }
    // rebar fibers
    m_reb_teq.clear();
    m_reb_angle.clear();
    m_reb_z.clear();
    for (const auto& r : rebars) {
        m_reb_teq.push_back(r.teq);
        m_reb_angle.push_back(r.angle);
        m_reb_z.push_back(r.z);
        m_fib_z.push_back(r.z);
        m_fib_w.push_back(r.teq);
        m_fib_mat.push_back(new PlateRebarMaterial(0, *r.mat, r.angle));
        m_fib_rebar.push_back(1);
    }
    // Shear_ElasticAuto: S1/S2 are NOT resolved here. See ensureAutoShear.
}

void ASDShellSection::ensureAutoShear()
{
    if (m_shear_mode != Shear_ElasticAuto || m_S_computed)
        return;
    m_S1 = m_S2 = 0.0;
    for (std::size_t i = 0; i < m_fib_mat.size(); ++i) {
        if (m_fib_rebar[i])
            continue;
        const Matrix& D0 = m_fib_mat[i]->getInitialTangent();
        m_S1 += m_k * m_fib_w[i] * D0(4, 4);
        m_S2 += m_k * m_fib_w[i] * D0(3, 3);
    }
    m_S_computed = true;
}

ASDShellSection::~ASDShellSection()
{
    for (NDMaterial* m : m_fib_mat)
        delete m;
}

SectionForceDeformation* ASDShellSection::getCopy()
{
    ASDShellSection* clone = new ASDShellSection();
    clone->setTag(getTag());
    clone->m_stack_t = m_stack_t;
    clone->m_stack_nip = m_stack_nip;
    clone->m_reb_teq = m_reb_teq;
    clone->m_reb_angle = m_reb_angle;
    clone->m_reb_z = m_reb_z;
    clone->m_fib_z = m_fib_z;
    clone->m_fib_w = m_fib_w;
    clone->m_fib_rebar = m_fib_rebar;
    clone->m_fib_mat.resize(m_fib_mat.size(), nullptr);
    for (std::size_t i = 0; i < m_fib_mat.size(); ++i)
        clone->m_fib_mat[i] = m_fib_mat[i]->getCopy();
    clone->m_h = m_h;
    clone->m_offset = m_offset;
    clone->m_k = m_k;
    clone->m_shear_mode = m_shear_mode;
    clone->m_S1 = m_S1;
    clone->m_S2 = m_S2;
    // Shear_ElasticAuto: force the copy to re-resolve on its own fibers,
    // which will carry the lch of the element that owns the copy
    clone->m_S_computed = (m_shear_mode == Shear_ElasticAuto) ? false : m_S_computed;
    clone->m_rho_extra = m_rho_extra;
    clone->strainResultant = strainResultant;
    return clone;
}

double ASDShellSection::getRho()
{
    double rhoH = m_rho_extra;
    for (std::size_t i = 0; i < m_fib_mat.size(); ++i)
        rhoH += m_fib_mat[i]->getRho() * m_fib_w[i];
    return rhoH;
}

int ASDShellSection::getOrder() const
{
    return 8;
}

const ID& ASDShellSection::getType()
{
    static bool initialized = false;
    if (!initialized) {
        array(0) = SECTION_RESPONSE_FXX;
        array(1) = SECTION_RESPONSE_FYY;
        array(2) = SECTION_RESPONSE_FXY;
        array(3) = SECTION_RESPONSE_MXX;
        array(4) = SECTION_RESPONSE_MYY;
        array(5) = SECTION_RESPONSE_MXY;
        array(6) = SECTION_RESPONSE_VXZ;
        array(7) = SECTION_RESPONSE_VYZ;
        initialized = true;
    }
    return array;
}

int ASDShellSection::commitState()
{
    int success = 0;
    for (NDMaterial* m : m_fib_mat)
        success += m->commitState();
    return success;
}

int ASDShellSection::revertToLastCommit()
{
    int success = 0;
    for (NDMaterial* m : m_fib_mat)
        success += m->revertToLastCommit();
    return success;
}

int ASDShellSection::revertToStart()
{
    int success = 0;
    for (NDMaterial* m : m_fib_mat)
        success += m->revertToStart();
    return success;
}

int ASDShellSection::setTrialSectionDeformation(const Vector& strain_from_element)
{
    ensureAutoShear();
    strainResultant = strain_from_element;

    // strain mapping identical to LayeredShellFiberSection: the fiber gets
    // [e1 - zeta k1, e2 - zeta k2, g12 - zeta k12, rk*sr(7), rk*sr(6)],
    // with zeta = z + offset. In the frozen-shear modes no transverse
    // shear strain reaches the fibers.
    static Vector strain(5);
    int success = 0;
    double rk = (m_shear_mode == Shear_Integrated) ? std::sqrt(m_k) : 0.0;
    for (std::size_t i = 0; i < m_fib_mat.size(); ++i) {
        double zeta = m_fib_z[i] + m_offset;
        strain(0) = strainResultant(0) - zeta * strainResultant(3);
        strain(1) = strainResultant(1) - zeta * strainResultant(4);
        strain(2) = strainResultant(2) - zeta * strainResultant(5);
        strain(4) = rk * strainResultant(6);
        strain(3) = rk * strainResultant(7);
        success += m_fib_mat[i]->setTrialStrain(strain);
    }
    return success;
}

const Vector& ASDShellSection::getSectionDeformation()
{
    return strainResultant;
}

const Vector& ASDShellSection::getStressResultant()
{
    ensureAutoShear();
    static Vector stress(5);
    stressResultant.Zero();
    double rk = (m_shear_mode == Shear_Integrated) ? std::sqrt(m_k) : 0.0;
    for (std::size_t i = 0; i < m_fib_mat.size(); ++i) {
        double zeta = m_fib_z[i] + m_offset;
        double weight = m_fib_w[i];
        stress = m_fib_mat[i]->getStress();
        // membrane
        stressResultant(0) += stress(0) * weight;
        stressResultant(1) += stress(1) * weight;
        stressResultant(2) += stress(2) * weight;
        // bending moments (the LayeredShell sign convention: m = +int(zeta*sigma))
        stressResultant(3) += (zeta * stress(0)) * weight;
        stressResultant(4) += (zeta * stress(1)) * weight;
        stressResultant(5) += (zeta * stress(2)) * weight;
        // shear
        stressResultant(6) += rk * stress(4) * weight;
        stressResultant(7) += rk * stress(3) * weight;
    }
    if (m_shear_mode != Shear_Integrated) {
        stressResultant(6) = m_S1 * strainResultant(6);
        stressResultant(7) = m_S2 * strainResultant(7);
    }
    return stressResultant;
}

const Matrix& ASDShellSection::getSectionTangent()
{
    ensureAutoShear();
    tangent.Zero();
    bool integrated = (m_shear_mode == Shear_Integrated);
    double rk = integrated ? std::sqrt(m_k) : 0.0;
    for (std::size_t i = 0; i < m_fib_mat.size(); ++i) {
        double zeta = m_fib_z[i] + m_offset;
        const Matrix& dd0 = m_fib_mat[i]->getTangent();
        static Matrix dd(5, 5);
        dd = dd0;
        dd *= m_fib_w[i];
        // resultant row map for stress component a: (row, factor) pairs;
        // strain column map for component b is the same with the LayeredShell
        // signs: membrane (b, 1) and bending (b+3, -zeta); shear channels
        // scaled by rk on both sides
        for (int a = 0; a < 5; ++a) {
            int rowA1 = -1, rowA2 = -1;
            double facA1 = 0.0, facA2 = 0.0;
            if (a < 3) {
                rowA1 = a;      facA1 = 1.0;
                rowA2 = a + 3;  facA2 = zeta;
            }
            else if (integrated) {
                rowA1 = (a == 4) ? 6 : 7;
                facA1 = rk;
            }
            else {
                continue;
            }
            for (int b = 0; b < 5; ++b) {
                int colB1 = -1, colB2 = -1;
                double facB1 = 0.0, facB2 = 0.0;
                if (b < 3) {
                    colB1 = b;      facB1 = 1.0;
                    colB2 = b + 3;  facB2 = -zeta;
                }
                else if (integrated) {
                    colB1 = (b == 4) ? 6 : 7;
                    facB1 = rk;
                }
                else {
                    continue;
                }
                double d = dd(a, b);
                tangent(rowA1, colB1) += facA1 * d * facB1;
                if (colB2 >= 0)
                    tangent(rowA1, colB2) += facA1 * d * facB2;
                if (rowA2 >= 0) {
                    tangent(rowA2, colB1) += facA2 * d * facB1;
                    if (colB2 >= 0)
                        tangent(rowA2, colB2) += facA2 * d * facB2;
                }
            }
        }
    }
    if (!integrated) {
        tangent(6, 6) = m_S1;
        tangent(7, 7) = m_S2;
    }
    return tangent;
}

const Matrix& ASDShellSection::getInitialTangent()
{
    // same assembly on the initial fiber tangents. The shear block of the
    // frozen modes is initial by construction.
    ensureAutoShear();
    tangent.Zero();
    bool integrated = (m_shear_mode == Shear_Integrated);
    double rk = integrated ? std::sqrt(m_k) : 0.0;
    for (std::size_t i = 0; i < m_fib_mat.size(); ++i) {
        double zeta = m_fib_z[i] + m_offset;
        const Matrix& dd0 = m_fib_mat[i]->getInitialTangent();
        double w = m_fib_w[i];
        for (int a = 0; a < 5; ++a) {
            int rowA1 = -1, rowA2 = -1;
            double facA1 = 0.0, facA2 = 0.0;
            if (a < 3) {
                rowA1 = a;      facA1 = 1.0;
                rowA2 = a + 3;  facA2 = zeta;
            }
            else if (integrated) {
                rowA1 = (a == 4) ? 6 : 7;
                facA1 = rk;
            }
            else {
                continue;
            }
            for (int b = 0; b < 5; ++b) {
                int colB1 = -1, colB2 = -1;
                double facB1 = 0.0, facB2 = 0.0;
                if (b < 3) {
                    colB1 = b;      facB1 = 1.0;
                    colB2 = b + 3;  facB2 = -zeta;
                }
                else if (integrated) {
                    colB1 = (b == 4) ? 6 : 7;
                    facB1 = rk;
                }
                else {
                    continue;
                }
                double d = w * dd0(a, b);
                tangent(rowA1, colB1) += facA1 * d * facB1;
                if (colB2 >= 0)
                    tangent(rowA1, colB2) += facA1 * d * facB2;
                if (rowA2 >= 0) {
                    tangent(rowA2, colB1) += facA2 * d * facB1;
                    if (colB2 >= 0)
                        tangent(rowA2, colB2) += facA2 * d * facB2;
                }
            }
        }
    }
    if (!integrated) {
        tangent(6, 6) = m_S1;
        tangent(7, 7) = m_S2;
    }
    return tangent;
}

void ASDShellSection::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ASDShellSection\", ";
        s << "\"totalHeight\": " << m_h << ", ";
        s << "\"offset\": " << m_offset << ", ";
        s << "\"nFibers\": " << static_cast<int>(m_fib_mat.size());
        s << "}";
        return;
    }
    s << "ASDShellSection tag: " << this->getTag() << endln;
    s << "  total stack height: " << m_h << ", offset: " << m_offset << endln;
    s << "  shear: " << (m_shear_mode == Shear_Integrated ?
        "integrated from the plies" : "frozen elastic")
        << ", k = " << m_k;
    if (m_shear_mode != Shear_Integrated)
        s << ", S = [" << m_S1 << ", " << m_S2 << "]";
    s << endln;
    double zb = -0.5 * m_h;
    int fib = 0;
    for (std::size_t i = 0; i < m_stack_t.size(); ++i) {
        if (m_stack_nip[i] < 0) {
            s << "  gap  t = " << m_stack_t[i] << endln;
        }
        else {
            s << "  ply  t = " << m_stack_t[i] << ", nip = " << m_stack_nip[i]
                << ", z in [" << zb << ", " << zb + m_stack_t[i] << "], material "
                << m_fib_mat[static_cast<std::size_t>(fib)]->getTag() << endln;
            fib += m_stack_nip[i];
        }
        zb += m_stack_t[i];
    }
    for (std::size_t i = 0; i < m_reb_teq.size(); ++i) {
        s << "  rebar  t_eq = " << m_reb_teq[i] << ", angle = " << m_reb_angle[i]
            << " deg, z = " << m_reb_z[i] << endln;
    }
}

int ASDShellSection::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;
    int dataTag = getDbTag();

    int nStack = static_cast<int>(m_stack_t.size());
    int nReb = static_cast<int>(m_reb_teq.size());
    int nFib = static_cast<int>(m_fib_mat.size());

    // INT data 1: fixed-size header
    static ID idData1(6);
    idData1(0) = getTag();
    idData1(1) = nStack;
    idData1(2) = nReb;
    idData1(3) = nFib;
    idData1(4) = m_shear_mode;
    idData1(5) = m_S_computed ? 1 : 0;
    res = theChannel.sendID(dataTag, commitTag, idData1);
    if (res < 0) {
        opserr << "WARNING ASDShellSection::sendSelf() - " << getTag() << " failed to send ID 1\n";
        return res;
    }

    // INT data 2: stack nips, fiber flags, fiber material (classTag, dbTag)
    ID idData2(nStack + 3 * nFib);
    int pos = 0;
    for (int i = 0; i < nStack; ++i)
        idData2(pos++) = m_stack_nip[static_cast<std::size_t>(i)];
    for (int i = 0; i < nFib; ++i) {
        NDMaterial* m = m_fib_mat[static_cast<std::size_t>(i)];
        int matDbTag = m->getDbTag();
        if (matDbTag == 0) {
            matDbTag = theChannel.getDbTag();
            if (matDbTag != 0)
                m->setDbTag(matDbTag);
        }
        idData2(pos++) = m->getClassTag();
        idData2(pos++) = matDbTag;
        idData2(pos++) = m_fib_rebar[static_cast<std::size_t>(i)];
    }
    res = theChannel.sendID(dataTag, commitTag, idData2);
    if (res < 0) {
        opserr << "WARNING ASDShellSection::sendSelf() - " << getTag() << " failed to send ID 2\n";
        return res;
    }

    // DOUBLE data: scalars, stack thicknesses, rebar data, fiber geometry
    Vector vectData(7 + nStack + 3 * nReb + 2 * nFib);
    pos = 0;
    vectData(pos++) = m_h;
    vectData(pos++) = m_offset;
    vectData(pos++) = m_k;
    vectData(pos++) = m_S1;
    vectData(pos++) = m_S2;
    vectData(pos++) = m_rho_extra;
    vectData(pos++) = 0.0; // reserved
    for (int i = 0; i < nStack; ++i)
        vectData(pos++) = m_stack_t[static_cast<std::size_t>(i)];
    for (int i = 0; i < nReb; ++i) {
        vectData(pos++) = m_reb_teq[static_cast<std::size_t>(i)];
        vectData(pos++) = m_reb_angle[static_cast<std::size_t>(i)];
        vectData(pos++) = m_reb_z[static_cast<std::size_t>(i)];
    }
    for (int i = 0; i < nFib; ++i) {
        vectData(pos++) = m_fib_z[static_cast<std::size_t>(i)];
        vectData(pos++) = m_fib_w[static_cast<std::size_t>(i)];
    }
    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDShellSection::sendSelf() - " << getTag() << " failed to send Vector\n";
        return res;
    }

    // the fiber materials
    for (int i = 0; i < nFib; ++i) {
        res = m_fib_mat[static_cast<std::size_t>(i)]->sendSelf(commitTag, theChannel);
        if (res < 0) {
            opserr << "WARNING ASDShellSection::sendSelf() - " << getTag()
                << " failed to send fiber material " << i + 1 << "\n";
            return res;
        }
    }

    return res;
}

int ASDShellSection::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;
    int dataTag = getDbTag();

    // INT data 1: header
    static ID idData1(6);
    res = theChannel.recvID(dataTag, commitTag, idData1);
    if (res < 0) {
        opserr << "WARNING ASDShellSection::recvSelf() - failed to receive ID 1\n";
        return res;
    }
    setTag(idData1(0));
    int nStack = idData1(1);
    int nReb = idData1(2);
    int nFib = idData1(3);
    m_shear_mode = idData1(4);
    m_S_computed = idData1(5) == 1;

    // INT data 2
    ID idData2(nStack + 3 * nFib);
    res = theChannel.recvID(dataTag, commitTag, idData2);
    if (res < 0) {
        opserr << "WARNING ASDShellSection::recvSelf() - failed to receive ID 2\n";
        return res;
    }
    int pos = 0;
    m_stack_nip.resize(static_cast<std::size_t>(nStack));
    for (int i = 0; i < nStack; ++i)
        m_stack_nip[static_cast<std::size_t>(i)] = idData2(pos++);
    std::vector<int> fibClassTags(static_cast<std::size_t>(nFib));
    std::vector<int> fibDbTags(static_cast<std::size_t>(nFib));
    m_fib_rebar.resize(static_cast<std::size_t>(nFib));
    for (int i = 0; i < nFib; ++i) {
        fibClassTags[static_cast<std::size_t>(i)] = idData2(pos++);
        fibDbTags[static_cast<std::size_t>(i)] = idData2(pos++);
        m_fib_rebar[static_cast<std::size_t>(i)] = static_cast<char>(idData2(pos++));
    }

    // DOUBLE data
    Vector vectData(7 + nStack + 3 * nReb + 2 * nFib);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDShellSection::recvSelf() - failed to receive Vector\n";
        return res;
    }
    pos = 0;
    m_h = vectData(pos++);
    m_offset = vectData(pos++);
    m_k = vectData(pos++);
    m_S1 = vectData(pos++);
    m_S2 = vectData(pos++);
    m_rho_extra = vectData(pos++);
    pos++; // reserved
    m_stack_t.resize(static_cast<std::size_t>(nStack));
    for (int i = 0; i < nStack; ++i)
        m_stack_t[static_cast<std::size_t>(i)] = vectData(pos++);
    m_reb_teq.resize(static_cast<std::size_t>(nReb));
    m_reb_angle.resize(static_cast<std::size_t>(nReb));
    m_reb_z.resize(static_cast<std::size_t>(nReb));
    for (int i = 0; i < nReb; ++i) {
        m_reb_teq[static_cast<std::size_t>(i)] = vectData(pos++);
        m_reb_angle[static_cast<std::size_t>(i)] = vectData(pos++);
        m_reb_z[static_cast<std::size_t>(i)] = vectData(pos++);
    }
    m_fib_z.resize(static_cast<std::size_t>(nFib));
    m_fib_w.resize(static_cast<std::size_t>(nFib));
    for (int i = 0; i < nFib; ++i) {
        m_fib_z[static_cast<std::size_t>(i)] = vectData(pos++);
        m_fib_w[static_cast<std::size_t>(i)] = vectData(pos++);
    }

    // the fiber materials, rebuilt through the broker when needed
    for (NDMaterial* m : m_fib_mat)
        delete m;
    m_fib_mat.assign(static_cast<std::size_t>(nFib), nullptr);
    for (int i = 0; i < nFib; ++i) {
        NDMaterial* m = theBroker.getNewNDMaterial(fibClassTags[static_cast<std::size_t>(i)]);
        if (m == nullptr) {
            opserr << "WARNING ASDShellSection::recvSelf() - could not get an NDMaterial with classTag "
                << fibClassTags[static_cast<std::size_t>(i)] << "\n";
            return -1;
        }
        m->setDbTag(fibDbTags[static_cast<std::size_t>(i)]);
        res = m->recvSelf(commitTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "WARNING ASDShellSection::recvSelf() - fiber material " << i + 1
                << " failed to recvSelf\n";
            return res;
        }
        m_fib_mat[static_cast<std::size_t>(i)] = m;
    }

    return res;
}

Response* ASDShellSection::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    Response* theResponse = 0;
    int nFib = static_cast<int>(m_fib_mat.size());

    if (strcmp(argv[0], "fiber") == 0 || strcmp(argv[0], "Fiber") == 0) {
        // fiber $i ... : 1-based over ALL fibers (ply points first,
        // bottom-to-top, then the rebar layers in input order)
        if (argc < 3) {
            opserr << "ASDShellSection::setResponse() - need to specify more data\n";
            return 0;
        }
        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= nFib) {
            std::size_t i = static_cast<std::size_t>(pointNum - 1);
            output.tag("FiberOutput");
            output.attr("number", pointNum);
            output.attr("zLoc", m_fib_z[i] + m_offset);
            output.attr("thickness", m_fib_w[i]);
            theResponse = m_fib_mat[i]->setResponse(&argv[2], argc - 2, output);
            output.endTag();
        }
    }
    else if (strcmp(argv[0], "rebar") == 0 || strcmp(argv[0], "Rebar") == 0) {
        // rebar $i ... : 1-based over the rebar layers only
        if (argc < 3) {
            opserr << "ASDShellSection::setResponse() - need to specify more data\n";
            return 0;
        }
        int nReb = static_cast<int>(m_reb_teq.size());
        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= nReb) {
            std::size_t i = static_cast<std::size_t>(nFib - nReb + pointNum - 1);
            output.tag("FiberOutput");
            output.attr("number", pointNum);
            output.attr("zLoc", m_fib_z[i] + m_offset);
            output.attr("thickness", m_fib_w[i]);
            theResponse = m_fib_mat[i]->setResponse(&argv[2], argc - 2, output);
            output.endTag();
        }
    }

    if (theResponse == 0)
        return SectionForceDeformation::setResponse(argv, argc, output);
    return theResponse;
}

int ASDShellSection::getResponse(int responseID, Information& secInfo)
{
    return SectionForceDeformation::getResponse(responseID, secInfo);
}

int ASDShellSection::setParameter(const char** argv, int argc, Parameter& param)
{
    // the same pass-through scheme as LayeredShellFiberSection: an explicit
    // "fiber id ..." targets one fiber, "fiber ..." broadcasts, anything
    // else is forwarded to every fiber
    int nFib = static_cast<int>(m_fib_mat.size());
    if (argc > 1) {
        if (strcmp(argv[0], "fiber") == 0 || strcmp(argv[0], "Fiber") == 0) {
            if (argc > 2) {
                int pointNum = atoi(argv[1]);
                if (pointNum > 0 && pointNum <= nFib)
                    return m_fib_mat[static_cast<std::size_t>(pointNum - 1)]->setParameter(&argv[2], argc - 2, param);
            }
            int mixed_result = -1;
            for (int i = 0; i < nFib; ++i) {
                if (m_fib_mat[static_cast<std::size_t>(i)]->setParameter(&argv[1], argc - 1, param) == 0)
                    mixed_result = 0;
            }
            return mixed_result;
        }
    }
    if (argc > 0) {
        int mixed_result = -1;
        for (int i = 0; i < nFib; ++i) {
            if (m_fib_mat[static_cast<std::size_t>(i)]->setParameter(argv, argc, param) == 0)
                mixed_result = 0;
        }
        return mixed_result;
    }
    return SectionForceDeformation::setParameter(argv, argc, param);
}

int ASDShellSection::updateParameter(int parameterID, Information& info)
{
    return SectionForceDeformation::updateParameter(parameterID, info);
}
