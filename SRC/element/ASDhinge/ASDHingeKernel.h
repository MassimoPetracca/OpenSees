/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2026/08/23 $

// Original implementation: Massimo Petracca (ASDEA)
//
// ASDHingeKernel: the constitutive core of ASDHinge.
//
// The element hands the kernel six local deformations and asks for six local
// forces and a 6x6 local tangent.  That is the ONLY contract, and it is the
// reason the class exists: today the tangent is diagonal (six independent
// uniaxial laws, exactly what a zeroLength does), tomorrow it can be full
// (N-M interaction, biaxial shear) without the element changing one line.
//
// The type id travels in sendSelf.  It costs one int now; adding it later
// would break every database and every parallel model already written.
//
// Slot layout is FIXED: 0..5 are the local Ux, Uy, Uz, Rx, Ry, Rz of the
// hinge, always all six, whatever the user assigned.  This is deliberate:
// with zeroLength the -dir/-mat lists have variable length, so "material 3"
// means Vz on one hinge and Mz on another, and a recorder column has no
// stable meaning across the model.

#ifndef ASDHingeKernel_h
#define ASDHingeKernel_h

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <OPS_Globals.h>

class ASDHingeKernel
{
public:
    enum KernelType {
        Uncoupled = 0
    };

    virtual ~ASDHingeKernel() {}

    virtual int getTypeId() const = 0;
    virtual ASDHingeKernel* getCopy() const = 0;

    // deformation -> force.  Returns 0, or the MOST NEGATIVE error code
    // returned by any slot.  See ASDHingeUncoupledKernel for why the codes
    // are not summed.
    virtual int setTrialDeformation(const Vector& e, const Vector& edot) = 0;

    // NOT const, and NOT cached at setTrialDeformation time: see the note in
    // ASDHingeUncoupledKernel::getForce.
    virtual const Vector& getForce() = 0;
    virtual const Matrix& getTangent() = 0;
    virtual const Matrix& getInitialTangent() = 0;
    virtual const Matrix& getDampTangent() = 0;

    virtual int commitState() = 0;
    virtual int revertToLastCommit() = 0;
    virtual int revertToStart() = 0;

    // slot introspection, for the element responses and for Print
    virtual UniaxialMaterial* getMaterial(int slot) const = 0;
    virtual double getSlotStiffness(int slot) const = 0;
    virtual int getSlotState(int slot) const = 0;

    virtual int sendSelf(int commitTag, Channel& theChannel) = 0;
    virtual int recvSelf(int commitTag, Channel& theChannel,
                         FEM_ObjectBroker& theBroker) = 0;

    virtual void Print(OPS_Stream& s) const = 0;

    // used by recvSelf: build an empty kernel of the received type
    static ASDHingeKernel* create(int typeId);
};


/**
The default kernel: six independent slots, each either

    FREE      no stiffness, no force (a released dof)
    LINEAR    a constant stiffness (this is the penalty tie of a beam end
              release: no uniaxialMaterial Elastic object needed for it)
    MATERIAL  an ASDHysteretic1DMaterial owned by the kernel

The element parser accepts nothing else in a MATERIAL slot: the IMPL-EX error
codes and the rich responses (damage, freeEnergy, limitStateRatio) are the
reason this element exists, and they are that material's.
*/
class ASDHingeUncoupledKernel : public ASDHingeKernel
{
public:
    enum SlotState {
        Free = 0,
        Linear = 1,
        Material = 2
    };

public:
    ASDHingeUncoupledKernel()
        : m_e(6), m_f(6), m_K(6, 6), m_K0(6, 6), m_C(6, 6)
    {
        m_e.Zero();
        for (int i = 0; i < 6; ++i) {
            m_state[i] = Free;
            m_k[i] = 0.0;
            m_mat[i] = 0;
        }
    }

    ~ASDHingeUncoupledKernel()
    {
        for (int i = 0; i < 6; ++i)
            if (m_mat[i]) delete m_mat[i];
    }

    // -- construction (used by the parser) -------------------------------

    void setFree(int slot)
    {
        clearSlot(slot);
        m_state[slot] = Free;
        m_k[slot] = 0.0;
    }

    void setLinear(int slot, double k)
    {
        clearSlot(slot);
        m_state[slot] = Linear;
        m_k[slot] = k;
    }

    // takes a COPY of the material
    void setMaterial(int slot, UniaxialMaterial* mat)
    {
        clearSlot(slot);
        m_state[slot] = Material;
        m_k[slot] = 0.0;
        m_mat[slot] = mat ? mat->getCopy() : 0;
        if (m_mat[slot] == 0)
            m_state[slot] = Free;
    }

    // must be called once the six slots are set
    void initialize()
    {
        m_K0.Zero();
        for (int i = 0; i < 6; ++i) {
            if (m_state[i] == Material)
                m_K0(i, i) = m_mat[i]->getInitialTangent();
            else if (m_state[i] == Linear)
                m_K0(i, i) = m_k[i];
        }
        m_K = m_K0;
        m_e.Zero();
        m_f.Zero();
        m_C.Zero();
    }

    // -- ASDHingeKernel ---------------------------------------------------

    int getTypeId() const { return ASDHingeKernel::Uncoupled; }

    ASDHingeKernel* getCopy() const
    {
        ASDHingeUncoupledKernel* other = new ASDHingeUncoupledKernel();
        for (int i = 0; i < 6; ++i) {
            other->m_state[i] = m_state[i];
            other->m_k[i] = m_k[i];
            other->m_mat[i] = m_mat[i] ? m_mat[i]->getCopy() : 0;
        }
        other->initialize();
        return other;
    }

    int setTrialDeformation(const Vector& e, const Vector& edot)
    {
        // Every slot is updated even after one of them fails: the materials
        // must all see the same trial state, otherwise a later revert leaves
        // them on different steps.  The codes are NOT summed the way
        // ZeroLength::update sums them -- the sum of two error codes is not an
        // error code, and EC_IMPLEX_Error_Control (-10) must reach the
        // analysis intact, because that is what makes it cut the step.
        int worst = 0;
        m_e = e;
        for (int i = 0; i < 6; ++i) {
            if (m_state[i] == Material) {
                int res = m_mat[i]->setTrialStrain(e(i), edot(i));
                if (res < worst) worst = res;
            }
        }
        return worst;
    }

    /**
    Read live from the materials, never cached at setTrialDeformation time.

    This is not a style choice, it is a defect that the comparison against
    zeroLength caught.  Under IMPL-EX, commitState RE-SOLVES IMPLICITLY and
    overwrites the material stress, so a value cached during the trial is the
    EXTRAPOLATED one and stops being what the material holds the moment the
    step is committed: measured 11.33 against zeroLength's 2.00 on the same
    hinge, with the deformation and the material's own getStress agreeing to
    the last bit.  A recorder would have shown Hinge.Force disagreeing with
    material $i stress on the same element -- exactly the confusion this
    element exists to remove.  A cache would also go stale on
    revertToLastCommit.  zeroLength reads getStress() live, and so do we.
    */
    const Vector& getForce()
    {
        for (int i = 0; i < 6; ++i) {
            if (m_state[i] == Material)
                m_f(i) = m_mat[i]->getStress();
            else if (m_state[i] == Linear)
                m_f(i) = m_k[i] * m_e(i);
            else
                m_f(i) = 0.0;
        }
        return m_f;
    }

    const Matrix& getTangent()
    {
        m_K.Zero();
        for (int i = 0; i < 6; ++i) {
            if (m_state[i] == Material)
                m_K(i, i) = m_mat[i]->getTangent();
            else if (m_state[i] == Linear)
                m_K(i, i) = m_k[i];
        }
        return m_K;
    }

    const Matrix& getInitialTangent() { return m_K0; }

    const Matrix& getDampTangent()
    {
        m_C.Zero();
        for (int i = 0; i < 6; ++i)
            if (m_state[i] == Material)
                m_C(i, i) = m_mat[i]->getDampTangent();
        return m_C;
    }

    int commitState()
    {
        int worst = 0;
        for (int i = 0; i < 6; ++i)
            if (m_state[i] == Material) {
                int res = m_mat[i]->commitState();
                if (res < worst) worst = res;
            }
        return worst;
    }

    int revertToLastCommit()
    {
        int worst = 0;
        for (int i = 0; i < 6; ++i)
            if (m_state[i] == Material) {
                int res = m_mat[i]->revertToLastCommit();
                if (res < worst) worst = res;
            }
        return worst;
    }

    int revertToStart()
    {
        int worst = 0;
        for (int i = 0; i < 6; ++i)
            if (m_state[i] == Material) {
                int res = m_mat[i]->revertToStart();
                if (res < worst) worst = res;
            }
        m_e.Zero();
        m_f.Zero();
        m_K = m_K0;
        return worst;
    }

    UniaxialMaterial* getMaterial(int slot) const
    {
        if (slot < 0 || slot > 5) return 0;
        return m_mat[slot];
    }

    double getSlotStiffness(int slot) const
    {
        if (slot < 0 || slot > 5) return 0.0;
        return m_k[slot];
    }

    int getSlotState(int slot) const
    {
        if (slot < 0 || slot > 5) return Free;
        return m_state[slot];
    }

    int sendSelf(int commitTag, Channel& theChannel)
    {
        // one ID and one Vector, then the materials: keeping the number and
        // the order of the messages fixed is what makes recvSelf possible
        static ID idata(19);
        idata(0) = getTypeId();
        int nmat = 0;
        for (int i = 0; i < 6; ++i) {
            idata(1 + i) = m_state[i];
            if (m_state[i] == Material) {
                int dbTag = m_mat[i]->getDbTag();
                if (dbTag == 0) {
                    dbTag = theChannel.getDbTag();
                    if (dbTag != 0)
                        m_mat[i]->setDbTag(dbTag);
                }
                idata(7 + i) = m_mat[i]->getClassTag();
                idata(13 + i) = dbTag;
                ++nmat;
            }
            else {
                idata(7 + i) = 0;
                idata(13 + i) = 0;
            }
        }
        if (theChannel.sendID(0, commitTag, idata) < 0) {
            opserr << "ASDHingeUncoupledKernel::sendSelf - failed to send ID\n";
            return -1;
        }
        static Vector vdata(6);
        for (int i = 0; i < 6; ++i)
            vdata(i) = m_k[i];
        if (theChannel.sendVector(0, commitTag, vdata) < 0) {
            opserr << "ASDHingeUncoupledKernel::sendSelf - failed to send Vector\n";
            return -1;
        }
        for (int i = 0; i < 6; ++i) {
            if (m_state[i] == Material) {
                if (m_mat[i]->sendSelf(commitTag, theChannel) < 0) {
                    opserr << "ASDHingeUncoupledKernel::sendSelf - material "
                           << i + 1 << " failed\n";
                    return -1;
                }
            }
        }
        return 0;
    }

    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
    {
        static ID idata(19);
        if (theChannel.recvID(0, commitTag, idata) < 0) {
            opserr << "ASDHingeUncoupledKernel::recvSelf - failed to recv ID\n";
            return -1;
        }
        static Vector vdata(6);
        if (theChannel.recvVector(0, commitTag, vdata) < 0) {
            opserr << "ASDHingeUncoupledKernel::recvSelf - failed to recv Vector\n";
            return -1;
        }
        for (int i = 0; i < 6; ++i) {
            clearSlot(i);
            m_state[i] = idata(1 + i);
            m_k[i] = vdata(i);
        }
        for (int i = 0; i < 6; ++i) {
            if (m_state[i] == Material) {
                int classTag = idata(7 + i);
                int dbTag = idata(13 + i);
                m_mat[i] = theBroker.getNewUniaxialMaterial(classTag);
                if (m_mat[i] == 0) {
                    opserr << "ASDHingeUncoupledKernel::recvSelf - no material "
                              "of class tag " << classTag << "\n";
                    return -1;
                }
                m_mat[i]->setDbTag(dbTag);
                if (m_mat[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
                    opserr << "ASDHingeUncoupledKernel::recvSelf - material "
                           << i + 1 << " failed\n";
                    return -1;
                }
            }
        }
        initialize();
        return 0;
    }

    void Print(OPS_Stream& s) const
    {
        static const char* names[6] = { "Ux", "Uy", "Uz", "Rx", "Ry", "Rz" };
        for (int i = 0; i < 6; ++i) {
            s << "    " << names[i] << ": ";
            if (m_state[i] == Material)
                s << "material " << m_mat[i]->getTag() << " ("
                  << m_mat[i]->getClassType() << ")\n";
            else if (m_state[i] == Linear)
                s << "linear, K = " << m_k[i] << "\n";
            else
                s << "free\n";
        }
    }

private:
    void clearSlot(int slot)
    {
        if (m_mat[slot]) {
            delete m_mat[slot];
            m_mat[slot] = 0;
        }
    }

private:
    int m_state[6];
    double m_k[6];
    UniaxialMaterial* m_mat[6];
    Vector m_e;   // the trial deformation, needed by the Linear slots
    Vector m_f;
    Matrix m_K;
    Matrix m_K0;
    Matrix m_C;
};


inline ASDHingeKernel* ASDHingeKernel::create(int typeId)
{
    switch (typeId) {
    case ASDHingeKernel::Uncoupled:
        return new ASDHingeUncoupledKernel();
    default:
        opserr << "ASDHingeKernel::create - unknown kernel type " << typeId << "\n";
        return 0;
    }
}

#endif // ASDHingeKernel_h
