/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: process-wide registry of the contiguous equation-ownership
// blocks produced by ParallelNumberer's Plain (by-subdomain) ordering in
// OpenSeesMP. Filled once per numbering by ParallelNumberer::numberDOF;
// consumed by distributed solvers that need row-block ownership of the
// equations (e.g. ClusterPardisoSOE) and, later, by the distributed-RHS /
// PARPACK layers. Invalid (empty) when a GraphNumberer (e.g. ParallelRCM)
// was used: that ordering is not block-contiguous.
//
// Block layout in equation order: [b_0 .. b_{N-1}, b_N, numEqs] where block
// j < N belongs to worker rank j+1 (channel order) and the LAST block
// belongs to rank 0 (P0's leftover vertices are numbered last).

#ifndef EquationPartition_h
#define EquationPartition_h

#include <vector>

class EquationPartition {
 public:
  static EquationPartition &instance();

  void clear();
  // boundaries has size numRanks+1 as described above; myRank is this
  // process' MP rank, numRanks the total number of ranks.
  void set(int myRank, int numRanks, const std::vector<int> &boundaries);

  bool isValid() const { return validFlag; }
  int getNumRanks() const { return nRanks; }
  int getNumEquations() const { return numEqs; }

  // owned half-open range [begin, end) of a rank
  int getRowBegin(int rank) const;
  int getRowEnd(int rank) const;
  int getMyRowBegin() const { return getRowBegin(myRankId); }
  int getMyRowEnd() const { return getRowEnd(myRankId); }

  // owning rank of a global equation number, O(log nRanks)
  int getOwner(int eq) const;

 private:
  EquationPartition() : validFlag(false), myRankId(0), nRanks(0), numEqs(0) {}

  bool validFlag;
  int myRankId;
  int nRanks;
  int numEqs;
  std::vector<int> blockBoundaries;  // size nRanks+1, ascending eq starts
  std::vector<int> blockRank;        // size nRanks, rank owning block j
};

#endif
