/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: see EquationPartition.h

#include <EquationPartition.h>
#include <OPS_Globals.h>
#include <stdlib.h>
#include <algorithm>

EquationPartition &
EquationPartition::instance()
{
  static EquationPartition theInstance;
  return theInstance;
}

void
EquationPartition::clear()
{
  validFlag = false;
  myRankId = 0;
  nRanks = 0;
  numEqs = 0;
  blockBoundaries.clear();
  blockRank.clear();
}

void
EquationPartition::set(int myRank, int numRanks, const std::vector<int> &boundaries)
{
  clear();

  if (numRanks < 1 || (int)boundaries.size() != numRanks + 1)
    return;
  for (int j = 0; j < numRanks; j++)
    if (boundaries[j] < 0 || boundaries[j] > boundaries[j + 1])
      return;

  myRankId = myRank;
  nRanks = numRanks;
  blockBoundaries = boundaries;
  numEqs = boundaries[numRanks];

  // block j < nRanks-1 -> worker rank j+1 (channel order); last block -> rank 0
  blockRank.resize(nRanks);
  for (int j = 0; j < nRanks - 1; j++)
    blockRank[j] = j + 1;
  blockRank[nRanks - 1] = 0;

  validFlag = true;

  if (getenv("OPS_PARTITION_DEBUG") != 0) {
    opserr << "EquationPartition rank " << myRankId << "/" << nRanks
	   << ": owns [" << getMyRowBegin() << ", " << getMyRowEnd()
	   << ") of " << numEqs << " eqs; blocks:";
    for (int j = 0; j < nRanks; j++)
      opserr << " r" << blockRank[j] << "=[" << blockBoundaries[j]
	     << "," << blockBoundaries[j + 1] << ")";
    opserr << "\n";
  }
}

int
EquationPartition::getRowBegin(int rank) const
{
  if (!validFlag)
    return -1;
  for (int j = 0; j < nRanks; j++)
    if (blockRank[j] == rank)
      return blockBoundaries[j];
  return -1;
}

int
EquationPartition::getRowEnd(int rank) const
{
  if (!validFlag)
    return -1;
  for (int j = 0; j < nRanks; j++)
    if (blockRank[j] == rank)
      return blockBoundaries[j + 1];
  return -1;
}

int
EquationPartition::getOwner(int eq) const
{
  if (!validFlag || eq < 0 || eq >= numEqs)
    return -1;
  // first boundary strictly greater than eq -> block index
  std::vector<int>::const_iterator it =
    std::upper_bound(blockBoundaries.begin() + 1, blockBoundaries.end(), eq);
  int j = (int)(it - blockBoundaries.begin()) - 1;
  return blockRank[j];
}
