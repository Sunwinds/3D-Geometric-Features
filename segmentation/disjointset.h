#include <vector>

/// <summary>
/// A Union-Find/Disjoint-Set data structure.
/// </summary>
class DisjointSet {
public:
	/// <summary>
	/// The number of elements in the universe.
	/// </summary>
	int Count;

	/// <summary>
	/// The parent of each element in the universe.
	/// </summary>
	std::vector<int> Parent;

	/// <summary>
	/// The rank of each element in the universe.
	/// </summary>
	std::vector<int> Rank;

	/// <summary>
	/// The size of each set.
	/// </summary>
	std::vector<int> SizeOfSet;

	/// <summary>
	/// The number of disjoint sets.
	/// </summary>
	int SetCount;

	/// <summary>
	/// Initializes a new Disjoint-Set data structure, with the specified amount of elements in the universe.
	/// </summary>
	/// <param name='count'>
	/// The number of elements in the universe.
	/// </param>
	DisjointSet(int count) {
		Count = count;
		SetCount = count;
		Parent.resize(Count, 0);
		Rank.resize(Count, 0);
		SizeOfSet.resize(Count, 0);

		for (int i = 0; i < Count; i++) {
			Parent[i] = i;
			Rank[i] = 0;
			SizeOfSet[i] = 1;
		}
	}

	/// <summary>
	/// Find the parent of the specified element.
	/// </summary>
	/// <param name='i'>
	/// The specified element.
	/// </param>
	/// <remarks>
	/// All elements with the same parent are in the same set.
	/// </remarks>
	int Find(int i) {
		if (Parent[i] == i) {
			return i;
		} else {
			// Recursively find the real parent of i, and then cache it for later lookups.
			Parent[i] = Find(Parent[i]);
			return Parent[i];
		}
	}

	/// <summary>
	/// Unite the sets that the specified elements belong to.
	/// </summary>
	/// <param name='i'>
	/// The first element.
	/// </param>
	/// <param name='j'>
	/// The second element.
	/// </param>
	void Union(int i, int j) {

		// Find the representatives (or the root nodes) for the set that includes i
		int irep = Find(i),
			// And do the same for the set that includes j
			jrep = Find(j),
			// Get the rank of i's tree
			irank = Rank[irep],
			// Get the rank of j's tree
			jrank = Rank[jrep];

		// Elements are in the same set, no need to unite anything.
		if (irep == jrep)
			return;

		SetCount--;

		// If i's rank is less than j's rank
		if (irank < jrank) {

			// Then move i under j
			Parent[irep] = jrep;
			SizeOfSet[jrep] += SizeOfSet[irep];

		} // Else if j's rank is less than i's rank
		else if (jrank < irank) {

			// Then move j under i
			Parent[jrep] = irep;
			SizeOfSet[irep] += SizeOfSet[jrep];

		} // Else if their ranks are the same
		else {

			// Then move i under j (doesn't matter which one goes where)
			Parent[irep] = jrep;
			SizeOfSet[jrep] += SizeOfSet[irep];

			// And increment the the result tree's rank by 1
			Rank[jrep]++;
		}
	}

	/// <summary>
	/// Return the element count of the set that the specified elements belong to.
	/// </summary>
	/// <param name='i'>
	/// The element.
	/// </param>
	int SetSize(int i) {
		return SizeOfSet[Find(i)];
	}
};
