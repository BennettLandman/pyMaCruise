package edu.jhu.ece.iacl.jist.structures.data;

// TODO: Auto-generated Javadoc
/**
 * Indexable object used for heap.
 * 
 * @param <E> Must be comparable
 * 
 * @author Blake Lucas
 */
public interface Indexable<E extends Comparable<? super E>> extends Comparable<E>{
	
	/**
	 * Set position in volume.
	 * 
	 * @param i row
	 * @param j column
	 * @param k slice
	 */
	public void setRefPosition(int i,int j,int k);
	
	/**
	 * Gets the row.
	 * 
	 * @return the row
	 */
	public int getRow();
	
	/**
	 * Gets the column.
	 * 
	 * @return the column
	 */
	public int getColumn();
	
	/**
	 * Gets the slice.
	 * 
	 * @return the slice
	 */
	public int getSlice();
	
	/**
	 * Gets the chain index.
	 * 
	 * @return the chain index
	 */
	public int getChainIndex();
	
	/**
	 * Sets the value.
	 * 
	 * @param obj the new value
	 */
	public void setValue(Comparable obj);
	
	/**
	 * Gets the value.
	 * 
	 * @return the value
	 */
	public Comparable getValue();
	
	/**
	 * Sets the chain index.
	 * 
	 * @param chainIndex the new chain index
	 */
	public void setChainIndex(int chainIndex);
	
	/**
	 * Set index into heap.
	 * 
	 * @param index the index
	 */
	public void setIndex(int index);
	
	/**
	 * Get heap index.
	 * 
	 * @return the index
	 */
	public int getIndex();
}
