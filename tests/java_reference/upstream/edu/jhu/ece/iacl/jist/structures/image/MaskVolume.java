package edu.jhu.ece.iacl.jist.structures.image;
// TODO: Auto-generated Javadoc

/**
 * Mask Interface.
 * 
 * @author Blake Lucas
 */
public interface MaskVolume {
	
	/**
	 * Gets the neighbors x.
	 * 
	 * @return the neighbors x
	 */
	public byte[] getNeighborsX();
	
	/**
	 * Gets the neighbors y.
	 * 
	 * @return the neighbors y
	 */
	public byte[] getNeighborsY();
	
	/**
	 * Gets the neighbors z.
	 * 
	 * @return the neighbors z
	 */
	public byte[] getNeighborsZ();	
}
