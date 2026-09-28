package edu.jhu.ece.iacl.jist.structures.image;
// TODO: Auto-generated Javadoc

/**
 * 26-way neighbors for volume
 * Note: the C implementations sometimes use a different ordering of these neighbors.
 * 
 * @author Blake Lucas
 */
public class MaskVolume26 implements MaskVolume{
	
	/** The xoff26. */
	static byte xoff26[] = {0,1,0, 0,-1, 0, 1, 1,-1,-1,1, 1,-1,-1,0, 0, 0, 0,-1,-1, 1, 1,-1,-1, 1, 1};
	
	/** The yoff26. */
	static byte yoff26[] = {1,0,0,-1, 0, 0, 1,-1, 1,-1,0, 0, 0, 0,1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1};
	
	/** The zoff26. */
	static byte zoff26[] = {0,0,1, 0, 0,-1, 0, 0, 0, 0,1,-1, 1,-1,1,-1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1};
	 
 	/** The Constant length. */
 	public static final int length=xoff26.length;
	
	/* (non-Javadoc)
	 * @see edu.jhu.ece.iacl.jist.structures.image.MaskVolume#getNeighborsX()
	 */
	public byte[] getNeighborsX() {
		return xoff26;
	}
	
	/* (non-Javadoc)
	 * @see edu.jhu.ece.iacl.jist.structures.image.MaskVolume#getNeighborsY()
	 */
	public byte[] getNeighborsY() {
		return yoff26;
	}
	
	/* (non-Javadoc)
	 * @see edu.jhu.ece.iacl.jist.structures.image.MaskVolume#getNeighborsZ()
	 */
	public byte[] getNeighborsZ() {
		return zoff26;
	}

}
