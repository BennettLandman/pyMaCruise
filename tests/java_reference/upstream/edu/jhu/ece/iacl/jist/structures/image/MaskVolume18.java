package edu.jhu.ece.iacl.jist.structures.image;
// TODO: Auto-generated Javadoc

/**
 * 18-way neighbors for volume
 * Note: the C implementations sometimes use a different ordering of these neighbors.
 * 
 * @author Blake Lucas
 */
public class MaskVolume18 implements MaskVolume{
	
	/** The xoff18. */
	private static byte xoff18[] = {1,0,0,-1,0,0, 1,-1, 1,0, 0, 0,1,1,-1,-1,0,-1};
	
	/** The yoff18. */
	private static byte yoff18[] = {0,1,0,0,-1,0, 1,1, -1,1, 1,-1,0,0, 0,-1,-1,0};
	
	/** The zoff18. */
	private static byte zoff18[] = {0,0,1,0, 0,-1,0,0, 0, 1,-1, 1,1,-1,1,0,-1,-1};
	 
 	/** The Constant length. */
 	public static final int length=xoff18.length;
	
	/* (non-Javadoc)
	 * @see edu.jhu.ece.iacl.jist.structures.image.MaskVolume#getNeighborsX()
	 */
	public byte[] getNeighborsX() {
		return xoff18;
	}
	
	/* (non-Javadoc)
	 * @see edu.jhu.ece.iacl.jist.structures.image.MaskVolume#getNeighborsY()
	 */
	public byte[] getNeighborsY() {
		return yoff18;
	}
	
	/* (non-Javadoc)
	 * @see edu.jhu.ece.iacl.jist.structures.image.MaskVolume#getNeighborsZ()
	 */
	public byte[] getNeighborsZ() {
		return zoff18;
	};
}
