package edu.jhu.ece.iacl.algorithms.ace;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.MaskVolume;

/**
 * Generic morphology routines. These are used by ACE, but Pilou's
 * implementation are much faster.
 * 
 * @author Blake Lucas
 * 
 */
public class Morphology {
	/**
	 * Dilate volume with mask
	 * 
	 * @param in
	 *            volume
	 * @param mask
	 *            mask
	 * @param iterations
	 *            iterations
	 * @return dilated volume
	 */
	public static ImageData dilate(ImageData in, MaskVolume mask, int iterations) {
		int i, j, k, index, ci, cj, ck, count;
		int XN = in.getRows();
		int YN = in.getCols();
		int ZN = in.getSlices();
		// ImageData tmpvol = in.clone();
		ImageDataInt tmpvol = new ImageDataInt(in);
		// ImageData out = in.clone();
		ImageDataInt out = new ImageDataInt(in);
		byte[] neighborsX = mask.getNeighborsX();
		byte[] neighborsY = mask.getNeighborsY();
		byte[] neighborsZ = mask.getNeighborsZ();
		int maskLength = neighborsX.length;
		for (count = 1; count <= iterations; count++) {
			for (i = 0; i < XN; i++) {
				for (j = 0; j < YN; j++) {
					for (k = 0; k < ZN; k++) {
						// if (tmpvol.get(i, j, k).byteValue() == 1) {
						if (tmpvol.getBoolean(i, j, k) == true) {
							// out.set(i, j, k, (byte) 1);
							out.set(i, j, k, (int) 1);
							continue;
						}
						out.set(i, j, k, 0);
						for (index = 0; index < maskLength; index++) {
							ci = i + neighborsX[index];
							cj = j + neighborsY[index];
							ck = k + neighborsZ[index];
							if (ci < 0 || ci >= XN || cj < 0 || cj >= YN
									|| ck < 0 || ck >= ZN)
								continue;
							// if (tmpvol.getUByte(ci, cj, ck) == 1) {
							if (tmpvol.getInt(ci, cj, ck) == 1) {
								out.set(i, j, k, (int) 1);
								break;
							}
						}
					}
				}
			}

			if (count < iterations) {
				for (i = 0; i < XN; i++) {
					for (j = 0; j < YN; j++) {
						for (k = 0; k < ZN; k++) {
							// tmpvol.set(i, j, k, out.getUByte(i, j, k));
							tmpvol.set(i, j, k, out.getInt(i, j, k));
						}
					}
				}
			}
		}
		return out;
	}

	/**
	 * Erode volume with mask
	 * 
	 * @param in
	 *            volume
	 * @param mask
	 *            mask
	 * @param iterations
	 *            iterations
	 * @return eroded volume
	 */
	public static ImageData erode(ImageData in, MaskVolume mask, int iterations) {
		int i, j, k, index, ci, cj, ck, count;
		int XN = in.getRows();
		int YN = in.getCols();
		int ZN = in.getSlices();
		// ImageData tmpvol = in.clone();
		// ImageData out = in.clone();
		ImageDataInt tmpvol = new ImageDataInt(in);
		ImageDataInt out = new ImageDataInt(in);
		byte[] neighborsX = mask.getNeighborsX();
		byte[] neighborsY = mask.getNeighborsY();
		byte[] neighborsZ = mask.getNeighborsZ();
		int maskLength = neighborsX.length;
		for (count = 1; count <= iterations; count++) {
			for (i = 0; i < XN; i++) {
				for (j = 0; j < YN; j++) {
					for (k = 0; k < ZN; k++) {
						// if (tmpvol.get(i, j, k).byteValue() == 0) {
						if (tmpvol.getInt(i, j, k) == 0) {
							// out.set(i, j, k, (byte) 0);
							out.set(i, j, k, (int) 0);
							continue;
						}
						out.set(i, j, k, 1);
						for (index = 0; index < maskLength; index++) {
							ci = i + neighborsX[index];
							cj = j + neighborsY[index];
							ck = k + neighborsZ[index];
							if (ci < 0 || ci >= XN || cj < 0 || cj >= YN
									|| ck < 0 || ck >= ZN)
								continue;
							// if (tmpvol.getUByte(ci, cj, ck) == 0) {
							if (tmpvol.getInt(ci, cj, ck) == 0) {
								out.set(i, j, k, (int) 0);
								break;
							}
						}
					}
				}
			}

			if (count < iterations) {
				for (i = 0; i < XN; i++) {
					for (j = 0; j < YN; j++) {
						for (k = 0; k < ZN; k++) {
							// tmpvol.set(i, j, k, out.getUByte(i, j, k));
							tmpvol.set(i, j, k, out.getInt(i, j, k));
						}
					}
				}
			}
		}
		return out;
	}

	/**
	 * Performs closing operation Erode(Dilate(A,Mask),Mask))
	 * 
	 * @param in
	 * @param mask
	 * @param iterations
	 * @return
	 */
	public static ImageData close(ImageData in, MaskVolume mask, int iterations) {
		ImageData out = Morphology.dilate(in, mask, iterations);
		out = Morphology.erode(out, mask, iterations);
		return out;
	}

	/**
	 * Performs closing operation Dilate(Erode(A,Mask),Mask))
	 * 
	 * @param in
	 * @param mask
	 * @param iterations
	 * @return
	 */
	public static ImageData open(ImageData in, MaskVolume mask, int iterations) {
		ImageData out = Morphology.erode(in, mask, iterations);
		out = Morphology.dilate(out, mask, iterations);
		return out;
	}
}
