package edu.jhu.ece.iacl.algorithms.ace;

import java.util.Vector;

import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.data.BinaryMinHeap;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.Voxel;
import edu.jhu.ece.iacl.jist.structures.image.VoxelFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelIndexed;

/**
 * Thins ACE skeleton because its generally too thick
 * 
 * @author Blake Lucas
 * 
 */
public class ThinACE extends AbstractCalculation {
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.1 $");
	}



	private static final int xoff6[] = { 1, 0, 0, -1, 0, 0 };
	private static final int yoff6[] = { 0, 1, 0, 0, -1, 0 };
	private static final int zoff6[] = { 0, 0, 1, 0, 0, -1 };

	private static final int xoff18[] = { 1, 0, 0, -1, 0, 0, 1, -1, 1, 0, 0, 0, 1, 1, -1, -1, 0, -1 };
	private static final int yoff18[] = { 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, 1, -1, 0, 0, 0, -1, -1, 0 };
	private static final int zoff18[] = { 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, -1, 1, 1, -1, 1, 0, -1, -1 };

	private static final int xoff26[] = { 0, 1, 0, 0, -1, 0, 1, 1, -1, -1, 1, 1, -1, -1, 0, 0, 0, 0, -1, -1, 1, 1, -1, -1, 1, 1 };
	private static final int yoff26[] = { 1, 0, 0, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
	private static final int zoff26[] = { 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1 };

	public ThinACE() {
		super();
		setLabel("Thin ACE");
	}

	public ThinACE(AbstractCalculation calc) {
		super(calc);
		setLabel("Thin ACE");
	}

	private int ConInCube18(int[][][] cubeMat, int label) {
		int i, j, k, index, index18;
		int number;
		Vector<VoxelIndexed<Voxel>> NeiQ = new Vector<VoxelIndexed<Voxel>>();
		int ci, cj, ck, ni, nj, nk;
		// int ioff, joff, koff;
		number = 0;
		for (index18 = 0; index18 < 18; index18++) {
			i = xoff18[index18];
			j = yoff18[index18];
			k = zoff18[index18];
			if (cubeMat[i + 1][j + 1][k + 1] != label)
				continue;
			VoxelIndexed<Voxel> cPt = new VoxelIndexed<Voxel>(null);
			cPt.setRefPosition(i, j, k);
			number++;
			if (number > 1) {
				return number;
			}
			cubeMat[i + 1][j + 1][k + 1] = 0;
			NeiQ.add(cPt);
			while (!NeiQ.isEmpty()) {
				cPt = NeiQ.remove(0);
				ci = cPt.getRow();
				cj = cPt.getColumn();
				ck = cPt.getSlice();
				for (index = 0; index < 18; index++) {
					ni = ci + xoff18[index];
					nj = cj + yoff18[index];
					nk = ck + zoff18[index];
					if (ni >= (-1) && ni <= 1 && nj >= (-1) && nj <= 1
							&& nk >= (-1) && nk <= 1) {
						if (cubeMat[ni + 1][nj + 1][nk + 1] == label) {
							VoxelIndexed<Voxel> nPt = new VoxelIndexed<Voxel>(
									null);
							nPt.setRefPosition(ni, nj, nk);
							cubeMat[ni + 1][nj + 1][nk + 1] = 0;
							NeiQ.add(nPt);
						}
					}
				}
			}
		}
		return number;
	}

	public ImageData solve(ImageData vol) {
		// int opt, errflg;
		ImageDataInt orivol = new ImageDataInt(vol);
		int[][][] orivolMat = orivol.toArray3d();
		int XN = orivol.getRows();
		int YN = orivol.getCols();
		int ZN = orivol.getSlices();

		int i, j, k, index, ci, cj, ck, flag, xoff, yoff, zoff;
		// int pointer;

		int fcon, bcon, obcount, count18, endpt;

		// int heapindex;

		// Part of Memory allocation
		ImageDataInt Label = new ImageDataInt(XN, YN, ZN);
		ImageDataFloat distvol;
		ImageDataFloat initphi = new ImageDataFloat(XN, YN, ZN);
		ImageDataInt cube = new ImageDataInt(3, 3, 3);
		int[][][] cubeMat = cube.toArray3d();
		int[][][] LabelMat = Label.toArray3d();
		float[][][] initphiMat = initphi.toArray3d();
		DistanceField distField=new DistanceField(this);
		// int count = 0;
		for (k = 0; k < ZN; k++) {
			for (j = 0; j < YN; j++) {
				for (i = 0; i < XN; i++) {
					// Clear boundary
					if (orivol.getInt(i, j, k) > 0) {
						orivolMat[i][j][k] = 1;
						initphiMat[i][j][k] = 1;
					} else {
						initphiMat[i][j][k] = -1;
					}
					LabelMat[i][j][k] = 0;
					if (k == 0 || k == (int) (ZN - 1) || j == 0
							|| j == (int) (YN - 1) || i == 0
							|| i == (int) (XN - 1)) {
						orivolMat[i][j][k] = 0;
					}
				}
			}
		}
		distvol = distField.solve(initphi, 4.0);
		float[][][] distvolMat = distvol.toArray3d();
		Vector<VoxelIndexed<Voxel>> List1 = new Vector<VoxelIndexed<Voxel>>();
		Vector<VoxelIndexed<Voxel>> List2 = new Vector<VoxelIndexed<Voxel>>();
		Vector<VoxelIndexed<Voxel>> tmpList;
		BinaryMinHeap heap = new BinaryMinHeap(XN * YN * ZN, XN, YN, ZN);

		/* Fjnd all border pojnts */
		for (k = 1; k < (ZN - 1); k++) {
			for (j = 1; j < (YN - 1); j++) {
				for (i = 1; i < (XN - 1); i++) {
//					if (orivol.getUByte(i, j, k) == 0)
					if (orivol.getInt(i, j, k) == 0)
						continue;
					/* Check wether jt js a border voxel */
					for (index = 0; index < 6; index++) {
						ck = k + zoff6[index];
						cj = j + yoff6[index];
						ci = i + xoff6[index];
//						if (orivol.getUByte(ci, cj, ck) == 0) {
						if (orivol.getInt(ci, cj, ck) == 0) {
							VoxelIndexed<Voxel> cPt = new VoxelIndexed<Voxel>(
									null);
							cPt.setRefPosition(i, j, k);
							List1.add(cPt);
							LabelMat[i][j][k] = 1;
							break;
						}
					}
				}
			}
		}
		flag = 1;

		while (flag != 0) {
			flag = 0;
			/* Find all the surface points in List1 */
			setTotalUnits(List1.size());
			while (!List1.isEmpty()) {
				VoxelIndexed<Voxel> cPt = List1.remove(0);
				i = cPt.getRow();
				j = cPt.getColumn();
				k = cPt.getSlice();
				// System.out.println("LIST ("+i+","+j+","+k+")");
				for (index = 0; index < 6; index++) {
					xoff = xoff26[index];
					yoff = yoff26[index];
					zoff = zoff26[index];
					ci = i + xoff;
					cj = j + yoff;
					ck = k + zoff;
					cubeMat[xoff + 1][yoff + 1][zoff + 1] = 0;
//					if (orivol.getUByte(ci, cj, ck) > 0) {
					if (orivol.getInt(ci, cj, ck) > 0) {
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = 1;
					} else {
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1; /*
																	 * background
																	 * points
																	 */
					}
				}
				for (index = 6; index < 18; index++) {
					xoff = xoff26[index];
					yoff = yoff26[index];
					zoff = zoff26[index];
					ci = i + xoff;
					cj = j + yoff;
					ck = k + zoff;
					// System.out.println("LIST ("+xoff+","+yoff+","+zoff+")
					// "+orivol.getDouble(ci, cj, ck)+" "+cube.getDouble(xoff+1,
					// 0, 0)+" "+cube.getByte(0, yoff + 1, 0)+"
					// "+cube.getDouble(0, 0, zoff+1));
					cubeMat[xoff + 1][yoff + 1][zoff + 1] = 0;
//					if (orivol.getUByte(ci, cj, ck) > 0) {
					if (orivol.getInt(ci, cj, ck) > 0) {
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = 1;/*
																	 * Object
																	 * points
																	 */
					} else { /* N_6^2 */
						if (zoff == 0) {
							if (cubeMat[1][yoff + 1][1] == -1
									|| cubeMat[xoff + 1][1][1] == -1)

								cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
						} else if (yoff == 0) {
							if (cubeMat[xoff + 1][1][1] == -1
									|| cubeMat[1][1][zoff + 1] == -1)

								cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
						} else { /* cj == j */
							if (cubeMat[1][1][zoff + 1] == -1
									|| cubeMat[1][yoff + 1][1] == -1)

								cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
						}
					}
				}
				for (index = 18; index < 26; index++) {
					xoff = xoff26[index];
					yoff = yoff26[index];
					zoff = zoff26[index];
					ci = i + xoff;
					cj = j + yoff;
					ck = k + zoff;
					cubeMat[xoff + 1][yoff + 1][zoff + 1] = 0;
//					if (orivol.getUByte(ci, cj, ck) > 0) {
					if (orivol.getInt(ci, cj, ck) > 0) {
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = 1;
					} else if (cubeMat[xoff + 1][yoff + 1][1] == -1
							|| cubeMat[xoff + 1][1][zoff + 1] == -1
							|| cubeMat[1][yoff + 1][zoff + 1] == -1) {
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
					}
				}
				fcon = ConInCube18(cubeMat, 1);
				bcon = ConInCube6(cubeMat, -1);
				if (fcon == 1 && bcon == 2) { /* surface point */
					LabelMat[i][j][k] = 4;
				} else if (fcon == 1 && bcon == 1) { /* Simple point */
					VoxelIndexed<VoxelFloat> v = new VoxelIndexed<VoxelFloat>(
							new VoxelFloat(distvolMat[i][j][k]));
					v.setRefPosition(i, j, k);
					heap.add(v);
					flag = 1; /* Simple point exists */
				} else { /* Put into List2 */
					List2.add(cPt);
					LabelMat[i][j][k] = 2; /* In List2; 1 will also do */
				}

			} /* End of While(List1) */
			// System.out.println("FIRST LIST "+count++);
			// orivol.printNonZeroString(100);
			ImageDataInt tmp = cube.mimic();
			int[][][] tmpMat = tmp.toArray3d();
			while (!heap.isEmpty()) {
				VoxelIndexed<Voxel> he = (VoxelIndexed<Voxel>) heap.remove();

				/* Get the top element, which has the smallest distance */
				i = he.getRow();
				j = he.getColumn();
				k = he.getSlice();
				Label.set(i, j, k, 0); /* No longer in list */
				// System.out.println("HEAP ("+i+","+j+","+k+") "+heap.size()+"
				// "+he.getFloat());
				/* Determine whether the point is simple */
				obcount = 0;
				/* Build N_6^3 */
				// System.out.println("CUBE 0");
				// cube.printNonZeroString(100);
				for (index = 0; index < 6; index++) {
					xoff = xoff26[index];
					yoff = yoff26[index];
					zoff = zoff26[index];
					ci = i + xoff;
					cj = j + yoff;
					ck = k + zoff;
					cubeMat[xoff + 1][yoff + 1][zoff + 1] = 0;
					tmpMat[xoff + 1][yoff + 1][zoff + 1] = 0;
					/*
					 * if(ci<0 || ci >= YN || cj<0 || cj >=XN || ck <0 || ck
					 * >=ZN) continue;
					 *//*
						 * Assume object doesn't touch image boundary
						 */
//					if (orivol.getUByte(ci, cj, ck) > 0) {
					if (orivol.getInt(ci, cj, ck) > 0) {
						tmpMat[xoff + 1][yoff + 1][zoff + 1] = 1;
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = 1; /*
																	 * Object
																	 * points
																	 */
						obcount++;
					} else {
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1; /*											 */
					}
				}
				// System.out.println("CUBE 1");
				// cube.printNonZeroString(100);
				for (index = 6; index < 18; index++) {
					xoff = xoff26[index];
					yoff = yoff26[index];
					zoff = zoff26[index];
					ci = i + xoff;
					cj = j + yoff;
					ck = k + zoff;
					cubeMat[xoff + 1][yoff + 1][zoff + 1] = 0;
					tmpMat[xoff + 1][yoff + 1][zoff + 1] = 0;
					/*
					 * if(ci<0 || ci >= YN || cj<0 || cj >=XN || ck <0 || ck
					 * >=ZN) continue;
					 */
					// System.out.println("("+ci+","+cj+","+ck+")"+orivol.getByte(ci,cj,ck));
//					if (orivol.getUByte(ci, cj, ck) > 0) {
					if (orivol.getInt(ci, cj, ck) > 0) {
						obcount++;
						tmpMat[xoff + 1][yoff + 1][zoff + 1] = 1;
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = 1; /*
																	 * Object
																	 * points
																	 */
					} else { /* N_6^2 */
						if (zoff == 0) {
							if (cubeMat[1][yoff + 1][1] == -1
									|| cubeMat[xoff + 1][1][1] == -1)

								cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
						} else if (yoff == 0) {
							if (cubeMat[xoff + 1][1][1] == -1
									|| cubeMat[1][1][zoff + 1] == -1)

								cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
						} else { /* cj == j */
							if (cubeMat[1][1][zoff + 1] == -1
									|| cubeMat[1][yoff + 1][1] == -1)

								cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
						}
					}
				}

				count18 = obcount;
				// System.out.println("CUBE 2");
				// cube.printNonZeroString(100);
				for (index = 18; index < 26; index++) {
					xoff = xoff26[index];
					yoff = yoff26[index];
					zoff = zoff26[index];
					ci = i + xoff;
					cj = j + yoff;
					ck = k + zoff;
					cubeMat[xoff + 1][yoff + 1][zoff + 1] = 0;
					tmpMat[xoff + 1][yoff + 1][zoff + 1] = 0;
					/*
					 * if(ci<0 || ci >= YN || cj<0 || cj >=XN || ck <0 || ck
					 * >=ZN) continue;
					 */
//					if (orivol.getUByte(ci, cj, ck) > 0) {
					if (orivol.getInt(ci, cj, ck) > 0) {
						obcount++;
						tmpMat[xoff + 1][yoff + 1][zoff + 1] = 1;
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = 1;
					} else if (cubeMat[xoff + 1][yoff + 1][1] == -1
							|| cubeMat[xoff + 1][1][zoff + 1] == -1
							|| cubeMat[1][yoff + 1][zoff + 1] == -1) {
						cubeMat[xoff + 1][yoff + 1][zoff + 1] = -1;
					}
				}

				fcon = ConInCube18(cubeMat, 1);
				bcon = ConInCube6(cubeMat, -1);

				if (!(fcon == 1 && bcon == 1)) {
					VoxelIndexed<Voxel> cPt = new VoxelIndexed<Voxel>(null);
					cPt.setRefPosition(i, j, k);
					List2.add(cPt);
					Label.set(i, j, k, 2);

					continue;
				}

				// tmp.printNonZeroString(100);
				if (count18 == 1)
					obcount = 1;

				endpt = 0;
				/* Determine whether the simple point is a surface-end */
				if (obcount == 2 || obcount == 3) { /* end points */
					endpt = 1;
				} else {
					if (obcount >= 4 && obcount <= 7) {
						endpt = 1;
						/*
						 * If more than 4 object points lie in one octant, the
						 * point is not an end point
						 */
						if ((tmp.getInt(0, 0, 0) + tmp.getInt(0, 0, 1)
								+ tmp.getInt(0, 1, 1) + tmp.getInt(0, 1, 0)
								+ tmp.getInt(1, 1, 0) + tmp.getInt(1, 0, 0) + tmp
								.getInt(1, 0, 1)) > 3)
							endpt = 0;
						if ((tmp.getInt(0, 0, 1) + tmp.getInt(0, 0, 2)
								+ tmp.getInt(0, 1, 2) + tmp.getInt(0, 1, 1)
								+ tmp.getInt(1, 1, 2) + tmp.getInt(1, 0, 1) + tmp
								.getInt(1, 0, 2)) > 3)
							endpt = 0;
						if ((tmp.getInt(0, 1, 0) + tmp.getInt(0, 1, 1)
								+ tmp.getInt(0, 2, 1) + tmp.getInt(0, 2, 0)
								+ tmp.getInt(1, 2, 0) + tmp.getInt(1, 1, 0) + tmp
								.getInt(1, 2, 1)) > 3)
							endpt = 0;
						if ((tmp.getInt(0, 1, 1) + tmp.getInt(0, 1, 2)
								+ tmp.getInt(0, 2, 2) + tmp.getInt(0, 2, 1)
								+ tmp.getInt(1, 2, 2) + tmp.getInt(1, 2, 1) + tmp
								.getInt(1, 1, 2)) > 3)
							endpt = 0;
						if ((tmp.getInt(1, 0, 0) + tmp.getInt(1, 0, 1)
								+ tmp.getInt(2, 1, 1) + tmp.getInt(1, 1, 0)
								+ tmp.getInt(2, 1, 0) + tmp.getInt(2, 0, 0) + tmp
								.getInt(2, 0, 1)) > 3)
							endpt = 0;
						if ((tmp.getInt(1, 0, 1) + tmp.getInt(1, 0, 2)
								+ tmp.getInt(1, 1, 2) + tmp.getInt(2, 1, 1)
								+ tmp.getInt(2, 1, 2) + tmp.getInt(2, 0, 1) + tmp
								.getInt(2, 0, 2)) > 3)
							endpt = 0;
						if ((tmp.getInt(1, 1, 0) + tmp.getInt(2, 1, 1)
								+ tmp.getInt(1, 2, 1) + tmp.getInt(1, 2, 0)
								+ tmp.getInt(2, 2, 0) + tmp.getInt(2, 1, 0) + tmp
								.getInt(2, 2, 1)) > 3)
							endpt = 0;
						if ((tmp.getInt(2, 1, 1) + tmp.getInt(1, 1, 2)
								+ tmp.getInt(1, 2, 2) + tmp.getInt(1, 2, 1)
								+ tmp.getInt(2, 2, 2) + tmp.getInt(2, 2, 1) + tmp
								.getInt(2, 1, 2)) > 3)
							endpt = 0;
					}
				}

				if (endpt == 1) { /* An end point */
					Label.set(i, j, k, 4);
				} else {
					/* Remove the current point */
					// System.out.println("REMOVE ("+i+","+j+","+k+")");
//					orivol.set(i, j, k, 0);
					orivol.set(i, j, k,(int) 0);

					/* Put its 6 object neighbor into List2 */
					for (index = 0; index < 6; index++) {
						ck = k + zoff26[index];
						cj = j + yoff26[index];
						ci = i + xoff26[index];
//						if (orivol.getUByte(ci, cj, ck) > 0
						if (orivol.getInt(ci, cj, ck) > 0
								&& Label.getInt(ci, cj, ck) == 0) {
							VoxelIndexed<Voxel> cPt = new VoxelIndexed<Voxel>(
									null);
							cPt.setRefPosition(ci, cj, ck);
							List2.add(cPt);
							Label.set(ci, cj, ck, 2);
						}
					}
				}

			} /* End of While(heap) */
			incrementCompletedUnits();
			tmpList = List2;
			List2 = List1;
			List1 = tmpList;
		} /* End of while(flag) */

		for (i = 0; i < XN; i++) {
			for (j = 0; j < YN; j++) {
				for (k = 0; k < ZN; k++) {
//					if (orivol.getUByte(i, j, k) > 0)
					if (orivol.getInt(i, j, k) > 0)
						orivol.set(i, j, k, (int) 254);
				}
			}
		}
		markCompleted();
		return orivol;
	}

	private int ConInCube6(int[][][] cubeMat, int label) {
		int i, j, k, index, index6;
		int number;
		Vector<VoxelIndexed<Voxel>> NeiQ = new Vector<VoxelIndexed<Voxel>>();
		int ci, cj, ck, ni, nj, nk;
		// int ioff, joff, koff;
		number = 0;
		for (index6 = 0; index6 < 6; index6++) {
			i = xoff6[index6];
			j = yoff6[index6];
			k = zoff6[index6];
			if (cubeMat[i + 1][j + 1][k + 1] != label)
				continue;
			VoxelIndexed<Voxel> cPt = new VoxelIndexed<Voxel>(null);
			cPt.setRefPosition(i, j, k);
			number++;
			if (number > 1) {
				return number;
			}
			cubeMat[i + 1][j + 1][k + 1] = 0;
			NeiQ.add(cPt);
			while (!NeiQ.isEmpty()) {
				cPt = NeiQ.remove(0);
				ci = cPt.getRow();
				cj = cPt.getColumn();
				ck = cPt.getSlice();
				for (index = 0; index < 6; index++) {
					ni = ci + xoff6[index];
					nj = cj + yoff6[index];
					nk = ck + zoff6[index];
					if (ni >= (-1) && ni <= 1 && nj >= (-1) && nj <= 1
							&& nk >= (-1) && nk <= 1) {
						if (cubeMat[ni + 1][nj + 1][nk + 1] == label) {
							VoxelIndexed<Voxel> nPt = new VoxelIndexed<Voxel>(
									null);
							nPt.setRefPosition(ni, nj, nk);
							cubeMat[ni + 1][nj + 1][nk + 1] = 0;
							NeiQ.add(nPt);
						}
					}
				}
			}
		}

		return number;
	}
}
