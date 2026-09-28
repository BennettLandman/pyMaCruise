package edu.jhu.ece.iacl.algorithms.volume;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.data.BinaryMinHeap;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.MaskVolume6;
import edu.jhu.ece.iacl.jist.structures.image.VoxelFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelIndexed;

/**
 * Calculate distance field, distance field characteristic and skeleton using
 * Fast-marching method
 * 
 * @author Blake Lucas
 */

public class DistanceFieldAndSkeleton extends AbstractCalculation {
	private static final byte ALIVE = 1;

	private static final byte NBAND = 2;

	private static final byte FARAWAY = 3;

	private static final float INFINITY = 10000;

	private static final DistanceFieldAndSkeleton sdf = new DistanceFieldAndSkeleton();

	public DistanceFieldAndSkeleton() {
		super();
		setLabel("Fast-Marching Distance");
	}

	public DistanceFieldAndSkeleton(AbstractCalculation parent) {
		super(parent);
		setLabel("Fast-Marching Distance");
	}

	public static ImageDataFloat doSolve(ImageDataFloat vol,
			ImageDataFloat csf) {
		return sdf.solve(vol, csf);
	}

	private ImageDataFloat distCharVol;

	private ImageDataUByte skel;

	public ImageData getDistanceFuncCharacteristic() {
		return distCharVol;
	}

	public ImageData getSkeleton() {
		return skel;
	}

	public ImageDataFloat solve(ImageDataFloat vol, ImageDataFloat speed) {
		String name = vol.getName();
		int XN = vol.getRows();
		int YN = vol.getCols();
		int ZN = vol.getSlices();
		int LX, HX, LY, HY, LZ, HZ;
		short NSFlag, WEFlag, FBFlag;
		int i, j, k, koff;
		int nj, nk, ni;
		double newvalue;
		double s = 0, t = 0, w = 0;
		double result;
		// double IFij;

		MaskVolume6 mask = new MaskVolume6();
		byte[] neighborsX = mask.getNeighborsX();
		byte[] neighborsY = mask.getNeighborsY();
		byte[] neighborsZ = mask.getNeighborsZ();
		VoxelIndexed<VoxelFloat> he;
		double Nv = 0, Sv = 0, Wv = 0, Ev = 0, Fv = 0, Bv = 0, Cv = 0;

		double Nd = 0, Sd = 0, Wd = 0, Ed = 0, Fd = 0, Bd = 0;
		double minT = 0, minD = 0;
		/*
		 * Value at six neighours of a pixel
		 */
		byte Nl, Sl, Wl, El, Fl, Bl; /* Label at six neighours of a pixel */
		setLabel("Fast-Marching Distance Field w/ Skeleton");
		setTotalUnits(ZN * 2);
		ImageDataFloat distVol = new ImageDataFloat(XN, YN, ZN);
		float[][][] distVolM = distVol.toArray3d();
		skel = new ImageDataUByte(XN, YN, ZN);
		byte[][][] skelM = skel.toArray3d();
		distCharVol = new ImageDataFloat(XN, YN, ZN);
		float[][][] distCharVolM = distCharVol.toArray3d();
		ImageDataUByte labelVol = new ImageDataUByte(XN, YN, ZN);
		byte[][][] labelVolM = labelVol.toArray3d();
		float[][][] csfM=speed.toArray3d();
		float[][][] volM = vol.toArray3d();
		for (i = 0; i < XN; i++) {
			for (j = 0; j < YN; j++) {
				for (k = 0; k < ZN; k++) {
					labelVolM[i][j][k]=( FARAWAY);
				}
			}
		}
		int countAlive = 0;
		for (k = 0; k < ZN; k++) {

			incrementCompletedUnits();
			for (j = 0; j < YN; j++) {
				for (i = 0; i < XN; i++) {

					if (volM[i][j][k] == 0) {
						distVolM[i][j][k] = (0);
						distCharVolM[i][j][k] = (0);
						labelVolM[i][j][k] = (ALIVE);
						countAlive++;
					} else {
						LX = (i == 0) ? 1 : 0;
						HX = (i == (XN - 1)) ? 1 : 0;

						LY = (j == 0) ? 1 : 0;
						HY = (j == (YN - 1)) ? 1 : 0;

						LZ = (k == 0) ? 1 : 0;
						HZ = (k == (ZN - 1)) ? 1 : 0;

						NSFlag = 0;
						WEFlag = 0;
						FBFlag = 0;

						Nv = volM[i][ j - 1 + LY][ k];
						Sv = volM[i][ j + 1 - HY][ k];
						Wv = volM[i - 1 + LX][ j][ k];
						Ev = volM[i + 1 - HX][ j][ k];
						Fv = volM[i][ j][ k + 1 - HZ];
						Bv = volM[i][ j][ k - 1 + LZ];
						Cv = volM[i][j][k];
						if (Nv * Cv < 0) {
							NSFlag = 1;
							s = Nv;
						}
						if (Sv * Cv < 0) {
							if (NSFlag == 0) {
								NSFlag = 1;
								s = Sv;
							} else {
								s = (Math.abs(Nv) > Math.abs(Sv)) ? Nv : Sv;
							}
						}
						if (Wv * Cv < 0) {
							WEFlag = 1;
							t = Wv;
						}
						if (Ev * Cv < 0) {
							if (WEFlag == 0) {
								WEFlag = 1;
								t = Ev;
							} else {
								t = (Math.abs(Ev) > Math.abs(Wv)) ? Ev : Wv;
							}
						}
						if (Fv * Cv < 0) {
							FBFlag = 1;
							w = Fv;
						}
						if (Bv * Cv < 0) {
							if (FBFlag == 0) {
								FBFlag = 1;
								w = Bv;
							} else {
								w = (Math.abs(Fv) > Math.abs(Bv)) ? Fv : Bv;
							}
						}

						result = 0;
						if (NSFlag != 0) {
							s = Cv / (Cv - s);
							result += 1.0 / (s * s);
						}
						if (WEFlag != 0) {
							t = Cv / (Cv - t);
							result += 1.0 / (t * t);
						}
						if (FBFlag != 0) {
							w = Cv / (Cv - w);
							result += 1.0 / (w * w);
						}
						if (result == 0)
							continue;
						/*
						 * if(count++<100){ System.out.println("ALIVE1
						 * ("+i+","+j+","+k+") "+result);
						 * System.out.println(Sv+" "+Wv+" "+Ev+" "+Fv+" "+Bv+"
						 * "+Cv);
						 *  }
						 */
						countAlive++;
						labelVolM[i][j][k] = (ALIVE);
						result = Math.sqrt(result);
						distCharVolM[i][j][k]=(float)( 1.0f / result);
						if (csfM[i][j][k] == 0) {
							distVolM[i][j][k]=( INFINITY);
						} else {
							/* IFij = 1.0/F[d][i][j]; */
							distVolM[i][j][k]=(float)( 1.0 / (result * csfM[
									i][ j][ k]));
						}
					}

				}
			}
		}
		BinaryMinHeap heap = new BinaryMinHeap(countAlive, XN, YN, ZN);
		/* Initialize NarrowBand Heap */
		for (k = 0; k < ZN; k++) {

			incrementCompletedUnits();
			for (j = 0; j < YN; j++) {
				for (i = 0; i < XN; i++) {

					if (labelVol.get(i, j, k).shortValue() != ALIVE)
						continue;
					/* Put its 6 neighbors into NarrowBand */
					for (koff = 0; koff < MaskVolume6.length; koff++) {/*
																		 * Find
																		 * six
																		 * neighbouring
																		 * points
																		 */
						ni = i + neighborsX[koff];
						nj = j + neighborsY[koff];
						nk = k + neighborsZ[koff];

						if (nj < 0 || nj >= YN || nk < 0 || nk >= ZN || ni < 0
								|| ni >= XN)
							continue; /* Out of computational Boundary */

						if (labelVolM[ni][nj][nk] != FARAWAY)
							continue;
						labelVolM[ni][nj][nk]=( NBAND);
						minT = 1E6f;
						/*
						 * Note: Only ALIVE points contribute to the distance
						 * computation
						 */
						/* Neighbour to the north */
						if (nj > 0) {
							Nv = distVolM[ni][ nj - 1][ nk];
							Nd = distCharVolM[ni][ nj - 1][ nk];
							Nl = labelVolM[ni][ nj - 1][ nk];
							if (Nl == ALIVE) {
								if (minT > Nv) {
									minT = Nv;
									minD = Nd;
								}
							}
						} else {
							Nl = 0;
						}
						/* Neighbour to the south */
						if (nj < YN - 1) {
							Sv = distVolM[ni][ nj + 1][ nk];
							Sd = distCharVolM[ni][ nj + 1][ nk];
							Sl = labelVolM[ni][ nj + 1][ nk];
							if (Sl == ALIVE) {
								if (minT > Sv) {
									minT = Sv;
									minD = Sd;
								}
							}
						} else {
							Sl = 0;
						}
						/* Neighbour to the east */
						if (nk < ZN - 1) {
							Ev = distVolM[ni][ nj][ nk + 1];
							Ed = distCharVolM[ni][ nj][ nk + 1];
							El = labelVolM[ni][ nj][ nk + 1];
							if (El == ALIVE) {
								if (minT > Ev) {
									minT = Ev;
									minD = Ed;
								}
							}
						} else {
							El = 0;
						}
						/* Neighbour to the west */
						if (nk > 0) {
							Wv = distVolM[ni][ nj][ nk - 1];
							Wd = distCharVolM[ni][ nj][ nk - 1];
							Wl = labelVolM[ni][ nj][ nk - 1];
							if (Wl == ALIVE) {
								if (minT > Wv) {
									minT = Wv;
									minD = Wd;
								}
							}
						} else {
							Wl = 0;
						}
						/* Neighbour to the front */
						if (ni < XN - 1) {
							Fv = distVolM[ni + 1][ nj][ nk];
							Fd = distCharVolM[ni + 1][ nj][ nk];
							Fl = labelVolM[ni + 1][ nj][ nk];
							if (Fl == ALIVE) {
								if (minT > Fv) {
									minT = Fv;
									minD = Fd;
								}
							}
						} else {
							Fl = 0;
						}
						/* Neighbour to the back */
						if (ni > 0) {
							Bv = distVolM[ni - 1][ nj][ nk];
							Bd = distCharVolM[ni - 1][ nj][ nk];
							Bl = labelVolM[ni - 1][ nj][ nk];
							if (Bl == ALIVE) {
								if (minT > Bv) {
									minT = Bv;
									minD = Bd;
								}
							}
						} else {
							Bl = 0;
						}
						/*
						 * Update the value of this to-be-updated NarrowBand
						 * point
						 */
						newvalue = march(Nv, Sv, Ev, Wv, Fv, Bv, Nd, Sd, Ed,
								Wd, Fd, Bd, minT, minD, Nl, Sl, El, Wl, Fl, Bl,
								csfM[ni][ nj][ nk], distCharVol, skel,
								ni, nj, nk);

						// newvalue = march(Nv, Sv, Ev, Wv, Fv, Bv, Nl, Sl,
						// El,Wl, Fl, Bl);

						distVolM[ni][nj][nk]=(float)( newvalue);

						VoxelIndexed<VoxelFloat> vox = new VoxelIndexed<VoxelFloat>(
								new VoxelFloat((float) newvalue));
						vox.setRefPosition(ni, nj, nk);
						heap.add(vox);
					}

				}
			}
		}
		/*
		 * Begin Fast Marching to get the unsigned distance function inwords and
		 * outwards simultaneously since points inside and outside the contour
		 * won't interfere with each other
		 */
		setLabel("Fast-Marching Distance Field w/ Skeleton");
		setTotalUnits(heap.size());
		while (!heap.isEmpty()) { /* There are still points not yet accepted */
			// if(heap.size()%50000==0)System.out.println(heap.size());
			incrementCompletedUnits();
			he = (VoxelIndexed<VoxelFloat>) heap.remove(); /*
															 * Label the point
															 * with smallest
															 * value among all
															 * NarrowBand points
															 * as ALIVE
															 */

			/* Put the smallest heap element to ALIVE */
			i = he.getRow();
			j = he.getColumn();
			k = he.getSlice();
			distVolM[i][j][k] = (he.getFloat());
			labelVolM[i][j][k] = (ALIVE);
			/* Update its neighbor */
			/*
			 * Put FARAWAY neighbour into NarrowBand, Recompute values at
			 * NarrowBand neighbours, Keep ALIVE (Accepted) neighbour unchanged
			 */
			for (koff = 0; koff < MaskVolume6.length; koff++) {
				ni = i + neighborsX[koff];
				nj = j + neighborsY[koff];
				nk = k + neighborsZ[koff];

				if (nj < 0 || nj >= YN || nk < 0 || nk >= ZN || ni < 0
						|| ni >= XN)
					continue; /* Out of boundary */
				if (labelVolM[ni][nj][nk] == ALIVE)
					continue; /* Don't change ALIVE neighbour */

				/* ReCompute the value at (nk, nj, ni) */
				/*
				 * Get the values and labels of six neighbours of the
				 * to-be-updated point. The labels are needed since only values
				 * at ALIVE neighbours will be used to update the value of the
				 * to-be-updated point
				 */

				/*
				 * Note: Only ALIVE points contribute to the distance
				 * computation
				 */
				/* Neighbour to the north */
				minT = 1E6f;
				if (nj > 0) {
					Nv = distVolM[ni][ nj - 1][ nk];
					Nd = distCharVolM[ni][ nj - 1][ nk];
					Nl = labelVolM[ni][ nj - 1][ nk];
					if (Nl == ALIVE) {
						if (minT > Nv) {
							minT = Nv;
							minD = Nd;
						}
					}
				} else
					Nl = 0;

				/* Neighbour to the south */
				if (nj < YN - 1) {
					Sv = distVolM[ni][ nj + 1][ nk];
					Sd = distCharVolM[ni][ nj + 1][ nk];
					Sl = labelVolM[ni][ nj + 1][ nk];
					if (Sl == ALIVE) {
						if (minT > Sv) {
							minT = Sv;
							minD = Sd;
						}
					}
				} else
					Sl = 0;

				/* Neighbour to the east */
				if (nk < ZN - 1) {
					Ev = distVolM[ni][ nj][ nk + 1];
					Ed = distCharVolM[ni][ nj][ nk + 1];
					El = labelVolM[ni][ nj][ nk + 1];
					if (El == ALIVE) {
						if (minT > Ev) {
							minT = Ev;
							minD = Ed;
						}
					}
				} else
					El = 0;

				/* Neighbour to the west */
				if (nk > 0) {
					Wv = distVolM[ni][ nj][ nk - 1];
					Wd = distCharVolM[ni][ nj][ nk - 1];
					Wl = labelVolM[ni][ nj][ nk - 1];
					if (Wl == ALIVE) {
						if (minT > Wv) {
							minT = Wv;
							minD = Wd;
						}
					}
				} else
					Wl = 0;

				/* Neighbour to the front */
				if (ni < XN - 1) {
					Fv = distVolM[ni + 1][ nj][ nk];
					Fd = distCharVolM[ni + 1][ nj][ nk];
					Fl = labelVolM[ni + 1][ nj][ nk];
					if (Fl == ALIVE) {
						if (minT > Fv) {
							minT = Fv;
							minD = Fd;
						}
					}
				} else
					Fl = 0;

				/* Neighbour to the back */
				if (ni > 0) {
					Bv = distVolM[ni - 1][ nj][ nk];
					Bd = distCharVolM[ni - 1][ nj][ nk];
					Bl = labelVolM[ni - 1][ nj][ nk];
					if (Bl == ALIVE) {
						if (minT > Bv) {
							minT = Bv;
							minD = Bd;
						}
					}
				} else
					Bl = 0;

				/* Update the value of this to-be-updated NarrowBand point */
				newvalue = march(Nv, Sv, Ev, Wv, Fv, Bv, Nd, Sd, Ed, Wd, Fd,
						Bd, minT, minD, Nl, Sl, El, Wl, Fl, Bl, csfM[
								ni][ nj][ nk], distCharVol, skel, ni, nj, nk);

				/*
				 * If it was a FARAWAY point, add it to the NarrowBand Heap;
				 * otherwise, just update its value using the backpointer
				 */
				VoxelIndexed<VoxelFloat> vox = new VoxelIndexed<VoxelFloat>(
						new VoxelFloat((float) newvalue));
				vox.setRefPosition(ni, nj, nk);
				if (labelVolM[ni][ nj][ nk] == NBAND) {
					heap.change(ni, nj, nk, vox);
				} else {
					decrementCompletedUnits();
					heap.add(vox);
					labelVolM[ni][nj][nk]=( NBAND);
				}
			} /* End of updating 6 neighbours */

		}/* End of marching loop */

		/* Add signs to the unsigned distance function */
		for (i = 0; i < XN; i++) {
			for (j = 0; j < YN; j++) {
				for (k = 0; k < ZN; k++) {
					if (volM[i][j][k] < 0) {
						distVolM[i][j][k] = (-distVolM[i][j][k]);

					}
				}
			}
		}
		markCompleted();
		distVol.setName(name);
		return distVol;
	}

	double march(double Nv, double Sv, double Ev, double Wv, double Fv,
			double Bv, double Nd, double Sd, double Ed, double Wd, double Fd,
			double Bd, double minT, double minD, int Nl, int Sl, int El,
			int Wl, int Fl, int Bl, double F, ImageDataFloat DV, ImageDataUByte skel,
			int i, int j, int k) {
		/*
		 * Compute the new value at a NarrowBand point using values of its
		 * ALIVE(Accepted) 6-connected neighbours Note: at least one of its 6
		 * neighbours must be ALIVE If ALIVE neighbours exist only in the
		 * north-south (or east-west, or front-back) neighbours, the updated
		 * value is equal to the value of the smaller one plus 1 (1/F(k,i,j) if
		 * F(k,i,j) ne 1). Otherwise, a quadratic equation need to be solved,
		 * and the bigger root is taken as the updated value
		 */
		/*
		 * The following program assumes F(k,i,j) = 1; if not, replace IFij and
		 * ISFij with the true values. Note, F(k,i,j) should be positive!
		 */

		/* Suppose a, b, and c are the minimum T value in three directions */
		double s, s2; /* s = a + b +c; s2 = a*a + b*b +c*c */
		double sD, sD2;
		double tmp;
		int count;
		double IFij, ISFij;

		if (F > 0)
			minT = minT + 1 / F;

		s = 0;
		s2 = 0;
		count = 0;
		sD = 0;
		sD2 = 0;
		float[][][] DVM=DV.toArray3d();
		byte[][][] skelM=skel.toArray3d();
		skelM[i][j][k]=( 0);

		if (Nl == ALIVE && Sl == ALIVE) {
			skelM[i][j][k]=( 1);
			tmp = Math.min(Nv, Sv); /* Take the smaller one if both ALIVE */
			if (tmp <= minT) {
				s += tmp;
				s2 += tmp * tmp;
				count++;
				if (tmp == Nv) {
					sD += Nd;
					sD2 += Nd * Nv;
				} else {
					sD += Sd;
					sD2 += Sd * Sv;
				}
			}
		} else if (Nl == ALIVE && (Nv <= minT)) {
			s += Nv; /* Else, take the ALIVE one */
			s2 += Nv * Nv;
			sD += Nd;
			sD2 += Nd * Nv;
			count++;
		} else if (Sl == ALIVE && (Sv <= minT)) {
			s += Sv;
			s2 += Sv * Sv;
			sD += Sd;
			sD2 += Sd * Sv;
			count++;
		}

		/*
		 * Similarly in the east-west direction to get correct approximation to
		 * the derivative in the x-direction
		 */
		if (El == ALIVE && Wl == ALIVE) {
			skelM[i][j][k]=( 1);
			tmp = Math.min(Ev, Wv); /* Take the smaller one if both ALIVE */
			if (tmp <= minT) {
				s += tmp;
				s2 += tmp * tmp;
				count++;
				if (tmp == Ev) {
					sD += Ed;
					sD2 += Ed * Ev;
				} else {
					sD += Wd;
					sD2 += Wd * Wv;
				}
			}
		} else if (El == ALIVE && (Ev <= minT)) {
			s += Ev; /* Else, take the ALIVE one */
			s2 += Ev * Ev;
			sD += Ed;
			sD2 += Ed * Ev;
			count++;
		} else if (Wl == ALIVE && (Wv <= minT)) {
			s += Wv;
			s2 += Wv * Wv;
			sD += Wd;
			sD2 += Wd * Wv;
			count++;
		}

		/*
		 * Similarly in the front-back direction to get correct approximation to
		 * the derivative in the z-direction
		 */
		if (Fl == ALIVE && Bl == ALIVE) {
			skelM[i][j][k]=( 1);
			tmp = Math.min(Fv, Bv); /* Take the smaller one if both ALIVE */
			if (tmp <= minT) {
				s += tmp;
				s2 += tmp * tmp;
				count++;
				if (tmp == Fv) {
					sD += Fd;
					sD2 += Fd * Fv;
				} else {
					sD += Bd;
					sD2 += Bd * Bv;
				}
			}
		} else if (Fl == ALIVE && (Fv <= minT)) {
			s += Fv; /* Else, take the ALIVE one */
			s2 += Fv * Fv;
			count++;
			sD += Fd;
			sD2 += Fd * Fv;
		} else if (Bl == ALIVE && (Bv <= minT)) {
			s += Bv;
			s2 += Bv * Bv;
			count++;
			sD += Bd;
			sD2 += Bd * Bv;
		}

		/*
		 * count must be greater than zero since there must be one ALIVE pt in
		 * the neighbors
		 */
		if (F == 0) {
			DVM[i][j][k]=(float)( sD / count + 1.0 / Math.sqrt(count));
			return INFINITY;
		}

		IFij = 1.0f / F;
		ISFij = IFij * IFij;

		tmp = s * s - count * (s2 - ISFij);
		if (tmp < 0) {
			DVM[i][j][k]=(float)( minD + 1.0);
			return minT;
		}

		tmp = (s + Math.sqrt((tmp))) / count;
		/* The larger root */
		DVM[i][j][k]=(float)( (IFij + sD * tmp - sD2) / (count * tmp - s));

		return tmp;
	}
}
