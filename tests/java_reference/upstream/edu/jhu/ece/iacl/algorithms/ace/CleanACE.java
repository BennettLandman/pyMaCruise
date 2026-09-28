package edu.jhu.ece.iacl.algorithms.ace;

import java.util.Vector;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.MaskVolume26;
import edu.jhu.ece.iacl.jist.structures.image.VoxelIndexed;
import edu.jhu.ece.iacl.jist.structures.image.VoxelInt;

/**
 * Clean ACE membership.
 * 
 * @author Blake Lucas
 * 
 */
public class CleanACE extends AbstractCalculation {
	public CleanACE() {
		super();
		setLabel("Clean ACE");
	}

	public CleanACE(AbstractCalculation calc) {
		super(calc);
		setLabel("Clean ACE");
	}

	public void solve(ImageData orivol, int thred) {

		int curlabel, curSize;
		int i, j, k, minX, maxX, minY, maxY, minZ, maxZ;
		int ci, cj, ck;
		int total, removed;
		int XN = orivol.getRows();

		VoxelIndexed<VoxelInt> cPt = new VoxelIndexed<VoxelInt>(new VoxelInt());
		int YN = orivol.getCols();
		int ZN = orivol.getSlices();
		ImageDataInt Label = new ImageDataInt(XN, YN, ZN);
		int[][][] LabelMat = Label.toArray3d();
		for (i = 0; i < XN; i++) {
			for (j = 0; j < YN; j++) {
				for (k = 0; k < ZN; k++) {
					LabelMat[i][j][k] = 0;
				}
			}
		}
		curlabel = 1;
		total = 0;
		removed = 0;
		Vector<VoxelIndexed<VoxelInt>> NeiQ = new Vector<VoxelIndexed<VoxelInt>>();
		int ni, nj, nk;
		int index = 0;
		MaskVolume26 mask = new MaskVolume26();
		byte neighborsX[] = mask.getNeighborsX();
		byte neighborsY[] = mask.getNeighborsY();
		byte neighborsZ[] = mask.getNeighborsZ();
		setTotalUnits(ZN);
		for (k = 1; k < (ZN - 1); k++) {

			incrementCompletedUnits();
			for (j = 1; j < (YN - 1); j++) {
				for (i = 1; i < (XN - 1); i++) {
//					if (!(orivol.getUByte(i, j, k) > 0 && LabelMat[i][j][k] == 0))
					if (!(orivol.getInt(i, j, k) > 0 && LabelMat[i][j][k] == 0))
						continue;
					VoxelIndexed<VoxelInt> Pt = new VoxelIndexed<VoxelInt>(
							new VoxelInt());
					Pt.setRefPosition(i, j, k);

					minX = i;
					maxX = i;
					minY = j;
					maxY = j;
					minZ = k;
					maxZ = k;

					LabelMat[i][j][k] = curlabel;
					NeiQ.add(Pt);

					curSize = 0;

					while (!NeiQ.isEmpty()) {
						cPt = NeiQ.remove(0);
						curSize++;
						ci = cPt.getRow();
						cj = cPt.getColumn();
						ck = cPt.getSlice();

						if (minX > ci)
							minX = ci;
						if (maxX < ci)
							maxX = ci;
						if (minY > cj)
							minY = cj;
						if (maxY < cj)
							maxY = cj;
						if (minZ > ck)
							minZ = ck;
						if (maxZ < ck)
							maxZ = ck;

						for (index = 0; index < 18; index++) {
							ni = ci + neighborsX[index];
							nj = cj + neighborsY[index];
							nk = ck + neighborsZ[index];
							if (ni >= 0 && ni < XN && nj >= 0 && nj < YN
									&& nk >= 0 && nk < ZN) {
								if (LabelMat[ni][nj][nk] == 0
//										&& orivol.getUByte(ni, nj, nk) > 0) {
										&& orivol.getInt(ni, nj, nk) > 0) {
									// Unvisited end point found
									VoxelIndexed<VoxelInt> nPt = new VoxelIndexed<VoxelInt>(
											new VoxelInt());

									nPt.setRefPosition(ni, nj, nk);
									Label.set(ni, nj, nk, curlabel);
									NeiQ.add(nPt);
								}
							}

						}
					}
					NeiQ.clear();
					if (curSize <= thred) { // Clear this ACE point CC
						total = total + curSize;
						removed++;
						for (ci = minX; ci <= maxX; ci++) {
							for (cj = minY; cj <= maxY; cj++) {
								for (ck = minZ; ck <= maxZ; ck++) {

									if (LabelMat[ci][cj][ck] == curlabel) {
										LabelMat[ci][cj][ck] = 0;
//										orivol.set(ci, cj, ck, 0);
										orivol.set(ci, cj, ck, (int) 0);
									}
								}
							}
						}
						curlabel--;
					}
					curlabel++;

				}
			}
		}
		markCompleted();
	}

}
