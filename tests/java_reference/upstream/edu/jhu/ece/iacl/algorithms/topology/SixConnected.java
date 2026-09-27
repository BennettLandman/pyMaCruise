package edu.jhu.ece.iacl.algorithms.topology;

import java.util.LinkedList;

import edu.jhu.ece.iacl.algorithms.tgdm.GenericTGDM.NBPoint;

public class SixConnected implements ConnectivityRule{
	
	public int check(int[][][] cube, int label) {
			/*
			 * All visited points will be assigned a label of 0, so the parameter
			 * 'label' cannot be 0.
			 */
			/* Assume cube is already in Geodesic Neighborhood */

			int i, j, k, index, index6;
			int number;
			NBPoint cPt, nPt;
			LinkedList<NBPoint> NeiQ = new LinkedList<NBPoint>();
			int ci, cj, ck, ni, nj, nk;

			number = 0;
			for (index6 = 0; index6 < 6; index6++) { /* 6-adjacent */
				k = xoff6[index6];
				i = yoff6[index6];
				j = zoff6[index6];
				if (cube[k + 1][i + 1][j + 1] != label)
					continue;
				cPt = new NBPoint();
				cPt.x = k;
				cPt.y = i;
				cPt.z = j;
				number++;
				if (number > 1) {
					NeiQ.clear();
					return number;
				}

				cube[k + 1][i + 1][j + 1] = 0;
				NeiQ.addLast(cPt);

				while (!NeiQ.isEmpty()) {
					cPt = (NBPoint) NeiQ.removeFirst();
					ci = cPt.y;
					cj = cPt.z;
					ck = cPt.x;
					for (index = 0; index < 6; index++) { /* 6 connectivity */
						ni = ci + yoff6[index];
						nj = cj + zoff6[index];
						nk = ck + xoff6[index];
						if (ni >= (-1) && ni <= 1 && nj >= (-1) && nj <= 1
								&& nk >= (-1) && nk <= 1) {
							if (cube[nk + 1][ni + 1][nj + 1] == label) {
								/* Unlabelled object point found */
								nPt = new NBPoint();
								nPt.x = nk;
								nPt.y = ni;
								nPt.z = nj;
								cube[nk + 1][ni + 1][nj + 1] = 0;
								NeiQ.addLast(nPt);

							}
						}
					}
				}
			}

			NeiQ.clear();
			return number;

		}// ConICube6

}
