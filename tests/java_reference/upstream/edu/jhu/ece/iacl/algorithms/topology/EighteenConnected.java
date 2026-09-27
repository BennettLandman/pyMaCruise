package edu.jhu.ece.iacl.algorithms.topology;

import java.util.LinkedList;

import edu.jhu.ece.iacl.algorithms.tgdm.GenericTGDM.NBPoint;

public class EighteenConnected implements ConnectivityRule{

	public int check(int[][][] cube, int label) {
			/*
			 * All visited points will be assigned a label of 0, so the parameter
			 * 'label' cannot be 0.
			 */

			int i, j, k, index, index18;
			int number;
			NBPoint cPt, nPt;
			LinkedList<NBPoint> NeiQ = new LinkedList<NBPoint>();
			int ci, cj, ck, ni, nj, nk;

			number = 0;
			for (index18 = 0; index18 < 18; index18++) { /* 18-adjacency */

				k = xoff18[index18];
				i = yoff18[index18];
				j = zoff18[index18];
				if (cube[k + 1][i + 1][j + 1] != label)
					continue;
				cPt = new NBPoint();
				cPt.x = k;
				cPt.y = i;
				cPt.z = j;
				number++;
				if (number > 1) {
					NeiQ.clear();
					return number;// isn't the returned # always 2?
				}
				cube[k + 1][i + 1][j + 1] = 0;
				NeiQ.addLast(cPt);

				while (!NeiQ.isEmpty()) {
					cPt = (NBPoint) NeiQ.removeFirst();
					ci = cPt.y;
					cj = cPt.z;
					ck = cPt.x;
					for (index = 0; index < 18; index++) { /* 18 connectivity */
						ni = ci + yoff18[index];
						nj = cj + zoff18[index];
						nk = ck + xoff18[index];
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
		}


}
