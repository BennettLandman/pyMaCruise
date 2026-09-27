package edu.jhu.ece.iacl.algorithms.topology;

public interface ConnectivityRule {
	public static final int CONNECT_6_26 = 0;

	public static final int CONNECT_6_18 = 6;

	public static final int CONNECT_18_6 = 18;

	public static final int CONNECT_26_6 = 26;

	final static int[] xoff6 = { 1, 0, 0, -1, 0, 0 };

	final static int[] yoff6 = { 0, 1, 0, 0, -1, 0 };

	final static int[] zoff6 = { 0, 0, 1, 0, 0, -1 };

	final static int[] xoff18 = { 1, 0, 0, -1, 0, 0, 1, -1, 1, 0, 0, 0, 1, 1,
			-1, -1, 0, -1 };

	final static int[] yoff18 = { 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, 1, -1, 0, 0,
			0, -1, -1, 0 };

	final static int[] zoff18 = { 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, -1, 1, 1, -1,
			1, 0, -1, -1 };
	
	final static int[] xoff26 = { 1, 0, 0, -1, 0, 0, 1, -1, 1, 0, 0, 0, 1, 1,
		-1, -1, 0, -1 ,   1,-1,1,-1,1,-1,1,-1 };

	final static int[] yoff26 = { 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, 1, -1, 0, 0,
			0, -1, -1, 0 , 1,1,-1,-1,1,1,-1,-1 };
	
	final static int[] zoff26 = { 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, -1, 1, 1, -1,
			1, 0, -1, -1, 1,1,1,1,-1,-1,-1,-1  };
/*
	final static int[] xoff26 = { 0, 1, 0, 0, -1, 0, 1, 1, -1, -1, 1, 1, -1,
			-1, 0, 0, 0, 0, -1, -1, 1, 1, -1, -1, 1, 1 };

	final static int[] yoff26 = { 1, 0, 0, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0,
			1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };

	final static int[] zoff26 = { 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 1, -1,
			1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1 };
*/
	public int check(int[][][] cube, int label);
	
}
