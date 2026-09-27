package edu.jhu.ece.iacl.algorithms.tgdm;

import javax.vecmath.Vector3d;

public abstract class ProfileTGDM {
	protected double MaxDist;
	protected double MaxSignDist;
	protected double finalMaxSignDist=20;
	protected double HowClose;
	protected double curvatureForce;
	protected double pressureForce;
	protected double externalForce;
	protected double augmentGradient=0;
	protected double externalForceThreshold=0;
	protected double deltaT;
	protected double finalDeltaT;
	protected int internalIters;
	protected int maxIters;
	protected int finalInternalIters;
	protected boolean skipLastIters=true;
	public enum Tracking{OFF,SPHERE,POINT,SHELL};
	public ProfileTGDM(
			double maxDist,
			double howClose,
			double maxSignDist,
			double finalMaxSignDist,
			int internalIters,
			int finalInternalIters,
			double deltaT,
			double finalDeltaT){
		this.MaxDist=maxDist;
		this.HowClose=howClose;
		this.MaxSignDist=maxSignDist;
		this.internalIters=internalIters;
		this.finalInternalIters=finalInternalIters;
		this.deltaT=deltaT;
		this.finalDeltaT=finalDeltaT;
	}
	public ProfileTGDM(
			double maxDist,
			double howClose,
			double maxSignDist,
			double finalMaxSignDist,
			int internalIters,
			int finalInternalIters,
			double deltaT,
			double finalDeltaT,
			double curvatureForce,
			double externalForce,
			double pressureForce,int maxiters){
		this(maxDist,howClose,maxSignDist,finalMaxSignDist,internalIters,finalInternalIters,deltaT,finalDeltaT);
		this.curvatureForce=curvatureForce;
		this.pressureForce=pressureForce;
		this.externalForce=externalForce;
		this.maxIters=maxiters;
	}
	/**
	 * 
	 * @param cu advection x
	 * @param cv advection y
	 * @param cw advection z
	 * @param cg 
	 * @param GPhi levelset gradient magnitude
	 * @param K laplacian curvature
	 * @param G mean curvature
	 * @param Vp pressure force
	 * @param Dmx level set diff -X
	 * @param Dpx level set diff +X
	 * @param Dmy level set diff -Y
	 * @param Dpy level set diff +Y
	 * @param Dmz level set diff -Z
	 * @param Dpz level set diff +Z
	 * @return
	 */
	//public abstract double getLevelSetChange(double phi,double cu,double cv,double cw,double cg, double GPhi, double K, double G,double Vp,double Dmx,double Dpx,double Dmy,double Dpy,double Dmz,double Dpz);
	public void setExternalForceThreshold(double thresh){
		this.externalForceThreshold=thresh;
	}
	public abstract double getEffectiveForce(double K,double G,Vector3d curvatureGrad,double regionForce,Vector3d regionGrad,Vector3d advectionForce,Vector3d advectionGrad);
	public abstract Vector3d getEffectiveVelocity(double K,double G,Vector3d curvatureGrad,double regionForce,Vector3d regionGrad,Vector3d advectionForce,Vector3d advectionGrad);
	
	public double getExternalForceThreshold(){
		return this.externalForceThreshold;
	}
	public void setGradientAugment(double aug){
		this.augmentGradient=aug;
	}
	public double getGradientAugment(){
		return augmentGradient;
	}
	public double getMaxSignDist(){
		return MaxSignDist;
	}
	public double getCurvatureForce(){
		return curvatureForce;
	}
	public double getPressureForce(){
		return pressureForce;
	}
	public double getMaxDist(){
		return MaxDist;
	}
	public double getHowClose(){
		return HowClose;
	}
	public double getDeltaT() {
		return deltaT;
	}
	public void setDeltaT(double deltaT) {
		this.deltaT = deltaT;
	}
	public double getFinalMaxSignDist(){
		return finalMaxSignDist;
	}
	public double getFinalDeltaT() {
		return finalDeltaT;
	}
	public void setFinalDeltaT(double finalDeltaT) {
		this.finalDeltaT = finalDeltaT;
	}
	public int getFinalInternalIters() {
		return finalInternalIters;
	}
	public void setFinalInternalIters(int finalInternalIters) {
		this.finalInternalIters = finalInternalIters;
	}
	public int getInternalIters() {
		return internalIters;
	}
	public void setInternalIters(int internalIters) {
		this.internalIters = internalIters;
	}
	public double getExternalForce() {
		return externalForce;
	}
	public void setExternalForce(double externalForce) {
		this.externalForce = externalForce;
	}
	public void setCurvatureForce(double curvatureForce) {
		this.curvatureForce = curvatureForce;
	}
	public void setPressureForce(double pressureForce) {
		this.pressureForce = pressureForce;
	}
	public void setSkipLastIters(boolean lastIters){
		this.skipLastIters=lastIters;
	}
	public boolean isSkipLastIters(){
		return skipLastIters;
	}
	public void setOuterIters(int iters){
		this.maxIters=iters;
	}
	public int getOuterIters(){
		return this.maxIters;
	}
}
