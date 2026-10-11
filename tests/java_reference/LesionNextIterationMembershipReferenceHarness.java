import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * Fixed-profile excerpt of VandyCRUISE LesionToadSegmentation.computeMemberships
 * for the 13-class MaCRUISEplus iteration-6 held-fixed comparison.
 *
 * The input is big-endian: nx,ny,nz; refreshed class centroids and standard
 * deviations; lesion centroid and standard deviation; channel weights;
 * shape and smoothing prior vectors; then voxel-major image, shape, direct
 * flags, class memberships, lesion/black-hole/outlier/lesion-relation values,
 * and ventricle/GM/brainstem distance values. Output is the dimensions,
 * source-ordered updated maps/energies, objective and explicit IEEE-NaN count.
 *
 * The caller supplies distance maps and direct relations as fixed inputs.
 * Those are separately materialized from the hash-verified source checkpoint
 * and source-checked distance refresh, not inferred from post-thinning labels.
 */
public final class LesionNextIterationMembershipReferenceHarness {
    private static final int NC = 2, C = 13, LESION = 10;
    private static final float PRIOR = 2.0f / 13.0f;
    private static final float PRIOR_SCALE = 0.1f;
    private static final float SMOOTH = (0.2f / 6.0f) / 14.0f;
    private static final float OUTLIER = 0.1f, ZERO = 1.0e-10f, INF = 1.0e10f, EXP = 1.2f;
    private int nx, ny, nz, n;
    private float[] mu, sd, lesionCentroid, lesionStd, weights, shapePriors, smoothingPriors;
    private float[] image, shape, memberships, classEnergy;
    private byte[] direct;
    private float[] lesion, black, outlier, lesionRelation, denominator;
    private int[] ventricleDistance, gmDistance, brainstemDistance;
    private float[] lesionEnergy, blackHoleEnergy, runningObjective;
    private float[] classSimilarity, lesionSimilarity;
    // computeMemberships declares these before the full x/y/z sweep. Keep
    // their values when a later voxel skips the corresponding class support.
    private final float[] energyState = new float[C];
    private float lesionEnergyCarry = 0.0f, blackHoleEnergyCarry = 0.0f;
    private float objective, previousObjective, energyRatio;

    private int voxel(int x, int y, int z) { return (x * ny + y) * nz + z; }
    private int classOffset(int v, int k) { return v * C + k; }
    private int imageOffset(int v, int c) { return v * NC + c; }

    private void read(String path) throws Exception {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(path)))) {
            nx = in.readInt(); ny = in.readInt(); nz = in.readInt(); n = nx * ny * nz;
            if (nx < 3 || ny < 3 || nz < 3) throw new IllegalArgumentException("invalid spatial dimensions");
            previousObjective = in.readFloat();
            mu = new float[NC * C]; sd = new float[NC * C];
            lesionCentroid = new float[NC]; lesionStd = new float[NC]; weights = new float[NC];
            shapePriors = new float[C]; smoothingPriors = new float[C];
            for (int i=0;i<mu.length;i++) mu[i]=in.readFloat();
            for (int i=0;i<sd.length;i++) sd[i]=in.readFloat();
            for (int i=0;i<NC;i++) lesionCentroid[i]=in.readFloat();
            for (int i=0;i<NC;i++) lesionStd[i]=in.readFloat();
            for (int i=0;i<NC;i++) weights[i]=in.readFloat();
            for (int k=0;k<C;k++) shapePriors[k]=in.readFloat();
            for (int k=0;k<C;k++) smoothingPriors[k]=in.readFloat();
            image = new float[n * NC]; shape = new float[n * C]; direct = new byte[n * C];
            memberships = new float[n * C]; classEnergy = new float[n * C];
            lesion = new float[n]; black = new float[n]; outlier = new float[n]; lesionRelation = new float[n];
            denominator = new float[n]; lesionEnergy = new float[n]; blackHoleEnergy = new float[n];
            ventricleDistance = new int[n]; gmDistance = new int[n]; brainstemDistance = new int[n];
            runningObjective = new float[n];
            for (int v=0;v<n;v++) {
                for (int c=0;c<NC;c++) image[imageOffset(v,c)]=in.readFloat();
                for (int k=0;k<C;k++) shape[classOffset(v,k)]=in.readFloat();
                for (int k=0;k<C;k++) direct[classOffset(v,k)]=in.readByte();
                for (int k=0;k<C;k++) memberships[classOffset(v,k)]=in.readFloat();
                lesion[v]=in.readFloat(); black[v]=in.readFloat(); outlier[v]=in.readFloat(); lesionRelation[v]=in.readFloat();
                ventricleDistance[v]=in.readInt(); gmDistance[v]=in.readInt(); brainstemDistance[v]=in.readInt();
            }
            if (in.read()!=-1) throw new IllegalArgumentException("trailing input bytes");
        }
        classSimilarity = new float[C*C]; lesionSimilarity = new float[C];
    }

    private static float similarity(float distance) {
        return 1.0f / (1.0f + ((1.0f - PRIOR_SCALE) / PRIOR_SCALE) * distance);
    }

    private double powerLookup(float value) {
        if (value < 0.0f || value > 4.0f) return Math.pow((double)value, (double)EXP);
        double step = (double)1.0e-6f;
        int index = (int)(((double)value / step) + 0.5);
        return Math.pow(index * step, (double)EXP);
    }

    private void computeSimilarities() {
        for (int k=0;k<C;k++) for (int m=0;m<C;m++) {
            if (k==m) { classSimilarity[k*C+m]=0.0f; continue; }
            float distance=0.0f;
            if (k==1 || m==1) {
                float q=mu[k]-mu[m], den=sd[k]+sd[m];
                distance += (q*q/den)/den;
            } else for (int c=0;c<NC;c++) {
                float q=mu[c*C+k]-mu[c*C+m];
                float weighted=(weights[c]*q)*q, den=sd[c*C+k]+sd[c*C+m];
                distance += (weighted/den)/den;
            }
            classSimilarity[k*C+m]=similarity(distance);
        }
        for (int k=0;k<C;k++) {
            float distance=0.0f;
            for (int c=0;c<NC;c++) {
                float q=lesionCentroid[c]-mu[c*C+k];
                float weighted=(weights[c]*q)*q, den=lesionStd[c]+sd[c*C+k];
                distance += (weighted/den)/den;
            }
            lesionSimilarity[k]=(k==LESION)?0.0f:similarity(distance);
        }
    }

    private static float shapeFactor(float value) { return 1.0f < 2.0f*value ? 1.0f : 2.0f*value; }
    private float face(float[] map, int v) { return map[v]; }

    private void voxelUpdate(int x, int y, int z) {
        int v=voxel(x,y,z);
        int[] face={voxel(x+1,y,z),voxel(x-1,y,z),voxel(x,y+1,z),voxel(x,y-1,z),voxel(x,y,z+1),voxel(x,y,z-1)};
        float den=0.0f;
        for (int k=0;k<C;k++) {
            int ck=classOffset(v,k);
            float sh=shape[ck];
            // Literal source guard skips the class body for nonpositive shape,
            // then its outer else clears only this class membership. Energy
            // slots and lesion auxiliary maps retain their previous state.
            if (!(sh>0.0f)) { memberships[ck]=0.0f; continue; }
            float raw=0.0f, rawLesion=0.0f, rawBlack=0.0f, rawOut=0.0f;
            if (k==0) {
                for (int c=0;c<NC;c++) raw += weights[c]*(image[imageOffset(v,c)]==0.0f?0.0f:1.0f);
            } else if (k==1 || k==2) {
                raw=OUTLIER;
            } else {
                for (int c=0;c<NC;c++) {
                    float d=image[imageOffset(v,c)]-mu[c*C+k];
                    float weighted=weights[c]*d;
                    raw += weighted*d;
                }
                if (k==LESION) {
                    float d=image[imageOffset(v,1)]-lesionCentroid[1];
                    rawLesion=(float)powerLookup(d*d);
                }
            }
            float other=0.0f, ls=0.0f, bs=0.0f, self=0.0f;
            for (int f=0;f<6;f++) {
                for (int m=0;m<C;m++) if (m!=k) {
                    float q=memberships[classOffset(face[f],m)]; other += q*q;
                }
                float l=face(lesion,face[f]), b=face(black,face[f]);
                ls += l*l; bs += b*b;
            }
            if (k==LESION) for (int f : new int[]{1,0,3,2,5,4}) {
                float q=memberships[classOffset(face[f],k)]; self += q*q;
            }
            float sum=other+ls+bs;
            float term=raw+(smoothingPriors[k]*SMOOTH)*sum;
            float e=raw+((0.5f*smoothingPriors[k])*SMOOTH)*sum;
            // Source assigns lesion_energy = numL before spatial smoothing
            // and prior terms, so the raw lesion data energy is part of the
            // carried energy state as well as the membership denominator.
            float termL=rawLesion, eL=rawLesion;
            if (k==LESION) {
                float lsFor=other+self+bs;
                termL+=(smoothingPriors[k]*SMOOTH)*lsFor;
                eL+=((0.5f*smoothingPriors[k])*SMOOTH)*lsFor;
            }
            float ip=0.0f, il=0.0f;
            for (int m=0;m<C;m++) {
                if (m!=k) {
                    float p=shapePriors[m]*classSimilarity[k*C+m];
                    p=p*shape[classOffset(v,m)]; p=p*shape[classOffset(v,m)]; ip+=p;
                }
                if (k==LESION) {
                    float p=shapePriors[m]*lesionSimilarity[m];
                    p=p*shape[classOffset(v,m)]; p=p*shape[classOffset(v,m)]; il+=p;
                }
            }
            float p=shapePriors[LESION]*lesionSimilarity[k];
            p=p*shape[classOffset(v,LESION)]; p=p*shape[classOffset(v,LESION)]; ip+=p;
            term+=PRIOR*ip; e+=((0.5f*PRIOR)*ip);
            if (k==LESION) { termL+=PRIOR*il; eL+=((0.5f*PRIOR)*il); }
            float cls=term>ZERO?1.0f/term:INF;
            if (direct[ck]==0) cls=shapeFactor(sh)*cls;
            memberships[ck]=cls; energyState[k]=e;
            if (k==LESION) {
                float l=termL>ZERO?1.0f/termL:INF;
                if (direct[ck]==0) l=shapeFactor(sh)*l;
                lesion[v]=l; lesionRelation[v]=(direct[ck]!=0)?1.0f:shapeFactor(sh); lesionEnergyCarry=eL;
                int vd=ventricleDistance[v];
                if (vd!=-1) {
                    if (vd==0) { lesion[v]=0.5f*lesion[v]/2.0f; lesionRelation[v]*=(0.5f/2.0f); }
                    else { lesion[v]=vd*lesion[v]/2.0f; lesionRelation[v]*=(vd/2); }
                } else if (gmDistance[v]!=-1) {
                    int gd=gmDistance[v];
                    if (gd==0) { lesion[v]=0.5f*lesion[v]/3.0f; lesionRelation[v]*=(0.5f/3.0f); }
                    else { lesion[v]=gd*lesion[v]/3.0f; lesionRelation[v]*=(gd/3); }
                }
                if (brainstemDistance[v]!=-1) {
                    lesion[v]=brainstemDistance[v]*lesion[v]/25.0f;
                    lesionRelation[v]*=(brainstemDistance[v]/25);
                }
                den+=lesion[v];
            }
            if (k==4 && ventricleDistance[v]!=-1) memberships[ck]=memberships[ck]*ventricleDistance[v]/2.0f;
            den+=memberships[ck];
        }
        for (int k=0;k<C;k++) classEnergy[classOffset(v,k)]=energyState[k];
        lesionEnergy[v]=lesionEnergyCarry;
        blackHoleEnergyCarry=0.0f; // This 13-class profile contains no BLACKHOLE class.
        blackHoleEnergy[v]=blackHoleEnergyCarry;
        denominator[v]=den;
        if (den>0.0f) {
            for (int k=0;k<C;k++) {
                int ck=classOffset(v,k);
                memberships[ck]=memberships[ck]/den;
                if (k==LESION) lesion[v]=lesion[v]/den;
            }
        }
        if (den>0.0f) for (int k=0;k<C;k++) {
            float q=classEnergy[classOffset(v,k)]*memberships[classOffset(v,k)];
            q=q*memberships[classOffset(v,k)]; objective+=q;
        }
        float q=lesionEnergy[v]*lesion[v]; q=q*lesion[v]; objective+=q;
        float b=blackHoleEnergy[v]*black[v]; b=b*black[v]; objective+=b;
        runningObjective[v]=objective;
    }

    private void run() {
        computeSimilarities(); objective=0.0f;
        for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) voxelUpdate(x,y,z);
        // AlgorithmLesionToads.calcSegmentation: Eprev=energy; energy=computeMemberships();
        // Eratio=Numerics.abs(2.0f*(Eprev-energy)/(Eprev+energy)). Numerics.abs(float)
        // is the selected source's comparison-and-negation expression, not Math.abs.
        float difference=previousObjective-objective;
        float numerator=2.0f*difference;
        float denominator=previousObjective+objective;
        float quotient=numerator/denominator;
        energyRatio=quotient>0.0f?quotient:-quotient;
    }
    private static int countNaN(float[] values) { int n=0; for(float v:values) if(Float.isNaN(v)) n++; return n; }
    private void write(String path) throws Exception {
        int nan=countNaN(memberships)+countNaN(classEnergy)+countNaN(lesion)+countNaN(black)+countNaN(outlier)
            +countNaN(lesionRelation)+countNaN(denominator)+countNaN(lesionEnergy)+countNaN(blackHoleEnergy)
            +countNaN(runningObjective)+(Float.isNaN(objective)?1:0);
        try (DataOutputStream out=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)))) {
            out.writeInt(nx); out.writeInt(ny); out.writeInt(nz);
            for (float v:memberships) out.writeFloat(v);
            for (float v:classEnergy) out.writeFloat(v);
            for (float v:lesion) out.writeFloat(v);
            for (float v:black) out.writeFloat(v);
            for (float v:outlier) out.writeFloat(v);
            for (float v:lesionRelation) out.writeFloat(v);
            for (float v:denominator) out.writeFloat(v);
            for (float v:lesionEnergy) out.writeFloat(v);
            for (float v:blackHoleEnergy) out.writeFloat(v);
            for (float v:runningObjective) out.writeFloat(v);
            out.writeFloat(objective); out.writeInt(nan); out.writeFloat(energyRatio);
        }
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=2) throw new IllegalArgumentException("input.bin output.bin");
        LesionNextIterationMembershipReferenceHarness h=new LesionNextIterationMembershipReferenceHarness();
        h.read(args[0]); h.run(); h.write(args[1]);
    }
}
