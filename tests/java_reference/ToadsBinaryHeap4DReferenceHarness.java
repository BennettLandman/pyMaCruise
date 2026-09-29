import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Deterministic ordering oracle for the selected source BinaryHeap4D. */
public final class ToadsBinaryHeap4DReferenceHarness {
    private static final class Op {
        final String type; final float value; final int x,y,z,k;
        Op(String t,float v,int a,int b,int c,int d){type=t;value=v;x=a;y=b;z=c;k=d;}
    }
    private static Op add(float v,int x,int y,int z,int k){return new Op("add",v,x,y,z,k);}
    private static Op op(String t){return new Op(t,0,0,0,0,0);}
    private static String bits(float v){return Integer.toString(Float.floatToIntBits(v));}
    private static void run(String name,int capacity,String initial,List<Op> ops,boolean first){
        if(!first)System.out.print(',');
        BinaryHeap4D heap=new BinaryHeap4D(capacity,initial.equals("min")?BinaryHeap4D.MINTREE:BinaryHeap4D.MAXTREE);
        System.out.print("{\"name\":\""+name+"\",\"capacity\":"+capacity+",\"initial_mode\":\""+initial+"\",\"operations\":[");
        for(int i=0;i<ops.size();i++){
            Op o=ops.get(i);if(i>0)System.out.print(',');
            if(o.type.equals("add")){
                System.out.print("[\"add\","+Float.toString(o.value)+","+o.x+","+o.y+","+o.z+","+o.k+"]");
                heap.addValue(o.value,o.x,o.y,o.z,o.k);
            } else {
                System.out.print("[\""+o.type+"\"]");
                if(o.type.equals("reset"))heap.reset();
                else if(o.type.equals("set_min"))heap.setMinTree();
                else if(o.type.equals("set_max"))heap.setMaxTree();
            }
        }
        System.out.print("],\"popped\":[");
        boolean firstItem=true;
        while(heap.isNotEmpty()){
            if(!firstItem)System.out.print(',');firstItem=false;
            System.out.print("["+bits(heap.getFirst())+","+heap.getFirstX()+","+heap.getFirstY()+","+heap.getFirstZ()+","+heap.getFirstK()+"]");
            heap.removeFirst();
        }
        System.out.print("]}");
    }
    private static List<Op> duplicates(){
        List<Op> a=new ArrayList<>();
        a.add(add(2.0f,1,0,0,0));a.add(add(1.0f,2,0,0,1));a.add(add(1.0f,3,0,0,2));
        a.add(add(-0.0f,4,0,0,3));a.add(add(0.0f,5,0,0,4));a.add(add(-2.5f,6,0,0,5));
        a.add(add(2.0f,7,0,0,6));a.add(add(1.0f,8,0,0,7));a.add(add(3.25f,32768,-32769,65536,65535));
        return a;
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE BinaryHeap4D\",\"source_sha256\":\"9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a\",\"cases\":[");
        run("min_growth_ties",2,"min",duplicates(),true);
        run("max_growth_ties",2,"max",duplicates(),false);
        List<Op> switched=new ArrayList<>();switched.add(add(1,1,1,1,0));switched.add(add(4,2,2,2,1));switched.add(add(2,3,3,3,2));switched.add(op("set_max"));switched.add(add(5,4,4,4,3));switched.add(add(3,5,5,5,4));
        run("change_mode_with_items",8,"min",switched,false);
        List<Op> reset=new ArrayList<>(duplicates());reset.add(op("reset"));reset.add(op("set_min"));reset.add(add(7,9,8,7,6));reset.add(add(7,10,9,8,7));reset.add(add(-1,11,10,9,8));
        run("reset_then_reuse",2,"max",reset,false);
        System.out.println("]}");
    }
}
